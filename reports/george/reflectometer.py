#!/usr/bin/env python3
"""
OFDR reflectometer processing.

Input : two channels of (wavelength, power) from the complementary outputs
        of the measurement combining coupler, sampled on the EXFO's
        per-trigger measured-wavelength table.
Output: reflectogram (reflection amplitude in dB vs fiber distance),
        as a plot and a CSV.

Pipeline (in order, each step matters):
  1. load both channels, align, sanity-check
  2. balanced subtraction in software: P = P1 - g(lambda)*P2
  3. convert wavelength -> optical frequency, enforce monotonicity
  4. resample onto a uniform frequency grid (PCHIP)
  5. remove residual baseline
  6. window
  7. FFT -> delay axis -> distance axis
  8. optional dispersion autofocus (quadratic spectral phase)

Usage:
  python reflectometer.py ch1.csv ch2.csv [options]
  python reflectometer.py ch1.csv            (single-ended: skips subtraction)

Options:
  --ng VALUE        group index (default 1.468 for SMF-28 near 1550 nm)
  --window NAME     hann | blackmanharris | kaiser  (default hann)
  --kaiser-beta B   beta for kaiser window (default 12)
  --autofocus       enable dispersion autofocus on the strongest peak
  --zmax METERS     limit plotted/saved range (default: full Nyquist range)
  --delimiter D     CSV delimiter (default: auto-sniff)
  --wl-unit UNIT    nm | m | um (default nm)
  --out PREFIX      output file prefix (default: reflectogram)

Expected file format: two columns, wavelength then power, header rows fine
(anything non-numeric is skipped). Wavelengths in nm unless --wl-unit says
otherwise. Power in linear units (W or mW, not dBm) -- if your export is in
dBm, convert first or the subtraction step is meaningless.
"""

import argparse
import csv
import io
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator

C = 299_792_458.0  # m/s


# ----------------------------------------------------------------------
# 1. loading
# ----------------------------------------------------------------------

def load_channel(path, delimiter=None, wl_unit="nm"):
    """Read a two-column (wavelength, power) file, skipping non-numeric rows."""
    with open(path, "r", errors="replace") as f:
        text = f.read()

    if delimiter is None:
        try:
            delimiter = csv.Sniffer().sniff(text[:4096], delimiters=",;\t ").delimiter
        except csv.Error:
            delimiter = ","

    wl, p = [], []
    for row in csv.reader(io.StringIO(text), delimiter=delimiter):
        row = [c for c in row if c.strip() != ""]
        if len(row) < 2:
            continue
        try:
            wl.append(float(row[0]))
            p.append(float(row[1]))
        except ValueError:
            continue  # header or comment line

    if len(wl) < 100:
        sys.exit(f"{path}: only {len(wl)} numeric rows parsed -- wrong "
                 f"delimiter or format? First 200 chars:\n{text[:200]}")

    wl = np.asarray(wl, dtype=float)
    p = np.asarray(p, dtype=float)

    scale = {"nm": 1e-9, "um": 1e-6, "m": 1.0}[wl_unit]
    return wl * scale, p


def align_channels(wl1, p1, wl2, p2):
    """Both channels come from the same trigger train and should share the
    wavelength table. Verify, and align by index if lengths differ by a
    sample or two (dropped first/last trigger)."""
    n = min(len(wl1), len(wl2))
    if abs(len(wl1) - len(wl2)) > 2:
        sys.exit(f"Channel lengths differ by {abs(len(wl1)-len(wl2))} samples "
                 f"({len(wl1)} vs {len(wl2)}) -- these are not the same sweep.")
    wl1, p1, wl2, p2 = wl1[:n], p1[:n], wl2[:n], p2[:n]

    dwl = np.max(np.abs(wl1 - wl2))
    if dwl > 0.5e-12:  # half a picometre
        print(f"  warning: wavelength tables differ by up to {dwl*1e12:.2f} pm "
              f"between channels; using channel 1's table.", file=sys.stderr)
    return wl1, p1, p2


# ----------------------------------------------------------------------
# 2. balanced subtraction in software
# ----------------------------------------------------------------------

def balanced_subtract(wl, p1, p2, nseg=32):
    """P = P1 - g(lambda) * P2 with g fitted piecewise so the pedestal
    (the slowly-varying common-mode part) cancels even though the two
    photodiode responsivities drift differently across 130 nm.

    g is estimated per segment as the ratio of smoothed powers, then
    interpolated smoothly. The fringes average out in the smoothing, so
    g tracks the pedestal ratio, not the interference term."""
    n = len(p1)
    seg = max(n // nseg, 256)
    centers, gains = [], []
    for i in range(0, n - seg + 1, seg):
        s1 = np.median(p1[i:i + seg])
        s2 = np.median(p2[i:i + seg])
        if s2 > 0:
            centers.append(0.5 * (wl[i] + wl[i + seg - 1]))
            gains.append(s1 / s2)
    if len(gains) < 4:
        g = np.full(n, np.median(p1) / np.median(p2))
    else:
        g = PchipInterpolator(centers, gains, extrapolate=True)(wl)
    return p1 - g * p2, g


# ----------------------------------------------------------------------
# 3-4. frequency axis and uniform resampling
# ----------------------------------------------------------------------

def to_uniform_frequency(wl, p):
    """Convert to optical frequency, enforce strict monotonicity, resample
    onto a uniform frequency grid with PCHIP (no overshoot at sharp
    features, unlike cubic splines)."""
    nu = C / wl

    order = np.argsort(nu)
    nu, p = nu[order], p[order]

    # collapse duplicate frequencies (quantized wavelength table entries)
    uniq, inv, counts = np.unique(nu, return_inverse=True, return_counts=True)
    if len(uniq) < len(nu):
        psum = np.zeros(len(uniq))
        np.add.at(psum, inv, p)
        p = psum / counts
        nu = uniq
        print(f"  note: {counts.sum() - len(uniq)} duplicate wavelength "
              f"entries averaged.", file=sys.stderr)

    n = len(nu)
    nu_u = np.linspace(nu[0], nu[-1], n)
    p_u = PchipInterpolator(nu, p)(nu_u)
    return nu_u, p_u


# ----------------------------------------------------------------------
# 5-7. baseline, window, transform
# ----------------------------------------------------------------------

def remove_baseline(nu, p, order=3):
    """Remove the slowly-varying residual pedestal (coupler ratio drift,
    imperfect subtraction). A low-order polynomial only touches the
    near-DC bins; reflector fringes are far above its bandwidth."""
    x = (nu - nu.mean()) / (np.ptp(nu) / 2)  # conditioning
    coef = np.polyfit(x, p, order)
    return p - np.polyval(coef, x)


def make_window(n, name, kaiser_beta):
    if name == "hann":
        return np.hanning(n)
    if name == "blackmanharris":
        from scipy.signal.windows import blackmanharris
        return blackmanharris(n)
    if name == "kaiser":
        return np.kaiser(n, kaiser_beta)
    sys.exit(f"unknown window: {name}")


def reflectogram(nu_u, p_u, ng):
    """FFT to delay domain; return one-sided distance axis and complex
    spectrum. Distance is one-way fiber distance: the factor 2 for the
    round trip is in the denominator."""
    n = len(p_u)
    dnu = nu_u[1] - nu_u[0]
    R = np.fft.rfft(p_u)
    tau = np.arange(len(R)) / (n * dnu)   # delay axis, seconds (round trip)
    z = C * tau / (2.0 * ng)              # one-way distance in fiber, metres
    return z, R


# ----------------------------------------------------------------------
# 8. dispersion autofocus
# ----------------------------------------------------------------------

def autofocus(nu_u, p_win, ng, verbose=True):
    """Find the quadratic spectral phase beta that sharpens the strongest
    peak: multiply the analytic signal by exp(-i*beta*(nu-nu0)^2), FFT,
    and maximize the peak amplitude (a sharper peak concentrates energy).
    Golden-section search over beta; returns corrected spectrum's beta."""
    n = len(p_win)
    x = nu_u - nu_u.mean()

    # analytic signal so the phase multiplication is well defined
    spec = np.fft.fft(p_win)
    spec[n // 2 + 1:] = 0.0
    spec[1:n // 2] *= 2.0
    analytic = np.fft.ifft(spec)

    def peak_height(beta):
        corrected = analytic * np.exp(-1j * beta * x**2)
        return np.max(np.abs(np.fft.fft(corrected)[: n // 2]))

    # bracket: |beta| up to ~ (few ps of GDD spread) / span^2 -- generous
    span = x[-1] - x[0]
    bmax = 50e-12 / span**2 * (2 * np.pi)   # very loose upper bound
    betas = np.linspace(-bmax, bmax, 41)
    heights = [peak_height(b) for b in betas]
    b0 = betas[int(np.argmax(heights))]

    # golden-section refine around b0
    lo, hi = b0 - bmax / 20, b0 + bmax / 20
    gr = (np.sqrt(5) - 1) / 2
    for _ in range(40):
        c1, c2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
        if peak_height(c1) > peak_height(c2):
            hi = c2
        else:
            lo = c1
    beta = 0.5 * (lo + hi)
    if verbose:
        gain = peak_height(beta) / peak_height(0.0)
        print(f"  autofocus: beta = {beta:.4e} rad/Hz^2, "
              f"peak amplitude gain = {gain:.3f}x")
    corrected = analytic * np.exp(-1j * beta * x**2)
    return corrected, beta


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ch1", help="CSV of (wavelength, power) for output 1")
    ap.add_argument("ch2", nargs="?", default=None,
                    help="CSV for complementary output 2 (optional)")
    ap.add_argument("--ng", type=float, default=1.468)
    ap.add_argument("--window", default="hann",
                    choices=["hann", "blackmanharris", "kaiser"])
    ap.add_argument("--kaiser-beta", type=float, default=12.0)
    ap.add_argument("--autofocus", action="store_true")
    ap.add_argument("--zmax", type=float, default=None)
    ap.add_argument("--delimiter", default=None)
    ap.add_argument("--wl-unit", default="nm", choices=["nm", "um", "m"])
    ap.add_argument("--out", default="reflectogram")
    args = ap.parse_args()

    print(f"loading {args.ch1} ...")
    wl, p1 = load_channel(args.ch1, args.delimiter, args.wl_unit)
    print(f"  {len(wl)} samples, {wl.min()*1e9:.3f} to {wl.max()*1e9:.3f} nm")

    if args.ch2:
        print(f"loading {args.ch2} ...")
        wl2, p2 = load_channel(args.ch2, args.delimiter, args.wl_unit)
        wl, p1, p2 = align_channels(wl, p1, wl2, p2)
        p, g = balanced_subtract(wl, p1, p2)
        print(f"  balanced subtraction: gain ratio spans "
              f"{g.min():.4f} to {g.max():.4f}")
    else:
        p = p1 - np.median(p1)
        print("  single channel: median-subtracted only "
              "(pedestal and intensity noise NOT cancelled)")

    span_nm = (wl.max() - wl.min()) * 1e9
    nu_u, p_u = to_uniform_frequency(wl, p)
    dnu = nu_u[1] - nu_u[0]
    dz = C / (2 * args.ng * (nu_u[-1] - nu_u[0]))
    zmax_nyq = C / (4 * args.ng * dnu)
    print(f"  span {span_nm:.1f} nm = {(nu_u[-1]-nu_u[0])/1e12:.2f} THz, "
          f"bin {dnu/1e6:.1f} MHz")
    print(f"  resolution {dz*1e6:.2f} um, Nyquist range {zmax_nyq:.3f} m")

    p_u = remove_baseline(nu_u, p_u)
    w = make_window(len(p_u), args.window, args.kaiser_beta)
    p_win = p_u * w

    if args.autofocus:
        analytic, beta = autofocus(nu_u, p_win, args.ng)
        n = len(analytic)
        R = np.fft.fft(analytic)[: n // 2 + 1]
        tau = np.arange(len(R)) / (n * dnu)
        z = C * tau / (2 * args.ng)
    else:
        z, R = reflectogram(nu_u, p_win, args.ng)

    mag = np.abs(R)
    mag_db = 20 * np.log10(mag / mag.max() + 1e-15)

    if args.zmax is not None:
        keep = z <= args.zmax
        z, mag_db = z[keep], mag_db[keep]

    # save CSV
    out_csv = f"{args.out}.csv"
    np.savetxt(out_csv, np.column_stack([z, mag_db]),
               delimiter=",", header="distance_m,amplitude_dB", comments="")
    print(f"wrote {out_csv}")

    # peak report: everything within 40 dB of the top, simple local-max scan
    thr = mag_db.max() - 40
    peaks = []
    for i in range(2, len(mag_db) - 2):
        if (mag_db[i] > thr and mag_db[i] >= mag_db[i-1]
                and mag_db[i] >= mag_db[i+1]
                and mag_db[i] > mag_db[i-2] and mag_db[i] > mag_db[i+2]):
            peaks.append((z[i], mag_db[i]))
    peaks.sort(key=lambda t: -t[1])
    print("strongest features (top 10 within 40 dB of max):")
    for zz, aa in peaks[:10]:
        print(f"    {zz*100:8.4f} cm   {aa:7.2f} dB")

    # plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(z * 100, mag_db, lw=0.6)
        ax.set_xlabel("one-way fiber distance (cm)")
        ax.set_ylabel("reflection amplitude (dB, rel. max)")
        ax.set_title(f"reflectogram  |  {span_nm:.0f} nm span, "
                     f"{dz*1e6:.1f} um resolution, {args.window} window")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(max(-120, mag_db.min() - 5), 5)
        fig.tight_layout()
        out_png = f"{args.out}.png"
        fig.savefig(out_png, dpi=150)
        print(f"wrote {out_png}")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
