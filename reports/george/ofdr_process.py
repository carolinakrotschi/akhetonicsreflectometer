#!/usr/bin/env python3
"""
OFDR processing for EXFO + CoreDAQ sweep data.

Tailored to the JSON export format:
    {"header": {...},
     "data": [{"Device Description": {...},
               "Wavelength [nm]": [...],   # per-trigger MEASURED table, 1 pm quantized
               "Ch1 [mW]": [...],          # complementary MZI/coupler output 1
               "Ch2 [mW]": [...]}]}        # complementary output 2

Pipeline (each step validated on real data in the 2026-08-12 session):
  1. load + cleanup      : drop retrace samples, enforce in-span monotonicity
  2. balanced subtract   : P = Ch1 - g(lambda)*Ch2, g fitted per segment
  3. lambda -> nu        : nominal 1 pm grid (table is quantized to the step,
                           so the nominal grid is as good; see HANDOVER.md)
  4. slow-bow correction : self-referenced on the dominant tone's Hilbert
                           phase, HEAVILY SMOOTHED so only the laser's slow
                           sweep nonlinearity enters the axis. This preserves
                           real reflections (diagnostic-grade). See --mode.
  5. window + FFT        : Kaiser beta=12 default (sidelobes < -80 dB)
  6. report + plot       : peak list, width check vs transform limit

Modes (--mode):
  diagnostic (default) : slow-only correction. Real features stay sharp and
                         at true positions. USE THIS for identifying peaks.
  cosmetic             : raw-phase self-reference. Sharpest possible main
                         peak, but ABSORBS/distorts weak additive features.
                         Use only to demonstrate transform-limited width.
  none                 : nominal grid only (shows the uncorrected smearing).

Usage:
  python ofdr_process.py scan.json [--mode diagnostic] [--out prefix]
                                   [--zmax m] [--window hann|blackmanharris|kaiser]
"""

import argparse
import json
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks
from scipy.signal.windows import kaiser, hann, blackmanharris

C = 299_792_458.0
NG = 1.468          # SMF-28 group index near 1550 nm
STEP_M = 1e-12      # nominal trigger step: 1 pm


# ----------------------------------------------------------------------
def load_exfo_json(path):
    """Load one sweep; return (wavelength_m, ch1, ch2, description)."""
    with open(path) as f:
        doc = json.load(f)
    if len(doc["data"]) != 1:
        print(f"note: file contains {len(doc['data'])} entries; using entry 0",
              file=sys.stderr)
    d = doc["data"][0]
    wl = np.asarray(d["Wavelength [nm]"], float) * 1e-9
    p1 = np.asarray(d["Ch1 [mW]"], float)
    p2 = np.asarray(d["Ch2 [mW]"], float)
    desc = d.get("Device Description", {})

    # cleanup: drop retrace / parked samples (laser returns to start after
    # sweep; observed as a final sample ~120 nm below the previous one).
    # Filter by value, not position: any sample breaking monotonicity by
    # more than 10 pm is discarded.
    keep = np.ones(len(wl), bool)
    keep[1:] = np.diff(wl) > -10e-12
    ndrop = int((~keep).sum())
    if ndrop:
        print(f"  dropped {ndrop} retrace/non-monotonic sample(s)",
              file=sys.stderr)
    return wl[keep], p1[keep], p2[keep], desc


# ----------------------------------------------------------------------
def balanced_subtract(p1, p2, nseg=32):
    """P = Ch1 - g*Ch2 with g estimated per segment from the median (the
    fringes average out; g tracks the pedestal ratio, which drifts across
    130 nm because the two photodiode responsivities differ)."""
    n = len(p1)
    seg = max(n // nseg, 256)
    g = np.empty(n)
    for i in range(0, n, seg):
        s = slice(i, min(i + seg, n))
        m2 = np.median(p2[s])
        g[s] = np.median(p1[s]) / m2 if m2 > 0 else 1.0
    g = uniform_filter1d(g, seg)   # smooth the staircase
    return p1 - g * p2, g


# ----------------------------------------------------------------------
def to_uniform_nu(wl, p):
    """Nominal-grid lambda -> nu resample. The measured table is quantized
    to the 1 pm step, so reconstructing the grid as start + k*1pm is as
    accurate as the table itself and avoids duplicate-value handling."""
    n = len(wl)
    lam = wl[0] + np.arange(n) * STEP_M
    nu = C / lam
    nu_u = np.linspace(nu[-1], nu[0], n)          # ascending
    p_u = PchipInterpolator(nu[::-1], p[::-1])(nu_u)
    return nu_u, p_u


# ----------------------------------------------------------------------
def analytic_signal(x):
    n = len(x)
    X = np.fft.fft(x)
    X[n // 2 + 1:] = 0.0
    X[1:n // 2] *= 2.0
    return np.fft.ifft(X)


def phase_correct(nu_u, p_u, mode, smooth_samples=2001):
    """Correct the frequency axis using the dominant tone's own phase.

    diagnostic: only components slower than `smooth_samples` (default
                ~200 pm periods) enter the correction -> laser slow bow
                removed, everything else untouched.
    cosmetic  : full raw phase -> absorbs weak additive features. Only for
                demonstrating transform-limited width of the main tone.
    none      : no correction.

    Returns (nu_grid, corrected_signal).
    """
    if mode == "none":
        return nu_u, p_u

    an = analytic_signal(p_u)
    phi = np.unwrap(np.angle(an))
    n = len(p_u)
    tau0 = (phi[-1] - phi[0]) / (2 * np.pi * (nu_u[-1] - nu_u[0]))

    if mode == "diagnostic":
        phi_used = uniform_filter1d(phi, smooth_samples)
        bow = phi_used - np.polyval(
            np.polyfit(np.arange(n), phi_used, 1), np.arange(n))
        nu_corr = nu_u + bow / (2 * np.pi * tau0)
        order = np.argsort(nu_corr)
        nu_g = np.linspace(nu_corr[order][0], nu_corr[order][-1], n)
        y = PchipInterpolator(nu_corr[order], p_u[order])(nu_g)
        return nu_g, y

    if mode == "cosmetic":
        phi_m = np.maximum.accumulate(phi) + np.arange(n) * 1e-7
        phi_u = np.linspace(phi_m[0], phi_m[-1], n)
        y = PchipInterpolator(phi_m, np.real(an))(phi_u)
        dnu_eff = (phi_u[1] - phi_u[0]) / (2 * np.pi * tau0)
        nu_g = nu_u[0] + np.arange(n) * dnu_eff
        return nu_g, y

    sys.exit(f"unknown mode: {mode}")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan", help="EXFO JSON sweep file")
    ap.add_argument("--mode", default="diagnostic",
                    choices=["diagnostic", "cosmetic", "none"])
    ap.add_argument("--window", default="kaiser",
                    choices=["hann", "blackmanharris", "kaiser"])
    ap.add_argument("--kaiser-beta", type=float, default=12.0)
    ap.add_argument("--zmax", type=float, default=None,
                    help="limit plot/report range in metres")
    ap.add_argument("--peak-floor-db", type=float, default=-45.0,
                    help="report peaks above this level rel. max")
    ap.add_argument("--out", default=None, help="output prefix")
    args = ap.parse_args()

    print(f"loading {args.scan} ...")
    wl, p1, p2, desc = load_exfo_json(args.scan)
    n = len(wl)
    print(f"  {n} samples, {wl[0]*1e9:.3f} to {wl[-1]*1e9:.3f} nm  "
          f"(info: {desc.get('info','-')})")

    p, g = balanced_subtract(p1, p2)
    print(f"  balanced subtraction, gain ratio {g.min():.3f}..{g.max():.3f}")

    nu_u, p_u = to_uniform_nu(wl, p)
    t = np.linspace(-1, 1, n)
    p_u = p_u - np.polyval(np.polyfit(t, p_u, 5), t)   # residual baseline

    nu_g, y = phase_correct(nu_u, p_u, args.mode)
    dnu = nu_g[1] - nu_g[0]
    span = nu_g[-1] - nu_g[0]

    dz_bin = C / (2 * NG * span)
    z_nyq = C / (4 * NG * dnu)
    print(f"  span {span/1e12:.2f} THz -> resolution {dz_bin*1e6:.2f} um, "
          f"Nyquist range {z_nyq:.3f} m")
    print(f"  !! any reflector beyond {z_nyq:.2f} m ALIASES into range as a "
          f"smeared band -- see HANDOVER.md")

    win = {"hann": hann(n), "blackmanharris": blackmanharris(n),
           "kaiser": kaiser(n, args.kaiser_beta)}[args.window]
    R = np.abs(np.fft.rfft(y * win))
    z = np.arange(len(R)) * C / (2 * NG * dnu * n)
    db = 20 * np.log10(R / R.max() + 1e-15)

    # main peak + width check
    i = np.argmax(R)
    half = R[i] / np.sqrt(2)
    l = r = i
    while l > 0 and R[l] > half:
        l -= 1
    while r < len(R) - 1 and R[r] > half:
        r += 1
    width = (r - l) * (z[1] - z[0])
    wlim = {"hann": 1.6, "blackmanharris": 2.7, "kaiser": 2.6}[args.window]
    print(f"\nmain peak: {z[i]*1000:.3f} mm, -3 dB width {width*1e6:.1f} um "
          f"(window limit ~{wlim*dz_bin*1e6:.1f} um)")
    if width > 2 * wlim * dz_bin:
        print("  WARNING: main peak more than 2x the window limit -- axis "
              "correction may have regressed; investigate before trusting "
              "positions.")
    print(f"  delay {2*NG*z[i]/C*1e12:.2f} ps | as transmission-MZI "
          f"imbalance: {2*z[i]*1000:.2f} mm fiber | as reflection: "
          f"{z[i]*1000:.2f} mm fiber one-way")

    # peak report
    pk, _ = find_peaks(db, height=args.peak_floor_db,
                       distance=max(3, int(200e-6 / (z[1] - z[0]))))
    print(f"\npeaks above {args.peak_floor_db:.0f} dB "
          f"(min separation 0.2 mm -- lower it in code for finer lists):")
    zmax = args.zmax if args.zmax else z[-1]
    for j in pk:
        if z[j] <= zmax:
            harm = z[j] / z[i]
            tag = (f"   <- {harm:.0f}x main: processing/detector harmonic, "
                   f"not a reflector"
                   if 1.5 < harm < 6 and abs(harm - round(harm)) < 0.03
                   else "")
            print(f"   {z[j]*1000:9.3f} mm  {db[j]:6.1f} dB{tag}")

    prefix = args.out or args.scan.rsplit(".", 1)[0]
    keep = z <= zmax
    np.savetxt(f"{prefix}_reflectogram.csv",
               np.column_stack([z[keep], db[keep]]), delimiter=",",
               header="distance_m,amplitude_dB", comments="")
    print(f"\nwrote {prefix}_reflectogram.csv")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(z[keep] * 1000, db[keep], lw=0.6)
        ax.set_xlabel("distance (mm, fiber one-way / reflection convention)")
        ax.set_ylabel("amplitude (dB rel. max)")
        ax.set_title(f"{args.scan}  |  mode={args.mode}, {args.window}, "
                     f"res {dz_bin*1e6:.1f} um, Nyquist {z_nyq:.2f} m")
        ax.grid(alpha=0.3)
        ax.set_ylim(max(-110, db[keep].min() - 5), 5)
        fig.tight_layout()
        fig.savefig(f"{prefix}_reflectogram.png", dpi=150)
        print(f"wrote {prefix}_reflectogram.png")
    except Exception as e:
        print(f"(plot skipped: {e})")


if __name__ == "__main__":
    main()
