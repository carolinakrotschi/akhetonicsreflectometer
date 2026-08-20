#!/usr/bin/env python3
"""
OFDR processing for FREE-RUNNING mode with an aux MZI (4 channels).

This is the evolution of process_reflectogram.py for the new setup.
The difference in one sentence: the frequency axis no longer comes from
the EXFO wavelength table (which no longer exists in single-trigger
mode), but from the phase of the aux MZI.

    Ch1 / Ch3 : aux MZI = ruler   (complementary, both channels weak)
    Ch2 / Ch4 : measurement interferometer, with the test fiber
                (complementary, both channels strong)

    (Corrected 2026-08-19: the pairing is NOT consecutive-numbered
    (Ch1/Ch2, Ch3/Ch4) as an earlier draft and the 2026-08-18 handover
    text assumed -- it interleaves. Ch1 and Ch3 are both weak and close
    to each other in median power; Ch2 and Ch4 are both strong and close
    to each other in median power. That is the real complementary
    pairing. Default channel selection here reflects this; override with
    --aux-a/--aux-b/--meas-a/--meas-b if your wiring differs.)

What is carried over UNCHANGED from process_reflectogram.py:
    balanced subtraction, windowing, FFT, peak list, width check
What is DROPPED:
    to_uniform_nu()   -- no wavelength table anymore
    phase_correct()   -- no self-referencing, so also no
                         diagnostic/cosmetic conflict anymore
What is NEW:
    resample_on_aux() -- one function, ~25 lines

IMPORTANT: the aux phase is NOT smoothed. It has to carry the fast
tuning ripple, otherwise it won't correct it away.

Calibration of tau_aux (once, in the OLD trigger mode with
wavelength table):
    tau_aux = total unwrapped aux phase / (2*pi * frequency span)
Pass the result via --tau-aux-ns. Without it, it is estimated from --dl
(tau_aux = n_g * dL / c), which is accurate to ~1% for positions.

Examples:
    python process_reflectogram_aux.py testdata.npz --dl 4
    python process_reflectogram_aux.py scan.json --tau-aux-ns 19.587 --zmax 2.5
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
NG = 1.468


# ---------------------------------------------------------------- loading
def load(path):
    """-> (ch1, ch2, ch3, ch4, meta). Accepts .npz and .json."""
    if path.endswith(".npz"):
        d = np.load(path)
        meta = {k: d[k].item() for k in
                ("step_us", "lam0_nm", "speed_nms") if k in d}
        if "truth_z" in d:
            meta["truth_z"] = d["truth_z"]
            meta["truth_db"] = d["truth_db"]
        return (d["ch1"], d["ch2"], d["ch3"], d["ch4"], meta)

    with open(path) as f:
        doc = json.load(f)
    e = doc["data"][0]
    want = ["Ch1 [mW]", "Ch2 [mW]", "Ch3 [mW]", "Ch4 [mW]"]
    missing = [k for k in want if k not in e]
    if missing:
        sys.exit(f"Missing channels in file: {missing}\n"
                 f"available keys: {list(e)}")
    ch = [np.asarray(e[k], float) for k in want]
    meta = {k: v for k, v in (e.get("Device Description") or {}).items()}
    return (*ch, meta)


# ------------------------------------------------- balanced subtraction
def balanced(a, b, nseg=32):
    """P = a - g*b, g computed segment-wise via median (identical to
    process_reflectogram.balanced_subtract)."""
    n = len(a)
    seg = max(n // nseg, 256)
    g = np.empty(n)
    for i in range(0, n, seg):
        s = slice(i, min(i + seg, n))
        m = np.median(b[s])
        g[s] = np.median(a[s]) / m if m > 0 else 1.0
    g = uniform_filter1d(g, seg)
    return a - g * b, g


def analytic(x):
    n = len(x)
    X = np.fft.fft(x)
    X[n // 2 + 1:] = 0.0
    X[1:n // 2] *= 2.0
    return np.fft.ifft(X)


# ------------------------------------------------------ the new core
def resample_on_aux(meas, aux, tau_aux, trim=0.01):
    """Resample the measurement signal onto uniform aux phase.

    Core idea: the aux phase is phi(t) = 2*pi*tau_aux*nu(t) + const.
    Uniform steps in phi are therefore uniform steps in nu --
    no matter how non-uniformly the laser was actually swept. That makes
    the axis ready and the FFT valid.

    Returns: (y, dnu, span_nu, diag)
    """
    n = len(aux)
    an = analytic(aux)
    phi = np.unwrap(np.angle(an))
    if phi[-1] < phi[0]:                 # nu decreases when lambda increases
        phi = -phi

    # Trim the edges: the laser is still ramping up there, and the
    # Hilbert transform has edge artifacts.
    k = max(1, int(trim * n))
    sl = slice(k, n - k)
    phi = phi[sl]
    meas = meas[sl]
    amp = np.abs(an)[sl]

    # Enforce monotonicity (noise can cause tiny backward steps).
    # The step must be large compared to the numerical precision of phi,
    # otherwise duplicates remain -- hence relative to the mean step.
    eps = 1e-6 * (phi[-1] - phi[0]) / len(phi)
    phi = np.maximum.accumulate(phi) + np.arange(len(phi)) * eps

    m = len(phi)
    phi_u = np.linspace(phi[0], phi[-1], m)
    y = PchipInterpolator(phi, meas)(phi_u)

    dphi = phi_u[1] - phi_u[0]
    dnu = dphi / (2 * np.pi * tau_aux)
    span_nu = (phi[-1] - phi[0]) / (2 * np.pi * tau_aux)
    fringes = (phi[-1] - phi[0]) / (2 * np.pi)
    diag = dict(fringes=fringes, pts_per_fringe=m / fringes,
                amp_min=amp.min() / amp.mean(), n_used=m)
    return y, dnu, span_nu, diag


# ---------------------------------------------------------------- main
def build_argparser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scan")
    p.add_argument("--tau-aux-ns", type=float, default=None,
                   help="calibrated aux delay in ns (preferred)")
    p.add_argument("--dl", type=float, default=None,
                   help="aux arm length difference in m (fallback, if tau_aux "
                        "has not been calibrated yet)")
    p.add_argument("--window", default="kaiser",
                   choices=["hann", "blackmanharris", "kaiser"])
    p.add_argument("--kaiser-beta", type=float, default=12.0)
    p.add_argument("--zmax", type=float, default=None)
    p.add_argument("--peak-floor-db", type=float, default=-45.0)
    p.add_argument("--trim", type=float, default=0.01,
                   help="fraction discarded at each edge")
    p.add_argument("--aux-a", type=int, default=1, choices=[1, 2, 3, 4],
                   help="first channel of the aux (calibration) pair "
                        "(default 1: aux = Ch1/Ch3, the weak pair -- "
                        "confirmed 2026-08-19, see HANDOVER.md)")
    p.add_argument("--aux-b", type=int, default=3, choices=[1, 2, 3, 4],
                   help="second channel of the aux (calibration) pair "
                        "(default 3)")
    p.add_argument("--meas-a", type=int, default=2, choices=[1, 2, 3, 4],
                   help="first channel of the measurement pair "
                        "(default 2: measurement = Ch2/Ch4, the strong pair)")
    p.add_argument("--meas-b", type=int, default=4, choices=[1, 2, 3, 4],
                   help="second channel of the measurement pair (default 4)")
    p.add_argument("--single", action="store_true",
                   help="skip balanced subtraction; use only the stronger "
                        "channel of each pair (by median power), high-pass "
                        "filtered on its own. Use this when a pair is very "
                        "unbalanced (one channel near the noise floor) and "
                        "balanced subtraction may be injecting more noise "
                        "than it cancels.")
    p.add_argument("--out", default=None)
    return p


def main():
    a = build_argparser().parse_args()
    process(a)


def process(a):
    warnings = []

    if a.tau_aux_ns is not None:
        tau_aux = a.tau_aux_ns * 1e-9
        src = "calibrated"
    elif a.dl is not None:
        tau_aux = NG * a.dl / C
        src = f"estimated from dL = {a.dl} m"
    else:
        sys.exit("specify --tau-aux-ns or --dl")

    ch1_, ch2_, ch3_, ch4_, meta = load(a.scan)
    chans = {1: ch1_, 2: ch2_, 3: ch3_, 4: ch4_}
    ch1, ch2 = chans[a.aux_a], chans[a.aux_b]
    ch3, ch4 = chans[a.meas_a], chans[a.meas_b]
    n = len(ch1)
    print(f"{a.scan}: {n:,} points x 4 channels")
    print(f"aux (calibration) pair: Ch{a.aux_a}/Ch{a.aux_b}   "
          f"measurement pair: Ch{a.meas_a}/Ch{a.meas_b}")
    print(f"tau_aux = {tau_aux*1e9:.4f} ns  ({src})"
          f"  -> aux appears at z = {C*tau_aux/(2*NG):.3f} m")

    if a.single:
        aux_ch, aux_name = ((ch2, f"Ch{a.aux_b}") if np.median(ch2) > np.median(ch1)
                             else (ch1, f"Ch{a.aux_a}"))
        meas_ch, meas_name = ((ch4, f"Ch{a.meas_b}") if np.median(ch4) > np.median(ch3)
                               else (ch3, f"Ch{a.meas_a}"))
        print(f"single-channel mode: aux={aux_name} "
              f"(median {np.median(aux_ch):.4g} mW), "
              f"meas={meas_name} (median {np.median(meas_ch):.4g} mW) "
              "-- stronger of each pair, high-pass filtered individually")
        seg_a = max(len(aux_ch) // 32, 256)
        seg_m = max(len(meas_ch) // 32, 256)
        aux = aux_ch - uniform_filter1d(aux_ch, seg_a)
        meas = meas_ch - uniform_filter1d(meas_ch, seg_m)
    else:
        aux, ga = balanced(ch1, ch2)
        meas, gm = balanced(ch3, ch4)
        print(f"balanced subtraction: g_meas {gm.min():.3f}..{gm.max():.3f}, "
              f"g_aux {ga.min():.3f}..{ga.max():.3f}")

    y, dnu, span_nu, diag = resample_on_aux(meas, aux, tau_aux, a.trim)
    m = len(y)
    print(f"Aux: {diag['fringes']:,.0f} fringes, "
          f"{diag['pts_per_fringe']:.1f} points/fringe, "
          f"amplitude minimum {diag['amp_min']*100:.0f} % of mean")
    if diag["pts_per_fringe"] < 4:
        print("  WARNING: fewer than 4 points per aux fringe -- sweep slower "
              "or shorten the aux path.")
        warnings.append("low_pts_per_fringe")
    if diag["amp_min"] < 0.2:
        print("  WARNING: aux contrast drops somewhere (polarization "
              "fading?). Phase unreliable there.")
        warnings.append("low_aux_contrast")

    # Remove residual baseline (as in process_reflectogram.py)
    t = np.linspace(-1, 1, m)
    y = y - np.polyval(np.polyfit(t, y, 5), t)

    dz_bin = C / (2 * NG * span_nu)
    z_nyq = C / (4 * NG * dnu)
    lam_mid = 1565e-9
    print(f"Span {span_nu/1e12:.3f} THz "
          f"(~{span_nu*lam_mid**2/C*1e9:.1f} nm), point spacing "
          f"{dnu/1e6:.2f} MHz")
    print(f"  RESOLUTION {dz_bin*1e6:.2f} um   RANGE {z_nyq:.3f} m   "
          f"cells {z_nyq/dz_bin:,.0f}")
    print(f"  !! no reflector beyond {z_nyq:.2f} m, otherwise it will alias")

    win = {"hann": hann(m), "blackmanharris": blackmanharris(m),
           "kaiser": kaiser(m, a.kaiser_beta)}[a.window]
    R = np.abs(np.fft.rfft(y * win))
    z = np.arange(len(R)) * C / (2 * NG * dnu * m)
    db = 20 * np.log10(R / R.max() + 1e-15)

    i = int(np.argmax(R))
    half = R[i] / np.sqrt(2)
    l = r = i
    while l > 0 and R[l] > half:
        l -= 1
    while r < len(R) - 1 and R[r] > half:
        r += 1
    width = (r - l) * (z[1] - z[0])
    wlim = {"hann": 1.6, "blackmanharris": 2.7, "kaiser": 2.6}[a.window]
    print(f"\nMain peak {z[i]*1000:.4f} mm, -3 dB width {width*1e6:.1f} um "
          f"(window limit ~{wlim*dz_bin*1e6:.1f} um)")
    if width > 2 * wlim * dz_bin:
        print("  WARNING: main peak > 2x window limit -- frequency axis "
              "suspicious. Investigate before trusting positions.")
        warnings.append("main_peak_too_wide")

    zmax = a.zmax if a.zmax else z[-1]
    pk, _ = find_peaks(db, height=a.peak_floor_db,
                       distance=max(3, int(200e-6 / (z[1] - z[0]))))
    print(f"\nPeaks above {a.peak_floor_db:.0f} dB:")
    for j in pk:
        if z[j] <= zmax:
            h = z[j] / z[i] if z[i] > 0 else 0
            tag = ("   <- harmonic, not a reflector"
                   if 1.5 < h < 6 and abs(h - round(h)) < 0.03 else "")
            print(f"   {z[j]*1000:10.4f} mm  {db[j]:6.1f} dB{tag}")

    # Self-test against known truth (only for synthetic data)
    comparison = []
    if "truth_z" in meta:
        print("\n--- Comparison with known truth ---")
        worst = 0.0
        any_unaccounted = False
        for zt, dbt in zip(meta["truth_z"], meta["truth_db"]):
            # A reflector beyond the Nyquist range doesn't vanish -- it
            # aliases (folds back into [0, z_nyq], HANDOVER.md 5). Compare
            # against the folded position in that case, instead of silently
            # skipping it (which used to let a genuinely unaccounted-for
            # reflector print PASSED).
            aliased = zt > z_nyq
            if aliased:
                period = 2 * z_nyq
                zt_expect = zt % period
                if zt_expect > z_nyq:
                    zt_expect = period - zt_expect
            else:
                zt_expect = zt

            w = (z > zt_expect - 20 * dz_bin) & (z < zt_expect + 20 * dz_bin)
            if not w.any():
                any_unaccounted = True
                print(f"   z_true {zt:.4f} m  -- unaccounted for "
                      f"(expected {'fold at ' + format(zt_expect, '.4f') + ' m' if aliased else 'on-axis'}, nothing found there)")
                comparison.append(dict(z_true=zt, db_true=dbt, z_found=None,
                                        db_found=None, err_um=None, err_cells=None,
                                        status="out_of_range_unaccounted"))
                continue
            jj = int(np.argmax(np.where(w, R, 0)))
            err = (z[jj] - zt_expect) * 1e6
            worst = max(worst, abs(err))
            status = "aliased_as_expected" if aliased else "ok"
            tag = f"  ({status}, folded from z_true)" if aliased else ""
            print(f"   z_true {zt:7.4f} m ({dbt:5.1f} dB) -> found "
                  f"{z[jj]:7.4f} m, error {err:+7.1f} um "
                  f"({err/(dz_bin*1e6):+.2f} cells), {db[jj]:6.1f} dB{tag}")
            comparison.append(dict(z_true=zt, db_true=dbt, z_found=z[jj],
                                    db_found=db[jj], err_um=err,
                                    err_cells=err / (dz_bin * 1e6), status=status))
        print(f"   largest error: {worst:.1f} um "
              f"(one cell = {dz_bin*1e6:.1f} um)")
        passed = worst < 2 * dz_bin * 1e6 and not any_unaccounted
        print("   -> PASSED" if passed
              else "   -> FAILED, check pipeline")

    prefix = a.out or a.scan.rsplit(".", 1)[0]
    keep = z <= zmax
    np.savetxt(f"{prefix}_reflectogram.csv",
               np.column_stack([z[keep], db[keep]]), delimiter=",",
               header="distance_m,amplitude_dB", comments="")
    print(f"\nwrote: {prefix}_reflectogram.csv")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(11, 5))
        ax.plot(z[keep], db[keep], lw=0.6)
        ax.set_xlabel("Distance (m, one-way / reflection convention)")
        ax.set_ylabel("Amplitude (dB rel. maximum)")
        ax.set_title(f"{a.scan} | aux-referenced, {a.window}, "
                     f"dz {dz_bin*1e6:.1f} um, Nyquist {z_nyq:.2f} m")
        ax.grid(alpha=0.3)
        ax.set_ylim(max(-110, db[keep].min() - 5), 5)
        fig.tight_layout()
        fig.savefig(f"{prefix}_reflectogram.png", dpi=150)
        print(f"wrote: {prefix}_reflectogram.png")
    except Exception as e:
        print(f"(plot skipped: {e})")

    return dict(z=z, db=db, R=R, dz_bin=dz_bin, z_nyq=z_nyq,
                main_peak_m=z[i], peak_width_um=width * 1e6,
                warnings=warnings, comparison=comparison)


if __name__ == "__main__":
    main()
