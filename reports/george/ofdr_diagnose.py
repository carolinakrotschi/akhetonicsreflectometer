#!/usr/bin/env python3
"""
OFDR diagnostics toolkit -- the discriminating tests from the 2026-08-12
artifact investigation, packaged. Run these BEFORE inventing explanations
for unexpected features; each test rejects a whole class of hypotheses.

Subcommands:

  pm-am SCAN.json --band-lo 3 --band-hi 6
      Additive vs multiplicative discriminator. Measures amplitude ripple
      and phase ripple of the dominant tone within a period band (in pm).
      ratio ~1  -> ADDITIVE optical path (real light via a real route:
                   reflection, multipath, aliased tone)
      ratio >>1 -> MULTIPLICATIVE phase noise (laser FM, acoustic pickup)
      History: this test killed the "laser FM" hypothesis (measured 1.04).

  compare SCAN_A.json SCAN_B.json [--z-lo 0.06 --z-hi 0.32]
      Two-scan envelope vs fine-structure correlation in a delay band.
      envelope ~+1, fine ~+1 -> FIXED artifact (etalon, connector pair):
                                calibratable, subtract a reference scan
      envelope ~+1, fine ~0  -> per-scan random realization (aliased tone
                                smeared by tuning noise; or SOP-dependent):
                                NOT calibratable, must be fixed optically
      History: measured +0.979 / +0.014 across the cap-removal pair.

  fold Z_TRUE [--step-pm 1.0 --span-nm 120 --lam-nm 1565]
      Where does a reflector at true distance Z_TRUE (metres, fiber) land
      after aliasing? Prints the folded apparent position and the fringe
      period vs the sampling step. Use to match a mystery band to a
      physical component beyond Nyquist.
      History: fold(1.046) -> 21 cm = the observed band. Case closed.

  tuning SCAN.json
      Extract the sweep's tuning error from the dominant tone: slow bow
      (rad p-p and pm-equivalent) and fast ripple (MHz rms, dominant
      period). Track these per session; they set how far the diagnostic
      correction can be trusted.

All subcommands expect the EXFO JSON format (see ofdr_process.py).
"""

import argparse
import json
import sys

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, filtfilt
from scipy.signal.windows import kaiser

C = 299_792_458.0
NG = 1.468
STEP_M = 1e-12


# ---------------------------------------------------------------- shared
def load(path):
    with open(path) as f:
        d = json.load(f)["data"][0]
    wl = np.asarray(d["Wavelength [nm]"], float) * 1e-9
    p1 = np.asarray(d["Ch1 [mW]"], float)
    p2 = np.asarray(d["Ch2 [mW]"], float)
    keep = np.ones(len(wl), bool)
    keep[1:] = np.diff(wl) > -10e-12
    return wl[keep], p1[keep], p2[keep]


def preprocess(wl, p1, p2):
    """Balanced subtraction + uniform-nu grid + baseline removal."""
    n = len(wl)
    seg = max(n // 32, 256)
    g = np.empty(n)
    for i in range(0, n, seg):
        s = slice(i, min(i + seg, n))
        m2 = np.median(p2[s])
        g[s] = np.median(p1[s]) / m2 if m2 > 0 else 1.0
    g = uniform_filter1d(g, seg)
    p = p1 - g * p2
    lam = wl[0] + np.arange(n) * STEP_M
    nu = C / lam
    nu_u = np.linspace(nu[-1], nu[0], n)
    p_u = PchipInterpolator(nu[::-1], p[::-1])(nu_u)
    t = np.linspace(-1, 1, n)
    return nu_u, p_u - np.polyval(np.polyfit(t, p_u, 5), t)


def analytic(x):
    n = len(x)
    X = np.fft.fft(x)
    X[n // 2 + 1:] = 0.0
    X[1:n // 2] *= 2.0
    return np.fft.ifft(X)


def spectrum_diag(nu_u, p_u, smooth=2001):
    """Diagnostic-grade (slow-only corrected) spectrum. Returns z, R."""
    n = len(p_u)
    an = analytic(p_u)
    phi = np.unwrap(np.angle(an))
    tau0 = (phi[-1] - phi[0]) / (2 * np.pi * (nu_u[-1] - nu_u[0]))
    phi_s = uniform_filter1d(phi, smooth)
    bow = phi_s - np.polyval(np.polyfit(np.arange(n), phi_s, 1), np.arange(n))
    nu_c = nu_u + bow / (2 * np.pi * tau0)
    o = np.argsort(nu_c)
    nu_g = np.linspace(nu_c[o][0], nu_c[o][-1], n)
    y = PchipInterpolator(nu_c[o], p_u[o])(nu_g)
    R = np.abs(np.fft.rfft(y * kaiser(n, 12)))
    z = np.arange(len(R)) * C / (2 * NG * (nu_g[1] - nu_g[0]) * n)
    return z, R


def bandpass_periods(x, per_lo, per_hi):
    """Bandpass keeping components with periods between per_lo and per_hi
    samples."""
    B, A = butter(3, [1 / per_hi / 0.5, 1 / per_lo / 0.5], btype="band")
    return filtfilt(B, A, x)


# ---------------------------------------------------------------- pm-am
def cmd_pm_am(args):
    wl, p1, p2 = load(args.scan)
    nu_u, p_u = preprocess(wl, p1, p2)
    n = len(p_u)
    an = analytic(p_u)
    A = np.abs(an)
    phi = np.unwrap(np.angle(an))
    am = bandpass_periods(A / uniform_filter1d(A, 2001),
                          args.band_lo, args.band_hi)
    pm = bandpass_periods(phi - uniform_filter1d(phi, 2001),
                          args.band_lo, args.band_hi)
    m = slice(n // 10, 9 * n // 10)
    ratio = pm[m].std() / am[m].std()
    print(f"band {args.band_lo}-{args.band_hi} pm periods:")
    print(f"  fractional AM ripple : {am[m].std()*100:.2f} %")
    print(f"  PM ripple            : {pm[m].std():.4f} rad")
    print(f"  PM/AM ratio          : {ratio:.2f}")
    verdict = ("ADDITIVE optical path(s): reflection, multipath, or an "
               "aliased tone. Hunt with `compare` and `fold`."
               if ratio < 3 else
               "MULTIPLICATIVE phase noise: laser tuning noise or acoustic "
               "path modulation. An aux-referenced axis will remove it.")
    print(f"  verdict: {verdict}")


# ---------------------------------------------------------------- compare
def cmd_compare(args):
    zA, RA = spectrum_diag(*preprocess(*load(args.scan_a)))
    zB, RB = spectrum_diag(*preprocess(*load(args.scan_b)))
    dbA = 20 * np.log10(RA / RA.max() + 1e-15)
    dbB = 20 * np.log10(RB / RB.max() + 1e-15)
    zg = np.linspace(args.z_lo, args.z_hi, 4000)
    bA = np.interp(zg, zA, dbA)
    bB = np.interp(zg, zB, dbB)
    envA, envB = uniform_filter1d(bA, 120), uniform_filter1d(bB, 120)
    fineA, fineB = bA - envA, bB - envB
    ce = np.corrcoef(envA, envB)[0, 1]
    cf = np.corrcoef(fineA, fineB)[0, 1]
    print(f"band {args.z_lo*100:.0f}-{args.z_hi*100:.0f} cm:")
    print(f"  A: median {np.median(bA):6.1f} dB, max {bA.max():6.1f} dB")
    print(f"  B: median {np.median(bB):6.1f} dB, max {bB.max():6.1f} dB")
    print(f"  envelope correlation      : {ce:+.3f}")
    print(f"  fine-structure correlation: {cf:+.3f}")
    if ce > 0.8 and cf > 0.8:
        print("  verdict: FIXED artifact -- calibratable; a reference scan "
              "subtracts it.")
    elif ce > 0.8:
        print("  verdict: stable envelope, random fine structure -- per-scan "
              "realization (aliased tone smeared by tuning noise, or "
              "SOP-dependent). NOT calibratable; fix optically. Use `fold` "
              "to test the aliasing candidate.")
    else:
        print("  verdict: band changed between scans -- something physical "
              "moved; re-run with an untouched setup.")


# ---------------------------------------------------------------- fold
def cmd_fold(args):
    lam = args.lam_nm * 1e-9
    dnu = C * args.step_pm * 1e-12 / lam**2
    z_nyq = C / (4 * NG * dnu)
    z = args.z_true
    period = 2 * z_nyq
    za = z % period
    if za > z_nyq:
        za = period - za
    fringe_pm = lam**2 / (2 * NG * z) * 1e12
    print(f"sampling: {args.step_pm} pm at {args.lam_nm} nm "
          f"-> dnu {dnu/1e6:.1f} MHz, Nyquist range {z_nyq:.3f} m")
    print(f"reflector at {z:.3f} m true:")
    if z <= z_nyq:
        print(f"  IN RANGE -- appears at its true position.")
    else:
        print(f"  BEYOND NYQUIST -- folds to {za:.3f} m "
              f"({za*100:.1f} cm) as a smeared band")
    print(f"  fringe period {fringe_pm:.2f} pm = "
          f"{fringe_pm/args.step_pm:.2f} points/fringe "
          f"({'OK' if fringe_pm/args.step_pm >= 2 else 'UNDERSAMPLED'})")


# ---------------------------------------------------------------- tuning
def cmd_tuning(args):
    wl, p1, p2 = load(args.scan)
    nu_u, p_u = preprocess(wl, p1, p2)
    n = len(p_u)
    an = analytic(p_u)
    phi = np.unwrap(np.angle(an))
    tau0 = (phi[-1] - phi[0]) / (2 * np.pi * (nu_u[-1] - nu_u[0]))
    resid = phi - np.polyval(np.polyfit(np.arange(n), phi, 1), np.arange(n))
    slow = uniform_filter1d(resid, 2001)
    fast = resid - uniform_filter1d(resid, 301)
    print(f"dominant tone: tau = {tau0*1e12:.2f} ps "
          f"(dL = {tau0*C/NG*1000:.2f} mm as MZI, "
          f"z = {tau0*C/NG/2*1000:.2f} mm as reflection)")
    lam_mid = 1.565e-6
    print(f"slow bow : {np.ptp(slow):.1f} rad p-p = "
          f"{np.ptp(slow)/(2*np.pi*tau0)/1e6:.0f} MHz p-p = "
          f"{np.ptp(slow)/(2*np.pi*tau0)*lam_mid**2/C*1e12:.1f} pm-equiv p-p")
    print(f"fast ripple (<~30 pm periods): "
          f"{fast.std()/(2*np.pi*tau0)/1e6:.0f} MHz rms")
    F = np.abs(np.fft.rfft(fast - fast.mean()))
    per = 1 / np.fft.rfftfreq(n, d=1.0)[1:]
    sel = (per > 2.5) & (per < 30)
    i = np.argmax(F[1:][sel])
    print(f"strongest fast component: period {per[sel][i]:.2f} pm")
    print("note: these numbers are measured THROUGH the dominant tone; "
          "additive interferers contaminate them. Trust trends, not the "
          "third digit.")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("pm-am")
    a.add_argument("scan")
    a.add_argument("--band-lo", type=float, default=3.0)
    a.add_argument("--band-hi", type=float, default=6.0)
    a.set_defaults(func=cmd_pm_am)

    a = sub.add_parser("compare")
    a.add_argument("scan_a")
    a.add_argument("scan_b")
    a.add_argument("--z-lo", type=float, default=0.06)
    a.add_argument("--z-hi", type=float, default=0.32)
    a.set_defaults(func=cmd_compare)

    a = sub.add_parser("fold")
    a.add_argument("z_true", type=float, help="true distance in metres")
    a.add_argument("--step-pm", type=float, default=1.0)
    a.add_argument("--lam-nm", type=float, default=1565.0)
    a.set_defaults(func=cmd_fold)

    a = sub.add_parser("tuning")
    a.add_argument("scan")
    a.set_defaults(func=cmd_tuning)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
