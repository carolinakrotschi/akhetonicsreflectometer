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

  alias SCAN.json [SCAN_B.json]
      Standalone aliasing check (ported from the old check_aliasing.py
      tool; self-contained, uses its own load/spectrum functions below).
      One file: is it even a peak (width vs. resolution), and what aliasing
      order does the fringe-frequency chirp across sub-bands imply.
      Two files (same optics, DIFFERENT step size): the definitive pair
      test -- a peak that stays put (in cm) is real, one that moves is
      aliased.

All subcommands expect the EXFO JSON format (see process_reflectogram.py).
"""

import argparse
import json
import sys
from pathlib import Path

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


# ---------------------------------------------------------------- alias
# Self-contained: deliberately keeps its own load/spectrum functions rather
# than reusing load()/spectrum_diag() above (different envelope removal,
# Hann window, k-domain FFT) so the numbers this test reports never change
# just because the rest of the toolkit's pipeline evolves.
K_BANDS = 8


def load_alias(path):
    raw = json.load(open(path, encoding="utf-8"))
    block = raw["data"][0]
    info = (block.get("Device Description") or {}).get("info", "")
    arr = {k: np.asarray(v, float) for k, v in block.items()
           if isinstance(v, list) and v and isinstance(v[0], (int, float))}
    wl_key = next(k for k in arr if k.lower().startswith("wavelength"))
    wl = arr.pop(wl_key)

    step = np.median(np.diff(wl))
    keep = np.ones(wl.size, bool)
    for i in np.where(np.abs(np.diff(wl) - step) > 50 * abs(step))[0]:
        keep[i + 1] = False
    wl, arr = wl[keep], {k: v[keep] for k, v in arr.items()}

    ports = [k for k in arr if k.lower().startswith("ch")]
    if len(ports) >= 2 and np.corrcoef(arr[ports[0]], arr[ports[1]])[0, 1] < -0.3:
        a, b = arr[ports[0]], arr[ports[1]]
        sig, chan = (a - b) / (a + b), f"balanced {ports[0]}/{ports[1]}"
    else:
        sig, chan = arr[ports[0]], ports[0]
    return wl, sig, chan, info


def spectrum_alias(sig, lam):
    """Reflectogram: amplitude vs. reflector distance [cm]."""
    w = 501 | 1
    env = np.convolve(np.pad(sig, (w // 2, w // 2), mode="edge"),
                      np.ones(w) / w, mode="valid")[: sig.size]
    v = (sig - env) * np.hanning(sig.size)
    k = 1.0 / lam
    o = np.argsort(k)
    ku = np.linspace(k[o][0], k[o][-1], k.size)
    vu = np.interp(ku, k[o], v[o])
    S = np.abs(np.fft.rfft(vu))
    f = np.fft.rfftfreq(vu.size, d=ku[1] - ku[0])
    return f / NG * 1e-9 * 50, S          # cm, amplitude


def band_freq_alias(sig, d_lam, f_min):
    w = 301 | 1
    env = np.convolve(np.pad(sig, (w // 2, w // 2), mode="edge"),
                      np.ones(w) / w, mode="valid")[: sig.size]
    v = (sig - env) * np.hanning(sig.size)
    S = np.abs(np.fft.rfft(v))
    f = np.fft.rfftfreq(v.size, d=d_lam)
    m = f > f_min
    return float((f[m] * S[m] ** 2).sum() / (S[m] ** 2).sum())


def analyse_alias(path, verbose=True):
    wl, sig, chan, info = load_alias(path)
    lam = np.linspace(wl[0], wl[-1], wl.size)
    d_lam = lam[1] - lam[0]
    lam0, span = float(lam.mean()), float(lam[-1] - lam[0])

    B = C * span * 1e-9 / (lam0 * 1e-9) ** 2
    d_res_cm = C / (2 * NG * B) * 100
    d_max = lam0 ** 2 / (4 * NG * d_lam) * 1e-9 * 100

    d, S = spectrum_alias(sig, lam)
    P = S ** 2
    m = d > 0.3
    dd, PP = d[m], P[m]

    # Noise floor from the topmost quarter, which the hump usually does not reach
    floor = np.median(PP[dd > 0.85 * d_max])
    sig_only = np.maximum(PP - floor, 0)
    pk = float(dd[np.argmax(sig_only)])

    # Width at half power (only the hump, noise floor subtracted)
    thr = sig_only.max() / 2
    above = dd[sig_only > thr]
    fwhm = float(above.max() - above.min()) if above.size else float("nan")
    broadening = fwhm / d_res_cm

    # Drift of the fringe frequency
    f_app = band_freq_alias(sig, d_lam, f_min=0.05 / d_lam)
    n = sig.size
    L, F = [], []
    for j in range(K_BANDS):
        s = slice(j * n // K_BANDS, (j + 1) * n // K_BANDS)
        L.append(float(lam[s].mean()))
        F.append(band_freq_alias(sig[s], d_lam, f_min=0.3 * f_app))
    L, F = np.array(L), np.array(F)
    x = 1.0 / L ** 2
    slope = float(np.polyfit(x, F, 1)[0])
    r = float(np.corrcoef(x, F)[0, 1])

    drift_meas = abs(F[-1] - F[0]) / F.mean()
    drift_theo = 1.0 - (lam[0] / lam[-1]) ** 2      # ~6.25 % at 1520..1570 nm
    drift_ratio = drift_meas / drift_theo
    opd_chirp_cm = abs(slope) / NG * 1e-9 * 100

    res = dict(name=Path(path).name, info=info, chan=chan, n=lam.size,
               d_lam_pm=d_lam * 1e3, span=span, B_THz=B / 1e12,
               d_res_cm=d_res_cm, d_max=d_max, peak=pk, fwhm=fwhm,
               broadening=broadening, slope=slope, r=r,
               drift_ratio=drift_ratio, opd_chirp_cm=opd_chirp_cm, L=L, F=F,
               floor_snr=sig_only.max() / floor)

    if verbose:
        report_alias(res)
    return res


def report_alias(res):
    print(f"File           : {res['name']}")
    print(f"Info field     : {res['info']}")
    print(f"Channel        : {res['chan']}")
    print(f"Points         : {res['n']}     delta-lambda = {res['d_lam_pm']:.4f} pm")
    print(f"Bandwidth      : {res['span']:.3f} nm = {res['B_THz']:.3f} THz")
    print(f"Resolution     : {res['d_res_cm']*1e4:.1f} um")
    print(f"Nyquist limit  : {res['d_max']:.2f} cm")
    print()
    print(f"Maximum at     : {res['peak']:.2f} cm")
    print(f"Width (FWHM)   : {res['fwhm']:.2f} cm  =  {res['broadening']:.0f} x resolution")
    print(f"Hump / noise   : {res['floor_snr']:.1f} x")
    print()

    # --- Test 1: is it a peak at all? ---
    if res["broadening"] > 20:
        print("[X] NOT A PEAK. The hump is "
              f"{res['broadening']:.0f}x wider than the resolution.")
        print("    A length reading from this is meaningless, whether aliased or not.")
    elif res["broadening"] > 5:
        print(f"[!] Peak clearly broadened ({res['broadening']:.0f}x). Length is roughly usable,")
        print("    but you won't reach the 16 um resolution this way.")
    else:
        print(f"[OK] Clean peak ({res['broadening']:.0f}x resolution).")
    print()

    # --- Test 2: drift evaluation, only if permissible ---
    print("Drift of the fringe frequency across the sweep:")
    for l, f in zip(res["L"], res["F"]):
        print(f"   {l:8.2f} nm : {f:9.2f} cycles/nm")
    print(f"   Fit against 1/lambda^2: slope {res['slope']:+.4g}   (r = {res['r']:+.4f})")
    print(f"   Drift magnitude     : {res['drift_ratio']:.2f} x theory "
          f"(a single target must give 1.00)")
    print()

    if abs(res["r"]) < 0.8 or not (0.6 < res["drift_ratio"] < 1.7):
        print("[?] DRIFT EVALUATION NOT PERMISSIBLE.")
        print("    The drift magnitude does not match a single reflector -- likely")
        print("    several superimposed (aliased) components. I cannot tell from THIS file")
        print("    alone whether, or how many times, it is aliased.")
        print("    -> Run the pair test (see below).")
        return

    parity = "EVEN (0, 2, 4, ...)" if res["slope"] > 0 else "ODD (1, 3, ...)"
    print(f"Aliasing order is {parity}")
    print(f"OPD from the drift magnitude : {res['opd_chirp_cm']:.1f} cm  ->  "
          f"distance {res['opd_chirp_cm']/2:.1f} cm")
    if res["opd_chirp_cm"] / 2 > res["d_max"]:
        print(f"    {res['opd_chirp_cm']/2:.1f} cm > Nyquist {res['d_max']:.2f} cm  "
              f"==>  ALIASED.")
    elif res["slope"] > 0:
        print("    within Nyquist and even order  ==>  order 0, measurement valid.")
    else:
        print("    odd order  ==>  ALIASED.")


def pair_alias(a, b):
    ra, rb = analyse_alias(a, verbose=False), analyse_alias(b, verbose=False)
    print("=" * 68)
    print("PAIR TEST")
    print("=" * 68)
    for r in (ra, rb):
        print(f"  {r['name'][:44]:44s} dl={r['d_lam_pm']:6.3f} pm  "
              f"Nyquist {r['d_max']:6.2f} cm  Maximum {r['peak']:6.2f} cm")
    tol = 3 * max(ra["d_res_cm"], rb["d_res_cm"]) + 0.05 * max(ra["fwhm"], rb["fwhm"])
    diff = abs(ra["peak"] - rb["peak"])
    print(f"\n  Shift: {diff:.2f} cm   (tolerance {tol:.2f} cm)")
    if abs(ra["d_lam_pm"] - rb["d_lam_pm"]) < 1e-3:
        print("\n  [?] Both measurements have the same delta-lambda. The pair test needs")
        print("      two DIFFERENT step sizes of the same optics.")
        return
    if diff <= tol:
        print("\n  [OK] Peak stays put  ==>  NOT aliased. The length is real.")
    else:
        print("\n  [X] Peak moves  ==>  at least one of the two is ALIASED.")

        def cands(r, kmax=6):
            P = 2 * r["d_max"]
            out = set()
            for k in range(kmax):
                out.add(round(k * P + r["peak"], 2))
                if k >= 1:
                    out.add(round(k * P - r["peak"], 2))
            return sorted(x for x in out if x > 0)
        ca, cb = cands(ra), cands(rb)
        both = [(x, y) for x in ca for y in cb if abs(x - y) < 2.0]
        print("\n      True distances consistent with both measurements:")
        for x, y in both[:12]:
            print(f"        {(x+y)/2:8.2f} cm")
        print("\n      Warning: if one delta-lambda is an integer multiple of the")
        print("      other, several candidates remain -- this is mathematically")
        print("      unavoidable. Use a pair like 1 pm / 1.5 pm, or go finer.")


def cmd_alias(args):
    if args.scan_b is None:
        analyse_alias(args.scan_a)
    else:
        pair_alias(args.scan_a, args.scan_b)


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

    a = sub.add_parser("alias")
    a.add_argument("scan_a")
    a.add_argument("scan_b", nargs="?", default=None,
                    help="second scan (different step size, same optics) "
                         "for the pair test; omit for a single-file check")
    a.set_defaults(func=cmd_alias)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
