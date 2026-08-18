#!/usr/bin/env python3
"""
Acceptance test for the aux MZI (the "ruler"), in FREE-RUNNING mode.

Run this BEFORE trusting anything about the measurement. It only looks
at Ch3/Ch4 (the aux) and answers four questions:

  1. CONTRAST      Does the aux fringe cleanly across the whole sweep, or
                    does it drop out somewhere? (Dropout = polarization
                    fading)
  2. MONOTONICITY   Does the laser always tune forward? If it tunes
                    backward even briefly, the phase becomes ambiguous and
                    NO algorithm can recover from that. This is the
                    go/no-go criterion.
  3. CALIBRATION    What is tau_aux actually? Via fringe counting:
                      tau_aux = number of fringes / frequency span
                    The frequency span comes from the start and stop
                    wavelengths you set on the laser. Accuracy ~7 ppm.
  4. LASER QUALITY  How strong is the slow bow, how strong the fast
                    ripple? The aux measures this CLEANLY -- unlike
                    diagnose_artifacts.py tuning, which measures through
                    the main tone and gets contaminated by interfering
                    signals.

IMPORTANT -- why free-running and not trigger mode:
    An aux with a 4 m arm difference appears at z = 2.0 m. The old trigger
    mode (1 pm step) only reaches up to 0.42 m. A 4 m aux is NOT
    measurable at all in trigger mode -- it aliases with itself. In
    trigger mode, only an aux with dL < 0.83 m could be checked at best.
    Therefore: always acceptance-test the aux in free-running mode.

Examples:
    python check_aux_interferometer.py scan.json --lam-start 1505 --lam-stop 1565
    python check_aux_interferometer.py sim.npz  --lam-start 1505 --lam-stop 1565
    python check_aux_interferometer.py scan.json --tau-aux-ns 19.587      # without calibration
"""

import argparse
import json
import sys

import numpy as np
from scipy.ndimage import uniform_filter1d

C = 299_792_458.0
NG = 1.468


def load(path):
    if path.endswith(".npz"):
        d = np.load(path)
        step = float(d["step_us"]) if "step_us" in d else 1.0
        return d["ch3"], d["ch4"], step
    with open(path) as f:
        e = json.load(f)["data"][0]
    need = ["Ch3 [mW]", "Ch4 [mW]"]
    if any(k not in e for k in need):
        sys.exit(f"Ch3/Ch4 missing. Available keys: {list(e)}")
    return (np.asarray(e["Ch3 [mW]"], float),
            np.asarray(e["Ch4 [mW]"], float), 1.0)


def balanced(a, b, nseg=32):
    n = len(a)
    seg = max(n // nseg, 256)
    g = np.empty(n)
    for i in range(0, n, seg):
        s = slice(i, min(i + seg, n))
        m = np.median(b[s])
        g[s] = np.median(a[s]) / m if m > 0 else 1.0
    return a - uniform_filter1d(g, seg) * b


def analytic(x):
    n = len(x)
    X = np.fft.fft(x)
    X[n // 2 + 1:] = 0.0
    X[1:n // 2] *= 2.0
    return np.fft.ifft(X)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scan")
    p.add_argument("--lam-start", type=float, default=None,
                   help="start wavelength set on the laser, in nm")
    p.add_argument("--lam-stop", type=float, default=None,
                   help="stop wavelength set on the laser, in nm")
    p.add_argument("--tau-aux-ns", type=float, default=None,
                   help="instead of calibrating: provide a known tau_aux")
    p.add_argument("--trim", type=float, default=0.02,
                   help="fraction discarded at each edge "
                        "(laser start-up, Hilbert edge artifacts)")
    p.add_argument("--out", default=None)
    a = p.parse_args()

    ch3, ch4, step_us = load(a.scan)
    n = len(ch3)
    print(f"{a.scan}: {n:,} points, sample interval {step_us} us "
          f"-> acquisition duration {n*step_us*1e-6:.3f} s\n")

    aux = balanced(ch3, ch4)
    k = max(1, int(a.trim * n))
    sl = slice(k, n - k)

    an = analytic(aux)
    amp = np.abs(an)[sl]
    phi = np.unwrap(np.angle(an))[sl]
    if phi[-1] < phi[0]:
        phi = -phi
    m = len(phi)

    # ---------------------------------------------------------- 1. Contrast
    print("1) CONTRAST")
    w = max(1000, m // 200)
    hi = uniform_filter1d(np.maximum.accumulate(ch3[sl] * 0), w)  # placeholder
    # Visibility, windowed, on the RAW channel
    nw = 100
    edges = np.linspace(0, m, nw + 1).astype(int)
    vis = []
    raw = ch3[sl]
    for i in range(nw):
        s = raw[edges[i]:edges[i + 1]]
        if len(s) < 10:
            continue
        lo, hi_ = np.percentile(s, 1), np.percentile(s, 99)
        vis.append((hi_ - lo) / (hi_ + lo) if (hi_ + lo) > 0 else 0.0)
    vis = np.array(vis)
    print(f"   Visibility across the sweep: min {vis.min():.3f}, "
          f"median {np.median(vis):.3f}, max {vis.max():.3f}")
    print(f"   Beat amplitude: min {amp.min()/amp.mean()*100:.0f} % "
          f"of mean")
    if vis.min() < 0.2:
        print("   >> PROBLEM: contrast drops out somewhere. Probably "
              "polarization fading.\n"
              "      Fix: add a polarization controller in one arm, or use PM fiber.")
    elif vis.min() < 0.4:
        print("   >> borderline, but usable. Keep an eye on it.")
    else:
        print("   >> good.")

    # -------------------------------------------------------- 2. Monotonicity
    print("\n2) MONOTONICITY  (Go/No-Go)")
    dphi = np.diff(phi)
    frac_back = float((dphi <= 0).mean())
    p01 = np.percentile(dphi, 0.1)
    print(f"   mean phase step {dphi.mean():.4f} rad")
    print(f"   0.1 percentile of steps {p01:+.4f} rad "
          f"(= {p01/dphi.mean()*100:.0f} % of the mean)")
    print(f"   absolute smallest step {dphi.min():+.4f} rad "
          f"(isolated outliers are usually detector noise)")
    print(f"   fraction of backward-running steps: {frac_back*100:.4f} %")
    if frac_back > 1e-4:
        print("   >> STOP: the laser systematically runs backward at "
              "times.\n"
              "      The phase is then ambiguous -- no aux and no "
              "software can fix that.\n"
              "      Fix: find the cause on the laser side, or use a "
              "shorter aux (makes the phase slower-varying).\n"
              "      Sweeping slower does NOT help: the limit does not "
              "depend on speed.")
    elif frac_back > 0:
        print("   >> isolated backward steps, probably noise. "
              "Usable, but keep an eye on it.")
    elif p01 / dphi.mean() < 0.2:
        print("   >> tight, but monotonic. OK.")
    else:
        print("   >> good, clear margin.")

    # ----------------------------------------------------- 3. Calibration
    print("\n3) CALIBRATION of tau_aux")
    fringes = (phi[-1] - phi[0]) / (2 * np.pi)
    print(f"   fringes counted across the sweep: {fringes:,.1f}")
    print(f"   points per fringe: {m/fringes:.2f}"
          + ("   >> too few, aim for at least 4"
             if m / fringes < 4 else "   >> ok"))
    if a.tau_aux_ns is not None:
        tau_aux = a.tau_aux_ns * 1e-9
        print(f"   tau_aux given: {tau_aux*1e9:.4f} ns")
    elif a.lam_start and a.lam_stop:
        # The trimmed range only covers (1-2*trim) of the sweep
        dnu_full = abs(C / (a.lam_start * 1e-9) - C / (a.lam_stop * 1e-9))
        dnu = dnu_full * (m / n)
        tau_aux = fringes / dnu
        print(f"   frequency span {a.lam_start}->{a.lam_stop} nm = "
              f"{dnu_full/1e12:.4f} THz "
              f"(used fraction {m/n*100:.0f} % = {dnu/1e12:.4f} THz)")
        print(f"   tau_aux = fringes / span = {tau_aux*1e9:.4f} ns")
        print(f"   accuracy for a +-1 fringe counting error: "
              f"{1/fringes*1e6:.1f} ppm")
    else:
        sys.exit("   For calibration, provide --lam-start and --lam-stop "
                 "(or --tau-aux-ns).")
    dl = C * tau_aux / NG
    print(f"   corresponds to arm length difference dL = {dl:.4f} m")
    print(f"   the aux appears in the reflectogram at z = {dl/2:.4f} m")
    print(f"\n   --> USE THIS NUMBER GOING FORWARD:  "
          f"process_reflectogram_aux.py ... --tau-aux-ns {tau_aux*1e9:.4f}")

    # ------------------------------------------------------- 4. Laser quality
    print("\n4) LASER QUALITY, cleanly measured via the aux")
    if not (a.lam_start and a.lam_stop):
        print("   (skipped -- needs --lam-start and --lam-stop)")
        lam_mid, to_pm, resid, slow, fast = 1565e-9, 0, None, None, None
    else:
        # nu(t) from the aux phase, made absolute using the configured
        # start wavelength. Then convert back to lambda(t): the laser is
        # SUPPOSED to tune linearly in lambda, so the deviation from the
        # straight line is the tuning error. (Measuring the deviation
        # directly in nu would pick up the curvature of nu = c/lambda --
        # here that is ~75 GHz and has nothing to do with the laser.)
        lam_a = a.lam_start * 1e-9
        lam_b = a.lam_stop * 1e-9
        lam_trim_start = lam_a + a.trim * (lam_b - lam_a)
        nu_rel = phi / (2 * np.pi * tau_aux)
        sgn = 1.0 if lam_b < lam_a else -1.0   # nu decreases when lambda increases
        nu_abs = C / lam_trim_start + sgn * (nu_rel - nu_rel[0])
        lam_meas = C / nu_abs
        idx = np.arange(m)
        resid = lam_meas - np.polyval(np.polyfit(idx, lam_meas, 1), idx)
        slow = uniform_filter1d(resid, max(101, m // 60))
        fast = resid - uniform_filter1d(resid, max(31, m // 4000))
        lam_mid = (lam_a + lam_b) / 2
        to_pm = 1e12                              # resid is already in meters
        to_mhz = C / lam_mid**2 / 1e6             # m -> MHz
        print(f"   slow bow        : {np.ptp(slow)*1e12:.1f} pm p-p "
              f"= {np.ptp(slow)*to_mhz*1e0/1e3:.2f} GHz p-p")
        print(f"     (handover spec quotes ~59 pm -- compare order of magnitude)")
        print(f"   fast ripple     : {fast.std()*to_mhz:.1f} MHz rms "
              f"= {fast.std()*1e12:.3f} pm rms")
        fast = fast * to_mhz * 1e6                # from here on in Hz, as expected below
    dlam_per_sample = (abs(a.lam_stop - a.lam_start) * 1e-9 / n
                       if a.lam_start and a.lam_stop else None)
    if dlam_per_sample and fast is not None:
        F = np.abs(np.fft.rfft(fast - fast.mean()))
        per_samples = np.divide(m, np.maximum(np.arange(len(F)), 1e-9))
        per_pm = per_samples * dlam_per_sample * 1e12
        sel = (per_pm > 2.0) & (per_pm < 50.0)
        if sel.any():
            j = int(np.argmax(np.where(sel, F, 0)))
            dom = per_pm[j]
            print(f"   strongest ripple period: {dom:.2f} pm")
            limit = C / lam_mid**2 * dom * 1e-12 / (2 * np.pi)
            print(f"   monotonicity limit at this period: "
                  f"~{limit/1e6:.0f} MHz  (measured: {fast.std()/1e6:.1f} MHz)")
            if fast.std() > 0.7 * limit:
                print("   >> WARNING: close to the limit. Sweeping slower "
                      "does not help (the limit does not depend on speed) "
                      "-- but a shorter aux makes the unwrapping more "
                      "robust.")
            else:
                print("   >> not critical.")

    # ------------------------------------------------------------- Plot
    prefix = a.out or a.scan.rsplit(".", 1)[0] + "_auxcheck"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 2, figsize=(13, 7))
        nshow = int(np.clip(6 * m / fringes, 30, 600))
        z0 = slice(0, min(nshow, m))
        ax[0, 0].plot(ch3[sl][z0], ".-", ms=3, lw=0.8, label="Ch3")
        ax[0, 0].plot(ch4[sl][z0], ".-", ms=3, lw=0.8, label="Ch4")
        ax[0, 0].set_title(f"Raw aux signal, first {nshow} points "
                           "(must be a clean sine wave)")
        ax[0, 0].set_xlabel("Point"); ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

        ax[0, 1].plot(np.linspace(0, 100, len(vis)), vis, lw=1.2)
        ax[0, 1].axhline(0.2, ls="--", color="r", lw=1)
        ax[0, 1].set_ylim(0, 1.05)
        ax[0, 1].set_title("Contrast across the sweep "
                           "(red = alarm threshold)")
        ax[0, 1].set_xlabel("% of sweep"); ax[0, 1].grid(alpha=.3)

        ax[1, 0].plot(np.linspace(0, 100, m - 1), dphi / dphi.mean(), lw=.4)
        ax[1, 0].axhline(0, ls="--", color="r", lw=1)
        ax[1, 0].set_title("Phase step / mean  "
                           "(must NEVER go below 0 = red)")
        ax[1, 0].set_xlabel("% of sweep"); ax[1, 0].grid(alpha=.3)

        if resid is not None:
            ax[1, 1].plot(np.linspace(0, 100, m), resid * to_pm, lw=.4,
                          label="total")
            ax[1, 1].plot(np.linspace(0, 100, m), slow * to_pm, lw=1.5,
                          label="slow bow")
            ax[1, 1].legend()
        ax[1, 1].set_title("Laser tuning error, measured by the aux")
        ax[1, 1].set_xlabel("% of sweep")
        ax[1, 1].set_ylabel("pm"); ax[1, 1].grid(alpha=.3)
        fig.tight_layout()
        fig.savefig(f"{prefix}.png", dpi=140)
        print(f"\nwritten: {prefix}.png")
    except Exception as e:
        print(f"(Plot skipped: {e})")


if __name__ == "__main__":
    main()
