#!/usr/bin/env python3
"""
Rebuild the two-panel diagnostic figure (originally saved as
real_data_diagnosis.png, referenced in OFDR_Erklaerung_DE.md and
OFDR_Rueckfragen_DE.md) for an arbitrary scan.

Uses the exact pipeline functions from ofdr_process.py (load, balanced
subtraction, nominal-grid resample, phase_correct) and the tuning-error
residual from ofdr_diagnose.py -- nothing reimplemented, just assembled
into the same layout:

  left  : main MZI peak, "EXFO grid (as recorded)" (mode=none) vs
          "self-referenced phase" (mode=cosmetic), zoomed +-1 mm around
          the peak.
  right : sweep tuning error (pm-equivalent) vs wavelength, from the
          dominant tone's unwrapped Hilbert phase residual.

Usage: python ofdr_realdata_diagnosis.py scan.json [--out prefix]
"""

import argparse

import numpy as np
from scipy.signal.windows import kaiser

import ofdr_process as op

C = op.C
NG = op.NG


def spectrum(nu_g, y):
    n = len(y)
    R = np.abs(np.fft.rfft(y * kaiser(n, 12)))
    dnu = nu_g[1] - nu_g[0]
    z = np.arange(len(R)) * C / (2 * NG * dnu * n)
    return z, 20 * np.log10(R / R.max() + 1e-15)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan", help="EXFO JSON sweep file")
    ap.add_argument("--out", default=None, help="output prefix")
    args = ap.parse_args()

    wl, p1, p2, desc = op.load_exfo_json(args.scan)
    n = len(wl)
    p, g = op.balanced_subtract(p1, p2)
    nu_u, p_u = op.to_uniform_nu(wl, p)
    t = np.linspace(-1, 1, n)
    p_u = p_u - np.polyval(np.polyfit(t, p_u, 5), t)

    # left panel: no correction vs full self-referenced phase
    z_none, db_none = spectrum(*op.phase_correct(nu_u, p_u, "none"))
    z_cos, db_cos = spectrum(*op.phase_correct(nu_u, p_u, "cosmetic"))

    i_peak = np.argmax(db_cos)
    z_peak_mm = z_cos[i_peak] * 1000
    xlo, xhi = z_peak_mm - 1.0, z_peak_mm + 1.0
    imbalance_cm = 2 * z_peak_mm / 10

    # right panel: dominant-tone phase residual -> pm-equivalent tuning error
    an = op.analytic_signal(p_u)
    phi = np.unwrap(np.angle(an))
    tau0 = (phi[-1] - phi[0]) / (2 * np.pi * (nu_u[-1] - nu_u[0]))
    resid = phi - np.polyval(np.polyfit(np.arange(n), phi, 1), np.arange(n))
    pm_equiv = resid / (2 * np.pi * tau0) * wl**2 / C * 1e12

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(14, 5))

    mL = (z_none * 1000 >= xlo) & (z_none * 1000 <= xhi)
    mC = (z_cos * 1000 >= xlo) & (z_cos * 1000 <= xhi)
    axL.plot(z_none[mL] * 1000, db_none[mL], label="EXFO grid (as recorded)")
    axL.plot(z_cos[mC] * 1000, db_cos[mC], label="self-referenced phase")
    axL.set_xlim(xlo, xhi)
    axL.set_ylim(-60, 3)
    axL.set_xlabel("apparent distance (mm)")
    axL.set_ylabel("dB rel. max")
    axL.set_title(f"MZI peak, {imbalance_cm:.2f} cm imbalance")
    axL.legend()
    axL.grid(alpha=0.3)

    axR.plot(wl * 1e9, pm_equiv, lw=0.5)
    axR.set_xlabel("wavelength (nm)")
    axR.set_ylabel("tuning error (pm equiv.)")
    axR.set_title("sweep error not captured by 1 pm grid")
    axR.grid(alpha=0.3)

    fig.tight_layout()
    prefix = args.out or args.scan.rsplit(".", 1)[0]
    out_png = f"{prefix}_diagnosis.png"
    fig.savefig(out_png, dpi=150)
    print(f"wrote {out_png}")
    print(f"  main peak: {z_peak_mm:.3f} mm -> {imbalance_cm:.2f} cm imbalance")
    print(f"  tuning error p-p: {np.ptp(pm_equiv):.1f} pm")


if __name__ == "__main__":
    main()
