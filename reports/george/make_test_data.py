#!/usr/bin/env python3
"""Generate synthetic two-channel OFDR data mimicking the EXFO + CoreDAQ
setup: 130k triggers at 1 pm spacing 1500-1630 nm, wavelength table with
small measured jitter, complementary coupler outputs with unequal detector
gains, chromatic dispersion, RIN, and shot-ish noise.

Ground truth reflectors (one-way fiber distance, amplitude rel. LO field):
    5.00 cm   -30 dB   (front connector)
   17.30 cm   -45 dB   (mid-span splice)
   17.65 cm   -50 dB   (second splice 3.5 mm away -- resolution test)
   30.00 cm   -55 dB   (far event)
"""
import numpy as np

rng = np.random.default_rng(7)
C = 299_792_458.0
NG = 1.468

# ---- wavelength table: nominal 1 pm grid + measured jitter (0.05 pm rms),
# reported to 0.1 pm resolution like a real table
n = 130_000
lam_nom = 1500e-9 + np.arange(n) * 1e-12
lam_true = lam_nom + rng.normal(0, 0.05e-12, n)
lam_table = np.round(lam_true / 0.1e-12) * 0.1e-12   # what the EXFO reports

nu = C / lam_true
nu0 = nu.mean()

# ---- reflectors: (one-way distance m, field amplitude)
events = [(0.0500, 10**(-30/20)),
          (0.1730, 10**(-45/20)),
          (0.1765, 10**(-50/20)),
          (0.3000, 10**(-55/20))]

# dispersion: beta2 of SMF-28 ~ -2.17e-26 s^2/m; quadratic spectral phase
# for round trip 2z: phi = beta2 * (2*pi*(nu-nu0))^2 * (2z) / 2
beta2 = -2.17e-26

P_LO = 1.0   # normalize LO power to 1
fringe = np.zeros(n)
for z, r in events:
    tau = 2 * z * NG / C
    phi_disp = 0.5 * beta2 * (2 * np.pi * (nu - nu0))**2 * (2 * z)
    fringe += 2 * np.sqrt(P_LO) * r * np.cos(2 * np.pi * tau * nu + phi_disp)

P_sig_total = sum(r**2 for _, r in events)

# RIN: common-mode multiplicative noise on the laser, -110 dBc/Hz-ish
rin = 1 + rng.normal(0, 3e-4, n)

# complementary outputs, distinct detector gains with wavelength slope
g1 = 1.00 + 0.03 * (lam_true - 1565e-9) / 65e-9
g2 = 0.93 - 0.02 * (lam_true - 1565e-9) / 65e-9

pedestal = (P_LO + P_sig_total) * rin
ch1 = g1 * 0.5 * (pedestal + fringe * rin)
ch2 = g2 * 0.5 * (pedestal - fringe * rin)

# additive detector noise
ch1 += rng.normal(0, 2e-6, n)
ch2 += rng.normal(0, 2e-6, n)

np.savetxt("ch1.csv", np.column_stack([lam_table * 1e9, ch1]),
           delimiter=",", header="wavelength_nm,power_mW", comments="")
np.savetxt("ch2.csv", np.column_stack([lam_table * 1e9, ch2]),
           delimiter=",", header="wavelength_nm,power_mW", comments="")
print("wrote ch1.csv, ch2.csv")
print("ground truth (one-way cm, dB):")
for z, r in events:
    print(f"   {z*100:7.3f}   {20*np.log10(r):6.1f}")
