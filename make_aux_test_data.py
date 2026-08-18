#!/usr/bin/env python3
"""
Synthetischer 4-Kanal-Freilauf-Datensatz mit BEKANNTER Wahrheit.

Simuliert genau das, was der schwarze CoreDAQ nach dem Single-Trigger
aufnimmt: 4 Kanaele, gleichmaessig in ZEIT abgetastet (1 us), keine
Wellenlaengentabelle.

    Ch1 / Ch2 : Mess-Interferometer (komplementaer)
    Ch3 / Ch4 : Aux-MZI, das Lineal (komplementaer)

Eingebaut, weil es in den echten Daten auch drin ist:
  - Sweep-Nichtlinearitaet: langsamer Bogen (~59 pm p-p, wie gemessen)
  - schnelle Abstimmwelligkeit (Standard 20 MHz rms, siehe --ripple-mhz)
  - unterschiedliche Photodioden-Empfindlichkeit (Gain-Drift 0.84 -> 1.14)
  - DC-Sockel und Detektorrauschen

Sinn: die Pipeline (ofdr_aux.py) VOR dem Aufbau gegen bekannte Wahrheit
pruefen. Wenn sie die Reflektoren nicht auf die eingegebenen Positionen
legt, ist die Software schuld -- nicht der Tisch.

Beispiel:
    python make_aux_test_data.py --points 1000000 --speed 60 --dl 4 \
        --reflectors 0.046:0 0.30:-25 1.05:-30 2.00:-35 --out testdata.npz
"""

import argparse
import numpy as np

C = 299_792_458.0
NG = 1.468


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--points", type=int, default=1_000_000)
    p.add_argument("--step-us", type=float, default=1.0)
    p.add_argument("--speed", type=float, default=60.0, help="nm/s")
    p.add_argument("--lam0", type=float, default=1535.0, help="Startwellenlaenge nm")
    p.add_argument("--dl", type=float, default=4.0, help="Aux-Armdifferenz m")
    p.add_argument("--reflectors", nargs="+",
                   default=["0.046:0", "0.30:-25", "1.05:-30", "2.00:-35"],
                   help="Liste z_in_m:dB")
    p.add_argument("--bow-pm", type=float, default=59.0,
                   help="langsamer Bogen p-p in pm")
    p.add_argument("--ripple-mhz", type=float, default=20.0,
                   help="schnelle Welligkeit rms in MHz. ACHTUNG: das "
                        "Handover nennt 90-110 MHz, gemessen DURCH den "
                        "Hauptton (also moeglicherweise ueberschaetzt). Bei "
                        "~90 MHz mit 4.5 pm Periode wird nu(t) nicht mehr "
                        "monoton -- dann kann kein Aux mehr helfen. Genau "
                        "das solltest du mit dem echten Aux nachmessen.")
    p.add_argument("--noise", type=float, default=2e-4, help="Detektorrauschen")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--out", default="testdata.npz")
    a = p.parse_args()

    rng = np.random.default_rng(a.seed)
    n = a.points
    dt = a.step_us * 1e-6
    t = np.arange(n) * dt

    # --- wahre Wellenlaenge vs Zeit: linear + Bogen + Welligkeit -----------
    lam_lin = a.lam0 * 1e-9 + a.speed * 1e-9 * t
    bow = a.bow_pm * 1e-12 * (np.sin(np.pi * t / t[-1]) - 0.5)
    lam_mid = lam_lin[n // 2]
    ripple_m = a.ripple_mhz * 1e6 * lam_mid**2 / C            # MHz -> m
    # schnelle Welligkeit: Bandbegrenzung auf Perioden 2.5 - 30 pm des
    # Sweeps, so wie in den echten Daten gemessen
    dlam_per_sample = a.speed * 1e-9 * dt                     # m pro Punkt
    per_lo_samp = 2.5e-12 / dlam_per_sample
    per_hi_samp = 30e-12 / dlam_per_sample
    spec = np.zeros(n // 2 + 1, complex)
    ks = np.arange(n // 2 + 1)
    per_samp = np.divide(n, np.maximum(ks, 1e-9))
    band = (per_samp > per_lo_samp) & (per_samp < per_hi_samp)
    spec[band] = np.exp(2j * np.pi * rng.random(band.sum()))
    ripple = np.fft.irfft(spec, n)
    ripple *= ripple_m / (ripple.std() + 1e-30)
    lam = lam_lin + bow + ripple
    nu = C / lam

    # --- Monotonie-Pruefung -----------------------------------------------
    # Wenn die Steigung der Welligkeit die nominelle Sweep-Steigung
    # uebersteigt, faehrt der Laser momentan RUECKWAERTS in nu. Dann ist die
    # Aux-Phase nicht mehr monoton und keine Software kann das retten.
    dnu = np.diff(nu)
    frac_back = float((dnu * np.sign(dnu[0]) < 0).mean())
    a_max = C / lam_mid**2 * 2.5e-12 / (2 * np.pi)     # kuerzeste Periode 2.5 pm
    print(f"  Welligkeit {a.ripple_mhz:.0f} MHz rms "
          f"(Monotonie-Grenze bei dieser Bandbreite: ~{a_max/1e6:.0f} MHz)")
    print(f"  Monotonie: {frac_back*100:.3f} % der Schritte laufen rueckwaerts"
          + ("  <-- PROBLEM: keine Software kann das retten"
             if frac_back > 0 else "  (ok)"))

    # --- Signale ----------------------------------------------------------
    refl = []
    for s in a.reflectors:
        z, db = s.split(":")
        refl.append((float(z), float(db)))

    meas = np.zeros(n)
    for z, db in refl:
        tau = 2 * NG * z / C                      # Hin- und Rueckweg
        amp = 10 ** (db / 20.0)
        ph = rng.random() * 2 * np.pi
        meas += amp * np.cos(2 * np.pi * tau * nu + ph)
    meas /= max(1.0, max(10 ** (d / 20.0) for _, d in refl))

    tau_aux = NG * a.dl / C                       # Transmission: nur ein Weg
    aux = np.cos(2 * np.pi * tau_aux * nu + rng.random() * 2 * np.pi)

    # --- Detektormodell: DC-Sockel, komplementaer, Gain-Drift, Rauschen ---
    gain = np.linspace(0.84, 1.14, n)
    vis_m, vis_a = 0.75, 0.9
    ch1 = 0.5 * (1 + vis_m * meas)
    ch2 = 0.5 * (1 - vis_m * meas) / gain
    ch3 = 0.5 * (1 + vis_a * aux)
    ch4 = 0.5 * (1 - vis_a * aux) / gain
    for arr in (ch1, ch2, ch3, ch4):
        arr += rng.normal(0, a.noise, n)

    span = a.speed * (n * dt)
    np.savez_compressed(
        a.out, ch1=ch1, ch2=ch2, ch3=ch3, ch4=ch4,
        step_us=a.step_us, lam0_nm=a.lam0, speed_nms=a.speed,
        truth_z=np.array([z for z, _ in refl]),
        truth_db=np.array([d for _, d in refl]),
        truth_tau_aux=tau_aux, truth_dl=a.dl)

    print(f"geschrieben: {a.out}")
    print(f"  {n:,} Punkte, {a.step_us} us Takt -> {n*dt:.3f} s")
    print(f"  Sweep {a.speed} nm/s ab {a.lam0} nm -> Span {span:.2f} nm")
    print(f"  Aux dL = {a.dl} m -> tau_aux = {tau_aux*1e9:.3f} ns "
          f"(erscheint bei z = {a.dl/2:.3f} m)")
    print("  WAHRHEIT (Reflektoren):")
    for z, db in refl:
        print(f"    z = {z:7.4f} m   {db:6.1f} dB   "
              f"tau = {2*NG*z/C*1e9:7.3f} ns")


if __name__ == "__main__":
    main()
