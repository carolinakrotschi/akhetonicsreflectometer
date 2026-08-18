#!/usr/bin/env python3
"""
Abnahmepruefung fuer das Aux-MZI (das "Lineal"), im FREILAUF-Modus.

Das hier laeufst du, BEVOR du der Messung irgendetwas glaubst. Es schaut nur
auf Ch3/Ch4 (das Aux) und beantwortet vier Fragen:

  1. KONTRAST     Flackert das Aux sauber ueber den ganzen Sweep, oder bricht
                  es irgendwo ein? (Einbruch = Polarisationsfading)
  2. MONOTONIE    Faehrt der Laser immer vorwaerts? Wenn er auch nur kurz
                  rueckwaerts faehrt, ist die Phase mehrdeutig und KEIN
                  Verfahren kann das retten. Das ist das Go/No-Go-Kriterium.
  3. KALIBRIERUNG Wie gross ist tau_aux wirklich? Ueber Fringe-Zaehlen:
                      tau_aux = Anzahl Fringes / Frequenzspanne
                  Die Frequenzspanne kommt aus Start- und Endwellenlaenge,
                  die du am Laser einstellst. Genauigkeit ~7 ppm.
  4. LASERGUETE   Wie stark ist der langsame Bogen, wie stark das schnelle
                  Zittern? Das Aux misst das SAUBER -- anders als
                  ofdr_diagnose.py tuning, das durch den Hauptton misst und
                  von Stoersignalen kontaminiert wird.

WICHTIG -- warum Freilauf und nicht Triggermodus:
    Ein Aux mit 4 m Armdifferenz erscheint bei z = 2,0 m. Der alte
    Triggermodus (1 pm Schritt) reicht nur bis 0,42 m. Ein 4-m-Aux ist im
    Triggermodus GAR NICHT messbar -- es aliasiert selbst. Im Triggermodus
    liesse sich hoechstens ein Aux mit dL < 0,83 m pruefen.
    Deshalb: Aux immer im Freilauf abnehmen.

Beispiele:
    python aux_check.py scan.json --lam-start 1505 --lam-stop 1565
    python aux_check.py sim.npz  --lam-start 1505 --lam-stop 1565
    python aux_check.py scan.json --tau-aux-ns 19.587      # ohne Kalibrierung
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
        sys.exit(f"Ch3/Ch4 fehlen. Vorhandene Schluessel: {list(e)}")
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
                   help="am Laser eingestellte Startwellenlaenge in nm")
    p.add_argument("--lam-stop", type=float, default=None,
                   help="am Laser eingestellte Endwellenlaenge in nm")
    p.add_argument("--tau-aux-ns", type=float, default=None,
                   help="statt zu kalibrieren: bekanntes tau_aux vorgeben")
    p.add_argument("--trim", type=float, default=0.02,
                   help="Anteil, der an jedem Rand verworfen wird "
                        "(Laser-Anlauf, Hilbert-Kantenartefakte)")
    p.add_argument("--out", default=None)
    a = p.parse_args()

    ch3, ch4, step_us = load(a.scan)
    n = len(ch3)
    print(f"{a.scan}: {n:,} Punkte, Takt {step_us} us "
          f"-> Aufnahmedauer {n*step_us*1e-6:.3f} s\n")

    aux = balanced(ch3, ch4)
    k = max(1, int(a.trim * n))
    sl = slice(k, n - k)

    an = analytic(aux)
    amp = np.abs(an)[sl]
    phi = np.unwrap(np.angle(an))[sl]
    if phi[-1] < phi[0]:
        phi = -phi
    m = len(phi)

    # ---------------------------------------------------------- 1. Kontrast
    print("1) KONTRAST")
    w = max(1000, m // 200)
    hi = uniform_filter1d(np.maximum.accumulate(ch3[sl] * 0), w)  # Platzhalter
    # Sichtbarkeit fensterweise auf dem ROHEN Kanal
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
    print(f"   Sichtbarkeit ueber den Sweep: min {vis.min():.3f}, "
          f"median {np.median(vis):.3f}, max {vis.max():.3f}")
    print(f"   Amplitude der Schwebung: min {amp.min()/amp.mean()*100:.0f} % "
          f"vom Mittel")
    if vis.min() < 0.2:
        print("   >> PROBLEM: Kontrast bricht irgendwo ein. Wahrscheinlich "
              "Polarisationsfading.\n"
              "      Abhilfe: Polarisationsregler in einen Arm, oder PM-Faser.")
    elif vis.min() < 0.4:
        print("   >> grenzwertig, aber brauchbar. Im Auge behalten.")
    else:
        print("   >> gut.")

    # -------------------------------------------------------- 2. Monotonie
    print("\n2) MONOTONIE  (Go/No-Go)")
    dphi = np.diff(phi)
    frac_back = float((dphi <= 0).mean())
    p01 = np.percentile(dphi, 0.1)
    print(f"   mittlerer Phasenschritt {dphi.mean():.4f} rad")
    print(f"   0,1-Perzentil der Schritte {p01:+.4f} rad "
          f"(= {p01/dphi.mean()*100:.0f} % des mittleren)")
    print(f"   absolut kleinster Schritt {dphi.min():+.4f} rad "
          f"(einzelne Ausreisser sind meist Detektorrauschen)")
    print(f"   Anteil rueckwaerts laufender Schritte: {frac_back*100:.4f} %")
    if frac_back > 1e-4:
        print("   >> STOPP: der Laser faehrt systematisch zwischendurch "
              "rueckwaerts.\n"
              "      Die Phase ist dann mehrdeutig -- kein Aux und keine "
              "Software kann das reparieren.\n"
              "      Abhilfe: Ursache am Laser suchen, oder kuerzeres Aux "
              "(macht die Phase traeger).\n"
              "      Langsamer sweepen hilft NICHT: die Grenze haengt nicht "
              "von der Geschwindigkeit ab.")
    elif frac_back > 0:
        print("   >> vereinzelte Rueckwaertsschritte, vermutlich Rauschen. "
              "Brauchbar, aber im Auge behalten.")
    elif p01 / dphi.mean() < 0.2:
        print("   >> knapp, aber monoton. Geht.")
    else:
        print("   >> gut, deutliche Reserve.")

    # ----------------------------------------------------- 3. Kalibrierung
    print("\n3) KALIBRIERUNG von tau_aux")
    fringes = (phi[-1] - phi[0]) / (2 * np.pi)
    print(f"   gezaehlte Fringes ueber den Sweep: {fringes:,.1f}")
    print(f"   Punkte pro Fringe: {m/fringes:.2f}"
          + ("   >> zu wenig, mindestens 4 anstreben"
             if m / fringes < 4 else "   >> ok"))
    if a.tau_aux_ns is not None:
        tau_aux = a.tau_aux_ns * 1e-9
        print(f"   tau_aux vorgegeben: {tau_aux*1e9:.4f} ns")
    elif a.lam_start and a.lam_stop:
        # Der getrimmte Bereich deckt nur (1-2*trim) des Sweeps ab
        dnu_full = abs(C / (a.lam_start * 1e-9) - C / (a.lam_stop * 1e-9))
        dnu = dnu_full * (m / n)
        tau_aux = fringes / dnu
        print(f"   Frequenzspanne {a.lam_start}->{a.lam_stop} nm = "
              f"{dnu_full/1e12:.4f} THz "
              f"(genutzter Anteil {m/n*100:.0f} % = {dnu/1e12:.4f} THz)")
        print(f"   tau_aux = Fringes / Spanne = {tau_aux*1e9:.4f} ns")
        print(f"   Genauigkeit bei +-1 Fringe Zaehlfehler: "
              f"{1/fringes*1e6:.1f} ppm")
    else:
        sys.exit("   Fuer die Kalibrierung --lam-start und --lam-stop angeben "
                 "(oder --tau-aux-ns).")
    dl = C * tau_aux / NG
    print(f"   entspricht Armdifferenz dL = {dl:.4f} m")
    print(f"   das Aux erscheint im Reflektogramm bei z = {dl/2:.4f} m")
    print(f"\n   --> DIESE ZAHL WEITERVERWENDEN:  "
          f"ofdr_aux.py ... --tau-aux-ns {tau_aux*1e9:.4f}")

    # ------------------------------------------------------- 4. Laserguete
    print("\n4) LASERGUETE, sauber gemessen durch das Aux")
    if not (a.lam_start and a.lam_stop):
        print("   (uebersprungen -- braucht --lam-start und --lam-stop)")
        lam_mid, to_pm, resid, slow, fast = 1565e-9, 0, None, None, None
    else:
        # nu(t) aus der Aux-Phase, absolut gemacht ueber die eingestellte
        # Startwellenlaenge. Dann zurueck nach lambda(t): der Laser SOLL
        # linear in lambda fahren, also ist die Abweichung von der Geraden
        # der Abstimmfehler. (Wuerde man die Abweichung direkt in nu messen,
        # bekaeme man die Kruemmung von nu = c/lambda mit -- das sind hier
        # ~75 GHz und hat mit dem Laser nichts zu tun.)
        lam_a = a.lam_start * 1e-9
        lam_b = a.lam_stop * 1e-9
        lam_trim_start = lam_a + a.trim * (lam_b - lam_a)
        nu_rel = phi / (2 * np.pi * tau_aux)
        sgn = 1.0 if lam_b < lam_a else -1.0   # nu faellt, wenn lambda steigt
        nu_abs = C / lam_trim_start + sgn * (nu_rel - nu_rel[0])
        lam_meas = C / nu_abs
        idx = np.arange(m)
        resid = lam_meas - np.polyval(np.polyfit(idx, lam_meas, 1), idx)
        slow = uniform_filter1d(resid, max(101, m // 60))
        fast = resid - uniform_filter1d(resid, max(31, m // 4000))
        lam_mid = (lam_a + lam_b) / 2
        to_pm = 1e12                              # resid ist schon in Metern
        to_mhz = C / lam_mid**2 / 1e6             # m -> MHz
        print(f"   langsamer Bogen : {np.ptp(slow)*1e12:.1f} pm p-p "
              f"= {np.ptp(slow)*to_mhz*1e0/1e3:.2f} GHz p-p")
        print(f"     (Handover nennt ~59 pm -- Groessenordnung vergleichen)")
        print(f"   schnelles Zittern: {fast.std()*to_mhz:.1f} MHz rms "
              f"= {fast.std()*1e12:.3f} pm rms")
        fast = fast * to_mhz * 1e6                # ab hier in Hz, wie unten erwartet
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
            print(f"   staerkste Zitterperiode: {dom:.2f} pm")
            limit = C / lam_mid**2 * dom * 1e-12 / (2 * np.pi)
            print(f"   Monotonie-Grenze bei dieser Periode: "
                  f"~{limit/1e6:.0f} MHz  (gemessen: {fast.std()/1e6:.1f} MHz)")
            if fast.std() > 0.7 * limit:
                print("   >> WARNUNG: nah an der Grenze. Langsamer sweepen "
                      "gibt keine Hilfe (die Grenze haengt nicht von der "
                      "Geschwindigkeit ab) -- aber ein kuerzeres Aux macht "
                      "das Unwrapping robuster.")
            else:
                print("   >> unkritisch.")

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
        ax[0, 0].set_title(f"Rohsignal Aux, erste {nshow} Punkte "
                           "(muss ein sauberer Sinus sein)")
        ax[0, 0].set_xlabel("Punkt"); ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

        ax[0, 1].plot(np.linspace(0, 100, len(vis)), vis, lw=1.2)
        ax[0, 1].axhline(0.2, ls="--", color="r", lw=1)
        ax[0, 1].set_ylim(0, 1.05)
        ax[0, 1].set_title("Kontrast ueber den Sweep "
                           "(rot = Alarmgrenze)")
        ax[0, 1].set_xlabel("% des Sweeps"); ax[0, 1].grid(alpha=.3)

        ax[1, 0].plot(np.linspace(0, 100, m - 1), dphi / dphi.mean(), lw=.4)
        ax[1, 0].axhline(0, ls="--", color="r", lw=1)
        ax[1, 0].set_title("Phasenschritt / Mittelwert  "
                           "(darf NIE unter 0 = rot)")
        ax[1, 0].set_xlabel("% des Sweeps"); ax[1, 0].grid(alpha=.3)

        if resid is not None:
            ax[1, 1].plot(np.linspace(0, 100, m), resid * to_pm, lw=.4,
                          label="gesamt")
            ax[1, 1].plot(np.linspace(0, 100, m), slow * to_pm, lw=1.5,
                          label="langsamer Bogen")
            ax[1, 1].legend()
        ax[1, 1].set_title("Abstimmfehler des Lasers, vom Aux gemessen")
        ax[1, 1].set_xlabel("% des Sweeps")
        ax[1, 1].set_ylabel("pm"); ax[1, 1].grid(alpha=.3)
        fig.tight_layout()
        fig.savefig(f"{prefix}.png", dpi=140)
        print(f"\ngeschrieben: {prefix}.png")
    except Exception as e:
        print(f"(Plot uebersprungen: {e})")


if __name__ == "__main__":
    main()
