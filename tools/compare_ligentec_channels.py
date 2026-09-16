#!/usr/bin/env python3
"""Mehrere Ligentec-Kanaele auf ihrer Facette A zentriert vergleichen.

Die Frage bei jedem neuen Kanal ist dieselbe: kommt hinter der Facette
etwas aus dem Chip zurueck? Ein Peak im Vorhersagefenster allein
beantwortet sie nicht -- in einem dichten Spektrum liegt fast ueberall
einer (Logbuch 2026-09-15). Dieses Werkzeug erzwingt deshalb zwei Dinge:

1. **Ausrichtung auf die eigene Facette A**, nicht auf die absolute
   Position. Die Arrayfasern sind unterschiedlich lang (fiber1 und
   fiber44 liegen 2.4 mm auseinander); wer absolut vergleicht, haelt
   diesen Faserversatz fuer ein Chipsignal.
2. **Nullverteilung.** Fuer JEDEN Scan wird berechnet, wie stark er die
   jeweils anderen in gleich breiten Fenstern uebertrifft. Ein Kandidat
   zaehlt nur, wenn sein Ueberschuss aus dieser Verteilung herausragt.
   Ein reiner Loopbackkanal wie fiber1, hinter dem nichts sein KANN,
   liefert den ehrlichen Massstab dafuer.

Aufruf
    python tools/compare_ligentec_channels.py \
        --scan 48 raw_data/2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48_reflectogram.csv \
        --scan 44 raw_data/2672_ligentechhi_2026-09-14-16-36_Ligentecfiber44_reflectogram.csv \
        --scan 1  raw_data/2672_ligentechhi_2026-09-14-16-02_Ligentecfiber1_reflectogram.csv \
        --control raw_data/2672_ligentechhi_2026-09-14-16-19_ligentecnofiber_reflectogram.csv
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ANCHOR_MM = 529.43       # interner Reflex des PM-Aufbaus, in allen Scans gleich
WIN_MM = 0.5             # Fensterbreite fuer die Nullverteilung
RANGE_MM = (0.2, 20.0)   # Suchbereich hinter Facette A


def load(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def facet_a(z, a, lo=800.0, hi=2000.0):
    m = (z > lo) & (z < hi)
    i = int(np.argmax(a[m]))
    return float(z[m][i]), float(a[m][i])


def anchor_level(z, a):
    m = (z > ANCHOR_MM - 0.05) & (z < ANCHOR_MM + 0.05)
    return float(a[m].max())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", nargs=2, action="append", required=True,
                    metavar=("FASER", "CSV"), help="Fasernummer und _reflectogram.csv")
    ap.add_argument("--control", default=None, help="Kontrollscan (nichts angeschlossen)")
    ap.add_argument("--predict", nargs=2, type=float, default=None,
                    metavar=("LO", "HI"), help="Vorhersagefenster hinter A in mm")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    S = {}
    for fib, path in a.scan:
        z, amp = load(path)
        za, aa = facet_a(z, amp)
        S["Faser %s" % fib] = dict(z=z, a=amp, zA=za, aA=aa,
                                   anch=anchor_level(z, amp))

    print("Facette A und Kopplung (Pegel relativ zum internen Anker %.2f mm):" % ANCHOR_MM)
    print("%-12s %12s %12s %10s" % ("Scan", "A [mm]", "A rel Anker", "Breite"))
    for k, s in S.items():
        m = (s["z"] > s["zA"] - 0.3) & (s["z"] < s["zA"] + 0.3)
        q, v = s["z"][m], s["a"][m]
        i = int(np.argmax(v)); h = v[i] - 3
        j = i
        while j > 0 and v[j] > h:
            j -= 1
        k2 = i
        while k2 < len(v) - 1 and v[k2] > h:
            k2 += 1
        print("%-12s %12.3f %12.2f %9.0f um" %
              (k, s["zA"], s["aA"] - s["anch"], (q[k2] - q[j]) * 1e3))

    # --- gemeinsames Raster relativ zu A, auf A normiert -----------------
    d = np.arange(RANGE_MM[0], RANGE_MM[1], 0.002)
    I = {k: np.interp(d, s["z"] - s["zA"], s["a"] - s["aA"]) for k, s in S.items()}

    if len(S) >= 2:
        print("\nNullverteilung: Ueberschuss jedes Scans gegenueber ALLEN anderen")
        print("(%.1f-mm-Fenster, auf die eigene Facette A normiert)" % WIN_MM)
        print("%-12s %8s %8s %8s %8s  %s" %
              ("Scan", "max", "95%", "Median", "min", "staerkste Stelle"))
        starts = np.arange(RANGE_MM[0] + 0.3, RANGE_MM[1] - WIN_MM, WIN_MM / 2)
        EX = {}
        for k in S:
            others = [o for o in S if o != k]
            ex = []
            for lo in starts:
                m = (d >= lo) & (d < lo + WIN_MM)
                ex.append(I[k][m].max() - max(I[o][m].max() for o in others))
            ex = np.array(ex)
            EX[k] = ex
            i = int(np.argmax(ex))
            print("%-12s %8.1f %8.1f %8.1f %8.1f  A+%.2f mm" %
                  (k, ex.max(), np.percentile(ex, 95), np.median(ex), ex.min(), starts[i]))
        print("\nLesart: ein Kandidat muss den max-Wert der ANDEREN Zeilen deutlich")
        print("uebertreffen. Tut er das nicht, ist er Instrumentenrock.")

    if a.predict:
        lo, hi = a.predict
        print("\nVorhersagefenster A+%.3f..%.3f mm (Pegel relativ zur eigenen Facette A):" % (lo, hi))
        for k in S:
            m = (d >= lo) & (d <= hi)
            print("  %-12s %7.1f dB" % (k, I[k][m].max()))

    # --- Plot -------------------------------------------------------------
    fig, ax = plt.subplots(2, 1, figsize=(11, 8), constrained_layout=True)
    for k in S:
        ax[0].plot(d, I[k], lw=0.7, label=k)
        ax[1].plot(d, I[k], lw=0.8, label=k)
    if a.control is not None:
        zc, ac = load(a.control)
        anc = anchor_level(zc, ac)
        ref = list(S.values())[0]
        off = ref["aA"] - ref["anch"]
        ic = np.interp(d, zc - 1582.0, ac - anc - off)
        for x in ax:
            x.plot(d, ic, lw=0.7, color="0.6", label="Kontrolle (nichts)")
    if a.predict:
        for x in ax:
            x.axvspan(a.predict[0], a.predict[1], color="tab:orange", alpha=0.18,
                      label="Vorhersage" if x is ax[0] else None)
    ax[0].set_xlim(*RANGE_MM); ax[0].set_ylim(-80, 5)
    ax[1].set_xlim(0.2, 3.0); ax[1].set_ylim(-60, 5)
    ax[0].set_title("Ligentec-Kanaele, auf die eigene Facette A zentriert und auf sie normiert")
    ax[1].set_title("Zoom: Nahbereich hinter A (Instrumenten-Seitenbaender)")
    for x in ax:
        x.set_xlabel("Abstand hinter Facette A [mm]")
        x.set_ylabel("Amplitude rel. A [dB]")
        x.grid(alpha=0.3); x.legend(fontsize=8, ncol=2)
    out = a.out or "results/2026-09-15/2672_ligentechhi_ligentec_kanalvergleich.png"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=140)
    print("\nPNG: %s" % out)


if __name__ == "__main__":
    main()
