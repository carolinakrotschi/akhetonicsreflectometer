#!/usr/bin/env python3
"""Alle Scans des PM-Aufbaus auf EINER Achse ab 0 mm.

Gemeinsam ist allen Messungen der Zirkulator. Danach folgt entweder das
Faserarray und der Chip (Ligentec-Scans) oder einfach die Faser unter
Test (Faserscans). Um beides nebeneinander lesen zu koennen, braucht es
einen Nullpunkt, den alle Scans teilen.

WAHL DES NULLPUNKTS: der interne Anker bei 529.43 mm, NICHT der
Steckerreflex. Begruendung aus den Daten:

    Anker      529.38 .. 529.43 mm  ueber acht Scans  -> Streuung 50 um
    Stecker    555.87 .. 557.85 mm  ueber acht Scans  -> Streuung 1.5 mm

Der Stecker wird bei jedem Faserwechsel neu gesteckt und wandert
deshalb; in den Ligentec-Scans ist er nicht einmal eindeutig als
einzelner Peak auffindbar. Der Anker sitzt fest im Zirkulatorteil und
ist das einzige Merkmal, das wirklich alle Scans teilen.

Der Steckerbereich wird als Band eingezeichnet: alles rechts davon ist
DUT, alles links davon Instrument.

Darstellung als Small Multiples (eine Zeile pro Scan, gemeinsame
x-Achse) statt als Ueberlagerung -- sechs dichte Reflektogramme
uebereinander ergeben nur einen Rauschblock.

Aufruf
    python tools/plot_overview_anchor_axis.py --out results/<datum>/uebersicht.png
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a85"
C_CHIP, C_FIBER = "#2a78d6", "#eb6834"      # validierte Referenzpalette, Slot 1 und 2
C_ARRAY = "#1baf7a"                         # Slot 3; Slot 1-3 sind all-pairs-validiert.
# Alle Arrays ohne Chip teilen sich diesen Slot: die Farbe traegt die Kategorie,
# den Hersteller traegt die Zeilenbeschriftung. Slot 4 (Gelb) faellt gegen
# C_FIBER (Orange) durch die All-pairs-Pruefung.
CONTROL = "#d6d5cf"
GRID = "#eceae5"

ANCHOR = 529.43                             # Pegelnormierung, in allen Scans vorhanden
CONN_ABS = (555.9, 557.9)                   # Steckerbereich, absolut

B = "raw_data/2026-09-15-08-19_ligentec2m07cmfiberwavelenght1520to1570nm"
M = "raw_data/2026-09-15-15-14meisusmfiberarrzchannel"
P = "raw_data/2026-09-15-15-14phixsmfiberarrzchannel"
SCANS = [
    ("Ligentec fiber 1  (channel 95, loopback)", C_CHIP,
     "raw_data/2672_ligentechhi_2026-09-14-16-02_Ligentecfiber1_reflectogram.csv", "Chip facet A"),
    ("Ligentec fiber 44  (channel 52)", C_CHIP,
     "raw_data/2672_ligentechhi_2026-09-14-16-36_Ligentecfiber44_reflectogram.csv", "Chip facet A"),
    ("Ligentec fiber 48  (channel 48)", C_CHIP,
     "raw_data/2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48_reflectogram.csv", "Chip facet A"),
    ("PM fiber a", C_FIBER, B + "_reflectogram.csv", "Fiber end"),
    ("PM fiber b", C_FIBER, B + "_b_reflectogram.csv", "Fiber end"),
    ("PM fiber c", C_FIBER, B + "_c_reflectogram.csv", "Fiber end"),
    ("Meisu array  (channel 4)", C_ARRAY, M + "4_reflectogram.csv", "Array facet"),
    ("Meisu array  (channel 4, remeasured)", C_ARRAY,
     M + "4remeasure_reflectogram.csv", "Array facet"),
    ("Meisu array  (channel 8)", C_ARRAY, M + "8_reflectogram.csv", "Array facet"),
    ("Meisu array  (channel 9)", C_ARRAY, M + "9_reflectogram.csv", "Array facet"),
    ("Meisu array  (channel 12)", C_ARRAY, M + "12_reflectogram.csv", "Array facet"),
    ("Meisu array  (channel 12, remeasured)", C_ARRAY,
     M + "12remeasure_reflectogram.csv", "Array facet"),
    ("PHIX array  (channel 1)", C_ARRAY, P + "1_reflectogram.csv", "Array facet"),
    ("PHIX array  (channel 8)", C_ARRAY, P + "8_reflectogram.csv", "Array facet"),
]
CONTROL_CSV = "raw_data/2672_ligentechhi_2026-09-14-16-19_ligentecnofiber_reflectogram.csv"
# Faserenden: erster Peak der Doublette (siehe logs/2026-09-15.md)
ENDS = {"PM fiber a": 2642.975, "PM fiber b": 2650.936, "PM fiber c": 2652.395}


def load(p):
    d = np.genfromtxt(p, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def on_anchor(p, zero):
    z, a = load(p)
    m = (z > ANCHOR - 0.12) & (z < ANCHOR + 0.12)
    return z - zero, a - a[m].max()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xmax", type=float, default=2800.0)
    ap.add_argument("--zero", type=float, default=ANCHOR,
                    help="absolute z-Position, die zu 0 mm wird "
                         "(Vorgabe: interner Anker 529.43; Steckerebene: 557.045)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    xc, yc = on_anchor(CONTROL_CSV, a.zero)
    n = len(SCANS)
    fig, ax = plt.subplots(n, 1, figsize=(13, 1.95 * n + 1.6), sharex=True,
                           constrained_layout=True, facecolor=SURFACE)

    print("%-42s %12s %9s" % ("Scan", "Feature [mm]", "dB"))
    for k, (label, color, path, what) in enumerate(SCANS):
        s = ax[k]
        x, y = on_anchor(path, a.zero)
        vis = (x > -10) & (x < a.xmax)
        top = float(y[vis].max())      # jede Zeile auf ihren EIGENEN Peak normiert;
        y = y - top                    # gemeinsam ist die x-Achse, nicht der Pegel
        s.set_facecolor(SURFACE)
        s.axvspan(CONN_ABS[0] - a.zero, CONN_ABS[1] - a.zero,
                  color=MUTED, alpha=0.18, zorder=0)
        mm = (xc > -10) & (xc < a.xmax)
        s.plot(xc[mm], yc[mm] - top, lw=0.7, color=CONTROL, zorder=1)
        mm = (x > -10) & (x < a.xmax)
        s.plot(x[mm], y[mm], lw=0.8, color=color, zorder=3)

        if label in ENDS:
            fx = ENDS[label] - a.zero
            w = (x > fx - 0.06) & (x < fx + 0.06)
        else:
            sel = (x > 900 - (a.zero - ANCHOR)) & (x < 1200 - (a.zero - ANCHOR))
            i = int(np.argmax(y[sel]))
            fx = float(x[sel][i])
            w = (x > fx - 0.06) & (x < fx + 0.06)
        fy = float(y[w].max())
        print("%-42s %12.2f %9.1f" % (label, fx, fy))
        s.plot([fx], [fy], "o", ms=7, mfc=color, mec=SURFACE, mew=1.8, zorder=5)
        s.annotate("%s  %.1f mm" % (what, fx), xy=(fx, fy), xytext=(fx + 70, fy + 3),
                   color=INK, fontsize=9.5, fontweight="bold", ha="left", va="center",
                   arrowprops=dict(arrowstyle="-", color=color, lw=1.1))

        s.set_xlim(-10, a.xmax)
        s.set_ylim(-82, 12)
        s.set_yticks([-75, -50, -25, 0])
        s.grid(color=GRID, lw=0.8, zorder=0)
        s.set_axisbelow(True)
        for sp in ("top", "right"):
            s.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            s.spines[sp].set_color(GRID)
        s.tick_params(colors=INK2, labelsize=8.5)
        s.set_ylabel("dB", color=INK2, fontsize=9)
        s.text(0.004, 0.94, label, transform=s.transAxes, color=INK,
               fontsize=10.5, fontweight="bold", va="top", ha="left")

    ax[0].text(CONN_ABS[1] - a.zero + 12, -76, "connector plane", color=INK2,
               fontsize=8.5, va="bottom", ha="left")
    ax[-1].set_xlabel("Distance from %s  [mm]        each row normalized to its own peak; gray = control with nothing connected"
                      % ("connector plane" if abs(a.zero - ANCHOR) > 1
                         else "internal anchor in the circulator section"),
                      color=INK2, fontsize=10)
    ax[0].set_title("Chip and fiber scans on one axis: top Ligentec array + chip, "
                    "middle bare fiber, bottom arrays without chip",
                    color=INK, fontsize=13.5, fontweight="bold", loc="left")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150, facecolor=SURFACE)
    print("PNG: %s" % a.out)


if __name__ == "__main__":
    main()
