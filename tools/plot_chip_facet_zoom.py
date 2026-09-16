#!/usr/bin/env python3
"""Die Chipfacette gross: alle Ligentec-Kanaele nebeneinander.

Zwei Ansichten derselben Daten, weil sie zwei verschiedene Fragen
beantworten:

  oben   auf der Steckerebene -- zeigt, WO die Facetten liegen. Die
         Unterschiede sind Laengenunterschiede der Arrayfasern, nicht
         des Chips.
  unten  jeder Kanal auf SEINE eigene Facette zentriert -- zeigt, was
         HINTER der Facette kommt. Das ist die eigentliche Frage: kommt
         aus dem Chip etwas zurueck? Wenn alle Kurven hier
         uebereinanderliegen, ist das Gesehene Instrument und nicht Chip.

Aufruf
    python tools/plot_chip_facet_zoom.py --out results/<datum>/facetten.png
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a85"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # validierte Referenzpalette
CONTROL = "#c9c8c2"
GRID = "#e6e5e0"

ANCHOR = 529.43
CONN = 557.045

SCANS = [
    ("Fiber 1  (channel 95, loopback)",
     "raw_data/2672_ligentechhi_2026-09-14-16-02_Ligentecfiber1_reflectogram.csv"),
    ("Fiber 44  (channel 52)",
     "raw_data/2672_ligentechhi_2026-09-14-16-36_Ligentecfiber44_reflectogram.csv"),
    ("Fiber 45  (channel 51)",
     "raw_data/2672_ligentechhi_2026-09-15-08-19_ligentecfiber45wavelenght1520to1570nm_reflectogram.csv"),
    ("Fiber 48  (channel 48)",
     "raw_data/2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48_reflectogram.csv"),
]
CONTROL_CSV = "raw_data/2672_ligentechhi_2026-09-14-16-19_ligentecnofiber_reflectogram.csv"


def load(p):
    d = np.genfromtxt(p, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--behind", type=float, default=14.0,
                    help="wie weit hinter die Facette gezeigt wird [mm]")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    fig, ax = plt.subplots(2, 1, figsize=(13, 9.5), constrained_layout=True,
                           facecolor=SURFACE)
    for s in ax:
        s.set_facecolor(SURFACE)

    print("%-32s %12s %10s %12s" % ("channel", "facet [mm]", "dyn. range", "width"))
    F = []
    for k, (label, path) in enumerate(SCANS):
        z, y = load(path)
        m = (z > 1560) & (z < 1600)
        i = int(np.argmax(y[m]))
        zf, yf = float(z[m][i]), float(y[m][i])
        y = y - yf                                  # auf die eigene Facette normiert
        F.append((label, z, y, zf, SERIES[k % len(SERIES)]))

        w = (z > zf - 0.3) & (z < zf + 0.3)
        q, v = z[w], y[w]
        j = int(np.argmax(v)); h = v[j] - 3
        lo = j
        while lo > 0 and v[lo] > h:
            lo -= 1
        hi = j
        while hi < len(v) - 1 and v[hi] > h:
            hi += 1
        far = (z > zf + 200) & (z < zf + 1200)
        fl = 20 * np.log10(np.sqrt(np.mean((10 ** (y[far] / 20)) ** 2)))
        print("%-32s %12.3f %9.1f dB %9.0f um"
              % (label, zf, -fl, (q[hi] - q[lo]) * 1e3))

        ax[0].plot(z - CONN, y, lw=1.0, color=F[-1][4],
                   label="%s  -  %.2f mm" % (label, zf - CONN))
        ax[1].plot(z - zf, y, lw=1.0, color=F[-1][4], label=label)

    zc, yc = load(CONTROL_CSV)
    mc = (zc > ANCHOR - 0.12) & (zc < ANCHOR + 0.12)
    # Kontrolle auf denselben Pegelbezug wie Kanal 48 bringen
    z48, y48, zf48 = F[3][1], F[3][2], F[3][3]
    ma = (z48 > ANCHOR - 0.12) & (z48 < ANCHOR + 0.12)
    ycn = (yc - yc[mc].max()) + y48[ma].max()
    ax[0].plot(zc - CONN, ycn, lw=0.8, color=CONTROL, zorder=1,
               label="control (nothing connected)")
    ax[1].plot(zc - 1582.0, ycn, lw=0.8, color=CONTROL, zorder=1,
               label="control (nothing connected)")

    lo0 = min(f[3] for f in F) - CONN - 4
    hi0 = max(f[3] for f in F) - CONN + a.behind
    ax[0].set_xlim(lo0, hi0)
    ax[1].set_xlim(-4, a.behind)
    for s in ax:
        s.set_ylim(-86, 10)
        s.grid(color=GRID, lw=0.8, zorder=0)
        s.set_axisbelow(True)
        for sp in ("top", "right"):
            s.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            s.spines[sp].set_color(GRID)
        s.tick_params(colors=INK2, labelsize=9)
        s.set_ylabel("amplitude rel. to own facet  [dB]", color=INK2, fontsize=10)
        s.legend(frameon=False, fontsize=8.5, labelcolor=INK2,
                 loc="upper left", ncol=1)

    ax[0].set_xlabel("distance from connector plane  [mm]", color=INK2, fontsize=10.5)
    ax[1].set_xlabel("distance behind own facet A  [mm]", color=INK2, fontsize=10.5)
    ax[0].set_title("Chip facet, all four Ligentec channels (MAP2672) -- where they sit",
                    color=INK, fontsize=13, fontweight="bold", loc="left")
    ax[1].set_title("Same data, each channel centred on its own facet -- "
                    "what comes behind it",
                    color=INK, fontsize=12, fontweight="bold", loc="left")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150, facecolor=SURFACE)
    print("PNG: %s" % a.out)


if __name__ == "__main__":
    main()
