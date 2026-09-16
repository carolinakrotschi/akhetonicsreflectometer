#!/usr/bin/env python3
"""Reflektogramme auf die Steckerebene beziehen, damit man Faserlaengen abliest.

Auf der absoluten Achse liegt ein Faserende irgendwo bei 2.65 m und der
aux-Geist bei ~2.10 m -- wer dort nach "meiner 2.07-m-Faser" sucht, findet
den aux. Bezieht man die Achse dagegen auf den Steckerreflex, steht jedes
Faserende bei seiner tatsaechlichen Laenge.

Mehrere Fasern koennen uebereinandergelegt werden (--scan mehrfach).

WARNUNG zum aux-Geist: seine Position ist ~ tau_aux und aendert sich
deshalb von Scan zu Scan, auch wenn sich physikalisch nichts bewegt.
Echte Reflektoren stehen still (der interne Anker bei 529.43 mm ist die
Probe aufs Exempel). Ein Peak, der zwischen zwei Scans mitwandert, ist
der aux und keine Faser.

Aufruf
    python tools/plot_fiber_length.py \
        --scan "Faser a" raw_data/..._reflectogram.csv \
        --scan "Faser b" raw_data/..._b_reflectogram.csv \
        --control raw_data/2672_ligentechhi_2026-09-14-16-19_ligentecnofiber_reflectogram.csv \
        --connector 557.045 --nominal 2070 --out results/<datum>/plot.png
"""

import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8a85"
# kategoriale Reihenfolge der validierten Referenzpalette, feste Slots
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
CONTROL = "#c9c8c2"        # Kontrolle = neutrale Referenzlinie, keine eigene Identitaet
GRID = "#e6e5e0"


def load(p):
    d = np.genfromtxt(p, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", nargs=2, action="append", required=True,
                    metavar=("LABEL", "CSV"))
    ap.add_argument("--control", default=None)
    ap.add_argument("--end", nargs=2, action="append", default=None,
                    metavar=("LABEL", "MM"),
                    help="Faserende explizit setzen (Abstand ab Steckerebene), statt einfach den staerksten Peak zu nehmen")
    ap.add_argument("--connector", type=float, required=True,
                    help="z des Steckerreflexes in mm (= neuer Nullpunkt)")
    ap.add_argument("--nominal", type=float, default=None,
                    help="Massband-Laenge in mm, wird als Marke eingezeichnet")
    ap.add_argument("--anchor", type=float, default=529.43)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    S = []
    for k, (label, path) in enumerate(a.scan):
        z, amp = load(path)
        m = (z > a.anchor - 0.08) & (z < a.anchor + 0.08)
        amp = amp - amp[m].max()               # gemeinsamer interner Anker
        x = z - a.connector                    # Steckerebene = 0
        sel = x > 1000.0
        i = int(np.argmax(amp[sel]))
        L, dB = float(x[sel][i]), float(amp[sel][i])
        for lab2, mm2 in (a.end or []):
            if lab2 == label:
                L = float(mm2)
                w = (x > L - 0.06) & (x < L + 0.06)
                dB = float(amp[w].max())
        S.append(dict(label=label, x=x, y=amp, color=SERIES[k % len(SERIES)],
                      L=L, dB=dB))

    if a.control:
        zc, ampc = load(a.control)
        mc = (zc > a.anchor - 0.08) & (zc < a.anchor + 0.08)
        ampc = ampc - ampc[mc].max()
        xc = zc - a.connector

    print("%-12s %12s %9s %12s" % ("Scan", "Laenge [mm]", "dB", "gg. Massband"))
    for s in S:
        d = ("%+.2f mm" % (s["L"] - a.nominal)) if a.nominal else "-"
        print("%-12s %12.2f %9.1f %12s" % (s["label"], s["L"], s["dB"], d))

    Ls = [s["L"] for s in S]
    lo2, hi2 = min(Ls) - 45, max(Ls) + 45
    if a.nominal:
        lo2 = min(lo2, a.nominal - 12)
    spans = [(0, max(max(Ls) * 1.10, 2300)), (lo2, hi2)]

    fig, ax = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True,
                           facecolor=SURFACE)
    for k, (xa, xb) in enumerate(spans):
        s = ax[k]
        s.set_facecolor(SURFACE)
        if a.control:
            mm = (xc > xa) & (xc < xb)
            s.plot(xc[mm], ampc[mm], lw=0.9, color=CONTROL, zorder=2,
                   label="Kontrolle - nichts angeschlossen")
        for j, d in enumerate(S):
            mm = (d["x"] > xa) & (d["x"] < xb)
            s.plot(d["x"][mm], d["y"][mm], lw=1.0, color=d["color"],
                   zorder=3 + j, alpha=0.9, label=d["label"])
        if a.nominal is not None:
            s.axvline(a.nominal, color=MUTED, ls=(0, (5, 4)), lw=1.2, zorder=1)
        for d in S:
            s.plot([d["L"]], [d["dB"]], "o", ms=8, mfc=d["color"],
                   mec=SURFACE, mew=2, zorder=9)
        s.set_xlim(xa, xb)
        s.set_ylim(-46, 14)
        s.grid(color=GRID, lw=0.8, zorder=0)
        s.set_axisbelow(True)
        for sp in ("top", "right"):
            s.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            s.spines[sp].set_color(GRID)
        s.tick_params(colors=INK2, labelsize=9)
        s.set_xlabel("Abstand ab Steckerebene  [mm]", color=INK2, fontsize=10)
        s.set_ylabel("Amplitude  [dB rel. Anker %.2f mm]" % a.anchor,
                     color=INK2, fontsize=10)

    # Direktbeschriftung nur an den Faserenden, gestaffelt damit nichts kollidiert
    order = sorted(range(len(S)), key=lambda j: S[j]["L"])
    for rank, j in enumerate(order):
        d = S[j]
        for s in (ax[0], ax[1]):
            xa, xb = s.get_xlim()
            ytxt = 11.5 - rank * 7.5
            s.annotate("%s   %.1f mm   %.1f dB" % (d["label"], d["L"], d["dB"]),
                       xy=(d["L"], d["dB"] + 0.8),
                       xytext=(d["L"] - 0.015 * (xb - xa), ytxt),
                       color=INK, fontsize=10, fontweight="bold",
                       ha="right", va="top",
                       arrowprops=dict(arrowstyle="->", color=d["color"], lw=1.3))
    if a.nominal is not None:
        xa, xb = ax[1].get_xlim()
        ax[1].annotate("Maßband %.0f mm" % a.nominal,
                       xy=(a.nominal, -26), xytext=(a.nominal - 0.015 * (xb - xa), -16),
                       color=INK2, fontsize=9.5, ha="right",
                       arrowprops=dict(arrowstyle="->", color=MUTED, lw=1))

    ax[0].legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper left")
    ax[0].set_title("Auf die Steckerebene bezogen liest man die Faserlänge direkt ab",
                    color=INK, fontsize=13, fontweight="bold", loc="left")
    ax[1].set_title("Zoom auf die Faserenden  -  jede Faser zeigt eine Doublette; markiert ist der ERSTE Peak  (Aufloesung 16.6 um)",
                    color=INK, fontsize=10, loc="left")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150, facecolor=SURFACE)
    print("PNG: %s" % a.out)


if __name__ == "__main__":
    main()
