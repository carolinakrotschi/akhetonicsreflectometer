#!/usr/bin/env python3
"""Zwei Reflektogramme desselben Aufbaus mit verschiedenem Sweepbereich
gegenueberstellen: was bringt der volle Span wirklich?

Die Rechnung sagt Faktor 2.3 (16.3 -> 6.9 um Binbreite). Was davon am
Ende ankommt, entscheidet nicht die Rechnung, sondern die gemessene
-3-dB-Breite des Facettenpeaks und der Rauschboden. Genau die beiden
Zahlen zieht dieses Skript aus den CSVs und stellt sie nebeneinander.

Gezeigt wird dreimal dasselbe Paar in wachsender Vergroesserung:
Uebersicht, Umgebung des Facettenpeaks, und der Peak selbst mit seinen
-3-dB-Flanken. Beide Kurven sind auf ihr eigenes Maximum normiert; die
absoluten Pegel sind nicht vergleichbar, die Formen schon.

Aufruf
    python tools/compare_span.py \
        --a raw_data/..._remeasure_reflectogram.csv "1520-1570 nm" \
        --b raw_data/..._und1505to1625_reflectogram.csv "1505-1625 nm" \
        --out results/2026-09-15/span_vergleich.png
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from process_reflectogram_aux import noise_floor        # noqa: E402

N_FIBER = 1.468


def load(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def envelope(z, db, nout):
    """Blockweises Maximum auf ~nout Punkte. Peaks bleiben erhalten, der
    Rauschboden wird zur oberen Huellkurve (deshalb liegt er im Plot etwas
    ueber der gestrichelten RMS-Linie -- das ist kein Widerspruch)."""
    k = max(len(z) // nout, 1)
    m = (len(z) // k) * k
    zz = z[:m].reshape(-1, k)
    dd = db[:m].reshape(-1, k)
    i = np.argmax(dd, axis=1)
    return zz[np.arange(len(i)), i], dd[np.arange(len(i)), i]


def load_marks(path):
    """Bauteile aus der GDS-Spur. Zurueck kommt (Name, s_ab_Facette_um) je
    Bauteil; Containerzellen ohne eigene Grenzflaeche fliegen raus, sonst
    steht im Plot dreimal derselbe Strich."""
    skip = ("XGM_2_Decoder", "XGMx3_SiN", "MZM_switch", "EdgeCoupler_BB",
            "MMI_1x2_BB", "_metal", "cross_array")
    out = []
    with open(path) as fh:
        head = fh.readline().strip().split(",")
        ic, isf = head.index("cell"), head.index("s_from_facet_um")
        for line in fh:
            f = line.strip().split(",")
            if not f or len(f) <= max(ic, isf) or any(s in f[ic] for s in skip):
                continue
            out.append((f[ic], float(f[isf])))
    out.sort(key=lambda t: t[1])
    return out


def panel_label(axk, text):
    """Panelueberschrift IN die Achse. Ueber der Achse kollidiert sie bei
    fuenf Panels mit den Ticklabels der daruberliegenden."""
    axk.text(0.5, 1.0, text, transform=axk.transAxes, ha="center",
             va="top", fontsize=10, bbox=dict(fc="white", ec="none",
                                              alpha=0.75, pad=1.5))


def short_name(cell):
    """Zellname -> Plotbeschriftung. Die GDS-Namen tragen Instanzsuffixe
    ($a5b0, _8721), die im Plot nur Platz kosten."""
    s = re.sub(r"(_\$[0-9a-f]+|_\d{3,})$", "", cell)
    s = s.replace("AN350BB_", "").replace("_symmetric_C", "")
    s = s.replace("EdgeCoupler_Lensed_C", "EdgeCoupler")
    s = s.replace("variable_pitch_splitting", "pitch split")
    return s.replace("_", " ").strip()


def peak_and_width(z, db, lo, hi):
    """Hauptpeak im Fenster [lo,hi] mm, -3 dB Breite durch lineare
    Interpolation auf den Flanken (nicht Bin-Zaehlen: bei 7 um Bins waeren
    das sonst Stufen von 17 um)."""
    m = (z >= lo) & (z <= hi)
    zz, aa = z[m], db[m]
    i = int(np.argmax(aa))
    top = aa[i]
    half = top - 3.0

    def cross(idx, step):
        j = idx
        while 0 < j < len(aa) - 1 and aa[j] > half:
            j += step
        if j == idx:
            return zz[idx]
        z1, z2 = zz[j - step], zz[j]
        a1, a2 = aa[j - step], aa[j]
        if a1 == a2:
            return z2
        return z1 + (half - a1) * (z2 - z1) / (a2 - a1)

    left, right = cross(i, -1), cross(i, +1)
    return float(zz[i]), float(top), float(right - left)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", nargs=2, required=True, metavar=("CSV", "LABEL"))
    ap.add_argument("--b", nargs=2, required=True, metavar=("CSV", "LABEL"))
    ap.add_argument("--c", nargs=2, default=None, metavar=("CSV", "LABEL"),
                    help="optionale dritte Spur (z.B. derselbe Span nach einem "
                         "Eingriff an der Facette)")
    ap.add_argument("--peak-window", nargs=2, type=float, default=[1600.0, 1700.0],
                    help="mm, Suchfenster fuer den Facettenpeak")
    ap.add_argument("--marks-csv", default=None,
                    help="CSV von gds_trace_ligentec.py: zeichnet die "
                         "Sollpositionen der Bauteile hinter die Facette")
    ap.add_argument("--ng", nargs=2, type=float, default=[1.80, 2.00],
                    help="Unsicherheitsbereich des Gruppenindex fuer die Marken")
    ap.add_argument("--marks-zmax", type=float, default=12.0,
                    help="mm hinter der Facette, bis wohin das Bauteilpanel geht")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    scans = []
    for csv, label in [a.a, a.b] + ([a.c] if a.c else []):
        z, db = load(csv)
        zp, top, w = peak_and_width(z, db, *a.peak_window)
        fl = noise_floor(db)
        scans.append(dict(z=z, db=db, label=label, zp=zp, w=w, floor=fl,
                          n=len(z), dz=float(np.median(np.diff(z))) * 1e3))
        print("%-14s  Peak %10.4f mm   -3 dB %7.1f um   Boden %6.1f dB   "
              "Bins %8d a %.2f um"
              % (label, zp, w * 1e3, fl, len(z), scans[-1]["dz"]))

    s0 = scans[0]
    for s in scans[1:]:
        print("\n%s -> %s:  Breite %.1f -> %.1f um (Faktor %.2f), "
              "Boden %.1f -> %.1f dB (%+.1f dB), Peakposition weicht um %.1f um ab"
              % (s0["label"], s["label"], s0["w"] * 1e3, s["w"] * 1e3,
                 s0["w"] / s["w"], s0["floor"], s["floor"],
                 s["floor"] - s0["floor"], abs(s0["zp"] - s["zp"]) * 1e3))

    marks = load_marks(a.marks_csv) if a.marks_csv else []

    nrow = 5 if marks else 3
    fig, ax = plt.subplots(nrow, 1, figsize=(11, 3.1 * nrow),
                           gridspec_kw=dict(hspace=0.32, top=0.965,
                                            bottom=0.045, left=0.075,
                                            right=0.985))
    col = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd")
    zc = float(np.mean([s["zp"] for s in scans]))
    spans = ((0.0, 3000.0), (zc - 2.0, zc + 2.0), (zc - 0.12, zc + 0.12))
    titles = ("overview", "facet peak, +/- 2 mm", "facet peak, +/- 120 um")

    for k, (axk, (lo, hi), ttl) in enumerate(zip(ax[:3], spans, titles)):
        # den dichter abgetasteten Scan nach hinten, sonst deckt sein
        # Rauschen den anderen komplett zu
        order = sorted(range(len(scans)), key=lambda i: -scans[i]["n"])
        for rank, i in enumerate(order):
            s, c = scans[i], col[i]
            m = (s["z"] >= lo) & (s["z"] <= hi)
            zz, dd = s["z"][m], s["db"][m]
            if k == 0:
                # Uebersicht: Huellkurve statt Rohpunkte. Zwei Scans mit
                # verschieden dichter Abtastung uebereinander sind sonst nur
                # zwei Rauschbaender, von denen eines das andere zudeckt.
                zz, dd = envelope(zz, dd, 3000)
            axk.plot(zz, dd, lw=0.8, color=c,
                     alpha=0.75 if rank == 0 and len(scans) > 1 else 1.0,
                     zorder=2 + rank,
                     label="%s  (%d bins, %.2f um)" % (s["label"], s["n"], s["dz"])
                     if k == 0 else None)
            axk.axhline(s["floor"], color=c, ls=":", lw=0.8)
            if k == 2:
                axk.axvline(s["zp"], color=c, ls="-", lw=0.6, alpha=0.6)
        axk.set_xlim(lo, hi)
        axk.set_ylabel("amplitude [dB]")
        panel_label(axk, ttl)
        axk.grid(alpha=0.3)
    ax[0].legend(loc="upper right", fontsize=8)
    ax[0].text(0.01, 0.03,
               "dotted: RMS noise floor  "
               + " / ".join("%.1f dB" % x["floor"] for x in scans),
               transform=ax[0].transAxes, fontsize=8)

    for s, c in zip(scans, col):
        ax[2].axhline(s["db"].max() - 3.0, color=c, ls="--", lw=0.6)
    ax[2].text(0.01, 0.04 + 0.075 * (1 + len(scans)),
               "-3 dB width  "
               + "  vs  ".join("%.1f um" % (s["w"] * 1e3) for s in scans)
               + "\ngain over %s:  " % s0["label"]
               + ",   ".join("%s x%.2f (bin limit x%.2f)"
                             % (s["label"], s0["w"] / s["w"], s0["dz"] / s["dz"])
                             for s in scans[1:])
               + "\npeak positions spread over %.1f um -- more than any single "
                 "resolution"
               % ((max(s["zp"] for s in scans)
                   - min(s["zp"] for s in scans)) * 1e3),
               transform=ax[2].transAxes, fontsize=9)
    ax[2].set_xlabel("distance [mm]")

    if marks:
        # Panel 4 und 5: jede Spur auf ihre EIGENE Facette genullt (die
        # beiden Facetten liegen zehner von um auseinander, und verglichen
        # wird, was DAHINTER passiert), darueber die Sollpositionen. Panel 5
        # zoomt auf den Bauteilbereich -- auf 12 mm Breite ist ein 200 um
        # breites MMI-Fenster sonst ein Strich.
        bands = [(short_name(nm), s_um * 1e-3 * a.ng[0] / N_FIBER,
                  s_um * 1e-3 * a.ng[1] / N_FIBER) for nm, s_um in marks]
        # Der Edge Coupler klebt an der Facette; wuerde er den Zoom
        # aufspannen, waere der Rest wieder zusammengequetscht.
        far = [b for b in bands if b[2] > 1.0] or bands
        pad = 0.15 * (max(b[2] for b in far) - min(b[1] for b in far) or 1.0)
        zoom = (min(b[1] for b in far) - pad, max(b[2] for b in far) + pad)

        for row, (lo_x, hi_x, ttl) in enumerate((
                (-0.3, a.marks_zmax,
                 "predicted component positions (grey: n_g %.2f-%.2f)" % tuple(a.ng)),
                (zoom[0], zoom[1], "component area, zoomed"))):
            axm = ax[3 + row]
            for s, c in zip(scans, col):
                x = s["z"] - s["zp"]
                m = (x >= lo_x) & (x <= hi_x)
                axm.plot(x[m], s["db"][m], lw=0.8, color=c, alpha=0.85,
                         label=s["label"] if row == 0 else None)
                axm.axhline(s["floor"], color=c, ls=":", lw=0.8)
            for i, (name, lo_, hi_) in enumerate(bands):
                axm.axvspan(lo_, hi_, color="0.5", alpha=0.22, lw=0)
                axm.annotate(name,
                             xy=(0.5 * (lo_ + hi_), 0.97 - 0.30 * (i % 3)),
                             xycoords=("data", "axes fraction"),
                             fontsize=7, rotation=90, ha="center", va="top",
                             color="0.2")
            axm.set_xlim(lo_x, hi_x)
            axm.set_ylabel("amplitude [dB]")
            panel_label(axm, ttl)
            axm.grid(alpha=0.3)
            if row == 0:
                axm.legend(loc="upper right", fontsize=8)
            else:
                axm.set_xlabel("distance behind facet A [mm]")
                over = []
                for s in scans:
                    x = s["z"] - s["zp"]
                    m = (x >= lo_x) & (x <= hi_x)
                    over.append("%s %+.1f dB"
                                % (s["label"], s["db"][m].max() - s["floor"]))
                axm.text(0.01, 0.06,
                         "peak in this window over own noise floor:   "
                         + ",   ".join(over),
                         transform=axm.transAxes, fontsize=9)

    fig.suptitle("Scan comparison: "
                 + "  vs  ".join(s["label"] for s in scans))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150)
    print("\nPlot: %s" % a.out)


if __name__ == "__main__":
    main()
