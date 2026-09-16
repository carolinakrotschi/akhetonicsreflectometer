#!/usr/bin/env python3
"""Beliebig viele Reflektogramme uebereinander legen -- fuer den Fall, wo
`compare_span.py`s Annahme (zwei Spans DERSELBEN Facette) nicht gilt, z.B.
"Faser 1 auf Chip A" gegen "Faser 1 auf Chip B": unterschiedliche Chips
heisst unterschiedliche physische Facettenposition, ein gemeinsames enges
Zoomfenster um eine gemittelte Position ergibt dort keinen Sinn.

Zwei Panels: Uebersicht (Huellkurve, alle Kurven auf eigenes Maximum
normiert) und ein Facettenfenster, das breit genug ist, alle uebergebenen
Peaks gemeinsam zu zeigen (nicht auf einen einzigen engen Bereich gezoomt).
Jede Kurve bekommt einen eigenen -3-dB-Peak-Eintrag in der Textausgabe und
eine eigene senkrechte Markierung im Plot.

Aufruf
    python tools/plot_chip_comparison.py \
        --scan "MAP2672 fiber1" raw_data/2672_ligentechhi_..._reflectogram.csv \
        --scan "MAP2680 fiber1 (1520-1570nm)" raw_data/..._fiber1_reflectogram.csv \
        --scan "MAP2680 fiber1 (voller Span)" raw_data/..._fiber1full_reflectogram.csv \
        --peak-window 1575 1590 \
        --out results/2026-09-16/fiber1_chip_vergleich.png
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
from process_reflectogram_aux import noise_floor        # noqa: E402
from compare_span import load_marks, short_name         # noqa: E402

COLORS = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b", "#17becf")
N_FIBER = 1.468


def load(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def envelope(z, db, nout):
    """Blockweises Maximum auf ~nout Punkte, s. compare_span.py."""
    k = max(len(z) // nout, 1)
    m = (len(z) // k) * k
    zz = z[:m].reshape(-1, k)
    dd = db[:m].reshape(-1, k)
    i = np.argmax(dd, axis=1)
    return zz[np.arange(len(i)), i], dd[np.arange(len(i)), i]


def peak_and_width(z, db, lo, hi):
    """Wie in compare_span.py: Hauptpeak im Fenster, -3 dB Breite per
    linearer Interpolation auf den Flanken."""
    m = (z >= lo) & (z <= hi)
    zz, aa = z[m], db[m]
    if len(aa) == 0:
        return None
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
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", nargs=2, action="append", required=True,
                    metavar=("LABEL", "CSV"),
                    help="mehrfach angebbar, mindestens zweimal")
    ap.add_argument("--peak-window", nargs=2, type=float, required=True,
                    metavar=("LO_MM", "HI_MM"),
                    help="mm, Suchfenster fuer den Hauptpeak jeder Kurve "
                         "-- muss alle Facetten enthalten")
    ap.add_argument("--zoom-pad", type=float, default=0.5,
                    help="mm Rand links/rechts der aeussersten Peaks im "
                         "Facettenpanel (default 0.5)")
    ap.add_argument("--align-pad", type=float, default=0.3,
                    help="mm links/rechts vom jeweils EIGENEN Peak im "
                         "zentrierten Panel -- zeigt die Struktur um die "
                         "Facette unabhaengig von deren absoluter Position "
                         "(default 0.3)")
    ap.add_argument("--marks-csv", default=None,
                    help="CSV von gds_trace_ligentec.py: zeichnet die "
                         "Sollpositionen der Bauteile hinter die Facette")
    ap.add_argument("--ng", nargs=2, type=float, default=[1.80, 2.00],
                    help="Unsicherheitsbereich des Gruppenindex fuer die Marken")
    ap.add_argument("--marks-zmax", type=float, default=13.0,
                    help="mm hinter der Facette, bis wohin das Bauteilpanel geht")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    scans = []
    for label, csv in a.scan:
        z, db = load(csv)
        pw = peak_and_width(z, db, *a.peak_window)
        if pw is None:
            raise SystemExit("kein Peak im Fenster %s fuer %r" % (a.peak_window, label))
        zp, top, w = pw
        fl = noise_floor(db)
        scans.append(dict(z=z, db=db, label=label, zp=zp, w=w, floor=fl,
                          n=len(z), dz=float(np.median(np.diff(z))) * 1e3))
        print("%-32s  Peak %10.4f mm   -3 dB %7.1f um   %+6.1f dB ueber Boden  "
              "Bins %8d a %.2f um"
              % (label, zp, w * 1e3, top - fl, len(z), scans[-1]["dz"]))

    zps = [s["zp"] for s in scans]
    print("\nPeakpositionen streuen ueber %.1f um (%.4f .. %.4f mm)"
          % ((max(zps) - min(zps)) * 1e3, min(zps), max(zps)))

    marks = load_marks(a.marks_csv) if a.marks_csv else []
    nrow = 3 + (2 if marks else 0)
    fig, ax = plt.subplots(nrow, 1, figsize=(11, 3.45 * nrow),
                           gridspec_kw=dict(hspace=0.36, top=0.965,
                                            bottom=0.045, left=0.08,
                                            right=0.98))

    # dichteste Abtastung nach hinten, sonst deckt ihr Rauschen die anderen zu
    order = sorted(range(len(scans)), key=lambda i: -scans[i]["n"])

    zoom = (min(zps) - a.zoom_pad, max(zps) + a.zoom_pad)
    spans = ((0.0, 3000.0), zoom, (-a.align_pad, a.align_pad))
    titles = ("overview", "facet region (+/- %.2f mm around outermost peaks)" % a.zoom_pad,
              "aligned to own facet -- same structure repeating?")

    for k, (axk, (lo, hi), ttl) in enumerate(zip(ax, spans, titles)):
        for rank, i in enumerate(order):
            s, c = scans[i], COLORS[i % len(COLORS)]
            x = s["z"] - s["zp"] if k == 2 else s["z"]
            m = (x >= lo) & (x <= hi)
            zz, dd = x[m], s["db"][m]
            if k == 0:
                zz, dd = envelope(zz, dd, 3000)
            axk.plot(zz, dd, lw=0.8, color=c,
                     alpha=0.75 if rank == 0 and len(scans) > 1 else 1.0,
                     zorder=2 + rank,
                     label="%s  (%d bins, %.2f um)" % (s["label"], s["n"], s["dz"]))
            axk.axhline(s["floor"], color=c, ls=":", lw=0.8)
            if k == 1:
                axk.axvline(s["zp"], color=c, ls="-", lw=0.7, alpha=0.6)
        axk.set_xlim(lo, hi)
        axk.set_ylabel("amplitude [dB]")
        axk.text(0.5, 1.0, ttl, transform=axk.transAxes, ha="center", va="top",
                 fontsize=10, bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
        axk.grid(alpha=0.3)
        axk.legend(loc="upper right", fontsize=8)
    ax[0].text(0.01, 0.03,
               "dotted: RMS noise floor  " + " / ".join("%.1f dB" % s["floor"] for s in scans),
               transform=ax[0].transAxes, fontsize=8)
    ax[1].set_xlabel("distance [mm]")
    ax[1].text(0.01, 0.07,
               "-3 dB width  " + "  /  ".join("%.1f um (%s)" % (s["w"] * 1e3, s["label"]) for s in scans),
               transform=ax[1].transAxes, fontsize=8,
               bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))
    ax[2].set_xlabel("distance from own facet peak [mm]")
    ax[2].axvline(0.0, color="0.3", ls="--", lw=0.6)

    if marks:
        # Bauteile liegen ON-CHIP; die z-Achse ist faseraequivalent
        # (n_g = 1.468). Deshalb die on-chip-Weglaenge mit n_g_chip/n_Faser
        # skalieren -- das graue Band ist die Unsicherheit in n_g_chip.
        bands = [(short_name(nm), s_um * 1e-3 * a.ng[0] / N_FIBER,
                  s_um * 1e-3 * a.ng[1] / N_FIBER) for nm, s_um in marks]
        # Der Edge Coupler klebt an der Facette; wuerde er den Zoom
        # aufspannen, waere der Rest wieder zusammengequetscht.
        far = [b for b in bands if b[2] > 1.0] or bands
        pad = 0.15 * (max(b[2] for b in far) - min(b[1] for b in far) or 1.0)
        zoom_m = (min(b[1] for b in far) - pad, max(b[2] for b in far) + pad)

        for row, (lo_x, hi_x, ttl) in enumerate((
                (-0.3, a.marks_zmax,
                 "predicted component positions (grey: n_g %.2f-%.2f)" % tuple(a.ng)),
                (zoom_m[0], zoom_m[1], "component area, zoomed"))):
            axm = ax[3 + row]
            for i, s_ in enumerate(scans):
                c = COLORS[i % len(COLORS)]
                x = s_["z"] - s_["zp"]
                m = (x >= lo_x) & (x <= hi_x)
                axm.plot(x[m], s_["db"][m], lw=0.8, color=c, alpha=0.9,
                         label=s_["label"] if row == 0 else None)
                axm.axhline(s_["floor"], color=c, ls=":", lw=0.8)
            for i, (name, lo_, hi_) in enumerate(bands):
                axm.axvspan(lo_, hi_, color="0.5", alpha=0.22, lw=0)
                axm.annotate(name,
                             xy=(0.5 * (lo_ + hi_), 0.97 - 0.30 * (i % 3)),
                             xycoords=("data", "axes fraction"),
                             fontsize=7, rotation=90, ha="center", va="top",
                             color="0.2")
            axm.set_xlim(lo_x, hi_x)
            axm.set_ylabel("amplitude [dB]")
            axm.text(0.5, 1.0, ttl, transform=axm.transAxes, ha="center",
                     va="top", fontsize=10,
                     bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
            axm.grid(alpha=0.3)
            if row == 0:
                axm.legend(loc="upper right", fontsize=8)
            else:
                axm.set_xlabel("distance behind facet A [mm]")
                over = []
                for s_ in scans:
                    x = s_["z"] - s_["zp"]
                    m = (x >= lo_x) & (x <= hi_x)
                    over.append("%s %+.1f dB"
                                % (s_["label"], s_["db"][m].max() - s_["floor"]))
                axm.text(0.01, 0.06,
                         "peak in this window over own noise floor:   "
                         + ",   ".join(over),
                         transform=axm.transAxes, fontsize=9,
                         bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))

    fig.suptitle("Facet comparison: " + "  vs  ".join(s["label"] for s in scans))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150)
    print("\nPlot: %s" % a.out)


if __name__ == "__main__":
    main()
