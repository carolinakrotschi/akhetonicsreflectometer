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
import csv
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


FIND_COLOR = {"Designgrenzflaeche": "#2ca02c", "Strecke gefuellt": "#2ca02c",
              "Seitenlinie": "0.45", "Doppelbounce": "0.45",
              "Fuss eines Stachels": "0.45", "unerklaert": "#ff7f0e"}
FIND_TAG = {"Seitenlinie": "SL", "Doppelbounce": "2x", "unerklaert": "?"}


def load_findings(paths, nscan):
    """Funde aus hhi_identify_peaks.py, eine Datei je Spur ('none' = keine)."""
    if len(paths) != nscan:
        raise SystemExit("--findings-csv %dx gegen --scan %dx"
                         % (len(paths), nscan))
    out = []
    for p in paths:
        if p.lower() == "none":
            out.append([])
            continue
        with open(p) as fh:
            out.append(list(csv.DictReader(fh)))
    return out


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


def draw_findings(axm, scans, finds, lo_x, hi_x, label_them):
    """Stacheln als Dreieck auf dem Peak, Buckel als Balken darunter.

    Gezeichnet wird, was in den DATEN gefunden wurde -- im Gegensatz zu den
    grauen Marken, die aus dem GDS kommen. Farbe = Deutung."""
    y0, y1 = axm.get_ylim()
    for i, (s_, rows) in enumerate(zip(scans, finds)):
        lane = y0 + (0.05 + 0.045 * i) * (y1 - y0)
        for r in rows:
            v = r["verdict"]
            c = FIND_COLOR.get(v, "0.45")
            if r["kind"] == "Stachel":
                zz = float(r["z_from_mm"])
                if not lo_x <= zz <= hi_x:
                    continue
                yy = float(r["db_over_floor"]) + s_["floor"]
                axm.plot([zz], [yy], marker="v", ms=5, color=c, mec="k",
                         mew=0.3, zorder=6)
                if label_them:
                    axm.annotate(FIND_TAG.get(v, r["candidate"].split()[0]),
                                 xy=(zz, yy), xytext=(0, 7),
                                 textcoords="offset points", fontsize=6.5,
                                 ha="center", color=c, zorder=6)
            else:
                z0, z1 = float(r["z_from_mm"]), float(r["z_to_mm"])
                if z1 < lo_x or z0 > hi_x:
                    continue
                axm.plot([max(z0, lo_x), min(z1, hi_x)], [lane, lane],
                         color=c, lw=3.2, solid_capstyle="butt", zorder=6)
    axm.set_ylim(y0, y1)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan", nargs=2, action="append", required=True,
                    metavar=("LABEL", "CSV"),
                    help="mehrfach angebbar, mindestens zweimal")
    ap.add_argument("--peak-window", nargs=2, type=float, default=None,
                    metavar=("LO_MM", "HI_MM"),
                    help="mm, Suchfenster fuer den Hauptpeak jeder Kurve "
                         "-- muss alle Facetten enthalten. Entfaellt, wenn "
                         "fuer JEDE Spur --facet gesetzt ist")
    ap.add_argument("--facet", action="append", default=None, metavar="MM",
                    help="Facette A dieser Spur in mm statt 'staerkster Peak "
                         "im Fenster' -- einmal je --scan, in derselben "
                         "Reihenfolge; 'auto' faellt auf die Peaksuche "
                         "zurueck. Noetig, wo der staerkste Peak NICHT die "
                         "Facette ist (HHI-Loopback: Facette B ist bis zu "
                         "12 dB staerker) oder wo die Facette gar nicht "
                         "sichtbar ist und aus der Pigtail-Familie kommt")
    ap.add_argument("--facet-unc", type=float, default=0.0, metavar="MM",
                    help="Unsicherheit der Facette-A-Position in mm. Die "
                         "Marken bekommen dann einen zweiten, helleren "
                         "Rand: um so viel kann das GANZE Bauteilmuster "
                         "verschoben sein (default 0 = nicht zeichnen)")
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
    ap.add_argument("--marks-zoom", nargs=2, type=float, default=None,
                    metavar=("LO_MM", "HI_MM"),
                    help="Fenster des untersten Panels von Hand setzen. Der "
                         "Automatik-Zoom spannt von der ersten bis zur "
                         "letzten Marke auf -- bei HHI-Schaltungskanaelen "
                         "sind das 33 Marken ueber 32 mm, da zoomt nichts "
                         "mehr")
    ap.add_argument("--only-marks", action="store_true",
                    help="nur die zwei Bauteilpanels zeichnen. Uebersicht, "
                         "Facettenfenster und zentriertes Panel sagen nichts "
                         "ueber die Bauteile aus und kosten zwei Drittel der "
                         "Bildhoehe")
    ap.add_argument("--findings-csv", action="append", default=None,
                    metavar="CSV",
                    help="Ausgabe von hhi_identify_peaks.py -- einmal je "
                         "--scan, in derselben Reihenfolge ('none' laesst "
                         "eine Spur aus). Zeichnet die GEFUNDENEN Stacheln "
                         "und Buckel mit ihrer Deutung ein")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    facets = a.facet if a.facet else ["auto"] * len(a.scan)
    if len(facets) != len(a.scan):
        raise SystemExit("--facet %dx gegen --scan %dx -- es muss zu jeder "
                         "Spur genau eine Angabe geben"
                         % (len(facets), len(a.scan)))
    if any(f == "auto" for f in facets) and not a.peak_window:
        raise SystemExit("--peak-window fehlt (wird fuer 'auto' gebraucht)")

    scans = []
    for (label, csv), fa in zip(a.scan, facets):
        z, db = load(csv)
        if fa == "auto":
            pw = peak_and_width(z, db, *a.peak_window)
            if pw is None:
                raise SystemExit("kein Peak im Fenster %s fuer %r"
                                 % (a.peak_window, label))
            zp, top, w = pw
            how = ""
        else:
            # Facette gesetzt: Hoehe/Breite trotzdem am naechstgelegenen
            # lokalen Peak ablesen, damit die Zahl vergleichbar bleibt --
            # aber die Achse wird auf den GESETZTEN Wert genullt.
            zp = float(fa)
            pw = peak_and_width(z, db, zp - 0.15, zp + 0.15)
            top, w = (pw[1], pw[2]) if pw else (float("nan"), float("nan"))
            how = " (gesetzt)"
        fl = noise_floor(db)
        scans.append(dict(z=z, db=db, label=label, zp=zp, w=w, floor=fl,
                          n=len(z), dz=float(np.median(np.diff(z))) * 1e3,
                          how=how))
        print("%-32s  Facette A %10.4f mm%-10s -3 dB %7.1f um   "
              "%+6.1f dB ueber Boden  Bins %8d a %.2f um"
              % (label, zp, how, w * 1e3, top - fl, len(z), scans[-1]["dz"]))

    zps = [s["zp"] for s in scans]
    print("\nPeakpositionen streuen ueber %.1f um (%.4f .. %.4f mm)"
          % ((max(zps) - min(zps)) * 1e3, min(zps), max(zps)))

    marks = load_marks(a.marks_csv) if a.marks_csv else []
    if a.only_marks and not marks:
        raise SystemExit("--only-marks ohne --marks-csv laesst nichts uebrig")
    finds = load_findings(a.findings_csv, len(scans)) if a.findings_csv else None
    base = 0 if a.only_marks else 3
    nrow = base + (2 if marks else 0)
    # Bei zwei Panels ist der 5-Panel-Rand viel zu knapp: 4.5 % von 6.9 Zoll
    # sind 0.3 Zoll und das Achsenlabel faellt aus dem Bild.
    geo = (dict(figsize=(13, 4.6 * nrow), top=0.935, bottom=0.085, hspace=0.24)
           if a.only_marks else
           dict(figsize=(11, 3.45 * nrow), top=0.965, bottom=0.045, hspace=0.36))
    fig, ax = plt.subplots(nrow, 1, figsize=geo["figsize"], squeeze=False,
                           gridspec_kw=dict(hspace=geo["hspace"],
                                            top=geo["top"],
                                            bottom=geo["bottom"], left=0.08,
                                            right=0.98))
    ax = ax[:, 0]

    # dichteste Abtastung nach hinten, sonst deckt ihr Rauschen die anderen zu
    order = sorted(range(len(scans)), key=lambda i: -scans[i]["n"])

    zoom = (min(zps) - a.zoom_pad, max(zps) + a.zoom_pad)
    spans = ((0.0, 3000.0), zoom, (-a.align_pad, a.align_pad))
    titles = ("overview", "facet region (+/- %.2f mm around outermost peaks)" % a.zoom_pad,
              "aligned to own facet -- same structure repeating?")

    for k, ((lo, hi), ttl) in enumerate(zip(spans, titles)) if not a.only_marks else []:
        axk = ax[k]
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
    if not a.only_marks:
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
        bands = [(short_name(nm), s0 * 1e-3 * a.ng[0] / N_FIBER,
                  s1 * 1e-3 * a.ng[1] / N_FIBER) for nm, s0, s1 in marks]
        # Der Edge Coupler klebt an der Facette; wuerde er den Zoom
        # aufspannen, waere der Rest wieder zusammengequetscht.
        if a.marks_zoom:
            zoom_m = tuple(a.marks_zoom)
        else:
            far = [b for b in bands if b[2] > 1.0] or bands
            pad = 0.15 * (max(b[2] for b in far) - min(b[1] for b in far) or 1.0)
            zoom_m = (min(b[1] for b in far) - pad, max(b[2] for b in far) + pad)

        for row, (lo_x, hi_x, ttl) in enumerate((
                (-0.3, a.marks_zmax,
                 "predicted component positions (grey: n_g %.4f-%.4f%s)"
                 % (a.ng[0], a.ng[1],
                    ", light: +/- %.2f mm facet A anchor" % a.facet_unc
                    if a.facet_unc else "")),
                (zoom_m[0], zoom_m[1], "component area, zoomed"))):
            axm = ax[base + row]
            for i, s_ in enumerate(scans):
                c = COLORS[i % len(COLORS)]
                x = s_["z"] - s_["zp"]
                m = (x >= lo_x) & (x <= hi_x)
                axm.plot(x[m], s_["db"][m], lw=0.8, color=c, alpha=0.9,
                         label=s_["label"] if row == 0 else None)
                axm.axhline(s_["floor"], color=c, ls=":", lw=0.8)
            for i, (name, lo_, hi_) in enumerate(bands):
                if a.facet_unc:
                    axm.axvspan(lo_ - a.facet_unc, hi_ + a.facet_unc,
                                color="0.5", alpha=0.10, lw=0)
                axm.axvspan(lo_, hi_, color="0.5", alpha=0.22, lw=0)
                # Bei HHI ist n_g auf 0.01 % bekannt und eine Einzelmarke
                # damit 4 um breit -- als Flaeche waere sie unsichtbar.
                if hi_ - lo_ < 0.01 * (hi_x - lo_x):
                    axm.axvline(0.5 * (lo_ + hi_), color="0.45", lw=0.8,
                                alpha=0.55, zorder=1)
                nrow_lab = 5 if len(bands) > 20 else (4 if len(bands) > 12 else 3)
                axm.annotate(name,
                             xy=(0.5 * (lo_ + hi_),
                                 0.98 - (0.94 / nrow_lab) * (i % nrow_lab)),
                             xycoords=("data", "axes fraction"),
                             fontsize=5.5 if len(bands) > 20 else
                             (6 if len(bands) > 12 else 7),
                             rotation=90, ha="center", va="top",
                             clip_on=True, color="0.2")
            if finds:
                draw_findings(axm, scans, finds, lo_x, hi_x, row == 1)
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
                txt = ("peak in this window over own noise floor:   "
                       + ",   ".join(over))
                if finds:
                    txt += ("\nfound in the DATA -- triangle = spike, bar = "
                            "stretch;  green = matches design,  grey = "
                            "sideband / double bounce,  orange = unexplained")
                axm.text(0.01, 0.06, txt,
                         transform=axm.transAxes, fontsize=9,
                         bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.5))

    fig.suptitle("Facet comparison: " + "  vs  ".join(s["label"] for s in scans))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150)
    print("\nPlot: %s" % a.out)


if __name__ == "__main__":
    main()
