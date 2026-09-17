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

Mit `--marks-csv` kommen zwei Bauteilpanels dazu. `--wide-panels K` setzt
K Kontextpanels davor: das erste zeigt die GANZE Achse, die weiteren zoomen
geometrisch bis zum Bauteilbereich hinein, so dass man von 3000 mm in vier
Schritten in den Zoom kommt und sieht, wo der Chip auf der Achse steht.
`--absolute-x` schreibt dabei die absolute Distanz vom Instrument auf die
Achse statt "hinter Facette A" -- die Bauteilmarken wandern dann um den
Anker mit.

Aufruf
    python tools/plot_chip_comparison.py \
        --scan "MAP2672 fiber1" raw_data/2672_ligentechhi_..._reflectogram.csv \
        --scan "MAP2680 fiber1 (1520-1570nm)" raw_data/..._fiber1_reflectogram.csv \
        --scan "MAP2680 fiber1 (voller Span)" raw_data/..._fiber1full_reflectogram.csv \
        --peak-window 1575 1590 \
        --out results/2026-09-16/fiber1_chip_vergleich.png

Fuer die HHI-1-Befunde-Plots ruft `tools/make_hhi_befunde_plots.py` dieses
Skript mit den Facette-A-Ankern aller Kanaele auf -- dort nachsehen statt
die Kommandozeile von Hand zusammenzusetzen.
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


# Abkuerzung -> Klartext. Gezeigt wird im Plot nur, was dort auch vorkommt.
GLOSS = [
    ("chip facet", "chip facet = Chipfacette (InP/Luft), der einzige harte Reflektor"),
    ("SSC end", "SSC = spot size converter; 'SSC end' = sein inneres Ende, Uebergang in den Waveguide"),
    ("loop", "loop = Wellenleiterschleife zwischen den zwei SSC (reflektiert nicht, streut nur)"),
    ("BJ", "BJ = butt joint, Stoss aktiv/passiv"),
    ("SOA", "SOA = Halbleiterverstaerker"),
    ("MMI2x2", "MMI2x2 = 2-auf-2-Koppler"),
    ("MMI1x2", "MMI1x2 = 1-auf-2-Leistungsteiler"),
    ("crossing", "crossing = Wellenleiterkreuzung"),
    ("PMTO", "PMTO = thermo-optischer Phasenschieber"),
    ("ISO", "ISO = elektrische Isolationssektion"),
    ("WGMETx", "WGMETx = Metallbruecke UEBER dem Waveguide (optisch fast nichts)"),
    ("EdgeCoupler", "EdgeCoupler = Kantenkoppler"),
    ("pitch split", "pitch split = Faecher, der den Kanalabstand aufweitet"),
    ("MZM", "MZM = Mach-Zehnder-Modulator"),
    ("Heater", "Heater = Heizelektrode (optisch unsichtbar)"),
]
SYMS = [("+", "'+n' = n weitere Bauteilsorten in derselben Markengruppe"),
        ("SL", "SL = Seitenlinie des Instruments"),
        ("2x", "2x = Doppelbounce"),
        ("?", "? = unerklaert")]


def glossary(labels):
    """Nur die Abkuerzungen erklaeren, die in diesem Plot auch auftauchen."""
    txt = " ".join(labels)
    out = [v for k, v in GLOSS if k in txt]
    out += [v for k, v in SYMS if k in txt]
    return out


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


def draw_findings(axm, scans, finds, lo_x, hi_x, label_them, foff=None):
    """Stacheln als Dreieck auf dem Peak, Buckel als Balken darunter.

    Gezeichnet wird, was in den DATEN gefunden wurde -- im Gegensatz zu den
    grauen Marken, die aus dem GDS kommen. Farbe = Deutung.

    `foff` verschiebt die Funde auf die Plotachse: die CSV von
    `hhi_identify_peaks.py` zaehlt ab Facette A, die Achse kann absolut
    sein. Ein Fund gehoert zu SEINEM Scan, also bekommt jede Spur ihren
    eigenen Offset (bei zwei Messtagen liegt Facette A nicht gleich)."""
    y0, y1 = axm.get_ylim()
    if foff is None:
        foff = [0.0] * len(scans)
    for i, (s_, rows) in enumerate(zip(scans, finds)):
        lane = y0 + (0.05 + 0.045 * i) * (y1 - y0)
        for r in rows:
            v = r["verdict"]
            c = FIND_COLOR.get(v, "0.45")
            if r["kind"] == "Stachel":
                zz = float(r["z_from_mm"]) + foff[i]
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
                z0 = float(r["z_from_mm"]) + foff[i]
                z1 = float(r["z_to_mm"]) + foff[i]
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
    ap.add_argument("--absolute-x", action="store_true",
                    help="x-Achse als ABSOLUTE Distanz vom Instrument statt "
                         "'hinter Facette A'. Die Marken werden dann um "
                         "Facette A der ERSTEN Spur verschoben, die Funde je "
                         "Spur um deren eigene Facette. Achtung: absolute "
                         "Lagen sind zwischen Messsitzungen nur auf ~1 mm "
                         "vergleichbar (die interne Referenz wandert), "
                         "Abstaende dagegen auf unter 1 um")
    ap.add_argument("--wide-panels", type=int, default=0, metavar="K",
                    help="K zusaetzliche Kontextpanels VOR den zwei "
                         "Bauteilpanels: das erste zeigt die ganze Achse, "
                         "die weiteren zoomen geometrisch bis zum "
                         "Bauteilbereich hinein. K=2 gibt vier Panels von "
                         "der ganzen Achse bis in den Zoom (default 0 = "
                         "Verhalten wie bisher)")
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
    # Das Glossar muss VOR der Figur stehen: es braucht Platz am unteren
    # Rand, und `fig.subplots_adjust` ist wirkungslos, sobald top/bottom im
    # gridspec stehen (matplotlib nimmt dann die dort gespeicherten Werte).
    # Vorher wurde der Platz deshalb nur scheinbar reserviert.
    gl = []
    if marks:
        used = [short_name(nm) for nm, _, _ in marks]
        if finds:
            used += [r["candidate"] for rr in finds for r in rr]
            used += [FIND_TAG.get(r["verdict"], "") for rr in finds for r in rr]
        gl = glossary(used)
    gl_txt = "   |   ".join(gl)
    gl_lin = (1 + len(gl_txt) // 190) if gl_txt else 0

    nwide = max(a.wide_panels, 0) if marks else 0
    base = 0 if a.only_marks else 3
    nrow = base + (2 + nwide if marks else 0)
    # Bei zwei Panels ist der 5-Panel-Rand viel zu knapp: 4.5 % von 6.9 Zoll
    # sind 0.3 Zoll und das Achsenlabel faellt aus dem Bild. Deshalb den
    # Rand in ZOLL festlegen und erst dann in Bruchteile umrechnen -- sonst
    # wird bei vier Panels aus demselben Bruchteil der doppelte Rand.
    bot_in = 0.50 + 0.15 * gl_lin      # Achsenlabel + Glossarzeilen
    if a.only_marks:
        h = 4.6 * nrow
        geo = dict(figsize=(13, h), top=1 - 0.6 / h, bottom=bot_in / h,
                   hspace=0.30)
    else:
        h = 3.45 * nrow
        geo = dict(figsize=(11, h), top=0.965, bottom=bot_in / h, hspace=0.36)
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
        #
        # `anchor` legt fest, was 0 auf der Achse heisst: bei --absolute-x
        # die Facette A der ersten Spur (auf der Achse steht dann die
        # absolute Distanz vom Instrument), sonst 0 (Distanz hinter der
        # Facette). Alles weiter unten rechnet in PLOTKOORDINATEN.
        anchor = scans[0]["zp"] if a.absolute_x else 0.0
        xoff = [(0.0 if a.absolute_x else -s_["zp"]) for s_ in scans]
        foff = [(s_["zp"] if a.absolute_x else 0.0) for s_ in scans]
        bands = [(short_name(nm), anchor + s0 * 1e-3 * a.ng[0] / N_FIBER,
                  anchor + s1 * 1e-3 * a.ng[1] / N_FIBER)
                 for nm, s0, s1 in marks]
        # Der Edge Coupler klebt an der Facette; wuerde er den Zoom
        # aufspannen, waere der Rest wieder zusammengequetscht.
        if a.marks_zoom:
            zoom_m = (anchor + a.marks_zoom[0], anchor + a.marks_zoom[1])
        else:
            far = [b for b in bands if b[2] - anchor > 1.0] or bands
            pad = 0.15 * (max(b[2] for b in far) - min(b[1] for b in far) or 1.0)
            zoom_m = (min(b[1] for b in far) - pad, max(b[2] for b in far) + pad)
            # Ein Kanal mit einer EINZIGEN Marke (b25 ist ein Stumpf) hat
            # keine Ausdehnung: "erste bis letzte Marke" waere dort 2 um
            # breit und das Panel bliebe leer. Dann stattdessen ein Drittel
            # des Bauteilbereichs um die Marke herum.
            span_det = a.marks_zmax + 0.3
            if zoom_m[1] - zoom_m[0] < 0.2 * span_det:
                w = span_det / 3.0
                c0 = 0.5 * (zoom_m[0] + zoom_m[1])
                lo_z = min(max(c0 - 0.5 * w, anchor - 0.3),
                           anchor + a.marks_zmax - w)
                zoom_m = (lo_z, lo_z + w)

        # Die zwei Bauteilpanels wie bisher ...
        det = [(anchor - 0.3, anchor + a.marks_zmax,
                "predicted component positions (grey: n_g %.4f-%.4f%s)"
                % (a.ng[0], a.ng[1], "")),
               (zoom_m[0], zoom_m[1], "component area, zoomed")]
        # ... und davor K Kontextpanels: das erste die ganze Achse, die
        # weiteren geometrisch bis zum Bauteilbereich hinein. Geometrisch,
        # nicht linear: von 3000 mm auf 36 mm ist der halbe Weg 330 mm; der
        # linear gemittelte waere 1518 mm und zeigte zweimal fast dasselbe.
        lo_full = min(float(s_["z"][0]) + o for s_, o in zip(scans, xoff))
        hi_full = max(float(s_["z"][-1]) + o for s_, o in zip(scans, xoff))
        s_det = det[0][1] - det[0][0]
        s_full = max(hi_full - lo_full, s_det)
        cen = 0.5 * (det[0][0] + det[0][1])
        wide = []
        for k in range(nwide):
            if k == 0:
                wide.append((lo_full, hi_full,
                             "whole axis %.0f-%.0f mm (envelope over all bins)"
                             % (lo_full, hi_full)))
                continue
            w = s_full * (s_det / s_full) ** (k / float(nwide))
            lo_w = min(max(cen - 0.5 * w, lo_full), hi_full - w)
            wide.append((lo_w, lo_w + w,
                         "zoomed in: %.0f mm window around the chip" % w))
        panels = wide + det
        print("\nPanels (%s):   %s"
              % ("absolute" if a.absolute_x else "hinter Facette A",
                 "  ->  ".join("%.2f-%.2f mm (%.1f breit)"
                               % (p0, p1, p1 - p0) for p0, p1, _ in panels)))

        for row, (lo_x, hi_x, ttl) in enumerate(panels):
            axm = ax[base + row]
            first, last = row == 0, row == len(panels) - 1
            # Auf einem 3000-mm-Panel sind 33 Marken ein einziger Strich und
            # die Funde liegen alle uebereinander -- dort nur der Bereich als
            # Flaeche; Einzelheiten erst, wo sie lesbar sind.
            detail = (hi_x - lo_x) <= 3.0 * s_det
            for i, s_ in enumerate(scans):
                c = COLORS[i % len(COLORS)]
                x = s_["z"] + xoff[i]
                m = (x >= lo_x) & (x <= hi_x)
                zz, dd = x[m], s_["db"][m]
                if not detail and len(zz) > 6000:
                    zz, dd = envelope(zz, dd, 3000)
                axm.plot(zz, dd, lw=0.8, color=c, alpha=0.9,
                         label=s_["label"] if first else None)
                axm.axhline(s_["floor"], color=c, ls=":", lw=0.8)
            if detail:
                for i, (name, lo_, hi_) in enumerate(bands):
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
            else:
                axm.axvspan(det[0][0], det[0][1], color="0.5", alpha=0.18,
                            lw=0, zorder=1)
                axm.axvline(anchor, color="0.35", lw=0.9, alpha=0.8, zorder=1)
                axm.text(0.012, 0.055,
                         "grey: chip / component region %.1f-%.1f mm "
                         "(%.1f mm wide, shown in the panels below)"
                         % (det[0][0], det[0][1], det[0][1] - det[0][0]),
                         transform=axm.transAxes, ha="left", va="bottom",
                         fontsize=8, color="0.25",
                         bbox=dict(fc="white", ec="none", alpha=0.8,
                                   pad=1.5))
            if finds and detail:
                draw_findings(axm, scans, finds, lo_x, hi_x, last, foff)
            if a.facet_unc and detail:
                # Als EIN Massstab, nicht als Hof um jede Marke: bei 24
                # Marken mit je +-0.5 mm ist sonst die ganze Achse grau und
                # es sieht so aus, als ueberlappten die Bauteile.
                x0 = lo_x + 0.02 * (hi_x - lo_x)
                axm.annotate("", xy=(x0, 0.93), xytext=(x0 + 2 * a.facet_unc, 0.93),
                             xycoords=("data", "axes fraction"),
                             arrowprops=dict(arrowstyle="<->", color="0.35",
                                             lw=1.1))
                axm.annotate("Anker Facette A +/- %.2f mm\n(alle Marken "
                             "wandern gemeinsam)" % a.facet_unc,
                             xy=(x0, 0.915), xycoords=("data", "axes fraction"),
                             ha="left", va="top", fontsize=7, color="0.25")
            axm.set_xlim(lo_x, hi_x)
            axm.set_ylabel("amplitude [dB]")
            axm.set_xlabel("absolute distance from instrument [mm]"
                           if a.absolute_x else "distance behind facet A [mm]")
            axm.text(0.5, 1.0, ttl, transform=axm.transAxes, ha="center",
                     va="top", fontsize=10,
                     bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
            axm.grid(alpha=0.3)
            if first:
                axm.legend(loc="upper right", fontsize=8)
                if a.absolute_x:
                    # Zwei Zeilen, nicht eine: einzeilig laeuft der Satz mit
                    # einem langen Spurnamen rechts aus dem Bild.
                    note = ("absolute axis: DISTANCES are good to <1 um, but "
                            "absolute positions shift by ~1 mm between "
                            "measurement sessions")
                    note += "\n(the internal reference drifts)."
                    if len(scans) > 1:
                        note += ("  The grey marks are anchored to facet A of "
                                 "%s." % scans[0]["label"])
                    axm.text(0.012, 0.115, note, transform=axm.transAxes,
                             fontsize=8, ha="left", va="bottom",
                             bbox=dict(fc="white", ec="none", alpha=0.8,
                                       pad=1.5))
            if last:
                over = []
                for i, s_ in enumerate(scans):
                    x = s_["z"] + xoff[i]
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

    if gl_txt:
        fig.text(0.008, 0.008, gl_txt, fontsize=7, color="0.25",
                 va="bottom", wrap=True)
    fig.suptitle("Facet comparison: " + "  vs  ".join(s["label"] for s in scans))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    fig.savefig(a.out, dpi=150)
    print("\nPlot: %s" % a.out)


if __name__ == "__main__":
    main()
