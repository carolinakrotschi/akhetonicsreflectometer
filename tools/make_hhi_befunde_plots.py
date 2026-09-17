#!/usr/bin/env python3
"""Die Befunde-Plots aller HHI-1-Kanaele aus einer Tabelle heraus bauen.

Bisher standen die Aufrufe von `plot_chip_comparison.py` nur in der
Kommandozeile, und mit ihnen die Facette-A-Anker -- die wichtigste Zahl im
ganzen Plot, denn sie legt fest, wo die Bauteilmarken sitzen. Hier stehen
sie an einer Stelle und sind nachlesbar.

Facette A kommt aus zwei verschiedenen Quellen, je nachdem was der Kanal
hergibt (s. `raw_data/hhi_1/CHIP_LENGTHS.md`):

  Loopback (b0, b28)   direkt gemessen -- zwei Facetten, bekannter Abstand,
                       die Topologie erzwingt die Lage
  einseitig (alle      eigene interne Referenz + 1096.31 mm Pigtail,
  anderen)             Unsicherheit +-0.5 mm

Beide Zahlen sind reproduzierbar: sie sind aus den Befunde-CSVs
zurueckgerechnet (der staerkste Stachel dort traegt eine bekannte
Amplitude an einer bekannten Position hinter Facette A) und stimmen mit
den im Log notierten Werten ueberein, z.B. f32 14.09. = 1668.6277 mm.

Aufruf
    python tools/make_hhi_befunde_plots.py --out-dir results/2026-09-17
    python tools/make_hhi_befunde_plots.py --out-dir ... --only f30
"""

import argparse
import csv
import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)

R11 = "results/2026-09-11/2026-09-11-12-16hhi1fiber%d_reflectogram.csv"
R16 = "results/2026-09-16"

# n_g des InP-Waveguides aus den beiden Loopbacks, auf 0.01 % bestimmt.
NG = ("3.4770", "3.4790")

# Facette-A-Anker in mm (absolut) je Scan. Die Unsicherheit steht daneben:
# bei den Loopbacks steckt sie nur in der Wahl des Subpeaks im
# Facettencluster (+-0.08 mm), bei den einseitigen Kanaelen im Pigtail
# (+-0.5 mm).
SCANS = {
    "f4_11":  dict(label="HHI-1 fiber 4 = b28 loopback (11.09.)",
                   csv=R11 % 4, facet=1670.2202, unc=0.08,
                   finds=R16 + "/hhi1_f4_11_befunde.csv"),
    "f4_14":  dict(label="HHI-1 fiber 4 = b28 loopback (14.09.)",
                   csv="results/2026-09-14/2026-09-14-09-03_fiber4_reflectogram.csv",
                   facet=1670.3037, unc=0.08,
                   finds=R16 + "/hhi1_f4_14_befunde.csv"),
    "f32_11": dict(label="HHI-1 fiber 32 = b0 loopback (11.09.)",
                   csv=R11 % 32, facet=1668.3832, unc=0.08,
                   finds=R16 + "/hhi1_f32_11_befunde.csv"),
    "f32_14": dict(label="HHI-1 fiber 32 = b0 loopback (14.09.)",
                   csv="results/2026-09-14/2026-09-14-09-35_fiber32_reflectogram.csv",
                   facet=1668.6277, unc=0.08,
                   finds=R16 + "/hhi1_f32_14_befunde.csv"),
    "f7":     dict(label="HHI-1 fiber 7 = b25 stub (11.09.)",
                   csv=R11 % 7, facet=1670.1890, unc=0.50,
                   finds=R16 + "/hhi1_f7_befunde.csv"),
    "f27":    dict(label="HHI-1 fiber 27 = b5 (11.09.)",
                   csv=R11 % 27, facet=1670.3510, unc=0.50,
                   finds=R16 + "/hhi1_f27_befunde.csv"),
    "f30":    dict(label="HHI-1 fiber 30 = b2 (11.09.)",
                   csv=R11 % 30, facet=1670.2380, unc=0.50,
                   finds=R16 + "/hhi1_f30_befunde.csv"),
    "f6":     dict(label="HHI-1 fiber 6 = b26 (11.09.)",
                   csv=R11 % 6, facet=1670.3870, unc=0.50,
                   finds=R16 + "/hhi1_f6_befunde.csv"),
    "f23":    dict(label="HHI-1 fiber 23 = b9 (14.09.)",
                   csv="results/2026-09-14/2026-09-14-09-58_fiber23_reflectogram.csv",
                   facet=1670.4710, unc=0.50,
                   finds=R16 + "/hhi1_f23_befunde.csv"),
}

# Ein Plot je Kanal; die zweimal gemessenen Loopbacks tragen beide Tage.
# Das Fenster des untersten Panels waehlt `zoom_for` je Kanal, s. dort.
PLOTS = [
    dict(name="hhi1_fiber4_b28",  scans=["f4_11", "f4_14"],
         marks=R16 + "/hhi1_b28_fiber4_bauteile.csv"),
    dict(name="hhi1_fiber32_b0",  scans=["f32_11", "f32_14"],
         marks=R16 + "/hhi1_b0_fiber32_bauteile.csv"),
    dict(name="hhi1_fiber7_b25",  scans=["f7"],
         marks=R16 + "/hhi1_b25_fiber7_bauteile.csv"),
    dict(name="hhi1_fiber27_b5",  scans=["f27"],
         marks=R16 + "/hhi1_b5_fiber27_bauteile.csv"),
    dict(name="hhi1_fiber30_b2",  scans=["f30"],
         marks=R16 + "/hhi1_b2_fiber30_bauteile.csv"),
    dict(name="hhi1_fiber6_b26",  scans=["f6"],
         marks=R16 + "/hhi1_b26_fiber6_bauteile.csv"),
    dict(name="hhi1_fiber23_b9",  scans=["f23"],
         marks=R16 + "/hhi1_b9_fiber23_bauteile.csv"),
]


def zmax_for(marks_csv, find_csvs, pad=1.0):
    """Wie weit hinter Facette A das Bauteilpanel reicht.

    Bis zum letzten Bauteil UND bis zum letzten Fund -- ein Fund hinter dem
    Modellende (f23/b9, f30 bei 35.7 mm) ist gerade der interessante Fall
    und darf nicht aus dem Bild fallen."""
    far = [float(r["dz_to_mm"]) for r in csv.DictReader(open(marks_csv))]
    for p in find_csvs:
        far += [float(r["z_to_mm"]) for r in csv.DictReader(open(p))]
    return round(max(far) + pad, 1)


def dense_window(marks_csv, width):
    """Das `width` mm breite Fenster mit den meisten Marken darin.

    Bei 33 Marken auf 36 mm ist "von der ersten bis zur letzten" kein Zoom.
    Interessant ist, wo die Bauteile dicht stehen -- dort sind die
    Beschriftungen im Panel darueber unlesbar und genau dort muss man
    naeher heran."""
    pos = sorted(0.5 * (float(r["dz_from_mm"]) + float(r["dz_to_mm"]))
                 for r in csv.DictReader(open(marks_csv)))
    best = (0, pos[0])
    for p0 in pos:
        n = sum(1 for q in pos if p0 <= q <= p0 + width)
        if n > best[0]:
            best = (n, p0)
    lo = best[1] - 0.08 * width
    return (lo, lo + width)


def zoom_for(marks_csv, zmax, min_factor=2.0):
    """Fenster des untersten Panels, oder None fuer die Toolautomatik.

    Die Automatik in `plot_chip_comparison.py` spannt von der ersten bis
    zur letzten Marke. Bei einem duennen Kanal (b28: drei Marken auf
    4.4 mm) ist das ein echter Zoom. Bei b2 und b26 stehen die Marken
    ueber die ganzen 36 mm, und dann zeigt das unterste Panel fast dasselbe
    wie das darueber -- dort stattdessen das dichteste Drittel.

    Die Schwelle ist "zoomt um weniger als `min_factor`", nicht eine Liste
    von Kanalnamen: so bleibt die Regel gueltig, wenn ein Kanal dazukommt."""
    pos = [(float(r["dz_from_mm"]), float(r["dz_to_mm"]))
           for r in csv.DictReader(open(marks_csv))]
    lo, hi = min(p[0] for p in pos), max(p[1] for p in pos)
    auto_w = (hi - lo) * 1.3          # 15 % Rand links und rechts
    if auto_w * min_factor > zmax + 0.3:
        return dense_window(marks_csv, (zmax + 0.3) / 3.0)
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--only", action="append", default=None, metavar="SCANKEY",
                    help="nur die Plots bauen, in denen dieser Scan vorkommt "
                         "(z.B. f30); mehrfach angebbar")
    ap.add_argument("--wide-panels", type=int, default=2,
                    help="Kontextpanels vor den zwei Bauteilpanels; 2 gibt "
                         "vier Panels von der ganzen Achse bis in den Zoom "
                         "(default 2)")
    ap.add_argument("--relative", action="store_true",
                    help="x-Achse wieder 'hinter Facette A' statt absolut")
    ap.add_argument("--suffix", default=None,
                    help="Dateinamenszusatz (default '_befunde_zoom4' bzw. "
                         "'_befunde' bei --wide-panels 0)")
    a = ap.parse_args()

    sfx = a.suffix
    if sfx is None:
        sfx = "_befunde_zoom%d" % (a.wide_panels + 2) if a.wide_panels else "_befunde"
    os.makedirs(a.out_dir, exist_ok=True)

    for pl in PLOTS:
        if a.only and not any(k in a.only for k in pl["scans"]):
            continue
        keys = pl["scans"]
        finds = [SCANS[k]["finds"] for k in keys]
        zmax = zmax_for(pl["marks"], finds)
        out = os.path.join(a.out_dir, pl["name"] + sfx + ".png")
        cmd = [sys.executable, os.path.join("tools", "plot_chip_comparison.py"),
               "--only-marks", "--marks-csv", pl["marks"],
               "--ng", NG[0], NG[1],
               "--marks-zmax", "%.1f" % zmax,
               "--facet-unc", "%.2f" % SCANS[keys[0]]["unc"],
               "--wide-panels", str(a.wide_panels),
               "--out", out]
        if not a.relative:
            cmd.append("--absolute-x")
        zm = zoom_for(pl["marks"], zmax)
        if zm:
            cmd += ["--marks-zoom", "%.2f" % zm[0], "%.2f" % zm[1]]
        for k in keys:
            cmd += ["--scan", SCANS[k]["label"], SCANS[k]["csv"],
                    "--facet", "%.4f" % SCANS[k]["facet"]]
        for p in finds:
            cmd += ["--findings-csv", p]
        print("\n=== %s ===" % pl["name"])
        subprocess.run(cmd, cwd=_ROOT, check=True)


if __name__ == "__main__":
    main()
