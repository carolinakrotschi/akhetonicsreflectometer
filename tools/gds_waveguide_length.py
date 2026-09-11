"""Exakte Waveguide-Laenge aus einer GDS-Datei (ohne KLayout-GUI).

Beispiel -- Schleife oben links in HHI_RUN_1.gds (Coupler-Port zu
Coupler-Port), Waveguide auf Layer 12/0:

    python tools/gds_waveguide_length.py "Akhelab Screenshots/HHI_RUN_1.gds" \
        --layer 12/0 --window 1180 3670 1470 3915 \
        --out results/2026-09-11/topleft_loop

METHODE
   Jedes Waveguide-Polygon ist ein Band konstanter Breite w mit
   Mittellinienlaenge L und stumpfen Enden, also exakt

       A = L*w              P = 2*L + 2*w

   -> L und w sind die Wurzeln von t^2 - (P/2)*t + A = 0 (groessere = L).
   Gilt auch fuer Boegen, weil sich der w/2-Ueberhang von Aussen- und
   Innenbogen weghebt. Die Breite muss also nicht bekannt sein, sie faellt
   als Kontrollwert mit heraus.

   Boegen liegen in GDS als Sekantenzug vor und sind damit systematisch zu
   kurz. Fuer Polygone mit vielen Stuetzpunkten wird deshalb ein Kreis
   gefittet und auf den echten Bogen hochkorrigiert:

       L_exakt = L_sekant * a/sin(a)    mit a = theta/(2*n_seg)
"""

import argparse
import math
import os

import numpy as np
import gdstk


def _perimeter(poly):
    q = poly.points
    return sum(math.dist(q[i], q[(i + 1) % len(q)]) for i in range(len(q)))


def ribbon(poly):
    """(Flaeche, Umfang, Mittellinienlaenge, Breite) eines Bandpolygons."""
    area = abs(poly.area())
    per = _perimeter(poly)
    half = per / 2.0
    disc = half * half - 4.0 * area
    if disc < 0:
        return area, per, float("nan"), float("nan")
    root = math.sqrt(disc)
    return area, per, (half + root) / 2.0, (half - root) / 2.0


def fit_circle(points):
    pts = np.asarray(points)
    x, y = pts[:, 0], pts[:, 1]
    mat = np.c_[2 * x, 2 * y, np.ones(len(x))]
    sol, *_ = np.linalg.lstsq(mat, x ** 2 + y ** 2, rcond=None)
    xc, yc = sol[0], sol[1]
    return xc, yc, math.sqrt(sol[2] + xc ** 2 + yc ** 2)


def analyse(poly, min_pts_for_arc=8):
    """Ein Polygon -> dict mit Laenge, Breite und (falls Bogen) R/theta."""
    area, per, length, width = ribbon(poly)
    (x0, y0), (x1, y1) = poly.bounding_box()
    row = dict(x0=x0, y0=y0, x1=x1, y1=y1, npoints=len(poly.points),
               area=area, perimeter=per, length_secant=length, width=width,
               kind="Gerade", radius=float("nan"), angle_deg=float("nan"),
               n_seg=0, length_exact=length)
    if len(poly.points) > min_pts_for_arc and length == length:
        _, _, radius = fit_circle(poly.points)
        theta = length / radius
        n_seg = max(len(poly.points) // 2 - 1, 1)
        a = theta / (2.0 * n_seg)
        row.update(kind="Bogen", radius=radius, angle_deg=math.degrees(theta),
                   n_seg=n_seg,
                   length_exact=length * a / math.sin(a) if a > 0 else length)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gds")
    ap.add_argument("--layer", default="12/0", help="Layer/Datatype, z.B. 12/0")
    ap.add_argument("--window", nargs=4, type=float, required=True,
                    metavar=("X0", "Y0", "X1", "Y1"),
                    help="nur Polygone komplett in diesem Fenster (um)")
    ap.add_argument("--out", help="Prefix fuer CSV + PNG (optional)")
    args = ap.parse_args()

    lay, dt = (int(v) for v in args.layer.split("/"))
    x0w, y0w, x1w, y1w = args.window

    lib = gdstk.read_gds(args.gds)
    top = lib.top_level()[0]
    flat = top.copy("_flat").flatten()
    sel = []
    for p in flat.polygons:
        if (p.layer, p.datatype) != (lay, dt):
            continue
        (a, b), (c, d) = p.bounding_box()
        if a > x0w and c < x1w and b > y0w and d < y1w:
            sel.append(p)

    rows = [analyse(p) for p in sel]
    rows.sort(key=lambda r: (-r["y1"], r["x0"]))

    hdr = ("%3s %10s %10s %6s %6s %12s %12s %9s %10s %9s %5s"
           % ("#", "x0", "y0", "dx", "dy", "L_sekant", "L_exakt", "Breite",
              "R", "Winkel", "nseg"))
    print("GDS: %s   Layer %d/%d   Fenster %s" % (args.gds, lay, dt, args.window))
    print("Polygone im Fenster: %d\n" % len(rows))
    print(hdr)
    print("-" * len(hdr))
    for i, r in enumerate(rows):
        print("%3d %10.3f %10.3f %6.1f %6.1f %12.4f %12.4f %9.4f %10.4f %9.3f %5d"
              % (i, r["x0"], r["y0"], r["x1"] - r["x0"], r["y1"] - r["y0"],
                 r["length_secant"], r["length_exact"], r["width"],
                 r["radius"], r["angle_deg"], r["n_seg"]))

    tot_s = sum(r["length_secant"] for r in rows)
    tot_e = sum(r["length_exact"] for r in rows)
    widths = [r["width"] for r in rows if r["width"] == r["width"]]
    print("-" * len(hdr))
    print("Summe Sekantenlaengen : %12.4f um" % tot_s)
    print("Summe bogenkorrigiert : %12.4f um   (+%.4f um = %.4f %%)"
          % (tot_e, tot_e - tot_s, 100 * (tot_e - tot_s) / tot_s))
    print("Breite min/max        : %.4f / %.4f um" % (min(widths), max(widths)))

    merged = gdstk.boolean(sel, [], "or", precision=1e-4)
    print("zusammenhaengende Komponenten: %d" % len(merged))
    for m in merged:
        print("   A=%10.3f um^2 -> L=A/w=%10.4f um  bbox=%s"
              % (abs(m.area()), abs(m.area()) / np.mean(widths),
                 tuple(round(v, 2) for bb in m.bounding_box() for v in bb)))

    if args.out:
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        import csv
        with open(args.out + "_segments.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print("\nCSV: %s_segments.csv" % args.out)

        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon as MPoly
        fig, ax = plt.subplots(figsize=(10, 8))
        for p in flat.polygons:
            if (p.layer, p.datatype) == (1021, 1):
                ax.add_patch(MPoly(p.points, closed=True, fc="#8fdc8f",
                                   ec="green", lw=.4, alpha=.55))
        for i, (p, r) in enumerate(zip(sorted(sel, key=lambda q: (-q.bounding_box()[1][1], q.bounding_box()[0][0])), rows)):
            ax.add_patch(MPoly(p.points, closed=True, fc="crimson", ec="k", lw=.3))
            ax.annotate("%d: %.2f um" % (i, r["length_exact"]),
                        ((r["x0"] + r["x1"]) / 2, (r["y0"] + r["y1"]) / 2),
                        fontsize=7, ha="center", va="center", color="navy")
        ax.set_xlim(x0w, x1w)
        ax.set_ylim(y0w, y1w)
        ax.set_aspect("equal")
        ax.grid(alpha=.3)
        ax.set_xlabel("x [um]")
        ax.set_ylabel("y [um]")
        ax.set_title("%s  Layer %d/%d\nGesamtlaenge Mittellinie = %.4f um (%d Segmente)"
                     % (os.path.basename(args.gds), lay, dt, tot_e, len(rows)))
        fig.savefig(args.out + ".png", dpi=130, bbox_inches="tight")
        print("PNG: %s.png" % args.out)


if __name__ == "__main__":
    main()
