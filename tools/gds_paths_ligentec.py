"""Weg zwischen zwei Facettenkanaelen im LIGENTEC-GDS, gegen Transmissionspeaks.

Fuer eine Transmissionsmessung Faser A -> Chip -> Faser B: welche Weglaenge
sagt das GDS fuer Kanal A -> Kanal B voraus, und wo muesste der Peak dann
relativ zu einer Loopback-Referenz im Reflektogramm stehen?

Verfahren
    1. Graph wie `gds_netlist_ligentec.build(bridge_boxes=True)`, zusaetzlich
       jede MMI-Zellinstanz durchlaessig (die MMI_2x2 der MZM-Zellen tragen
       keine Blackbox auf 13/2, ohne das endet jeder MZM am Kombinierer).
    2. Dijkstra ab Kanal A, Pfad zu Kanal B, entlang des Pfads die
       Zellinstanzen ausgeben.
    3. Spruenge zwischen zwei Ports auf DERSELBEN Seite eines MMI
       (Abstand < 10 um, gleiche y) markieren: die Bruecke verbindet alle
       Ports einer Box paarweise, physikalisch laeuft Licht aber nur von
       einer Seite zur anderen. Ein Pfad mit solchen Spruengen ist kein
       Entwurfsweg, sondern nur ueber Rueckreflex/Streuung moeglich.
    4. Erwartete Lage gegen die Referenz (Transmission, Reflexionskonvention):
           dz = (L_AB + 2*COUPLER - L_ref) * n_g / (2 * n_Faser)
       mit L_ref = Facette-zu-Facette-Weg der Referenz (Loop 1<->2 bzw.
       94<->95: 1039.5 um). Arrayfaserlaengen streuen um einige mm, das
       verschiebt alles um +-~3 mm z (Loop A gegen B: 2.70 mm).

Grenzen: nur der SiN-Teil (y < SIN_YMAX). Wege durch den aufgeklebten
InP-Chip sind nicht im Graphen. Laengen in Blackboxen sind Luftlinie.

Aufruf
    python tools/gds_paths_ligentec.py --from 48 --to 47 49 \
        --peaks-dz 9.09 14.79 23.22 46.66
"""

import argparse
import csv
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gdstk                                                      # noqa: E402
import gds_netlist_ligentec as G                                  # noqa: E402
from gds_trace_ligentec import cell_instances, dijkstra_pred      # noqa: E402

L_REF_UM = 1039.5      # Loop 1<->2 / 94<->95, Facette zu Facette (logs/2026-09-15.md)
SAME_SIDE_UM = 10.0


def build_graph(gds):
    nodes, edges, _ = G.build(gds, bridge_boxes=True)
    lib = gdstk.read_gds(gds)
    top = lib.top_level()[0]
    inst = cell_instances(top, G.SIN_YMAX)
    pts = np.array([n[:2] for n in nodes])
    seen = set()
    for name, _, b in inst:
        if not re.search(r"MMI", name, re.I) or math.hypot(b[2] - b[0], b[3] - b[1]) > 600:
            continue
        key = tuple(round(v, 1) for v in b)
        if key in seen:
            continue
        seen.add(key)
        sel = np.where((pts[:, 0] >= b[0] - 2) & (pts[:, 0] <= b[2] + 2)
                       & (pts[:, 1] >= b[1] - 2) & (pts[:, 1] <= b[3] + 2))[0]
        for ii in range(len(sel)):
            for jj in range(ii + 1, len(sel)):
                i, j = int(sel[ii]), int(sel[jj])
                d = math.dist(pts[i], pts[j])
                if d > 1.0:
                    edges.setdefault(i, []).append((j, d))
                    edges.setdefault(j, []).append((i, d))
    small = [(n, dep, b) for (n, dep, b) in inst
             if math.hypot(b[2] - b[0], b[3] - b[1]) < 3000]
    return nodes, edges, small


def cell_at(inst, p):
    hits = sorted((dep, n) for n, dep, b in inst
                  if b[0] - 1 <= p[0] <= b[2] + 1 and b[1] - 1 <= p[1] <= b[3] + 1)
    return hits[-1][1] if hits else "-"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gds", default=G.GDS)
    ap.add_argument("--from", dest="src", type=int, required=True, help="Chipkanal (Faser = 96 - Kanal)")
    ap.add_argument("--to", type=int, nargs="+", required=True)
    ap.add_argument("--ng", type=float, nargs=3, default=[1.80, 1.90, 2.00])
    ap.add_argument("--peaks-dz", type=float, nargs="*", default=[],
                    help="gemessene Peaks, mm hinter der Referenz (Transmission)")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    nodes, edges, inst = build_graph(a.gds)
    _, dist, pred = dijkstra_pred(nodes, edges, (G.channel_x(a.src), G.WG_START_Y), 2e5)
    rows = []
    for dst in a.to:
        t = min(range(len(nodes)),
                key=lambda i: math.dist(nodes[i][:2], (G.channel_x(dst), G.WG_START_Y)))
        print("\n=== Kanal %d (Faser %d) -> Kanal %d (Faser %d)"
              % (a.src, 96 - a.src, dst, 96 - dst))
        if t not in dist:
            print("  nicht verbunden")
            continue
        path, i = [], t
        while i is not None:
            path.append(i)
            i = pred[i]
        path.reverse()
        jumps, last = [], None
        for k in range(1, len(path)):
            p0, p1 = nodes[path[k - 1]][:2], nodes[path[k]][:2]
            c = cell_at(inst, p1)
            if 1.0 < math.dist(p0, p1) < SAME_SIDE_UM and abs(p0[1] - p1[1]) < 1.0:
                # nicht cell_at: eine tiefer verschachtelte Container-Zelle
                # (XGMx3_SiN) kann die MMI-Box ueberdecken
                mmi = [n for n, _, b in inst if re.search("MMI", n)
                       and all(b[0] - 1 <= p[0] <= b[2] + 1 and b[1] - 1 <= p[1] <= b[3] + 1
                               for p in (p0, p1))]
                if mmi:
                    jumps.append((dist[path[k]], mmi[0]))
            if c != last:
                print("  s = %8.1f um  (%7.1f, %7.1f)  %s" % (dist[path[k]], p1[0], p1[1], c))
                last = c
        L = dist[t] + 2 * G.COUPLER_UM
        dz = [(L - L_REF_UM) * ng / (2 * G.N_FIBER) * 1e-3 for ng in a.ng]
        print("  Weg Facette-Facette %.1f um; dz gegen Loop %.2f / %.2f / %.2f mm (n_g %s)"
              % (L, *dz, " / ".join("%.2f" % g for g in a.ng)))
        if jumps:
            print("  ACHTUNG: %d Sprung/Spruenge zwischen Ports derselben MMI-Seite -> "
                  "kein Entwurfsweg:" % len(jumps))
            for s, c in jumps:
                print("    bei s = %.1f um in %s" % (s, c))
        else:
            print("  Entwurfsweg (keine Spruenge auf derselben MMI-Seite)")
        for pk in a.peaks_dz:
            ng = pk * 2 * G.N_FIBER / ((L - L_REF_UM) * 1e-3)
            print("  Peak dz %.2f mm -> implizites n_g %.3f" % (pk, ng))
        rows.append([a.src, dst, round(L, 1), *[round(v, 3) for v in dz], len(jumps)])

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["from_ch", "to_ch", "L_facet_um",
                        *["dz_mm_ng%.2f" % g for g in a.ng], "same_side_jumps"])
            w.writerows(rows)
        print("CSV: %s" % a.csv)


if __name__ == "__main__":
    main()
