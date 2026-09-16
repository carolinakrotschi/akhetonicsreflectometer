"""Was liegt bei welcher Weglaenge? Strukturprofil eines LIGENTEC-Kanals.

`gds_netlist_ligentec.py` liefert nur die Bauteile, die Ports auf Layer 1002
tragen -- die AN350-Blackboxen und die meisten Layoutzellen tun das nicht,
die Tabelle bleibt fuer die meisten Kanaele leer. Dieses Skript geht
stattdessen den Pfad selbst ab und fragt an jeder Stelle, in welcher
Zellinstanz er gerade steckt.

Verfahren
    1. Graph wie im Netzlistenskript (Wellenleiterpolygone auf 2/0,
       Ribbon-Bogenlaenge als Kantengewicht), Dijkstra ab dem
       Wellenleiterport des Edge Couplers bei (x_Kanal, y=420).
    2. Pfad zum entferntesten erreichbaren Knoten rekonstruieren.
    3. Jedes Pfadsegment gegen die Bounding Box JEDER Zellinstanz
       schneiden (Liang-Barsky). Das ist noetig, weil Knoten nur an
       Polygonenden sitzen: zwischen y=2798 und y=8255 liegt auf Kanal 52
       ein einziges 5.5 mm langes Polygon, und alles, was es unterwegs
       durchquert, wuerde eine reine Knotenbetrachtung verschlucken.
    4. Eintritts- und Austrittslaenge ueber den Segmentparameter auf die
       Bogenlaenge abbilden.

Eine Box-Ueberschneidung ist nicht dasselbe wie "durchlaeuft das Bauteil":
bei einer L-foermigen oder sehr grossen Zelle kann der Pfad die Box
streifen, ohne die Struktur zu beruehren. Deshalb wird die Boxgroesse
mitgedruckt und tief verschachtelte (= kleine, spezifische) Instanzen
zuerst genannt.

Aufruf
    python tools/gds_trace_ligentec.py --channel 48
    python tools/gds_trace_ligentec.py --channel 52 --bridge --csv out.csv
"""

import argparse
import csv
import heapq
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gdstk                                                      # noqa: E402
import gds_netlist_ligentec as G                                  # noqa: E402
from gds_netlist_reflectors import _transform                     # noqa: E402

# Zellnamen, die als optische Struktur berichtet werden. Alles andere
# (Bends, Routing, Metall, Pads) ist Verdrahtung und nur Rauschen in der
# Tabelle.
INTERESTING = re.compile(
    r"cross|MMI|MZM|MZI|coupler|spiral|trench|heater|threshold|"
    r"splitting|combiner|decoder|LUT|XGM|SOA|loop|taper|termination|RR",
    re.I)
MAX_DIAG = 6000.0      # um, groessere Instanzen gelten als Container
PAD = 0.5              # um, Toleranz um die Box


def cell_instances(top, ymax):
    """(name, tiefe, (x0,y0,x1,y1)) jeder Referenz unterhalb `top`."""
    out = []
    bb_cache = {}

    def walk(cell, xf, depth):
        for ref in cell.references:
            def nxf(pt, ref=ref, xf=xf):
                return xf(_transform(pt, ref.origin, ref.rotation,
                                     ref.x_reflection, ref.magnification))
            child = ref.cell
            if child.name not in bb_cache:
                bb_cache[child.name] = child.bounding_box()
            bb = bb_cache[child.name]
            if bb is not None:
                (ax, ay), (bx, by) = bb
                cs = [nxf((ax, ay)), nxf((bx, ay)), nxf((bx, by)), nxf((ax, by))]
                xs = [c[0] for c in cs]
                ys = [c[1] for c in cs]
                if min(ys) < ymax:
                    out.append((child.name, depth,
                                (min(xs), min(ys), max(xs), max(ys))))
            walk(child, nxf, depth + 1)

    walk(top, lambda p: p, 0)
    return out


def dijkstra_pred(nodes, edges, start_xy, max_um):
    start = min(range(len(nodes)), key=lambda i: math.dist(nodes[i][:2], start_xy))
    dist = {start: 0.0}
    pred = {start: None}
    pq = [(0.0, start)]
    while pq:
        d, i = heapq.heappop(pq)
        if d > dist.get(i, 1e18) or d > max_um:
            continue
        for j, w in edges.get(i, ()):
            nd = d + w
            if nd < dist.get(j, 1e18) - 1e-9:
                dist[j] = nd
                pred[j] = i
                heapq.heappush(pq, (nd, j))
    return start, dist, pred


def clip(p0, p1, box):
    """Liang-Barsky: Parameterintervall [t0,t1], in dem p0->p1 in `box` liegt."""
    x0, y0 = p0
    dx, dy = p1[0] - x0, p1[1] - y0
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x0 - (box[0] - PAD)), (dx, (box[2] + PAD) - x0),
                 (-dy, y0 - (box[1] - PAD)), (dy, (box[3] + PAD) - y0)):
        if abs(p) < 1e-12:
            if q < 0:
                return None
            continue
        t = q / p
        if p < 0:
            if t > t1:
                return None
            t0 = max(t0, t)
        else:
            if t < t0:
                return None
            t1 = min(t1, t)
    return (t0, t1) if t1 > t0 else None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gds", default=G.GDS)
    ap.add_argument("--channel", type=int, default=None, help="Chipkanal 1..95")
    ap.add_argument("--fiber", type=int, default=None, help="Faser N = Kanal 96-N")
    ap.add_argument("--bridge", action="store_true",
                    help="Blackboxen durchlaessig machen (Laenge dahinter genaehert)")
    ap.add_argument("--ng", type=float, default=1.90)
    ap.add_argument("--max-um", type=float, default=200000.0)
    ap.add_argument("--min-len", type=float, default=0.0,
                    help="Strukturen kuerzer als dies nicht auflisten")
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    ch = a.channel if a.channel else (96 - a.fiber if a.fiber else None)
    if ch is None:
        ap.error("--channel oder --fiber angeben")

    lib = gdstk.read_gds(a.gds)
    top = lib.top_level()[0]
    inst = [(n, d, b) for (n, d, b) in cell_instances(top, G.SIN_YMAX)
            if INTERESTING.search(n)
            and math.hypot(b[2] - b[0], b[3] - b[1]) < MAX_DIAG]
    boxes = np.array([b for (_, _, b) in inst])

    nodes, edges, skipped = G.build(a.gds, bridge_boxes=a.bridge)
    x0 = G.channel_x(ch)
    start, dist, pred = dijkstra_pred(nodes, edges, (x0, G.WG_START_Y), a.max_um)
    end = max(dist, key=lambda i: dist[i])

    path = []
    i = end
    while i is not None:
        path.append(i)
        i = pred[i]
    path.reverse()

    hits = {}      # (name, box) -> [s_ein, s_aus]
    for k in range(1, len(path)):
        i, j = path[k - 1], path[k]
        p0, p1 = nodes[i][:2], nodes[j][:2]
        s0, s1 = dist[i], dist[j]
        chord = math.dist(p0, p1)
        if chord < 1e-9:
            continue
        lo = np.minimum(p0, p1) - PAD
        hi = np.maximum(p0, p1) + PAD
        sel = np.where((boxes[:, 0] <= hi[0]) & (boxes[:, 2] >= lo[0]) &
                       (boxes[:, 1] <= hi[1]) & (boxes[:, 3] >= lo[1]))[0]
        for ii in sel:
            name, depth, box = inst[ii]
            tt = clip(p0, p1, box)
            if tt is None:
                continue
            key = (name, box, depth)
            sa, sb = s0 + tt[0] * (s1 - s0), s0 + tt[1] * (s1 - s0)
            if key in hits:
                hits[key][0] = min(hits[key][0], sa)
                hits[key][1] = max(hits[key][1], sb)
            else:
                hits[key] = [sa, sb]

    rows = sorted(((v[0], v[1], k[0], k[2], k[1]) for k, v in hits.items()))

    print("Kanal %d  (Faser %d)   Facette x = %.1f um, y = 0"
          % (ch, 96 - ch, x0))
    print("Blackboxen %s   Knoten %d   Instanzen gepruft %d"
          % ("durchlaessig" if a.bridge else "undurchlaessig", len(nodes), len(inst)))
    print("Pfadende bei (%.1f, %.1f), %.1f um Wellenleiter ab Coupler-Port, "
          "%.1f um ab Facette" % (nodes[end][0], nodes[end][1], dist[end],
                                  dist[end] + G.COUPLER_UM))
    print("\n%9s %9s %8s  %-46s %s"
          % ("ein [um]", "aus [um]", "dz [mm]", "Struktur", "Box [um]"))
    out = []
    for sa, sb, name, depth, box in rows:
        if sb - sa < a.min_len:
            continue
        s_facet = G.COUPLER_UM + sa
        dz = s_facet * 1e-3 * a.ng / G.N_FIBER
        print("%9.1f %9.1f %8.3f  %-46s %.0fx%.0f"
              % (sa, sb, dz, ("  " * min(depth, 6)) + name[:44],
                 box[2] - box[0], box[3] - box[1]))
        out.append([ch, round(sa, 2), round(sb, 2), round(s_facet, 2),
                    round(dz, 4), name, depth,
                    round(box[2] - box[0], 2), round(box[3] - box[1], 2)])

    print("\n(dz = Position hinter Facette A im Reflektogramm, n_g = %.3f, "
          "n_Faser = %.3f)" % (a.ng, G.N_FIBER))

    # Sackgassen: erreichbare Knoten, von denen aus es nicht weiter geht.
    # Bei Verzweigungen (Splitter) endet der Kanal nicht an einer Stelle.
    leaves = []
    for i, d in dist.items():
        if all(dist.get(j, -1) <= d + 1e-9 for j, _ in edges.get(i, ())):
            leaves.append((d, i))
    leaves.sort(reverse=True)
    shown, seen = 0, []
    print("\nEndpunkte des erreichbaren Netzes (Sackgassen, weiteste zuerst):")
    for d, i in leaves:
        p = nodes[i][:2]
        if any(math.dist(p, q) < 30.0 for q in seen):
            continue
        seen.append(p)
        here = sorted((dep, nm) for (nm, dep, b) in inst
                      if b[0] - PAD <= p[0] <= b[2] + PAD
                      and b[1] - PAD <= p[1] <= b[3] + PAD)
        print("  %9.1f um  (%9.1f,%9.1f)   %s"
              % (d, p[0], p[1], here[-1][1] if here else "-"))
        shown += 1
        if shown >= 8:
            break

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["channel", "s_in_um", "s_out_um", "s_from_facet_um",
                        "dz_mm", "cell", "depth", "box_w_um", "box_h_um"])
            w.writerows(out)
        print("CSV: %s" % a.csv)


if __name__ == "__main__":
    main()
