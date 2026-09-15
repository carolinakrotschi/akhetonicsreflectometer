"""Optische Netzliste des LIGENTEC-SiN-Chips: Weglaenge Coupler -> Bauteil.

Gegenstueck zu `gds_netlist_reflectors.py` (InP/HHI). Der SiN-Chip hat
andere Konventionen, deshalb ein eigenes Skript statt Sonderfaelle im
anderen:

    Wellenleiter   Layer 2/0 ("X1"), Regelbreite 1.000 um
    Bauteilports   Layer 1002, Namen a0,a1,... (eine Seite) und
                   b0,b1,... (andere Seite) -- nicht o1/o2 wie bei HHI
    Facette        Layer 101/2, 95 Stueck bei y=0, Pitch exakt 127.000 um
    Edge Coupler   AN350BB_EdgeCoupler_Lensed_C, ohne Ports: die Facette
                   liegt bei y=0, der Wellenleiter beginnt bei y=420

Kanal n (1..95) sitzt bei x = 500 + (n-1)*127. Patchcord-Faser N liegt an
Kanal 96-N, Faser 44 also an Kanal 52.

Die Geometrie (Ribbon-Formel) kommt aus dem HHI-Skript.

Aufruf
    python tools/gds_netlist_ligentec.py --channel 52
    python tools/gds_netlist_ligentec.py --fiber 44 --ng 1.9 --top 40
"""

import argparse
import csv
import heapq
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gdstk                                                      # noqa: E402
from gds_netlist_reflectors import (ribbon_length, ribbon_ends,   # noqa: E402
                                    _transform)

GDS = "raw_data/hhi_ligentec/Ligentec_HHI_combined 4.gds"
WG_LAYER = (2, 0)
PORT_LAYER = 1002
SIN_YMAX = 19313.0     # darueber beginnt der aufgeklebte InP-Chip
SNAP = 1.2             # um, Ports gelten als verbunden
N_FIBER = 1.468
CH1_X = 500.0
PITCH = 127.0
WG_START_Y = 420.0     # Wellenleiterport des Edge Couplers
COUPLER_UM = 420.0     # Facette (y=0) -> Wellenleiterport


def channel_x(ch):
    return CH1_X + (ch - 1) * PITCH


def collect_ports(lib, top):
    """Ports pro INSTANZ, absolut. Die Gruppierung kommt aus der
    Hierarchie, nicht aus raeumlicher Naehe -- die Bauteile stehen dicht
    an dicht, der naechstgelegene b-Port gehoert oft zum Nachbarn."""
    out = []

    def walk(cell, xf):
        for ref in cell.references:
            def nxf(pt, ref=ref, xf=xf):
                return xf(_transform(pt, ref.origin, ref.rotation,
                                     ref.x_reflection, ref.magnification))
            child = ref.cell
            if not child.references:
                seen, mem = set(), []
                for lab in child.labels:
                    if lab.layer != PORT_LAYER:
                        continue
                    p = nxf(lab.origin)
                    key = (lab.text, round(p[0], 3), round(p[1], 3))
                    if key in seen:        # trench_train zeichnet jeden doppelt
                        continue
                    seen.add(key)
                    mem.append((lab.text, p))
                if len(mem) >= 2:
                    out.append((child.name, mem))
            walk(child, nxf)

    walk(top, lambda p: p)
    return out


def optical_pairs(cellname, members):
    """Welche Ports eines Bauteils sind optisch verbunden.

    a* ist die eine, b* die andere Seite. `trench_train` fuehrt ZWEI
    parallele, unabhaengige Wellenleiter (a0(0.3,0)->b0(109.7,0) und
    a1(0.3,24.3)->b1(109.7,24.3)); dort darf nur nach Index gepaart
    werden, sonst springt der Pfad quer auf den Nachbarwellenleiter --
    derselbe Fehler, der beim HHI-Chip das Phantom-Loopback erzeugt hat.
    Richtkoppler und MMI koppeln dagegen wirklich jeden Eingang auf jeden
    Ausgang.
    """
    a = [i for i, m in enumerate(members) if m[0].startswith("a")]
    b = [i for i, m in enumerate(members) if m[0].startswith("b")]
    if not a or not b:
        return []
    if "trench" in cellname.lower():
        bi = {members[j][0][1:]: j for j in b}
        return [(i, bi[members[i][0][1:]]) for i in a
                if members[i][0][1:] in bi]
    return [(i, j) for i in a for j in b]


def build(gds):
    lib = gdstk.read_gds(gds)
    top = lib.top_level()[0]
    flat = top.copy("_flat").flatten()
    wgs = [p for p in flat.polygons
           if (p.layer, p.datatype) == WG_LAYER
           and p.points[:, 1].mean() < SIN_YMAX]

    nodes, edges = [], {}

    def add(x, y, tag):
        nodes.append((x, y, tag))
        return len(nodes) - 1

    def link(i, j, w):
        edges.setdefault(i, []).append((j, w))
        edges.setdefault(j, []).append((i, w))

    skipped = 0
    for poly in wgs:
        L, w = ribbon_length(poly)
        if L != L or w != w or not (0.05 < w < 20.0):
            skipped += 1
            continue
        ends = ribbon_ends(poly, w)
        if ends is None:
            skipped += 1
            continue
        link(add(ends[0][0], ends[0][1], "wg"),
             add(ends[1][0], ends[1][1], "wg"), L)

    for cellname, members in collect_ports(lib, top):
        if np.mean([m[1][1] for m in members]) > SIN_YMAX:
            continue
        idx = [add(p[0], p[1], "dev:%s:%s" % (cellname, nm)) for nm, p in members]
        for i, j in optical_pairs(cellname, members):
            link(idx[i], idx[j],
                 math.dist(nodes[idx[i]][:2], nodes[idx[j]][:2]))

    pts = np.array([[n[0], n[1]] for n in nodes])
    order = np.argsort(pts[:, 0])
    for ii in range(len(order)):
        i = order[ii]
        for jj in range(ii + 1, len(order)):
            j = order[jj]
            if pts[j, 0] - pts[i, 0] > SNAP:
                break
            if abs(pts[j, 1] - pts[i, 1]) <= SNAP and \
                    math.dist(pts[i], pts[j]) <= SNAP:
                link(i, j, 0.0)
    return nodes, edges, skipped


def dijkstra(nodes, edges, start_xy, max_um):
    start = min(range(len(nodes)),
                key=lambda i: math.dist(nodes[i][:2], start_xy))
    off = math.dist(nodes[start][:2], start_xy)
    dist = {start: 0.0}
    pq = [(0.0, start)]
    while pq:
        d, i = heapq.heappop(pq)
        if d > dist.get(i, 1e18) or d > max_um:
            continue
        for j, w in edges.get(i, ()):
            nd = d + w
            if nd < dist.get(j, 1e18):
                dist[j] = nd
                heapq.heappush(pq, (nd, j))
    return dist, off


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gds", default=GDS)
    ap.add_argument("--channel", type=int, default=None, help="Chipkanal 1..95")
    ap.add_argument("--fiber", type=int, default=None, help="Faser N -> Kanal 96-N")
    ap.add_argument("--ng", type=float, default=1.90, help="Gruppenindex SiN")
    ap.add_argument("--max-um", type=float, default=60000.0)
    ap.add_argument("--top", type=int, default=25)
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()

    ch = a.channel if a.channel else (96 - a.fiber if a.fiber else None)
    if ch is None:
        ap.error("--channel oder --fiber angeben")
    fib = 96 - ch

    nodes, edges, skipped = build(a.gds)
    print("Knoten %d   Wellenleiterpolygone uebersprungen %d" % (len(nodes), skipped))

    x0 = channel_x(ch)
    print("\n=== Kanal %d  (Faser %d)   Facette bei x = %.1f um, y = 0"
          % (ch, fib, x0))
    dist, off = dijkstra(nodes, edges, (x0, WG_START_Y), a.max_um)
    print("    Startknoten %.2f um neben dem erwarteten Wellenleiterport" % off)

    refl = sorted((d, nodes[i][2]) for i, d in dist.items()
                  if nodes[i][2].startswith("dev:"))
    uniq = []
    for d, tag in refl:
        if uniq and abs(d - uniq[-1][0]) < 1.0:
            continue
        uniq.append((d, tag))

    print("\n%10s %12s %11s  %s"
          % ("Weg [um]", "ab Facette", "dz [mm]", "Bauteil : Port"))
    rows = []
    for d, tag in uniq[:a.top]:
        s = COUPLER_UM + d
        dz = s * 1e-3 * a.ng / N_FIBER
        print("%10.1f %12.1f %11.4f  %s" % (d, s, dz, tag[4:]))
        rows.append([ch, fib, round(d, 2), round(s, 2), round(dz, 4), tag[4:]])
    if len(uniq) > a.top:
        print("   ... %d weitere bis %.0f um" % (len(uniq) - a.top, a.max_um))
    print("\n(dz = Abstand hinter Facette A im Reflektogramm, n_g = %.3f, "
          "n_Faser = %.3f)" % (a.ng, N_FIBER))

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["channel", "fiber", "waveguide_path_um",
                        "path_from_facet_um", "expected_dz_mm", "device_port"])
            w.writerows(rows)
        print("CSV: %s" % a.csv)


if __name__ == "__main__":
    main()
