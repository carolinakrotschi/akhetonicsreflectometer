"""Optische Netzliste aus einem GDS: alle Reflektorabstaende ab jedem Coupler.

Motivation: `gds_waveguide_length.py` verfolgt nur den Waveguide bis zum
ersten Bauteil. Das Licht laeuft aber DURCH die Bauteile weiter, und jede
Bauteilgrenzflaeche ist ein potenzieller Reflektor. Fuer die Deutung eines
Reflektogramms braucht man deshalb die vollstaendige Liste der
Weglaengen vom Coupler zu jeder Grenzflaeche.

Modell
   Knoten = optische Ports:
     - die zwei Enden jedes Waveguide-Polygons (Layer 12/0). Sie werden
       ueber die Stirnkanten gefunden: bei einem Band konstanter Breite w
       sind genau zwei Polygonkanten ~w lang, ihre Mittelpunkte sind die
       Mittellinien-Endpunkte.
     - die `o1`, `o2`, ... Ports jeder HHI-Primitivzelle (Labels auf
       Layer 235/0 in Zellen ohne eigene Referenzen).
   Kanten:
     - innerhalb eines Waveguide-Polygons: Ribbon-Laenge (A/P-Formel)
     - innerhalb eines Bauteils: geometrischer Abstand der Portpaare
     - zwischen zusammenfallenden Ports (< SNAP um): Laenge 0

   Dijkstra ab dem Port eines Couplers liefert dann die Weglaenge zu jeder
   Grenzflaeche im Chip.

Aufruf
   python tools/gds_netlist_reflectors.py <gds> [--ng 3.4803] [--max-um 12000]
                                          [--channel b2] [--csv out.csv]
"""

import argparse
import csv
import heapq
import math
import re
from collections import defaultdict

import gdstk
import numpy as np

WG_LAYER = (12, 0)
PORT_LAYER = 235
SNAP = 1.6          # um, Ports gelten als verbunden
N_FIBER = 1.468


# ---------------------------------------------------------------- Geometrie
def _transform(pt, origin, rotation, x_reflection, magnification):
    x, y = pt
    if x_reflection:
        y = -y
    m = magnification if magnification else 1.0
    x, y = x * m, y * m
    if rotation:
        c, s = math.cos(rotation), math.sin(rotation)
        x, y = x * c - y * s, x * s + y * c
    return x + origin[0], y + origin[1]


def collect_primitive_ports(lib, top):
    """Ports PRO INSTANZ: [(zellname, [(portname,(x,y)), ...]), ...], absolut.

    Wichtig: die Gruppierung muss aus der Hierarchie kommen, nicht aus
    raeumlicher Naehe -- die Bauteile stehen in Ketten nur 30 um
    auseinander, und der naechstgelegene `o2` zu einem `o1` gehoert dann
    zum NACHBARbauteil (0.6 um entfernt) statt zur eigenen Instanz
    (29.4 um entfernt).
    """
    out = []

    def walk(cell, xf):
        for ref in cell.references:
            def nxf(pt, ref=ref, xf=xf):
                return xf(_transform(pt, ref.origin, ref.rotation,
                                     ref.x_reflection, ref.magnification))
            child = ref.cell
            if not child.references:            # Primitivzelle = eine Instanz
                members = [(lab.text, nxf(lab.origin)) for lab in child.labels
                           if lab.layer == PORT_LAYER
                           and re.fullmatch(r"o\d+", lab.text)]
                if len(members) >= 2:
                    out.append((child.name, members))
            walk(child, nxf)

    walk(top, lambda pt: pt)
    return out


def _perimeter(poly):
    q = poly.points
    return sum(math.dist(q[i], q[(i + 1) % len(q)]) for i in range(len(q)))


def ribbon_length(poly):
    area = abs(poly.area())
    per = _perimeter(poly)
    half = per / 2.0
    disc = half * half - 4.0 * area
    if disc < 0:
        return float("nan"), float("nan")
    root = math.sqrt(disc)
    return (half + root) / 2.0, (half - root) / 2.0


def ribbon_ends(poly, width):
    """Mittelpunkte der beiden Stirnkanten (= Mittellinien-Endpunkte)."""
    q = poly.points
    n = len(q)
    cand = []
    for i in range(n):
        a, b = q[i], q[(i + 1) % n]
        d = math.dist(a, b)
        if 0.45 * width < d < 2.2 * width:
            cand.append((d, ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)))
    if len(cand) < 2:
        return None
    # die zwei weitest voneinander entfernten Kandidaten
    best, pair = -1.0, None
    for i in range(len(cand)):
        for j in range(i + 1, len(cand)):
            d = math.dist(cand[i][1], cand[j][1])
            if d > best:
                best, pair = d, (cand[i][1], cand[j][1])
    return pair


# ---------------------------------------------------------------- Netzliste
def build(gds):
    lib = gdstk.read_gds(gds)
    top = lib.top_level()[0]
    flat = top.copy("_flat").flatten()
    wgs = [p for p in flat.polygons if (p.layer, p.datatype) == WG_LAYER]

    nodes = []          # (x, y, tag)
    edges = defaultdict(list)   # node index -> [(other, weight)]

    def add_node(x, y, tag):
        nodes.append((x, y, tag))
        return len(nodes) - 1

    skipped = 0
    for poly in wgs:
        L, w = ribbon_length(poly)
        if L != L or w != w:
            skipped += 1
            continue
        ends = ribbon_ends(poly, w)
        if ends is None:
            skipped += 1
            continue
        a = add_node(ends[0][0], ends[0][1], "wg")
        b = add_node(ends[1][0], ends[1][1], "wg")
        edges[a].append((b, L))
        edges[b].append((a, L))

    comp_nodes = collect_primitive_ports(lib, top)

    for cellname, members in comp_nodes:
        idx = [add_node(p[0], p[1], "dev:%s:%s" % (cellname, nm))
               for nm, p in members]
        for i in range(len(idx)):
            for j in range(i + 1, len(idx)):
                d = math.dist(nodes[idx[i]][:2], nodes[idx[j]][:2])
                edges[idx[i]].append((idx[j], d))
                edges[idx[j]].append((idx[i], d))

    # koinzidente Knoten verschmelzen (Laenge 0)
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
                edges[i].append((j, 0.0))
                edges[j].append((i, 0.0))
    return lib, top, nodes, edges, skipped


def coupler_ports(lib):
    eca = {c.name: c for c in lib.cells}["edge_couplers_array_$9eca"]
    return {lab.text: (lab.origin[0], lab.origin[1] + 112.0)
            for lab in eca.labels if lab.layer == PORT_LAYER
            and re.fullmatch(r"b\d+", lab.text)}


def reflectors(nodes, edges, start_xy, max_um):
    start = min(range(len(nodes)),
                key=lambda i: math.dist(nodes[i][:2], start_xy))
    dist = {start: 0.0}
    pq = [(0.0, start)]
    while pq:
        d, i = heapq.heappop(pq)
        if d > dist.get(i, 1e18) or d > max_um:
            continue
        for j, w in edges[i]:
            nd = d + w
            if nd < dist.get(j, 1e18):
                dist[j] = nd
                heapq.heappush(pq, (nd, j))
    out = [(d, nodes[i][2]) for i, d in dist.items() if nodes[i][2].startswith("dev:")]
    return sorted(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gds")
    ap.add_argument("--ng", type=float, default=3.4803)
    ap.add_argument("--max-um", type=float, default=12000.0)
    ap.add_argument("--channel", default=None, help="nur dieser Kanal, z.B. b2")
    ap.add_argument("--csv", default=None)
    ap.add_argument("--top", type=int, default=12, help="wie viele Reflektoren je Kanal")
    a = ap.parse_args()

    lib, top, nodes, edges, skipped = build(a.gds)
    print("Knoten: %d   Waveguide-Polygone uebersprungen: %d" % (len(nodes), skipped))
    cps = coupler_ports(lib)

    rows = []
    chans = [a.channel] if a.channel else sorted(cps, key=lambda s: int(s[1:]))
    for ch in chans:
        refl = reflectors(nodes, edges, cps[ch], a.max_um)
        # gleiche Distanz zusammenfassen
        uniq = []
        for d, tag in refl:
            if uniq and abs(d - uniq[-1][0]) < 1.0:
                continue
            uniq.append((d, tag))
        print("\n=== %s : %d Grenzflaechen bis %.0f um" % (ch, len(uniq), a.max_um))
        print("%12s %14s  %s" % ("Weg [um]", "dz [mm]", "Bauteil / Port"))
        for d, tag in uniq[:a.top]:
            dz = (1200.0 + d) * 1e-3 * a.ng / N_FIBER
            print("%12.1f %14.4f  %s" % (d, dz, tag[4:]))
            rows.append([ch, round(d, 2), round(dz, 4), tag[4:]])
        if len(uniq) > a.top:
            print("   ... %d weitere" % (len(uniq) - a.top))

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["channel", "waveguide_path_um", "expected_delta_z_mm", "device_port"])
            w.writerows(rows)
        print("\nCSV: %s" % a.csv)


if __name__ == "__main__":
    main()
