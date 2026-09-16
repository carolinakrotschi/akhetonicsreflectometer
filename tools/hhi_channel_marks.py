#!/usr/bin/env python3
"""Bauteil-Sollpositionen eines HHI-Kanals als Marken-CSV fuer die Plots.

Die optische Netzliste (`tools/gds_netlist_reflectors.py`) gibt pro Kanal
jede Bauteilgrenzflaeche als Weglaenge `waveguide_path_um` ab dem
INNEREN Ende des eigenen SSC aus. Die Reflektogramm-Achse zaehlt dagegen
ab Facette A, also ab dem AEUSSEREN Ende desselben SSC:

    s_ab_Facette = SSC + waveguide_path       mit SSC = 1199.4 um

Das ist dieselbe Umrechnung, die `identify_chip_scan.py` benutzt
(`dz = (SSC + d) * n_g / n_Faser`), hier nur als Tabelle statt als
Einzelvergleich.

ZWEI DINGE, DIE DIE NETZLISTE ROH NICHT LIEFERT

1. **Die eigene Facette taucht als Scheinreflektor auf.** Dijkstra startet
   am inneren SSC-Ende und erreicht das aeussere (= Facette A) nach
   1199.4 um -- in der Tabelle steht sie damit bei s = 2398.8 um. Physisch
   liegt sie bei s = 0: das ist der Weg nach DRAUSSEN, kein Reflektor im
   Chip. Diese eine Zeile wird verworfen.

2. **Bis zu 95 Grenzflaechen pro Kanal sind als Einzelstriche unlesbar.**
   Die Bauteile stehen in Ketten 25-30 um auseinander, bei 7.0 um
   On-Chip-Aufloesung ist das nicht trennbar. Grenzflaechen mit weniger
   als `--gap` um Abstand werden deshalb zu einer Gruppe zusammengefasst;
   die Marke im Plot ist dann der Bereich von der ersten bis zur letzten
   Grenzflaeche der Gruppe.

Aufruf
    python tools/hhi_channel_marks.py --fiber 4 \
        --out results/2026-09-16/hhi1_b28_bauteile.csv
"""

import argparse
import csv
import os

SSC = 1199.4        # um, Port-zu-Port eines SSC laut Netzliste
N_FIBER = 1.468

# GDS-Zellname -> kurze Plotbeschriftung
SHORT = [
    ("HHI_SSCLATE1700", "SSC"),
    ("HHI_WGMETxE1700twin", "WGMETx twin"),
    ("HHI_WGMETxE1700single", "WGMETx"),
    ("HHI_MMI1x2E1700", "MMI1x2"),
    ("HHI_MMI2x2E1700", "MMI2x2"),
    ("HHI_ISOsectionSingle", "ISO"),
    ("HHI_PMTOE1700", "PMTO"),
    ("HHI_BJsingle", "BJ"),
    ("HHI_SOAsection", "SOA"),
    ("HHI_SOA_", "SOA"),
    ("crossing_e1700", "crossing"),
]


def short(cell, port=""):
    """Zellname (+ Port) -> kurze Plotbeschriftung.

    Der SSC wird nach Port unterschieden: `o1` ist sein inneres Ende
    (SSC -> Waveguide), `o2` das aeussere, also eine Chipfacette. Beim
    Loopback ist genau dieses `o2` des Partner-SSC die Facette B, der
    einzige zweite HARTE Reflektor im Chip -- als "SSC" beschriftet waere
    das im Plot nicht von der Zwischengrenzflaeche zu unterscheiden."""
    for key, name in SHORT:
        if cell.startswith(key):
            if name == "SSC":
                return "chip facet" if port == "o2" else "SSC end"
            return name
    return cell.split("_$")[0]


def load(path, channel):
    """[(s_ab_Facette_um, kurzname), ...] eines Kanals, aufsteigend."""
    out, dropped = [], 0
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["channel"] != channel:
                continue
            d = float(r["waveguide_path_um"])
            cell, _, port = r["device_port"].partition(":")
            # Punkt 1 im Docstring: das aeussere Ende des EIGENEN SSC ist
            # Facette A selbst, erreicht auf dem Weg nach DRAUSSEN.
            if abs(d - SSC) < 1.0 and short(cell, port) == "chip facet":
                dropped += 1
                continue
            out.append((SSC + d, short(cell, port)))
    out.sort()
    return out, dropped


def group(items, gap):
    """Grenzflaechen mit < gap um Abstand zu einer Marke zusammenfassen."""
    groups = []
    for s, name in items:
        if groups and s - groups[-1][1] < gap:
            groups[-1][1] = s
            groups[-1][2].append(name)
        else:
            groups.append([s, s, [name]])
    return groups


# Wie hart reflektiert die Grenzflaeche? Nur die Reihenfolge zaehlt: bei
# einer Gruppe aus mehreren Bauteilen wird die vorderste im Plot genannt.
# Facette und Butt Joint sind echte Indexspruenge (Luft/InP bzw.
# aktiv/passiv), MMI- und Crossing-Grenzflaechen leben von
# Modenfehlanpassung, WGMETx und ISO sind Metallbruecke bzw.
# Implantation UEBER einem durchlaufenden Waveguide -- optisch fast nichts.
RANK = ["chip facet", "BJ", "SOA", "MMI2x2", "MMI1x2", "crossing", "SSC end",
        "PMTO", "ISO", "WGMETx twin", "WGMETx"]


def label(names, maxparts=2):
    """['BJ','WGMETx','WGMETx'] -> 'BJ + 2x WGMETx'.

    Bei mehr als `maxparts` verschiedenen Bauteilen werden die optisch
    haertesten genannt und der Rest als '+n' angehaengt -- sonst reicht die
    Beschriftung einer einzigen Marke im Plot ueber die halbe Achse."""
    seen = []
    for n in names:
        if n not in seen:
            seen.append(n)
    seen.sort(key=lambda n: RANK.index(n) if n in RANK else len(RANK))
    parts = []
    for n in seen[:maxparts]:
        k = names.count(n)
        parts.append("%dx %s" % (k, n) if k > 1 else n)
    out = " + ".join(parts)
    if len(seen) > maxparts:
        out += " +%d" % (len(seen) - maxparts)
    return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--netlist",
                    default="results/2026-09-16/hhi1_netlist_reflectors.csv")
    ap.add_argument("--fiber", type=int, default=None,
                    help="Fasernummer, umgerechnet mit b = 32 - Faser")
    ap.add_argument("--channel", default=None, help="Kanal direkt, z.B. b28")
    ap.add_argument("--gap", type=float, default=60.0,
                    help="um On-Chip-Abstand, ab dem eine neue Marke beginnt "
                         "(default 60, = 8.5 Aufloesungszellen)")
    ap.add_argument("--loop-span", action="store_true",
                    help="Loopback-Kanaele: die Schleife zwischen den beiden "
                         "SSC-Enden als eigene, beschriftete Marke eintragen. "
                         "Sie ist KEIN Reflektor -- der Waveguide laeuft "
                         "stufenlos in den SSC -- aber ohne sie steht im Plot "
                         "nur 'SSC end, SSC end' und man sieht nicht, dass "
                         "der Abstand dazwischen die Schleife IST")
    ap.add_argument("--ng", type=float, default=3.47794,
                    help="Gruppenindex des Chips fuer die dz-Spalte")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if not (a.fiber is None) ^ (a.channel is None):
        raise SystemExit("genau eines von --fiber / --channel angeben")
    ch = a.channel if a.channel else "b%d" % (32 - a.fiber)

    items, dropped = load(a.netlist, ch)
    if not items:
        raise SystemExit("Kanal %s steht nicht in %s" % (ch, a.netlist))
    groups = group(items, a.gap)

    if a.loop_span:
        ends = [g for g in groups if g[2] == ["SSC end"]]
        if len(ends) < 2:
            raise SystemExit("--loop-span: %s hat keine zwei SSC-Enden, ist "
                             "also kein Loopback" % ch)
        s0, s1 = ends[0][0], ends[1][0]
        groups.append([s0, s1, ["loop %.1f um (kein Reflex)" % (s1 - s0)]])
        groups.sort()

    print("Kanal %s: %d Grenzflaechen -> %d Marken "
          "(%d als eigene Facette verworfen)"
          % (ch, len(items), len(groups), dropped))
    print("%10s %10s  %8s %8s  %s"
          % ("s_von[um]", "s_bis[um]", "dz_von", "dz_bis", "Bauteile"))
    rows = []
    for s0, s1, names in groups:
        lab = names[0] if names[0].startswith("loop ") else label(names)
        dz0, dz1 = (s0 * 1e-3 * a.ng / N_FIBER, s1 * 1e-3 * a.ng / N_FIBER)
        print("%10.1f %10.1f  %8.4f %8.4f  %s" % (s0, s1, dz0, dz1, lab))
        rows.append([ch, lab, round(s0, 1), round(s1, 1), len(names),
                     round(dz0, 4), round(dz1, 4)])

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["channel", "cell", "s_from_facet_um", "s_to_um",
                    "n_interfaces", "dz_from_mm", "dz_to_mm"])
        w.writerows(rows)
    print("\nCSV: %s" % a.out)


if __name__ == "__main__":
    main()
