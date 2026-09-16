#!/usr/bin/env python3
"""Aus den DATEN heraus fragen, was ein Reflektogramm zeigt -- und was davon
Zufall ist.

Gegenstueck zu `hhi_channel_marks.py`: dort wird aus dem GDS gerechnet, wo
etwas sein MUESSTE, und ein Strich gezogen. Hier andersherum -- erst finden,
was in den Daten steht, dann fuer jeden Fund fragen, welche
Designgrenzflaeche in Frage kaeme und wie wahrscheinlich diese Zuordnung
reiner Zufall ist.

ZWEI SORTEN FUND, WEIL ES ZWEI SORTEN REFLEXION GIBT

1. **Stachel** -- eine einzelne Grenzflaeche (Facette, Butt Joint,
   MMI-Kante). Erkannt als Peak ueber dem OERTLICHEN Untergrund
   (gleitender Median, `--bg-mm`), nicht ueber dem Rauschboden: der
   Multipath-Wald steht 30 dB ueber dem Boden und gaebe sonst hunderte
   Scheintreffer.

2. **Buckel** -- eine erhoehte STRECKE. Ein Waveguide reflektiert nicht an
   einer Stelle, er streut ueber seine ganze Laenge zurueck. So ist die
   Schleife der Loopback-Kanaele sichtbar (b0: Soll 1.318 mm, gemessen
   1.47 mm breit), und so sehen die breiten Erhebungen in den
   Schaltungskanaelen aus. Erkannt als zusammenhaengender Bereich, in dem
   die Huellkurve (gleitendes 75-%-Quantil) ueber ihrem eigenen breiten
   Untergrund liegt.

WARUM ZU JEDEM FUND EINE ZUFALLSZAHL GEHOERT

Ein Schaltungskanal hat bis zu 81 Grenzflaechen auf 35 mm. Bei 100 um
Toleranz ist dann ein Drittel der Achse "in der Naehe von irgendetwas",
und eine Zuordnung ohne diese Zahl daneben sagt nichts. Ausgegeben wird

    p = Anteil der Achse, der genauso nah oder naeher an einer
        Grenzflaeche liegt wie dieser Fund

p = 0.05 heisst: ein an zufaelliger Stelle gewuerfelter Peak haette in 5 %
der Faelle einen mindestens so guten Treffer. Klein ist gut. Fuer Buckel
wird zusaetzlich der ueberdeckte Achsenanteil ausgegeben -- sind 50 % der
Achse "erhoeht", ist ein einzelner Buckel keine Aussage.

WAS VORHER RAUSFLIEGT
   Die Instrument-Satelliten bei +-0.70 und +-1.53 mm um jeden starken
   Peak (`drop_satellites` aus `identify_chip_scan.py`) und alles, was als
   Doppelbounce 2 x dz eines staerkeren Peaks erklaerbar ist -- ein
   Fabry-Perot zwischen zwei Reflektoren erzeugt den zwangslaeufig.

Aufruf
    python tools/hhi_identify_peaks.py \
        --csv results/2026-09-14/2026-09-14-09-35_fiber32_reflectogram.csv \
        --fiber 32 --facet 1668.6277 --zmax 16 \
        --out results/2026-09-16/hhi1_f32_peaks.csv
"""

import argparse
import csv
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
from process_reflectogram_aux import noise_floor              # noqa: E402
from identify_chip_scan import SSC, N_FIBER, drop_satellites  # noqa: E402
from hhi_channel_marks import short                            # noqa: E402


def load(path):
    d = np.genfromtxt(path, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def interfaces(netlist, channel, ng):
    """[(dz_mm, Bauteilname), ...] -- jede EINZELNE Grenzflaeche, ungruppiert.

    Fuer die Zuordnung darf nicht gruppiert werden: die Gruppen aus
    `hhi_channel_marks.py` sind eine Lesehilfe fuer den Plot, geprueft wird
    hier gegen die echten Positionen."""
    out = []
    with open(netlist) as fh:
        for r in csv.DictReader(fh):
            if r["channel"] != channel:
                continue
            d = float(r["waveguide_path_um"])
            cell, _, port = r["device_port"].partition(":")
            if abs(d - SSC) < 1.0 and "SSCLATE" in cell:
                continue          # eigene Facette, s. hhi_channel_marks.py
            out.append(((SSC + d) * 1e-3 * ng / N_FIBER, short(cell, port)))
    return sorted(out)


def spikes(x, db, lo, hi, bg_mm, min_prom):
    """Peaks ueber dem OERTLICHEN Untergrund (gleitender Median)."""
    from scipy.signal import find_peaks
    from scipy.ndimage import median_filter
    m = (x >= lo) & (x <= hi)
    xx, dd = x[m], db[m]
    k = max(int(bg_mm / float(np.median(np.diff(xx)))) | 1, 3)
    bg = median_filter(dd, size=k, mode="nearest")
    prom = dd - bg
    idx, _ = find_peaks(dd, distance=3)
    return [(float(xx[i]), float(dd[i]), float(prom[i]))
            for i in idx if prom[i] >= min_prom]


def humps(x, db, lo, hi, env_mm, base_mm, thr_db, min_w):
    """Erhoehte STRECKEN: Huellkurve ueber ihrem eigenen breiten Untergrund.

    Huellkurve = gleitendes 75-%-Quantil ueber `env_mm` (folgt der Oberkante
    des Specklemusters, ohne den einzelnen Ausreisser mitzunehmen);
    Untergrund = gleitender Median ueber `base_mm`.
    -> (z_von, z_bis, Breite, Maximum_dB, Flaeche_ueber_Untergrund)"""
    from scipy.ndimage import median_filter, percentile_filter
    m = (x >= lo) & (x <= hi)
    xx, dd = x[m], db[m]
    step = float(np.median(np.diff(xx)))
    env = percentile_filter(dd, 75, size=max(int(env_mm / step) | 1, 3),
                            mode="nearest")
    base = median_filter(dd, size=max(int(base_mm / step) | 1, 3),
                         mode="nearest")
    over = env - base >= thr_db
    out, i = [], 0
    while i < len(over):
        if not over[i]:
            i += 1
            continue
        j = i
        while j < len(over) and over[j]:
            j += 1
        w = xx[j - 1] - xx[i]
        if w >= min_w:
            sl = slice(i, j)
            out.append((float(xx[i]), float(xx[j - 1]), float(w),
                        float(dd[sl].max()),
                        float(np.sum(env[sl] - base[sl]) * step)))
        i = j
    return out, float(np.mean(over))


def sidebands(pk, tol_mm, min_gap, min_lead=10.0):
    """Peaks, die als SPIEGELPAAR um einen staerkeren Peak auftreten.

    `drop_satellites` kennt nur die zwei fest bekannten Instrumentabstaende
    (+-0.70 / +-1.53 mm). Um Facette B stehen aber weitere Paare, z.B. bei
    +-1.15 mm. Ein Reflektor im Chip kann kein Spiegelpaar erzeugen -- dass
    links UND rechts im selben Abstand etwas gleich Starkes steht, ist die
    Zwei Bedingungen, sonst erklaert sich der Wald selbst: der Traeger muss
    mindestens `min_lead` dB STAERKER sein als die Seitenlinie, und das
    Spiegelbild muss auf 6 dB genauso hoch sein wie sie.
    -> {z: (z_traeger, abstand)}"""
    out = {}
    for zc, dc, _ in sorted(pk, key=lambda t: -t[1]):
        for zi, di, _ in pk:
            d = zi - zc
            if abs(d) < min_gap or dc - di < min_lead or zi in out:
                continue
            if any(abs((zj - zc) + d) < tol_mm and abs(dj - di) <= 6.0
                   for zj, dj, _ in pk):
                out[zi] = (zc, abs(d))
    return out


def p_chance(delta_mm, ifaces, lo, hi, step=0.001):
    """Anteil der Achse [lo,hi], der <= delta von einer Grenzflaeche weg ist."""
    if not ifaces:
        return 1.0
    grid = np.arange(lo, hi, step)
    pos = np.array([d for d, _ in ifaces])
    near = np.min(np.abs(grid[:, None] - pos[None, :]), axis=1)
    return float(np.mean(near <= delta_mm))


def hump_verdict(z0, z1, w, segs, ifaces):
    """(Urteil, Kandidatentext, Abweichung_um) fuer einen Buckel."""
    mid = 0.5 * (z0 + z1)
    # Nur Grenzflaechen ECHT innerhalb zaehlen. Ein Buckel, der eine
    # Strecke ausfuellt, beruehrt zwangslaeufig deren beide Enden -- die
    # duerfen ihn nicht als "ueberspannt 2 Grenzflaechen" disqualifizieren.
    eps = max(0.1 * (z1 - z0), 0.02)
    n_in = sum(1 for d, _ in ifaces if z0 + eps < d < z1 - eps)
    if n_in > 1:
        return ("ueberspannt %d Grenzfl." % n_in,
                "kein einzelner Waveguide", None)
    # Die Designstrecke, in der die Buckelmitte LIEGT -- nicht die mit der
    # aehnlichsten Laenge. Ein Buckel 2 mm neben der Strecke ist keine
    # Zuordnung, auch wenn die Laenge zufaellig passt.
    seg = next((t for t in segs if t[0] <= mid <= t[1]), None)
    if seg is None:
        return ("unerklaert", "ausserhalb aller Designstrecken", None)
    s0, s1, sl = seg
    txt = "Strecke %.2f-%.2f mm (%.2f mm lang)" % (s0, s1, sl)
    # Fuellt der Buckel die ganze Strecke (= Rueckstreuung DIESES
    # Waveguides) oder nur ein Stueck davon?
    good = abs(w - sl) < 0.25 * sl
    return ("Strecke gefuellt" if good else "Teil der Strecke",
            txt, (w - sl) * 1e3)


def p_hump(w, dlt_um, segs, lo, hi, n=8000, seed=0):
    """Anteil zufaellig platzierter Buckel gleicher Breite, deren Treffer
    MINDESTENS SO GUT waere wie der beobachtete.

    Gleiche Definition wie `p_chance` bei den Stacheln: nicht "passt
    irgendwie", sondern "passt mindestens so genau". Die Kontrolle ist
    noetig, weil bei 81 Grenzflaechen auf 35 mm fast jede Breite
    irgendwohin passt."""
    rng = np.random.default_rng(seed)
    tol = abs(dlt_um) * 1e-3
    z0 = rng.uniform(lo, max(hi - w, lo + 1e-6), size=n)
    hit = 0
    for a in z0:
        mid = a + 0.5 * w
        seg = next((t for t in segs if t[0] <= mid <= t[1]), None)
        if seg is not None and abs(w - seg[2]) <= tol:
            hit += 1
    return hit / float(n)


def segments(ifaces):
    """Design-Waveguidestrecken: (von, bis, Laenge) zwischen aufeinander
    folgenden Grenzflaechen. Ein Buckel sollte auf eine davon passen."""
    return [(ifaces[i][0], ifaces[i + 1][0], ifaces[i + 1][0] - ifaces[i][0])
            for i in range(len(ifaces) - 1)]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="Reflektogramm-CSV")
    ap.add_argument("--facet", type=float, required=True,
                    help="Facette A in mm (absolut)")
    ap.add_argument("--fiber", type=int, default=None)
    ap.add_argument("--channel", default=None)
    ap.add_argument("--netlist",
                    default="results/2026-09-16/hhi1_netlist_reflectors.csv")
    ap.add_argument("--ng", type=float, default=3.47794)
    ap.add_argument("--zmin", type=float, default=0.5,
                    help="mm hinter Facette A, ab wo gesucht wird (default 0.5 "
                         "-- davor ist alles Flanke der Facette selbst)")
    ap.add_argument("--zmax", type=float, default=40.0)
    ap.add_argument("--bg-mm", type=float, default=1.0,
                    help="Fensterbreite des gleitenden Medians fuer Stacheln")
    ap.add_argument("--min-prom", type=float, default=10.0,
                    help="dB, die ein Stachel ueber seinem oertlichen "
                         "Untergrund stehen muss")
    ap.add_argument("--env-mm", type=float, default=0.30,
                    help="Fensterbreite der Huellkurve fuer Buckel")
    ap.add_argument("--base-mm", type=float, default=4.0,
                    help="Fensterbreite des Untergrunds fuer Buckel")
    ap.add_argument("--hump-db", type=float, default=6.0,
                    help="dB, die ein Buckel ueber seinem Untergrund steht")
    ap.add_argument("--hump-min-mm", type=float, default=0.40,
                    help="Mindestbreite eines Buckels")
    ap.add_argument("--tol-um", type=float, default=100.0,
                    help="bis zu welchem Abstand eine Grenzflaeche noch als "
                         "Kandidat gilt")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if not (a.fiber is None) ^ (a.channel is None):
        raise SystemExit("genau eines von --fiber / --channel angeben")
    ch = a.channel if a.channel else "b%d" % (32 - a.fiber)

    z, db = load(a.csv)
    fl = noise_floor(db)
    x = z - a.facet
    ifs = interfaces(a.netlist, ch, a.ng)
    segs = segments(ifs)

    # ------------------------------------------------------------ Stacheln
    pk = spikes(x, db, a.zmin, a.zmax, a.bg_mm, a.min_prom)
    n_raw = len(pk)
    keep = {p[0] for p in drop_satellites([(p[0], p[1]) for p in pk])}
    pk = [p for p in pk if p[0] in keep]
    n_sat = n_raw - len(pk)

    side = sidebands(pk, 2.0 * a.tol_um * 1e-3, 0.25)
    pk.sort(key=lambda t: -t[1])
    echo = {}
    for i, (zi, _, _) in enumerate(pk):
        if zi in side:
            continue
        for zj, _, _ in pk[:i]:
            if zj in side:
                continue
            if abs(zi - 2.0 * zj) < 2.0 * a.tol_um * 1e-3:
                echo[zi] = zj
                break

    rows, spike_verdict = [], {}
    for zi, di, pi in sorted(pk):
        if zi in side:
            zc, d = side[zi]
            rows.append(["Stachel", round(zi, 4), round(zi, 4), 0.0,
                         round(di - fl, 1), round(pi, 1), "Seitenlinie",
                         "+-%.2f mm um %.3f mm" % (d, zc), "", ""])
            spike_verdict[zi] = (di, "Seitenlinie von %.3f mm" % zc)
            continue
        if zi in echo:
            rows.append(["Stachel", round(zi, 4), round(zi, 4), 0.0,
                         round(di - fl, 1), round(pi, 1), "Doppelbounce",
                         "2 x %.3f mm" % echo[zi], round((zi - 2 * echo[zi]) * 1e3, 1), ""])
            spike_verdict[zi] = (di, "Doppelbounce von %.3f mm" % echo[zi])
            continue
        dzm, name = min(ifs, key=lambda t: abs(t[0] - zi))
        dlt = (zi - dzm) * 1e3
        if abs(dlt) <= a.tol_um:
            rows.append(["Stachel", round(zi, 4), round(zi, 4), 0.0,
                         round(di - fl, 1), round(pi, 1), "Designgrenzflaeche",
                         name, round(dlt, 1),
                         round(p_chance(abs(zi - dzm), ifs, a.zmin, a.zmax), 4)])
            spike_verdict[zi] = (di, name)
        else:
            rows.append(["Stachel", round(zi, 4), round(zi, 4), 0.0,
                         round(di - fl, 1), round(pi, 1), "unerklaert",
                         "naechste: %s %+.2f mm" % (name, -dlt * 1e-3), "", ""])

    # -------------------------------------------------------------- Buckel
    hp, cover = humps(x, db, a.zmin, a.zmax, a.env_mm, a.base_mm,
                      a.hump_db, a.hump_min_mm)
    for z0, z1, w, mx, area in hp:
        mid = 0.5 * (z0 + z1)
        # Die Designstrecke, in der die Buckelmitte LIEGT -- nicht die mit
        # der aehnlichsten Laenge. Ein Buckel, der 2 mm neben der Strecke
        # sitzt, ist keine Zuordnung, auch wenn die Laenge zufaellig passt.
        # Ein Buckel um einen bereits erklaerten Stachel herum ist dessen
        # Fuss, kein eigener Befund -- sonst taucht Facette B dreimal auf.
        ins = [(d, v) for zz, (d, v) in spike_verdict.items() if z0 <= zz <= z1]
        if ins:
            d, v = max(ins)
            rows.append(["Buckel", round(z0, 4), round(z1, 4), round(w, 4),
                         round(mx - fl, 1), round(mx, 1), "Fuss eines Stachels",
                         v, "", ""])
            continue
        verdict, cand, dlt = hump_verdict(z0, z1, w, segs, ifs)
        ph = (p_hump(w, dlt, segs, a.zmin, a.zmax)
              if verdict == "Strecke gefuellt" else None)
        rows.append(["Buckel", round(z0, 4), round(z1, 4), round(w, 4),
                     round(mx - fl, 1), round(mx, 1), verdict, cand,
                     round(dlt, 1) if dlt is not None else "",
                     round(ph, 4) if ph is not None else ""])

    # ------------------------------------------------------------- Ausgabe
    print("Kanal %s, %d Grenzflaechen im Design, Fenster %.1f-%.1f mm"
          % (ch, len(ifs), a.zmin, a.zmax))
    print("Stacheln: %d ueber %+.0f dB oertlichem Untergrund, %d als bekannter "
          "Satellit verworfen, %d als Seitenlinie erkannt"
          % (n_raw, a.min_prom, n_sat, len(side)))
    print("Buckel:   %d, sie bedecken %.0f %% der Achse\n"
          % (len(hp), 100.0 * cover))
    print("%-8s %9s %9s %7s %8s  %-20s %-34s %9s %7s"
          % ("Art", "z von", "z bis", "Breite", "ueb.Bod", "Deutung",
             "Kandidat", "Abw [um]", "p"))
    for r in rows:
        print("%-8s %9.3f %9.3f %7.3f %8.1f  %-20s %-34s %9s %7s"
              % (r[0], r[1], r[2], r[3], r[4], r[6], r[7], r[8], r[9]))

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["kind", "z_from_mm", "z_to_mm", "width_mm",
                    "db_over_floor", "db_over_background", "verdict",
                    "candidate", "delta_um", "p_chance"])
        w.writerows(rows)
    print("\nCSV: %s" % a.out)


if __name__ == "__main__":
    main()
