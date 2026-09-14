"""Rohscan rein -> Reflektogramm, Anker, und Abgleich gegen einen Kanal.

Nimmt eine rohe HHI-Chipmessung und bestimmt daraus:
  1. tau_aux (aus dem aux-Interferometer des Scans selbst)
  2. das Reflektogramm
  3. die interne Referenzreflexion und die Chip-Eingangsfacette
  4. den Abgleich der Peaks mit der optischen Netzliste aus dem GDS
     (`gds_netlist_reflectors.py`), je Bauteil den implizierten n_g

WICHTIGE GRENZE -- bitte lesen, bevor man dem Ergebnis traut:
   Den Kanal AUS DEN DATEN zu erkennen funktioniert NICHT. Die Netzliste
   hat 690 Grenzflaechen auf 29 Kanaelen; bei einer Toleranz von zwei
   Aufloesungszellen findet sich zu praktisch jeder Vorhersage ein Peak.
   Nachweis: variiert man n_g von 3.0 bis 4.0 (`--scan-ng`), bleibt der
   beste Score praktisch konstant (7.3 .. 19.3, ohne Maximum bei 3.48).
   Ein flacher Score heisst: die Daten enthalten die Information nicht.

   Deshalb ist `--channel bNN` der vorgesehene Betrieb: man sagt, welcher
   Kanal gesteckt war, und das Tool PRUEFT die Vorhersage. Ohne
   `--channel` wird eine Rangliste ausgegeben, die ausdruecklich NICHT
   als Identifikation zu lesen ist.

   Belastbar ist der Abgleich nur dort, wo die Topologie die Zuordnung
   erzwingt: bei den Loopback-Paaren (b0/b1 und b27/b28) gibt es genau
   zwei Facetten und keinen Interpretationsspielraum.

Aufruf -- ein Befehl, alles landet im richtigen Ergebnisordner:

    python tools/identify_chip_scan.py raw_data/2026-09-11-12-16hhi1chipsignal.json --channel b27

Das legt automatisch an (Datum kommt aus dem Dateinamen, sonst heute):

    results/2026-09-11/<Scanname>_reflectogram.csv    Reflektogramm als Daten
    results/2026-09-11/<Scanname>_reflectogram.png    Reflektogramm-Plot
    results/2026-09-11/<Scanname>_identified.csv      Peak -> Bauteil, je n_g
    results/2026-09-11/<Scanname>_identified.png      annotierter Zoom

Ein bereits gerechnetes Reflektogramm geht genauso (CSV wird an der
Endung erkannt, --from-csv ist nicht noetig):

    python tools/identify_chip_scan.py results/2026-09-11/x_reflectogram.csv --channel b0

Optionen
    --channel   bekannten Kanal pruefen (empfohlen, siehe Grenze oben)
    --netlist   Reflektortabelle (Default results/2026-09-11/netlist_reflectors.csv)
    --ng        Gruppenindex des Chips fuer die Vorhersage (Default 3.4791)
    --tol       Toleranz in Aufloesungszellen (Default 2.0)
    --scan-ng   Bestimmbarkeitstest ueber n_g = 3.0 .. 4.0
    --out       eigener Prefix statt der Automatik
"""

import argparse
import csv
import datetime
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import process_reflectogram_aux as P   # noqa: E402  (load/balanced/analytic/process)

N_FIBER = 1.468
SSC = 1199.4          # um, Port-zu-Port eines SSC laut Netzliste
REF_WIN = (500.0, 650.0)      # mm, wo die interne Referenzreflexion liegt
PIGTAIL = (1090.0, 1103.0)    # mm, gemessene Pigtail-Familie des Faserarrays


# --------------------------------------------------------------- Hilfen
def auto_tau_aux(scan, lam_start, lam_stop, aux=(2, 4), trim=0.02):
    """tau_aux aus dem aux-Paar des Scans selbst (wie check_aux_interferometer)."""
    ch = dict(zip((1, 2, 3, 4), P.load(scan)[:4]))
    a, b = ch[aux[0]], ch[aux[1]]
    s, _ = P.balanced(a, b)
    n = len(s)
    k0, k1 = int(trim * n), int((1 - trim) * n)
    phi = np.unwrap(np.angle(P.analytic(s[k0:k1])))
    fringes = abs(phi[-1] - phi[0]) / (2 * np.pi)
    c = 299_792_458.0
    span = c * (1.0 / (lam_start * 1e-9) - 1.0 / (lam_stop * 1e-9))
    span *= (k1 - k0) / n
    return fringes / span, fringes


def clusters(z, db, lo, hi, floor, tol_mm=0.25):
    from scipy.signal import find_peaks
    m = (z >= lo) & (z <= hi)
    zz, dd = z[m], db[m]
    if not len(zz):
        return []
    pk, _ = find_peaks(dd, height=floor, distance=3)
    zp, dp = zz[pk], dd[pk]
    used = np.zeros(len(zp), bool)
    out = []
    for i in np.argsort(dp)[::-1]:
        if used[i]:
            continue
        used |= (np.abs(zp - zp[i]) < tol_mm) & (~used)
        out.append((zp[i], dp[i]))
    return sorted(out)


def load_netlist(path):
    net = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            net.setdefault(r["channel"], []).append(
                (float(r["waveguide_path_um"]), r["device_port"]))
    for k in net:
        net[k] = sorted(set(net[k]))
    return net


# --------------------------------------------------------- Identifikation
SAT_MM = (0.70, 1.53)     # bekannte Instrument-Satelliten um jeden starken Peak
SAT_TOL = 0.09


def drop_satellites(zpk):
    """Peaks entfernen, die als +-0.70/+-1.53-mm-Satellit eines staerkeren
    Peaks erklaerbar sind. Diese Seitenlinien treten um JEDEN starken
    Reflektor auf (auch um die interne Referenz) und sind kein Chipsignal --
    ohne sie wegzuwerfen erzeugen sie Scheintreffer."""
    keep = []
    for zz, dd in zpk:
        sat = any(abs(abs(zz - z2) - s) < SAT_TOL
                  for z2, d2 in zpk if d2 > dd + 3.0 for s in SAT_MM)
        if not sat:
            keep.append((zz, dd))
    return keep


def identify(zpk, facet_candidates, net, ng, res_mm, tol_cells):
    """-> nach Score sortierte Liste (score, -mean_err, channel, facetA, hits).

    Score = Summe der amplitudengewichteten Treffer, abgestraft mit der
    Streuung des implizierten Gruppenindex. Letzteres ist der eigentliche
    Diskriminator: ein zufaellig passender Peaksatz ergibt widerspruechliche
    n_g-Werte, der richtige Kanal ergibt denselben Wert an jeder
    Grenzflaeche.
    """
    tol = tol_cells * res_mm
    ranked = []
    for fa, fa_db in facet_candidates:
        obs = [(zz - fa, dd) for zz, dd in zpk if zz - fa > 0.5]
        if not obs:
            continue
        for ch, ifaces in net.items():
            pred = [((SSC + d) * 1e-3 * ng / N_FIBER, dev, d)
                    for d, dev in ifaces]
            hits, errs, ngs = [], [], []
            for dzo, dbo in obs:
                best = min(pred, key=lambda t: abs(t[0] - dzo))
                if abs(best[0] - dzo) < tol:
                    hits.append((dzo, dbo, best[0], best[1],
                                 (dzo - best[0]) * 1e3))
                    errs.append(abs(dzo - best[0]) * 1e3)
                    ngs.append(N_FIBER * dzo / ((SSC + best[2]) * 1e-3))
            if not hits:
                continue
            # Amplitudengewicht: 0 dB zaehlt 10x so viel wie -40 dB
            w = sum(1.0 + 9.0 * 10 ** (d / 20.0) for _, d, _, _, _ in hits)
            # Konsistenzstrafe: streuen die implizierten n_g, war es Zufall
            scat = (max(ngs) - min(ngs)) / float(np.mean(ngs)) if len(ngs) > 1 else 0.0
            score = w / (1.0 + 60.0 * scat)
            ranked.append((round(score, 2), -float(np.mean(errs)), ch, fa,
                           fa_db, hits))
    ranked.sort(reverse=True)
    return ranked


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan", help=".json Rohscan, oder ein *_reflectogram.csv mit --from-csv")
    ap.add_argument("--from-csv", action="store_true")
    ap.add_argument("--netlist", default="results/2026-09-11/netlist_reflectors.csv")
    ap.add_argument("--ng", type=float, default=3.4791)
    ap.add_argument("--tol", type=float, default=2.0, help="Toleranz in Aufloesungszellen")
    ap.add_argument("--lam-start", type=float, default=1520.0)
    ap.add_argument("--lam-stop", type=float, default=1570.0)
    ap.add_argument("--zmax", type=float, default=3.0)
    ap.add_argument("--peak-floor-db", type=float, default=-30.0)
    ap.add_argument("--channel", default=None,
                    help="bekannten Kanal pruefen statt raten, z.B. b27")
    ap.add_argument("--scan-ng", action="store_true",
                    help="Score ueber n_g=3.0..4.0 auftragen (Bestimmbarkeitstest)")
    ap.add_argument("--out", default=None,
                    help="Prefix fuer die Ausgaben. Ohne Angabe automatisch "
                         "results/<Datum>/<Scanname> -- Datum aus dem "
                         "Dateinamen, sonst heute.")
    a = ap.parse_args()

    # CSV wird am Suffix erkannt, damit --from-csv nicht noetig ist
    if a.scan.lower().endswith(".csv"):
        a.from_csv = True

    # Ausgabeordner automatisch: results/<Datum>/<Scanname ohne Endung>
    if a.out is None:
        base = os.path.splitext(os.path.basename(a.scan))[0]
        for suf in ("_reflectogram", "_trimmed1530", "_full"):
            if base.endswith(suf):
                base = base[:-len(suf)]
        m = re.match(r"(\d{4}-\d{2}-\d{2})", base)
        day = m.group(1) if m else datetime.date.today().isoformat()
        a.out = os.path.join("results", day, base)
        os.makedirs(os.path.dirname(a.out), exist_ok=True)
        print("Ausgabeordner: %s/  (Prefix %s)"
              % (os.path.dirname(a.out), os.path.basename(a.out)))

    # ---------------------------------------------------------- 1) Daten
    if a.from_csv:
        d = np.genfromtxt(a.scan, delimiter=",", names=True)
        z, db = d["distance_m"] * 1e3, d["amplitude_dB"]
        res_mm = float(z[1] - z[0])
        print("Reflektogramm aus CSV: %s" % a.scan)
    else:
        tau, fringes = auto_tau_aux(a.scan, a.lam_start, a.lam_stop)
        print("Scan: %s" % a.scan)
        print("  tau_aux automatisch aus dem aux-Paar: %.4f ns  "
              "(%.0f Fringes)  ->  aux dL = %.4f m"
              % (tau * 1e9, fringes, 299_792_458.0 * tau / 1.468))
        ns = P.build_argparser().parse_args(
            [a.scan, "--tau-aux-ns", "%.6f" % (tau * 1e9),
             "--zmax", str(a.zmax), "--peak-floor-db", "0",
             "--out", a.out] if a.out else
            [a.scan, "--tau-aux-ns", "%.6f" % (tau * 1e9),
             "--zmax", str(a.zmax), "--peak-floor-db", "0"])
        res = P.process(ns)
        z, db = res["z"] * 1e3, res["db"]
        res_mm = float(res["dz_bin"]) * 1e3
    print("  Auflaesung: %.2f um faseraequivalent = %.2f um on-chip"
          % (res_mm * 1e3, res_mm * 1e3 * N_FIBER / a.ng))

    # ------------------------------------------- 2) Anker: Referenz + Facette
    refs = clusters(z, db, REF_WIN[0], REF_WIN[1], -45)
    if not refs:
        sys.exit("keine interne Referenzreflexion in %s mm gefunden" % (REF_WIN,))
    ref = max(refs, key=lambda t: t[1])
    print("  interne Referenzreflexion: %.4f mm (%.2f dB)" % ref)

    fa_c = clusters(z, db, ref[0] + PIGTAIL[0], ref[0] + PIGTAIL[1], -34)
    if not fa_c:
        sys.exit("keine Facette-A-Kandidaten im Pigtail-Fenster")
    print("  Facette-A-Kandidaten (Pigtail %s mm): %s"
          % (PIGTAIL, ", ".join("%.4f (%.1f dB, Pigtail %.2f mm)"
                                % (p[0], p[1], p[0] - ref[0]) for p in fa_c)))

    # ------------------------------------------------------- 3) Abgleich
    zpk_all = clusters(z, db, ref[0] + PIGTAIL[0], ref[0] + PIGTAIL[1] + 40,
                       a.peak_floor_db)
    zpk = drop_satellites(zpk_all)
    print("  Peaks im Chipbereich: %d, davon %d als Instrument-Satelliten "
          "verworfen -> %d verwendet"
          % (len(zpk_all), len(zpk_all) - len(zpk), len(zpk)))
    net = load_netlist(a.netlist)
    if a.channel:
        if a.channel not in net:
            sys.exit("Kanal %s nicht in der Netzliste" % a.channel)
        net = {a.channel: net[a.channel]}
    ranked = identify(zpk, fa_c, net, a.ng, res_mm, a.tol)

    if a.scan_ng:
        full = load_netlist(a.netlist)
        sc = []
        for ngx in np.arange(3.0, 4.02, 0.02):
            r = identify(zpk, fa_c, full, ngx, res_mm, a.tol)
            sc.append((ngx, r[0][0] if r else 0.0))
        arr = np.array([t[1] for t in sc])
        print("\n  Bestimmbarkeitstest: bester Score ueber n_g = 3.0 .. 4.0")
        print("     Spanne %.2f .. %.2f, Maximum bei n_g = %.2f"
              % (arr.min(), arr.max(), sc[int(np.argmax(arr))][0]))
        print("     %s" % ("-> FLACH: n_g/Kanal sind aus diesen Daten NICHT bestimmbar"
                           if arr.max() < 1.5 * np.median(arr[arr > 0])
                           else "-> ausgepraegtes Maximum: bestimmbar"))
    if not ranked:
        print("\nKEINE Zuordnung gefunden -- kein Kanal erklaert die Peaks.")
        print("Das ist bei Stumpf-Kanaelen normal (die haben keine Bauteile).")
        return

    print("\n%s" % ("=" * 78))
    print("ABGLEICH GEGEN %s" % (a.channel if a.channel else
                                 "ALLE KANAELE -- Rangliste, KEINE Identifikation"))
    print("=" * 78)
    if not a.channel and len(ranked) > 1 and ranked[1][0] > 0.8 * ranked[0][0]:
        print("  WARNUNG: Rang 1 (%.2f) und Rang 2 (%.2f) liegen dicht beieinander"
              % (ranked[0][0], ranked[1][0]))
        print("  -> die Zuordnung ist NICHT eindeutig. Mit --channel arbeiten.")
    print("%6s %9s %8s %10s %s" % ("Rang", "Kanal", "Score", "Fehler", "Facette A / Pigtail"))
    for i, (sc, negerr, ch, fa, fadb, hits) in enumerate(ranked[:5]):
        print("%6d %9s %8.2f %8.1f um %10.4f mm / %.2f mm%s"
              % (i + 1, ch, sc, -negerr, fa, fa - ref[0],
                 "   <== bester Fit (nicht: bewiesen)" if i == 0 else ""))

    sc, negerr, ch, fa, fadb, hits = ranked[0]
    print("\nINNERE STRUKTUR von %s  (Facette A bei %.4f mm, %.2f dB)" % (ch, fa, fadb))
    print("%10s %7s %12s %12s %9s %10s  %s"
          % ("dz [mm]", "dB", "on-chip um", "Design um", "Fehler", "-> n_g", "Bauteil / Port"))
    rows = []
    for dzo, dbo, dzp, dev, err in sorted(hits):
        onchip = dzo * 1e3 * N_FIBER / a.ng
        design = dzp * 1e3 * N_FIBER / a.ng
        ng_i = N_FIBER * dzo / (design * 1e-3)
        print("%10.4f %7.1f %12.1f %12.1f %+8.1f %10.5f  %s"
              % (dzo, dbo, onchip, design, err, ng_i, dev))
        rows.append([ch, round(dzo, 4), round(dbo, 2), round(onchip, 1),
                     round(design, 1), round(err, 1), round(ng_i, 5), dev])
    ngs = [r[6] for r in rows]
    print("\n  %d Bauteilgrenzflaechen zugeordnet, Toleranz %.1f Zellen (%.1f um)"
          % (len(rows), a.tol, a.tol * res_mm * 1e3))
    if len(ngs) > 1:
        print("  implizierter Gruppenindex: %.5f .. %.5f  (Streuung %.2f %%)"
              % (min(ngs), max(ngs), 100 * (max(ngs) - min(ngs)) / np.mean(ngs)))

    # --------------------------------------------------------- 4) Ausgabe
    if a.out:
        with open(a.out + "_identified.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["channel", "delta_z_mm", "amplitude_dB", "onchip_um",
                        "design_um", "error_um", "implied_n_g", "device_port"])
            w.writerows(rows)
        print("\nCSV: %s_identified.csv" % a.out)
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(13, 6))
            m = (z > fa - 3) & (z < fa + 22)
            ax.plot(z[m] - fa, db[m], lw=.9, color="navy")
            ax.axvline(0, color="crimson", ls="--", lw=1.3)
            ax.annotate("facet A\n(fibre -> chip)", (0, -48), fontsize=9,
                        ha="center", color="crimson",
                        bbox=dict(fc="w", ec="crimson", lw=.8, alpha=.9))
            for k, (dzo, dbo, dzp, dev, err) in enumerate(sorted(hits)):
                ax.axvline(dzo, color="seagreen", ls=":", lw=1.3)
                ax.annotate("%s\n%.4f mm (%+.1f um)"
                            % (dev.split(":")[0].replace("HHI_", ""), dzo, err),
                            (dzo, -8 - 9 * (k % 3)), fontsize=8, ha="center",
                            color="darkgreen",
                            bbox=dict(fc="w", ec="seagreen", lw=.7, alpha=.93))
            ax.set_xlabel("distance behind facet A  [mm]")
            ax.set_ylabel("amplitude  [dB]")
            ax.set_ylim(-55, 4)
            ax.grid(alpha=.3)
            ax.set_title("%s  ->  identified as channel %s : %d component "
                         "interfaces matched (n_g = %.4f)"
                         % (os.path.basename(a.scan), ch, len(rows), a.ng))
            fig.tight_layout()
            fig.savefig(a.out + "_identified.png", dpi=130)
            print("PNG: %s_identified.png" % a.out)
        except Exception as exc:
            print("(Plot uebersprungen: %s)" % exc)


if __name__ == "__main__":
    main()
