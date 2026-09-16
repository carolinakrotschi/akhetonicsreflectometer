"""Vorhersage aus dem Ligentec-GDS gegen einen gemessenen Scan halten.

Rechnet mit `gds_netlist_ligentec.py` die Weglaenge Facette -> erstes
Bauteil eines Kanals und sieht im Reflektogramm nach, ob dort ein Peak
ist. Entscheidend ist dabei NICHT, ob im Vorhersagefenster ein Peak
liegt, sondern ob dort MEHR ist als

  a) in einem beliebigen anderen Fenster desselben Scans, und
  b) an derselben Stelle in einem Kanal, der dieses Bauteil gar nicht hat.

Ohne diese zwei Kontrollen findet man in einem dichten Spektrum immer
einen Treffer -- derselbe Fehler, der beim HHI-Chip einmal ein
"unabhaengig bestimmtes n_g" vorgetaeuscht hat.

Aufruf
    python tools/check_ligentec_channel.py --fiber 44 \
        --scan raw_data/2672_ligentechhi_2026-09-14-16-36_Ligentecfiber44_reflectogram.csv \
        --ref  raw_data/2672_ligentechhi_2026-09-14-16-02_Ligentecfiber1_reflectogram.csv
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from process_reflectogram_aux import noise_floor          # noqa: E402
import gds_netlist_ligentec as G                          # noqa: E402

N_FIBER = 1.468
NG_LO, NG_HI = 1.80, 2.00     # plausibler Bereich fuer AN350-SiN


def load(p):
    d = np.genfromtxt(p, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]


def first_device(ch):
    """(Weg ab Facette in um, Bauteilname) des ersten erreichten Bauteils."""
    nodes, edges, _ = G.build(G.GDS)
    dist, off = G.dijkstra(nodes, edges, (G.channel_x(ch), G.WG_START_Y), 1e9)
    dev = sorted((d, nodes[i][2]) for i, d in dist.items()
                 if nodes[i][2].startswith("dev:"))
    reach = max(dist.values()) if dist else 0.0
    if dev:
        return G.COUPLER_UM + dev[0][0], dev[0][1][4:], reach, off
    return G.COUPLER_UM + reach, None, reach, off


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fiber", type=int, required=True)
    ap.add_argument("--scan", required=True, help="_reflectogram.csv des Kanals")
    ap.add_argument("--ref", default=None,
                    help="_reflectogram.csv eines Kanals OHNE dieses Bauteil")
    ap.add_argument("--path-um", type=float, default=None,
                    help="Weglaenge ab Facette selbst vorgeben (statt GDS)")
    ap.add_argument("--label", default="Bauteil")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    ch = 96 - a.fiber
    if a.path_um:
        s_um, dev = a.path_um, a.label
    else:
        s_um, dev, reach, off = first_device(ch)
        print("Kanal %d (Faser %d): Startversatz %.3f um, Pfad reicht %.1f um"
              % (ch, a.fiber, off, reach))
        if dev is None:
            print("  KEIN Bauteil erreicht -- der Wellenleiter endet an einer "
                  "Black Box, deren Inneres nicht auf Layer 2/0 gezeichnet ist.")
            dev = a.label
    print("  Weg Facette -> %s : %.1f um on-chip" % (dev, s_um))

    dz_lo = s_um * NG_LO / N_FIBER * 1e-3
    dz_hi = s_um * NG_HI / N_FIBER * 1e-3
    print("  erwartet %.3f .. %.3f mm hinter Facette A (n_g %.2f..%.2f)"
          % (dz_lo, dz_hi, NG_LO, NG_HI))

    z, db = load(a.scan)
    nf = noise_floor(db)
    za = float(z[np.argmax(db)])
    m = (z >= za + dz_lo) & (z <= za + dz_hi)
    snr = db[m].max() - nf if m.any() else float("nan")
    print("\n  Facette A %.4f mm, RMS-Floor %.1f dB" % (za, nf))
    print("  im Vorhersagefenster: max %.1f dB  -> SNR %.1f dB"
          % (db[m].max(), snr))

    # Kontrolle a): wie typisch ist dieser SNR im selben Scan?
    width = dz_hi - dz_lo
    offs = np.arange(1.0, 40.0, width)
    best = []
    for o in offs:
        w = (z >= za + o) & (z < za + o + width)
        if w.any():
            best.append(db[w].max() - nf)
    best = np.array(best)
    rank = int((best >= snr).sum())
    print("  Kontrolle a) gleich breite Fenster +1..+40 mm: %d von %d "
          "erreichen >= diesen SNR  -> Rang %d"
          % (rank, len(best), rank))

    # Kontrolle b): derselbe Ort in einem Kanal ohne dieses Bauteil
    ref = None
    if a.ref:
        zr, dbr = load(a.ref)
        nfr = noise_floor(dbr)
        zar = float(zr[np.argmax(dbr)])
        mr = (zr >= zar + dz_lo) & (zr <= zar + dz_hi)
        snrr = dbr[mr].max() - nfr
        ref = (zr, dbr, nfr, zar, snrr)
        print("  Kontrolle b) Referenzkanal an derselben Stelle: SNR %.1f dB"
              % snrr)
        print("\n  URTEIL: %s" % (
            "kein Nachweis -- der Referenzkanal ist dort genauso stark"
            if snrr >= snr - 3.0 else
            "im Vorhersagefenster deutlich mehr als im Referenzkanal"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(13, 8), height_ratios=[1, 1.3])

    w0 = (z >= za - 5) & (z <= za + 45)
    ax[0].plot(z[w0] - za, db[w0], lw=.7, color="navy", label="Faser %d" % a.fiber)
    ax[0].axhline(nf, color="darkviolet", ls="--", lw=1.1)
    ax[0].annotate("RMS-Floor %.1f dB" % nf, (44, nf + 1), fontsize=8,
                   ha="right", color="darkviolet")
    ax[0].axvspan(dz_lo, dz_hi, color="gold", alpha=.5, zorder=0)
    ax[0].set_title("1)  Faser %d (Kanal %d), ganze Chipregion hinter Facette A"
                    % (a.fiber, ch), fontsize=11, loc="left")
    ax[0].legend(fontsize=8, loc="upper right")

    w1 = (z >= za + dz_lo - 6) & (z <= za + dz_hi + 6)
    ax[1].plot(z[w1] - za, db[w1], lw=1.0, color="navy",
               label="Faser %d (hat das Bauteil), Floor %.1f dB" % (a.fiber, nf))
    ax[1].axhline(nf, color="navy", ls="--", lw=.9, alpha=.6)
    if ref:
        zr, dbr, nfr, zar, snrr = ref
        wr = (zr >= zar + dz_lo - 6) & (zr <= zar + dz_hi + 6)
        ax[1].plot(zr[wr] - zar, dbr[wr], lw=1.0, color="darkorange", alpha=.85,
                   label="Referenzkanal (hat es NICHT), Floor %.1f dB" % nfr)
        ax[1].axhline(nfr, color="darkorange", ls="--", lw=.9, alpha=.6)
    ax[1].axvspan(dz_lo, dz_hi, color="gold", alpha=.5, zorder=0)
    ax[1].annotate("Vorhersage %s\n%.0f um on-chip\nn_g %.2f..%.2f"
                   % (dev, s_um, NG_LO, NG_HI),
                   ((dz_lo + dz_hi) / 2, ax[1].get_ylim()[1]), fontsize=8.5,
                   ha="center", va="top",
                   bbox=dict(fc="white", ec="darkgoldenrod", alpha=.9))
    ax[1].set_title("2)  Zoom auf das Vorhersagefenster -- entscheidend ist der "
                    "Vergleich mit dem Referenzkanal", fontsize=11, loc="left")
    ax[1].legend(fontsize=8, loc="upper right")

    for x in ax:
        x.grid(alpha=.3)
        x.set_xlabel("Abstand hinter Facette A  [mm]")
        x.set_ylabel("Amplitude  [dB]")

    out = a.out or ("results/%s/ligentec_fiber%d_vorhersage.png"
                    % (os.path.basename(a.scan)[:10], a.fiber))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=125)
    print("\nPNG: %s" % out)


if __name__ == "__main__":
    main()
