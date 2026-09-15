"""Ligentec-SiN-Chip: Loopback-Kanaele aus dem GDS und Abgleich mit einem Scan.

Der kombinierte GDS `Ligentec_HHI_combined 4.gds` enthaelt ZWEI Chips
uebereinander:

    y =     0 .. 19313 um   LIGENTEC SiN   (Layer 2/x, 30..33, 101/2)
    y = 19372 .. 35372 um   HHI InP        (Layer 12/0, 18/0, 51..54)

Optische Schnittstellen:
    untere Kante des SiN-Chips (y=0)   95 Ports, Pitch 127 um -> FASERARRAY
    obere  Kante des SiN-Chips          132 Ports, Pitch  92 um -> Chip-zu-Chip
                                        zum InP-Chip, KEIN Faserarray

Nummerierung: Patchcord-Faser N liegt an Chipkanal 96 - N.

Die vier aeussersten Kanaele sind paarweise durch eine kurze U-Schleife
verbunden (Kanal 1<->2 und Kanal 94<->95). Das sind die einzigen Kanaele,
bei denen die Topologie die Zuordnung erzwingt -- dieselbe Situation wie
b0/b28 beim InP-Chip.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

N_FIBER = 1.468

# --- Designgeometrie, aus dem GDS gemessen -------------------------------
COUPLER_UM = 420.0      # Facette -> Wellenleiterport, AN350BB_EdgeCoupler_Lensed_C
LOOP_UM = 209.5         # Mittellinie der U-Schleife (Flaeche/Breite, w = 1.000 um)
L_LOOPBACK_UM = 2 * COUPLER_UM + LOOP_UM        # = 1049.5 um Facette zu Facette
PITCH_UM = 127.0
X0_UM = 500.0           # Kanal 1 (bzw. Slot 0) liegt bei x = 500 um


def expected_dz_mm(n_chip, L_um=L_LOOPBACK_UM):
    """Facette-A -> Facette-B im Reflektogramm, bei Chip-Gruppenindex n_chip."""
    return L_um * n_chip / N_FIBER * 1e-3


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scan", help=".json Rohscan")
    ap.add_argument("--fiber", type=int, default=None,
                    help="Patchcord-Fasernummer (Kanal = 96 - Faser)")
    ap.add_argument("--lam-start", type=float, default=1520.0)
    ap.add_argument("--lam-stop", type=float, default=1570.0)
    ap.add_argument("--zmax", type=float, default=4.0)
    ap.add_argument("--zoom", nargs=2, type=float, default=(1575.0, 1592.0),
                    metavar=("LO", "HI"))
    ap.add_argument("--control", default=None,
                    help="Kontrollscan (nichts angeschlossen), gleicher Aufbau")
    ap.add_argument("--anchor", type=float, default=529.43,
                    help="interner Reflex, auf den beide Scans normiert werden")
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import process_reflectogram_aux as P
    from identify_chip_scan import auto_tau_aux, clusters

    name = os.path.splitext(os.path.basename(a.scan))[0]
    date = name[:10] if name[:4].isdigit() else "unsortiert"
    outdir = a.outdir or os.path.join("results", date)
    os.makedirs(outdir, exist_ok=True)

    tau, fringes = auto_tau_aux(a.scan, a.lam_start, a.lam_stop)
    dl_aux = 299_792_458.0 * tau / N_FIBER
    print("Scan: %s" % a.scan)
    print("  tau_aux (aus diesem Scan) = %.4f ns  ->  aux dL = %.3f m"
          % (tau * 1e9, dl_aux))

    ns = P.build_argparser().parse_args(
        [a.scan, "--tau-aux-ns", "%.6f" % (tau * 1e9),
         "--zmax", str(a.zmax), "--peak-floor-db", "0"])
    res = P.process(ns)
    z, db = res["z"] * 1e3, res["db"]
    res_mm = float(res["dz_bin"]) * 1e3

    chan = (96 - a.fiber) if a.fiber else None
    loop = chan in (1, 2, 94, 95) if chan else None

    zmain = float(z[np.argmax(db)])
    print("  staerkster Reflex bei %.3f mm" % zmain)
    if chan:
        print("  Faser %d  ->  Chipkanal %d  (%s)"
              % (a.fiber, chan, "LOOPBACK" if loop else "kein Loopback"))

    # -------------------------------------------------- Kontrollscan
    zc = dbc = None
    if a.control:
        tau_c, _ = auto_tau_aux(a.control, a.lam_start, a.lam_stop)
        nsc = P.build_argparser().parse_args(
            [a.control, "--tau-aux-ns", "%.6f" % (tau_c * 1e9),
             "--zmax", str(a.zmax), "--peak-floor-db", "0"])
        rc = P.process(nsc)
        zc, dbc = rc["z"] * 1e3, rc["db"]

        def anchor_level(zz, dd):
            m = (zz > a.anchor - 0.15) & (zz < a.anchor + 0.15)
            return float(dd[m].max())
        off, off_c = anchor_level(z, db), anchor_level(zc, dbc)
        db, dbc = db - off, dbc - off_c     # beide relativ zum internen Anker
        print("  Kontrolle: tau_aux %.4f ns, Normierung auf %.2f mm "
              "(Versatz %.2f dB)" % (tau_c * 1e9, a.anchor, off - off_c))

    # -------------------------------------------------- Autokorrelation
    def acf(zz, dd, lo, hi):
        m = (zz >= lo) & (zz <= hi)
        v = 10 ** (dd[m] / 20.0)
        v = v - v.mean()
        r = np.correlate(v, v, "full")[len(v) - 1:]
        return np.arange(len(r)) * float(zz[1] - zz[0]), r / r[0]

    lag, r = acf(z, db, a.zoom[0] + 3.0, a.zoom[1] - 2.0)
    sel = (lag > 0.25) & (lag < 3.0)
    lag_best = float(lag[sel][np.argmax(r[sel])])
    print("  Autokorrelation der Chip-Region: staerkster Lag %.4f mm (r=%.3f)"
          % (lag_best, r[sel].max()))

    # -------------------------------------------------- Plot
    nrow = 4 if a.control else 3
    fig, ax = plt.subplots(nrow, 1, figsize=(14, 3.6 * nrow))

    ax[0].plot(z, db, lw=0.4, color="navy")
    ax[0].set_xlim(0, a.zmax * 1e3)
    ax[0].set_title("1)  ganzes Reflektogramm   %s   (tau_aux = %.4f ns, aux dL = %.2f m)"
                    % (name, tau * 1e9, dl_aux))
    zaux = 299_792_458.0 * tau / (2 * N_FIBER) * 1e3
    for zz, lab, col in [(zmain, "staerkster Reflex %.3f mm" % zmain, "r"),
                         (zaux, "aux %.0f mm" % zaux, "g"),
                         (zmain + zaux, "aux-Geist %.0f mm" % (zmain + zaux), "g")]:
        ax[0].axvline(zz, color=col, ls=":", lw=0.9)
        ax[0].text(zz, 2, lab, rotation=90, fontsize=7, color=col, va="bottom")

    lo, hi = a.zoom
    m = (z >= lo) & (z <= hi)
    ax[1].plot(z[m], db[m], lw=1.0, color="navy")
    ax[1].set_xlim(lo, hi)
    ax[1].set_title("2)  Zoom auf die Chip-Region")
    for zz, dd in sorted(clusters(z, db, lo, hi, -40.0, tol_mm=0.30),
                         key=lambda t: -t[1])[:6]:
        ax[1].annotate("%.3f mm\n%.1f dB\n(%+.3f)" % (zz, dd, zz - zmain),
                       (zz, dd), textcoords="offset points", xytext=(0, 8),
                       ha="center", fontsize=6.5, color="darkred")

    if a.control:
        ax[1].plot(zc[(zc >= lo) & (zc <= hi)], dbc[(zc >= lo) & (zc <= hi)],
                   lw=0.8, color="gray", label="Kontrolle, nichts angeschlossen")
        ax[1].legend(fontsize=7, loc="upper right")

    ax[2].plot(z[m], db[m], lw=1.0, color="navy")
    ax[2].set_xlim(lo, hi)
    ax[2].axvline(zmain, color="r", lw=1.2)
    ax[2].text(zmain, 4, "Facette A (angenommen)", rotation=90, fontsize=7,
               color="r", va="bottom")
    ax[2].set_title("3)  wo Facette B beim LOOPBACK liegen muesste "
                    "(Design %.1f um = 2x%.0f um Koppler + %.1f um Schleife)"
                    % (L_LOOPBACK_UM, COUPLER_UM, LOOP_UM))
    for n_chip, col in [(1.80, "tab:orange"), (1.90, "tab:green"), (2.00, "tab:purple")]:
        zb = zmain + expected_dz_mm(n_chip)
        ax[2].axvline(zb, color=col, ls="--", lw=1.1)
        ax[2].text(zb, 4, "n=%.2f  ->  %.3f mm" % (n_chip, zb), rotation=90,
                   fontsize=7, color=col, va="bottom")

    nf = P.noise_floor(db)
    print("  RMS-Rauschpegel %.1f dB" % nf)
    top = float(db.max()) + 10.0
    for x in ax[:3]:
        x.grid(alpha=0.3)
        x.set_xlabel("z [mm]")
        x.set_ylabel("Amplitude [dB]" + (" (rel. interner Reflex)" if a.control else ""))
        x.axhline(nf, color="darkviolet", ls="--", lw=1.1, zorder=4)
        x.set_ylim(min(top - 90.0, nf - 5.0), top)
    ax[0].annotate("RMS-Rauschpegel %.1f dB" % nf, (a.zmax * 1e3, nf + 0.8),
                   fontsize=8, ha="right", va="bottom", color="darkviolet",
                   zorder=5, bbox=dict(fc="white", ec="none", alpha=.75, pad=1.0))

    if a.control:
        sel2 = (lag > 0.20) & (lag < 2.2)
        ax[3].plot(lag[sel2], r[sel2], lw=1.0, color="navy")
        ax[3].axvline(lag_best, color="r", ls="-", lw=1.1)
        ax[3].text(lag_best, r[sel].max(), "  staerkster Lag %.4f mm" % lag_best,
                   fontsize=8, color="r", va="top")
        for n_chip, col in [(1.80, "tab:orange"), (1.90, "tab:green"),
                            (2.00, "tab:purple")]:
            ax[3].axvline(expected_dz_mm(n_chip), color=col, ls="--", lw=1.1)
        ax[3].axvspan(expected_dz_mm(1.80), expected_dz_mm(2.00),
                      color="tab:green", alpha=0.12)
        ax[3].text(expected_dz_mm(1.90), r[sel].max() * 0.6,
                   " Loopback-Fenster\n n = 1.8 .. 2.0", fontsize=8,
                   color="tab:green")
        ax[3].set_title("4)  Autokorrelation der Chip-Region — ein Facettenpaar "
                        "im Abstand dz gibt hier einen Peak bei dz")
        ax[3].set_xlabel("Lag [mm]")
        ax[3].set_ylabel("Autokorrelation")
        ax[3].grid(alpha=0.3)

    plt.tight_layout()
    png = os.path.join(outdir, name + "_ligentec.png")
    plt.savefig(png, dpi=110)
    print("  geschrieben: %s" % png)
    print("  Aufloesung %.2f um faseraequivalent" % (res_mm * 1e3))
    print("\n  Loopback-Erwartung (Design %.1f um):" % L_LOOPBACK_UM)
    for n_chip in (1.80, 1.90, 2.00):
        print("     n_chip = %.2f  ->  Facette B bei A + %.3f mm"
              % (n_chip, expected_dz_mm(n_chip)))


if __name__ == "__main__":
    main()
