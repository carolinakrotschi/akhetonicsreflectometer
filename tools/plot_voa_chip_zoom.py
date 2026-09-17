#!/usr/bin/env python3
"""Zeigt eine VOA-Serie im CHIPBEREICH -- kommen die Chipreflexionen heraus?

Der VOA sitzt im Referenzarm des Messinterferometers (hinter dem 20-dB-
Koppler geht ein Teil ueber den Zirkulator auf den Chip, der andere ueber
den VOA; beide treffen sich im 3-dB-Koppler wieder). Er stellt also die
**LO-Leistung**, nicht das Licht auf dem Chip. Fuer die Frage "sehe ich
die Chipreflexionen?" heisst das:

  Schwebungsamplitude  ~  sqrt(P_LO * P_sig)

Mehr LO hebt JEDE Reflexion um denselben Faktor -- und mit ihr das
Schrotrauschen, aber nur mit der Wurzel. Solange das Verstaerkerrauschen
dominiert, gewinnt man deshalb SNR; sobald das Schrotrauschen dominiert,
gewinnt man nichts mehr. Genau das entscheidet die Serie.

Deshalb ist die y-Achse hier **dB ueber dem eigenen Rauschboden**, nicht
dB auf eine gemeinsame Referenz. Der absolute Pegel ist die falsche
Frage: sichtbar ist, was aus dem eigenen Rauschen herausragt.
`tools/plot_voa_series.py` macht die andere Ansicht.

Der Rauschboden wird in einem **Kontrollfenster hinter dem Chip** bestimmt
(`--floor-mm`), nicht im Chipbereich selbst -- dort steckt ja das Signal,
das gemessen werden soll.

ZWEI KONTROLLEN, OHNE DIE DAS NICHTS HEISST

In einem dichten Spektrum findet man immer einen Peak. Verglichen wird
deshalb gegen zwei Fenster, in denen nichts sein kann:

1. **Faser weit vor der Facette** (`--control-mm`). Ruhige Faser, misst
   den nackten Boden. Steigt die Zahl der Funde dort genauso, ist der
   Zugewinn kein Chipsignal, sondern nur ein tiefer gelegter Boden.

2. **Das Spiegelfenster** dA -(hi..2.2) mm, also derselbe Bereich vor der
   Facette statt dahinter. Das ist die schaerfere Kontrolle: ein starker
   Peak traegt einen Rauschsaum (Laserphasenrauschen, spektrale Leckage),
   und der ist **symmetrisch**. Was im Chipfenster steht, aber im
   Spiegelfenster genauso, ist Saum der Facette und kein Bauteil. Nur der
   Ueberschuss ueber das Spiegelfenster ist ein Kandidat.

Dasselbe gilt fuer das GDS-Fenster: ausgegeben wird auch sein Spiegelbild
A-8.81..-7.93 mm.

Ebenfalls eingezeichnet und NICHT als Befund zu lesen:
  - die Instrumenten-Satelliten +-0.70 / +-1.53 mm um jeden starken Peak
    (`identify_chip_scan.SAT_MM`), hier um die Facette
  - der Nahbereich dA < 2.2 mm, in dem am 15.09. Faser 44 und 48
    pegelgleich waren -- Seitenbaender, kein Chip
  - das GDS-Fenster A+7.93..8.81 mm (Kanal 48, 6469.8 um ab Facette,
    n_g 1.80-2.00), siehe logs/2026-09-15.md

Aufruf
    python tools/plot_voa_chip_zoom.py \
        --series 2680_20260917-095757_fiber48_fullnmrange \
        --out-dir results/2026-09-17
"""

import argparse
import datetime
import os
import sys

import numpy as np
from scipy.signal import find_peaks

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import process_reflectogram_aux as pra      # noqa: E402
import plot_voa_series as pv                # noqa: E402
from identify_chip_scan import SAT_MM       # noqa: E402

ROOT = os.path.dirname(_HERE)

# Vorhersage aus dem Ligentec-GDS fuer Kanal 48, logs/2026-09-15.md
GDS_WINDOW_MM = (7.93, 8.81)
GDS_LABEL = "GDS Kanal 48: erste Black Box"
# Nahbereich, am 15.09. in Faser 44 und 48 pegelgleich -> Instrument
NEARFIELD_MM = 2.2


def facet_of(r, lo, hi):
    """Facette A = staerkster Peak im Suchfenster [lo, hi] in m."""
    m = (r["z"] >= lo) & (r["z"] <= hi)
    idx = np.flatnonzero(m)
    return r["z"][idx[int(np.argmax(r["R"][m]))]]


def count_visible(da, snr, win, thr, min_sep_mm=0.2):
    """(n Peaks ueber thr, max SNR, Sockelhoehe) im Fenster win [mm].

    Die Sockelhoehe ist der Median -- fuer eine breite Erhebung (Saum,
    Multipath-Wald) sagt sie mehr als das Maximum, das eine einzelne
    Rauschspitze sein kann."""
    m = (da >= win[0]) & (da <= win[1])
    if not m.any():
        return 0, float("nan"), float("nan")
    d, s = da[m], snr[m]
    step = d[1] - d[0]
    pk, _ = find_peaks(s, height=thr, distance=max(3, int(min_sep_mm / step)))
    return len(pk), float(s.max()), float(np.median(s))


def build_argparser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--series", default=None,
                   help="series_id -- alle Laeufe einer VOA-Serie")
    p.add_argument("--scan", nargs=2, action="append", default=None,
                   metavar=("LABEL", "PFAD"),
                   help="einzelner Scan mit freier Beschriftung, statt "
                        "--series und mehrfach angebbar. Damit laeuft "
                        "dieselbe Spiegelfenster-Kontrolle auf Scans aus "
                        "verschiedenen Aufbauzustaenden")
    p.add_argument("--tag", default=None,
                   help="Praefix der Ausgabedateien (default: --series)")
    p.add_argument("--raw-dir", default=os.path.join(ROOT, "raw_data"))
    p.add_argument("--out-dir", default=None)
    p.add_argument("--facet-search-mm", nargs=2, type=float,
                   default=[1570.0, 1590.0],
                   help="Suchfenster fuer Facette A (Default passt zu "
                        "MAP2680 Faser 48, s. logs/2026-09-16.md)")
    p.add_argument("--chip-mm", nargs=2, type=float, default=[-40.0, 40.0],
                   help="dargestellter Bereich relativ zu Facette A. Die "
                        "negative Seite ist das Spiegelfenster und gehoert "
                        "mit ins Bild -- ohne sie sieht man nicht, dass der "
                        "Facettensaum symmetrisch ist")
    p.add_argument("--floor-mm", nargs=2, type=float, default=[60.0, 160.0],
                   help="Kontrollfenster HINTER dem Chip, in dem der "
                        "Rauschboden bestimmt wird")
    p.add_argument("--control-mm", nargs=2, type=float,
                   default=[-200.0, -40.0],
                   help="Kontrollfenster VOR der Facette (Faser, dort kann "
                        "nichts sein) fuer die Nullverteilung")
    p.add_argument("--snr-thresh", type=float, default=10.0,
                   help="ab wie viel dB ueber dem Boden ein Peak als "
                        "sichtbar gezaehlt wird")
    p.add_argument("--window", default="kaiser",
                   choices=["hann", "blackmanharris", "kaiser"])
    p.add_argument("--kaiser-beta", type=float, default=12.0)
    p.add_argument("--trim", type=float, default=0.01)
    p.add_argument("--tau-aux-ns", type=float, default=None)
    p.add_argument("--aux-a", type=int, default=2, choices=[1, 2, 3, 4])
    p.add_argument("--aux-b", type=int, default=4, choices=[1, 2, 3, 4])
    p.add_argument("--meas-a", type=int, default=1, choices=[1, 2, 3, 4])
    p.add_argument("--meas-b", type=int, default=3, choices=[1, 2, 3, 4])
    return p


def plot(runs, a, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    nrun = len(runs)
    fig, axes = plt.subplots(nrun, 1, figsize=(13.5, 1.5 * nrun + 1.3),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    lo, hi = a.chip_mm

    tops = [r["snr"][(r["da"] >= lo) & (r["da"] <= hi)].max() for r in runs]
    ymax = max(tops) + 0.22 * max(tops)

    for k, (ax, r) in enumerate(zip(axes, runs)):
        m = (r["da"] >= lo) & (r["da"] <= hi)
        col = cmap(0.08 + 0.82 * k / max(nrun - 1, 1))

        # Zonen, die KEIN Befund sind
        ax.axvspan(lo, -NEARFIELD_MM, color="#dbe7f3", zorder=0)   # Spiegel
        ax.axvspan(-NEARFIELD_MM, NEARFIELD_MM, color="0.86", zorder=0)
        ax.axvspan(*GDS_WINDOW_MM, color="#ffe08a", alpha=0.6, zorder=0)
        ax.axvspan(-GDS_WINDOW_MM[1], -GDS_WINDOW_MM[0], color="#ffe08a",
                   alpha=0.3, zorder=0)
        for s in SAT_MM:
            for sgn in (1, -1):
                ax.axvline(sgn * s, color="0.55", ls=":", lw=0.9, zorder=1)

        ax.plot(r["da"][m], r["snr"][m], lw=0.55, color=col, zorder=3)
        ax.axhline(0, color="darkviolet", ls="--", lw=0.9, zorder=4)
        ax.axhline(a.snr_thresh, color="crimson", ls="-.", lw=0.9, zorder=4)

        if r.get("label"):
            head = r["label"]
        else:
            tag = "VOA %s" % r["voa"] if r["voa"] is not None else "VOA ?"
            head = "Lauf %d  |  %s" % (r["run"], tag)
        ax.text(0.006, 0.94, head,
                transform=ax.transAxes, ha="left", va="top", fontsize=9,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec="gray", alpha=0.9))
        ax.text(0.994, 0.94,
                "Facette %.3f mm, SNR %.1f dB   |   Chip dA %.1f..%.0f mm: "
                "max %.1f, Sockel %.1f dB   |   Spiegel: max %.1f, "
                "Sockel %.1f dB   |   Ueberschuss %+.1f / %+.1f dB"
                % (r["facet"] * 1e3, r["snr_facet"], NEARFIELD_MM, hi,
                   r["max_chip"], r["med_chip"], r["max_mir"], r["med_mir"],
                   r["max_chip"] - r["max_mir"],
                   r["med_chip"] - r["med_mir"]),
                transform=ax.transAxes, ha="right", va="top", fontsize=7.2,
                color="0.2",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                          alpha=0.8))
        ax.grid(alpha=0.28)
        ax.set_ylim(-12, ymax)
        ax.set_ylabel("dB ueber\nBoden", fontsize=7.5)
        ax.tick_params(labelsize=8)

    axes[-1].set_xlabel("dA -- Abstand hinter Chipfacette A (mm, "
                        "einfacher Weg; Achse mit n_g = 1.468 der Faser)")
    axes[0].set_xlim(lo, hi)

    handles = [
        plt.Line2D([], [], color="darkviolet", ls="--", lw=1.2,
                   label="eigener Rauschboden (0 dB)"),
        plt.Line2D([], [], color="crimson", ls="-.", lw=1.2,
                   label="Sichtbarkeitsschwelle %.0f dB" % a.snr_thresh),
        plt.Line2D([], [], color="0.55", ls=":", lw=1.2,
                   label="Instrument-Satelliten %s mm"
                         % "/".join("%.2f" % s for s in SAT_MM)),
        plt.Rectangle((0, 0), 1, 1, fc="0.86",
                      label="Nahbereich +-%.1f mm: Seitenbaender, kein Chip"
                            % NEARFIELD_MM),
        plt.Rectangle((0, 0), 1, 1, fc="#dbe7f3",
                      label="Spiegelfenster (Faser -- hier kann nichts sein)"),
        plt.Rectangle((0, 0), 1, 1, fc="#ffe08a", alpha=0.6,
                      label="%s (A+%.2f..%.2f mm) und sein Spiegel"
                            % (GDS_LABEL, *GDS_WINDOW_MM)),
    ]
    axes[0].legend(handles=handles, fontsize=7, ncol=3, loc="lower left",
                   bbox_to_anchor=(0.0, 1.02), frameon=False)

    fig.suptitle(
        "%s  --  Chipbereich, VOA im Referenzarm (LO)\n"
        "Rauschboden je Lauf aus dA %.0f..%.0f mm (hinter dem Chip); "
        "Kontrollfenster Faser dA %.0f..%.0f mm"
        % (a.tag, *a.floor_mm, *a.control_mm), fontsize=10.5)
    # Raender in ZOLL reservieren, nicht als fester Bruchteil: bei drei
    # Panels ist die Figur 5.8 statt 16.3 Zoll hoch, und 0.05 davon reichen
    # fuer das Achsenlabel nicht mehr aus (es wurde abgeschnitten).
    h = 1.5 * nrun + 1.3
    fig.subplots_adjust(top=1 - 1.45 / h, bottom=max(0.05, 0.5 / h),
                        left=0.065, right=0.99, hspace=0.12)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote: %s" % path)


def main():
    a = build_argparser().parse_args()
    out_dir = a.out_dir or os.path.join(ROOT, "results",
                                        datetime.date.today().isoformat())
    os.makedirs(out_dir, exist_ok=True)

    if bool(a.series) == bool(a.scan):
        sys.exit("entweder --series oder --scan angeben, nicht beides")
    if a.series:
        files = [(run, None, path)
                 for run, path in pv.series_files(a.raw_dir, a.series)]
        if not files:
            sys.exit("keine Dateien fuer series_id %r" % a.series)
    else:
        files = [(i + 1, lab, path) for i, (lab, path) in enumerate(a.scan)]
    a.tag = a.tag or a.series or "chipzoom"

    runs = []
    for run, label, path in files:
        r = pv.reflectogram(path, a)
        r["run"], r["voa"] = run, pv.voa_label(r["meta"])
        r["label"] = label
        r["facet"] = facet_of(r, a.facet_search_mm[0] * 1e-3,
                              a.facet_search_mm[1] * 1e-3)
        r["da"] = (r["z"] - r["facet"]) * 1e3          # mm hinter Facette A

        db = 20 * np.log10(r["R"] / r["R"].max() + 1e-15)
        fw = (r["da"] >= a.floor_mm[0]) & (r["da"] <= a.floor_mm[1])
        r["nf"] = pra.noise_floor(db[fw])
        r["snr"] = db - r["nf"]

        i = int(np.argmin(np.abs(r["da"])))
        r["snr_facet"] = float(r["snr"][i])
        chip_win = (NEARFIELD_MM, a.chip_mm[1])
        mirror_win = (-a.chip_mm[1], -NEARFIELD_MM)
        r["n_chip"], r["max_chip"], r["med_chip"] = count_visible(
            r["da"], r["snr"], chip_win, a.snr_thresh)
        r["n_mir"], r["max_mir"], r["med_mir"] = count_visible(
            r["da"], r["snr"], mirror_win, a.snr_thresh)
        r["n_ctrl"], r["max_ctrl"], r["med_ctrl"] = count_visible(
            r["da"], r["snr"], tuple(a.control_mm), a.snr_thresh)
        r["n_gds"], r["max_gds"], _ = count_visible(
            r["da"], r["snr"], GDS_WINDOW_MM, a.snr_thresh)
        r["n_gdsm"], r["max_gdsm"], _ = count_visible(
            r["da"], r["snr"], (-GDS_WINDOW_MM[1], -GDS_WINDOW_MM[0]),
            a.snr_thresh)
        runs.append(r)
        print("  %s verarbeitet" % (label or "Lauf %d" % run), flush=True)

    plot(runs, a, os.path.join(out_dir, "%s_chipzoom.png" % a.tag))

    # Die Kontrollen sind pro Millimeter zu rechnen, die Fenster sind
    # unterschiedlich breit. Chip- und Spiegelfenster sind gleich breit.
    w_chip = a.chip_mm[1] - NEARFIELD_MM
    w_ctrl = a.control_mm[1] - a.control_mm[0]
    csv = os.path.join(out_dir, "%s_chipzoom.csv" % a.tag)
    with open(csv, "w", newline="") as f:
        f.write("run,voa_setting,facet_mm,noise_floor_dB_rel_own_max,"
                "snr_facet_dB,n_chip,max_chip_dB,median_chip_dB,n_mirror,"
                "max_mirror_dB,median_mirror_dB,excess_max_dB,"
                "excess_median_dB,n_control,peaks_per_mm_control,"
                "max_control_dB,max_gds_dB,max_gds_mirror_dB\n")
        print("\n| Lauf | VOA | SNR Facette | Chip n / max / Sockel | "
              "Spiegel n / max / Sockel | Ueberschuss max / Sockel | "
              "Faser je mm | GDS / Spiegel |")
        print("|" + "---|" * 8)
        for r in runs:
            exc = r["max_chip"] - r["max_mir"]
            excm = r["med_chip"] - r["med_mir"]
            f.write("%d,%s,%.4f,%.2f,%.2f,%d,%.2f,%.2f,%d,%.2f,%.2f,%.2f,"
                    "%.2f,%d,%.4f,%.2f,%.2f,%.2f\n"
                    % (r["run"], (r.get("label") or "").replace(",", ";") or r["voa"],
                       r["facet"] * 1e3, r["nf"],
                       r["snr_facet"], r["n_chip"], r["max_chip"],
                       r["med_chip"], r["n_mir"], r["max_mir"], r["med_mir"],
                       exc, excm, r["n_ctrl"], r["n_ctrl"] / w_ctrl,
                       r["max_ctrl"], r["max_gds"], r["max_gdsm"]))
            print("| %d | %s | %.1f | %d / %.1f / %.1f | %d / %.1f / %.1f | "
                  "**%+.1f / %+.1f** | %.2f | %.1f / %.1f |"
                  % (r["run"], r.get("label") or r["voa"], r["snr_facet"],
                     r["n_chip"],
                     r["max_chip"], r["med_chip"], r["n_mir"], r["max_mir"],
                     r["med_mir"], exc, excm, r["n_ctrl"] / w_ctrl,
                     r["max_gds"], r["max_gdsm"]))
    print("\nChipfenster dA %.1f..%.0f mm, Spiegelfenster dA %.0f..%.1f mm "
          "(gleich breit, %.1f mm)" % (NEARFIELD_MM, a.chip_mm[1],
                                       -a.chip_mm[1], -NEARFIELD_MM, w_chip))
    print("wrote: %s" % csv)


if __name__ == "__main__":
    main()
