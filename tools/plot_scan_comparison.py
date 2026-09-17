#!/usr/bin/env python3
"""Zwei oder mehr Scans aus VERSCHIEDENEN Aufbauzustaenden vergleichen --
z.B. "vorher 10-dB-Splitter" gegen "nachher 20-dB-Splitter".

Warum nicht `plot_voa_series.py`: das haengt an einer `series_id` und
beschriftet die Spuren mit Laufnummer und VOA-Stellung. Hier kommen die
Dateien einzeln und mit freier Beschriftung herein, duerfen aus
verschiedenen Aufnahmen stammen und trotzdem NICHT jede auf ihr eigenes
Maximum normiert werden -- sonst steht der staerkste Peak jeder Spur per
Definition auf 0 dB und der gesuchte Unterschied ist weg. Die Verarbeitung
(Crop, aux-Resampling, fenstersummen-normierte FFT) ist dieselbe und wird
aus `plot_voa_series.py` importiert.

Drei Dinge macht dieses Skript zusaetzlich:

1. **Zoomleiter auf absoluter mm-Achse.** Vier Panels: ganze Achse, dann
   geometrisch bis in den Bauteilbereich hinein. Geschrieben wird die
   absolute Distanz vom Instrument, nicht "hinter Facette A" -- die Facette
   liegt bei zwei Aufnahmen nicht an derselben Stelle, und wo der Chip auf
   der Achse steht, ist beim Aufbauvergleich mit die Frage.

2. **Facette je Spur einzeln gesucht** (`--facet-window`), und die
   Bauteilmarken wandern je Spur um deren eigene Facette mit.

3. **Fenstertest** (`--marks-csv`): Maximum im Vorhersagefenster eines
   Bauteils gegen gleich breite Kontrollfenster im selben Abstandsbereich.
   Das ist die Zahl, die "sehe ich die Struktur?" beantwortet -- ein Peak,
   den auch jedes zweite Kontrollfenster erreicht, ist kein Bauteil,
   sondern der Multipath-Wald.

Aufruf
    python tools/plot_scan_comparison.py \\
        --scan "vorher (10 dB)"  raw_data/..._fullnmrange_10.npz \\
        --scan "nachher (20 dB)" raw_data/2026-09-17-10-33_..._20dbsplitter.npz \\
        --marks-csv results/2026-09-16/ligentec_kanal48_struktur_bridged.csv \\
        --tag fiber48_splitter --out-dir results/2026-09-17
"""

import argparse
import datetime
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import process_reflectogram_aux as pra          # noqa: E402
import plot_voa_series as pvs                   # noqa: E402
from compare_span import load_marks, short_name  # noqa: E402

ROOT = os.path.dirname(_HERE)
N_FIBER = 1.468
COLORS = ("#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#8c564b", "#17becf")


# ------------------------------------------------------- Facette, Marken
def find_facet(z, db, lo_mm, hi_mm):
    """Staerkster Peak im Fenster -> (Position mm, Pegel dB)."""
    m = (z * 1e3 >= lo_mm) & (z * 1e3 <= hi_mm)
    if not m.any():
        sys.exit("Facettenfenster %g..%g mm liegt ausserhalb der Achse"
                 % (lo_mm, hi_mm))
    j = int(np.flatnonzero(m)[np.argmax(db[m])])
    return z[j] * 1e3, db[j]


def mark_windows(marks, ng):
    """Bauteil -> (Name, lo_mm, hi_mm) hinter Facette A, fuer n_g = ng[0..1].

    Der on-chip-Weg wird auf die faseraequivalente z-Achse umgerechnet;
    die n_g-Unsicherheit des Chips (1.80-2.00) macht die Fensterbreite."""
    out = []
    for name, s_from, s_to in marks:
        lo = s_from * ng[0] / N_FIBER / 1e3
        hi = s_to * ng[1] / N_FIBER / 1e3
        out.append((short_name(name), lo, hi))
    return out


# ------------------------------------------------------------ Fenstertest
def window_test(z, db, facet_mm, win, all_wins, ctrl_lo, ctrl_hi):
    """Ueberschuss ueber den Boden im Fenster `win`, und wie viele gleich
    breite Kontrollfenster ihn erreichen.

    Die Kontrollfenster kacheln `ctrl_lo..ctrl_hi` mm hinter der Facette
    und lassen alle Bauteilfenster aus -- sonst pruefte man das Bauteil
    gegen sich selbst. Der Boden wird ueber denselben Kontrollbereich
    genommen, nicht ueber die ganze Achse: der Pigtail und der Anstieg ab
    2.6 m haetten dort sonst mit im Mittel gestanden."""
    zm = z * 1e3 - facet_mm
    cm = (zm >= ctrl_lo) & (zm <= ctrl_hi)
    nf = pra.noise_floor(db[cm])

    def excess(lo, hi):
        m = (zm >= lo) & (zm <= hi)
        return (db[m].max() - nf) if m.any() else None

    val = excess(*win)
    if val is None:
        return None
    w = win[1] - win[0]
    hits = n = 0
    x = ctrl_lo
    while x + w <= ctrl_hi:
        if not any(x < b and a < x + w for _, a, b in all_wins):
            e = excess(x, x + w)
            if e is not None:
                n += 1
                hits += e >= val
        x += w
    return dict(excess=val, hits=hits, n=n, floor=nf,
                frac=(hits / n if n else float("nan")))


# ------------------------------------------------------------------- Plot
def zoom_ladder(traces, comp_lo, comp_hi, z_end, n_panels=4):
    """Fenster in mm (absolut): ganze Achse, dann geometrisch bis in den
    Bauteilbereich. Linear gemittelt zeigte Panel 2 und 3 fast dasselbe
    Bild -- von 3000 mm auf 30 mm ist der halbe Weg geometrisch 300 mm."""
    full = (0.0, z_end * 1e3)
    tgt = (comp_lo, comp_hi)
    wins = [full]
    w0, w1 = full[1] - full[0], tgt[1] - tgt[0]
    c = 0.5 * (tgt[0] + tgt[1])
    for k in range(1, n_panels - 1):
        w = w0 * (w1 / w0) ** (k / (n_panels - 1.0))
        lo = max(c - w / 2, 0.0)
        wins.append((lo, min(lo + w, z_end * 1e3)))
    wins.append(tgt)
    return wins


def stack_plot(traces, wins, comp, path, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nw = len(wins)
    fig, axes = plt.subplots(nw, 1, figsize=(13.5, 2.6 * nw + 1.4))
    axes = np.atleast_1d(axes)

    for k, (ax, win) in enumerate(zip(axes, wins)):
        span = win[1] - win[0]
        detail = comp and span < 3.0 * (comp[-1][2] - comp[0][1])
        for t, col in zip(traces, COLORS):
            zm = t["z"] * 1e3
            m = (zm >= win[0]) & (zm <= win[1])
            if not m.any():
                continue
            ax.plot(zm[m], t["db"][m], lw=0.6, color=col, alpha=0.9,
                    label="%s   Facette %.3f mm, Boden %.1f dB"
                          % (t["label"], t["facet_mm"], t["nf_ctrl"]))
            if comp:
                if detail:
                    for name, lo, hi in comp:
                        a, b = t["facet_mm"] + lo, t["facet_mm"] + hi
                        if b < win[0] or a > win[1]:
                            continue
                        ax.axvspan(a, b, color=col, alpha=0.10, lw=0)
                else:
                    a = t["facet_mm"] + comp[0][1]
                    b = t["facet_mm"] + comp[-1][2]
                    ax.axvspan(a, b, color=col, alpha=0.12, lw=0)
                ax.axvline(t["facet_mm"], color=col, ls=":", lw=1.0)

        if comp and detail:
            # Namen nur einmal, an der Facette der ersten Spur
            t0 = traces[0]
            ylim = ax.get_ylim()
            for i, (name, lo, hi) in enumerate(comp):
                a, b = t0["facet_mm"] + lo, t0["facet_mm"] + hi
                if b < win[0] or a > win[1]:
                    continue
                ax.text(0.5 * (a + b), ylim[1] - 0.04 * (ylim[1] - ylim[0])
                        - 0.07 * (ylim[1] - ylim[0]) * (i % 3),
                        name, ha="center", va="top", fontsize=7,
                        rotation=0, color="0.25",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white",
                                  ec="none", alpha=0.75))

        ax.set_xlim(*win)
        ax.grid(alpha=0.28)
        ax.set_ylabel("dB (gemeinsame Referenz)", fontsize=8)
        ax.tick_params(labelsize=8)
        ax.text(0.004, 0.96, "Panel %d  |  Fenster %.3f .. %.3f mm  "
                "(Breite %.3f mm)" % (k + 1, win[0], win[1], span),
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray",
                          alpha=0.9))
        if k == 0:
            ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
            ax.text(0.004, 0.06,
                    "x = absolute Distanz vom Instrument. Abstaende sind auf "
                    "< 1 um genau, absolute Lagen wandern zwischen "
                    "Messsitzungen um ~1 mm (interne Referenz).",
                    transform=ax.transAxes, fontsize=7, color="0.35")

    axes[-1].set_xlabel("Distanz vom Instrument "
                        "(mm, einfacher Weg / Reflexionskonvention)")
    fig.suptitle(title, fontsize=10.5)
    fig.subplots_adjust(top=1 - 0.75 / (2.6 * nw + 1.4), bottom=0.05,
                        left=0.06, right=0.99, hspace=0.22)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote: %s" % path)


# ------------------------------------------------------------------- main
def build_argparser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scan", nargs=2, action="append", required=True,
                   metavar=("LABEL", "PFAD"),
                   help="Beschriftung und .npz, mehrfach angebbar. Rohe "
                        "JSON-Scans vorher durch tools/json_to_npz.py")
    p.add_argument("--out-dir", default=None, help="default results/<heute>")
    p.add_argument("--tag", default="scan_vergleich",
                   help="Praefix der Ausgabedateien")
    p.add_argument("--facet-window", nargs=2, type=float, default=[1500, 1700],
                   metavar="MM", help="Suchfenster fuer Facette A je Spur")
    p.add_argument("--marks-csv", default=None,
                   help="Bauteile aus gds_trace_ligentec.py")
    p.add_argument("--ng", nargs=2, type=float, default=[1.80, 2.00],
                   help="Gruppenindex des Chips (Fensterbreite)")
    p.add_argument("--ctrl-range", nargs=2, type=float, default=[1.0, 40.0],
                   metavar="MM", help="Kontrollbereich hinter Facette A")
    p.add_argument("--zoom", nargs="*", default=None, metavar="LO:HI",
                   help="Fenster in mm absolut statt der Zoomleiter")
    p.add_argument("--relative", action="store_true",
                   help="jede Spur auf ihr eigenes Maximum (versteckt genau "
                        "den gesuchten Unterschied)")
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


def main():
    a = build_argparser().parse_args()
    out_dir = a.out_dir or os.path.join(ROOT, "results",
                                        datetime.date.today().isoformat())
    os.makedirs(out_dir, exist_ok=True)

    traces = []
    for label, path in a.scan:
        t = pvs.reflectogram(path, a)
        t["label"] = label
        traces.append(t)
        print("  %-26s tau_aux %.4f ns  %9.0f Fringes  "
              "Aux-Kontrast min %3.0f %%  P_meas %.4f mW  P_aux %.4f mW"
              % (label, t["tau_aux"] * 1e9, t["fringes"],
                 t["amp_min"] * 100, t["p_meas"], t["p_aux"]))

    # gemeinsame Referenz: staerkster Peak ueber alle Spuren
    ref = max(traces, key=lambda t: t["R"].max())
    for t in traces:
        base = t["R"].max() if a.relative else ref["R"].max()
        t["db"] = 20 * np.log10(t["R"] / base + 1e-15)

    for t in traces:
        t["facet_mm"], t["facet_db"] = find_facet(
            t["z"], t["db"], a.facet_window[0], a.facet_window[1])
        zm = t["z"] * 1e3 - t["facet_mm"]
        cm = (zm >= a.ctrl_range[0]) & (zm <= a.ctrl_range[1])
        t["nf_ctrl"] = pra.noise_floor(t["db"][cm])
        t["nf_full"] = pra.noise_floor(t["db"])

    comp = []
    if a.marks_csv:
        comp = mark_windows(load_marks(a.marks_csv), a.ng)
        comp.sort(key=lambda c: c[1])

    z_end = min(t["z"][-1] for t in traces)
    if a.zoom:
        wins = [tuple(float(v) for v in s.split(":")) for s in a.zoom]
    elif comp:
        # Der Edge Coupler sitzt 0.5 mm hinter der Facette und ist von ihr
        # nicht zu trennen -- er wuerde den Zoom bis an die Facette
        # aufspannen und den Rest zusammenquetschen.
        inner = [c for c in comp if c[1] > 1.0] or comp
        f0 = traces[0]["facet_mm"]
        pad = 0.15 * (inner[-1][2] - inner[0][1]) + 0.2
        wins = zoom_ladder(traces, f0 + inner[0][1] - pad,
                           f0 + inner[-1][2] + pad, z_end)
    else:
        f0 = traces[0]["facet_mm"]
        wins = zoom_ladder(traces, f0 - 5, f0 + 5, z_end)

    norm = ("jede Spur auf ihr EIGENES Maximum" if a.relative else
            "alle Spuren auf dieselbe Referenz (staerkster Peak von %r)"
            % ref["label"])
    title = ("%s -- %d Scans, absolute Achse\naux-referenziert, %s, "
             "dz %.1f um, Nyquist %.2f m  |  0 dB = %s"
             % (a.tag, len(traces), a.window, traces[0]["dz_bin"] * 1e6,
                traces[0]["z_nyq"], norm))
    stack_plot(traces, wins, comp,
               os.path.join(out_dir, "%s_zoom4.png" % a.tag), title)

    # ------------------------------------------------- Kennzahlen
    csv = os.path.join(out_dir, "%s_kennzahlen.csv" % a.tag)
    with open(csv, "w", newline="") as f:
        f.write("label,tau_aux_ns,fringes,aux_contrast_min,P_meas_mW,"
                "P_aux_mW,facet_mm,facet_dB,floor_ctrl_dB,floor_full_dB,"
                "dynamic_range_dB\n")
        for t in traces:
            f.write("%s,%.4f,%.0f,%.4f,%.6f,%.6f,%.4f,%.2f,%.2f,%.2f,%.2f\n"
                    % (t["label"], t["tau_aux"] * 1e9, t["fringes"],
                       t["amp_min"], t["p_meas"], t["p_aux"], t["facet_mm"],
                       t["facet_db"], t["nf_ctrl"], t["nf_full"],
                       t["facet_db"] - t["nf_ctrl"]))
    print("wrote: %s" % csv)

    print("\n| Groesse | %s |" % " | ".join(t["label"] for t in traces))
    print("|---|%s" % ("---|" * len(traces)))
    rows = [("P_aux [mW]", "%.4f", "p_aux"),
            ("P_meas [mW]", "%.4f", "p_meas"),
            ("Aux-Kontrast min [%]", "%.0f", None),
            ("Fringes", "%.0f", "fringes"),
            ("Facette A [mm]", "%.4f", "facet_mm"),
            ("Facettenpeak [dB]", "%.2f", "facet_db"),
            ("Boden 1-40 mm hinter Facette [dB]", "%.2f", "nf_ctrl"),
            ("Boden ganze Achse [dB]", "%.2f", "nf_full")]
    for name, fmt, key in rows:
        if key is None:
            vals = [fmt % (t["amp_min"] * 100) for t in traces]
        else:
            vals = [fmt % t[key] for t in traces]
        print("| %s | %s |" % (name, " | ".join(vals)))
    print("| **Dynamik Facette ueber Boden [dB]** | %s |"
          % " | ".join("**%.2f**" % (t["facet_db"] - t["nf_ctrl"])
                       for t in traces))

    # ------------------------------------------------- Fenstertest
    if comp:
        wcsv = os.path.join(out_dir, "%s_fenstertest.csv" % a.tag)
        with open(wcsv, "w", newline="") as f:
            f.write("label,bauteil,lo_mm,hi_mm,excess_dB,hits,n_ctrl,"
                    "frac_percent\n")
            print("\n| Bauteil | Fenster hinter Facette | %s |"
                  % " | ".join(t["label"] for t in traces))
            print("|---|---|%s" % ("---|" * len(traces)))
            for name, lo, hi in comp:
                if lo < 1.0:
                    continue          # Edge Coupler: nicht von der Facette trennbar
                cells = []
                for t in traces:
                    r = window_test(t["z"], t["db"], t["facet_mm"], (lo, hi),
                                    comp, a.ctrl_range[0], a.ctrl_range[1])
                    if r is None:
                        cells.append("-")
                        continue
                    cells.append("+%.1f dB, %.0f %%"
                                 % (r["excess"], 100 * r["frac"]))
                    f.write("%s,%s,%.4f,%.4f,%.2f,%d,%d,%.1f\n"
                            % (t["label"], name, lo, hi, r["excess"],
                               r["hits"], r["n"], 100 * r["frac"]))
                print("| %s | %.2f-%.2f mm | %s |"
                      % (name, lo, hi, " | ".join(cells)))
        print("\n(Prozent = Anteil gleich breiter Kontrollfenster in "
              "%.0f-%.0f mm hinter der Facette, die den Wert erreichen. "
              "Unter 5 %% waere eine Detektion.)"
              % (a.ctrl_range[0], a.ctrl_range[1]))
        print("wrote: %s" % wcsv)


if __name__ == "__main__":
    main()
