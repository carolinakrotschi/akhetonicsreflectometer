#!/usr/bin/env python3
"""VOA-Serie: alle Laeufe einer Messreihe untereinander als Reflektogramm.

Gedacht fuer die Frage "macht es einen Unterschied, wenn ich den VOA
verdrehe?". Dafuer reicht ein Stapel Reflektogramme nicht -- sie muessen
auch VERGLEICHBAR sein. Zwei Dinge macht dieses Skript deshalb anders als
`plot_reflectogram_only.py`:

1. **Gemeinsame Normierung.** Jede Spur wird auf DIESELBE Referenz bezogen
   (den staerksten Peak der ganzen Serie), nicht auf ihr eigenes Maximum.
   Sonst steht jeder Peak per Definition auf 0 dB und genau der Unterschied,
   nach dem gefragt wird, ist wegnormiert. Damit unterschiedlich lange
   Aufnahmen vergleichbar bleiben, wird die FFT vorher durch die
   Fenstersumme geteilt.  `--relative` gibt die alte Eigen-Normierung.

2. **Crop aus den Metadaten.** Am Ende jeder Aufnahme steht ein geparkter
   Laser (Fringes hoeren auf, Licht bleibt). Das sind hier ~4.7 % des
   Records -- mehr als `--trim` wegnimmt -- und es wuerde die
   Aux-Phasen-Resampling-Achse verderben. `meta["crop"]` schneidet es weg.

tau_aux kommt aus `meta["wavelength_axis_aux"]["tau_aux_implied_s"]` der
jeweiligen Aufnahme (jede zaehlt ihre eigenen Fringes), `--tau-aux-ns`
ueberschreibt das fuer alle.

Markiert werden die `--n-peaks` staerksten Peaks je Spur, die zugleich
ueber `--peak-floor-db` liegen. Ein Prominenz- oder Rauschboden-Kriterium
taugt hier nicht: der Multipath-Wald steht 10-15 dB ueber dem Boden und
liefert ~1900 "Peaks" je Spur.

Aufruf
    python tools/plot_voa_series.py --series 2680_20260917-095757_fiber48_fullnmrange
    python tools/plot_voa_series.py --series <id> --out-dir results/2026-09-17
"""

import argparse
import datetime
import json
import os
import re
import sys

import numpy as np
from scipy.signal import find_peaks
from scipy.signal.windows import blackmanharris, hann, kaiser

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import process_reflectogram_aux as pra  # noqa: E402

ROOT = os.path.dirname(_HERE)
MIN_PEAK_SPACING_M = 200e-6      # wie in process_reflectogram_aux.py


# ------------------------------------------------------------- laden
def series_files(raw_dir, series_id):
    """[(run_index, pfad), ...], nach run_index sortiert."""
    pat = re.compile(re.escape(series_id) + r"_(\d+)\.npz$")
    out = []
    for name in os.listdir(raw_dir):
        m = pat.match(name)
        if m:
            out.append((int(m.group(1)), os.path.join(raw_dir, name)))
    return sorted(out)


def reflectogram(path, a):
    """Ein Lauf -> dict mit z [m], linearer Amplitude (fenstersummen-
    normiert), Metadaten und Diagnose. Die dB-Skala kommt erst spaeter,
    wenn die gemeinsame Referenz ueber alle Laeufe bekannt ist."""
    d = np.load(path, allow_pickle=True)
    meta = json.loads(str(d["meta"])) if "meta" in d.files else {}

    crop = meta.get("crop")
    sl = slice(int(crop[0]), int(crop[1])) if crop else slice(None)
    ch = {i: d["ch%d" % i][sl] for i in (1, 2, 3, 4)}

    if a.tau_aux_ns is not None:
        tau_aux = a.tau_aux_ns * 1e-9
    else:
        tau_aux = meta["wavelength_axis_aux"]["tau_aux_implied_s"]

    aux, _ = pra.balanced(ch[a.aux_a], ch[a.aux_b])
    meas, _ = pra.balanced(ch[a.meas_a], ch[a.meas_b])

    y, dnu, span_nu, diag = pra.resample_on_aux(meas, aux, tau_aux, a.trim)
    m = len(y)
    t = np.linspace(-1, 1, m)
    y = y - np.polyval(np.polyfit(t, y, 5), t)

    win = {"hann": hann(m), "blackmanharris": blackmanharris(m),
           "kaiser": kaiser(m, a.kaiser_beta)}[a.window]
    R = np.abs(np.fft.rfft(y * win)) / win.sum()   # /sum: laengenunabhaengig
    z = np.arange(len(R)) * pra.C / (2 * pra.NG * dnu * m)

    return dict(
        path=path, meta=meta, z=z, R=R, tau_aux=tau_aux,
        dz_bin=pra.C / (2 * pra.NG * span_nu),
        z_nyq=pra.C / (4 * pra.NG * dnu),
        fringes=diag["fringes"], amp_min=diag["amp_min"],
        pts_per_fringe=diag["pts_per_fringe"],
        p_meas=float(np.mean(ch[a.meas_a]) + np.mean(ch[a.meas_b])),
        p_aux=float(np.mean(ch[a.aux_a]) + np.mean(ch[a.aux_b])),
    )


def voa_label(meta):
    """meta['voa_setting'] ist in Lauf 10 dieser Serie versehentlich mit der
    Kommandozeile ueberschrieben -- nur kurze Werte sind eine Stellung."""
    v = str(meta.get("voa_setting", "")).strip()
    return v if 0 < len(v) <= 12 and " " not in v else None


# ------------------------------------------------------------- peaks
def pick_peaks(z, db, n, floor_db, win):
    """Die n staerksten Peaks ueber floor_db im Fenster win = (lo, hi) [m]."""
    dz = z[1] - z[0]
    pk, _ = find_peaks(db, height=floor_db,
                       distance=max(3, int(MIN_PEAK_SPACING_M / dz)))
    pk = pk[(z[pk] >= win[0]) & (z[pk] <= win[1])]
    if len(pk) > n:
        pk = pk[np.argsort(db[pk])[::-1][:n]]
    return np.sort(pk)


def label_set(z, db, pk, n, span):
    """Aus pk die n staerksten auswaehlen, die untereinander mindestens 6 %
    der Fensterbreite auseinanderliegen -- sonst schreiben sich die
    Beschriftungen gegenseitig zu."""
    out = []
    for j in pk[np.argsort(db[pk])[::-1]]:
        if len(out) >= n:
            break
        if all(abs(z[j] - z[o]) > 0.06 * span for o in out):
            out.append(j)
    return out


# -------------------------------------------------------------- plot
def axis_unit(win):
    """mm, sobald das Fenster schmaler als 50 mm ist -- darueber m."""
    return ("mm", 1e3) if (win[1] - win[0]) < 0.05 else ("m", 1.0)


def stack_plot(runs, ref, win, a, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    nrun = len(runs)
    unit, scale = axis_unit(win)
    span = win[1] - win[0]

    fig, axes = plt.subplots(nrun, 1, figsize=(13, 1.5 * nrun + 1.1),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)

    inside = [(r["z"] >= win[0]) & (r["z"] <= win[1]) for r in runs]
    lo = max(min(r["db"][m].min() for r, m in zip(runs, inside)) - 3, a.ymin)
    hi = max(r["db"][m].max() for r, m in zip(runs, inside))
    # Kopfraum fuer die Beschriftungen, in Bruchteilen der Panelhoehe
    top = hi + 0.30 * (hi - lo)

    for k, (ax, r, m) in enumerate(zip(axes, runs, inside)):
        col = cmap(0.08 + 0.82 * k / max(nrun - 1, 1))
        ax.plot(r["z"][m] * scale, r["db"][m], lw=0.55, color=col)
        ax.axhline(r["nf"], color="darkviolet", ls="--", lw=0.9, zorder=4)

        # Schwelle je Spur: bei gemeinsamer Normierung stehen die
        # gedaempften Laeufe 50 dB tiefer, eine absolute Schwelle liesse
        # sie ohne einen einzigen Marker. Relativ zum eigenen Maximum
        # bekommt jede Spur ihre eigenen staerksten Peaks markiert.
        floor = (a.peak_floor_db if a.peak_floor_db is not None
                 else max(r["db"][m].max() - a.peak_span_db, r["nf"] + 6.0))
        pk = pick_peaks(r["z"], r["db"], a.n_peaks, floor, win)
        if len(pk):
            ax.plot(r["z"][pk] * scale, r["db"][pk], "v", ms=4.5,
                    color="crimson", zorder=5)
            for j in label_set(r["z"], r["db"], pk, a.n_labels, span):
                # an den Fensterraendern nach innen ausrichten, sonst steht
                # die Beschriftung halb ausserhalb des Panels
                f = (r["z"][j] - win[0]) / span
                ha = "left" if f < 0.08 else "right" if f > 0.92 else "center"
                dx = 4 if ha == "left" else -4 if ha == "right" else 0
                ax.annotate("%.3f mm\n%.1f dB" % (r["z"][j] * 1e3, r["db"][j]),
                            (r["z"][j] * scale, r["db"][j]),
                            textcoords="offset points", xytext=(dx, 6),
                            ha=ha, va="bottom", fontsize=6.5,
                            color="crimson", annotation_clip=True,
                            clip_on=True)

        tag = "VOA %s" % r["voa"] if r["voa"] is not None else "VOA ?"
        ax.text(0.006, 0.94, "Lauf %d  |  %s" % (r["run"], tag),
                transform=ax.transAxes, ha="left", va="top", fontsize=9,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec="gray", alpha=0.9))
        ax.text(0.994, 0.94,
                "max %.1f dB   Boden %.1f dB   Dynamik %.1f dB"
                % (r["db"][m].max(), r["nf"], r["db"][m].max() - r["nf"]),
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                color="0.25",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                          alpha=0.75))
        ax.grid(alpha=0.28)
        ax.set_ylim(lo, top)
        ax.set_ylabel("dB", fontsize=8)
        ax.tick_params(labelsize=8)

    axes[-1].set_xlabel("Distanz (%s, einfacher Weg / Reflexionskonvention)"
                        % unit)
    axes[0].set_xlim(win[0] * scale, win[1] * scale)

    r0 = runs[0]
    norm = ("jede Spur auf ihr EIGENES Maximum" if a.relative else
            "alle Spuren auf dieselbe Referenz (staerkster Peak der Serie, "
            "Lauf %d)" % ref["run"])
    fig.suptitle(
        "%s  --  %d Laeufe untereinander, VOA verdreht   "
        "[Fenster %.3f .. %.3f m]\n"
        "aux-referenziert, %s, dz %.1f um, Nyquist %.2f m  |  0 dB = %s"
        % (a.series, nrun, win[0], win[1], a.window,
           r0["dz_bin"] * 1e6, r0["z_nyq"], norm),
        fontsize=10.5)
    fig.subplots_adjust(top=1 - 0.62 / (1.5 * nrun + 1.1), bottom=0.05,
                        left=0.055, right=0.99, hspace=0.12)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote: %s" % path)


def overlay_plot(runs, win, a, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("viridis")
    unit, scale = axis_unit(win)
    fig, ax = plt.subplots(figsize=(13, 6))
    floors, tops = [], []
    for k, r in enumerate(runs):
        m = (r["z"] >= win[0]) & (r["z"] <= win[1])
        tag = "VOA %s" % r["voa"] if r["voa"] is not None else "VOA ?"
        ax.plot(r["z"][m] * scale, r["db"][m], lw=0.7,
                color=cmap(0.08 + 0.82 * k / max(len(runs) - 1, 1)),
                label="Lauf %d (%s)" % (r["run"], tag))
        floors.append(pra.noise_floor(r["db"][m]))
        tops.append(r["db"][m].max())
    # Die Interferenznullen reichen 40 dB unter den Boden und wuerden die
    # Achse so stauchen, dass die Bodentreppe -- der eigentliche Inhalt --
    # nicht mehr zu sehen ist.
    ax.set_ylim(min(floors) - 15, max(tops) + 12)
    ax.set_xlabel("Distanz (%s, einfacher Weg / Reflexionskonvention)" % unit)
    ax.set_ylabel("Amplitude (dB, gemeinsame Referenz)")
    ax.set_title("%s -- alle %d Laeufe uebereinander  [Fenster %.3f .. %.3f m]"
                 % (a.series, len(runs), win[0], win[1]))
    ax.grid(alpha=0.3)
    ax.set_xlim(win[0] * scale, win[1] * scale)
    ax.legend(fontsize=7.5, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote: %s" % path)


# -------------------------------------------------------------- main
def build_argparser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--series", required=True,
                   help="series_id, z.B. 2680_20260917-095757_fiber48_fullnmrange")
    p.add_argument("--raw-dir", default=os.path.join(ROOT, "raw_data"))
    p.add_argument("--out-dir", default=None, help="default results/<heute>")
    p.add_argument("--zoom", nargs="*", default=None, metavar="LO:HI",
                   help="Fenster in mm, z.B. 1500:1700 . Default: ganze "
                        "Achse, dann zweimal um den staerksten Peak der "
                        "Serie herum (+-100 mm, +-5 mm), plus das "
                        "Nahfeld 0:10 mm")
    p.add_argument("--n-peaks", type=int, default=10,
                   help="markierte Peaks je Spur (staerkste zuerst)")
    p.add_argument("--n-labels", type=int, default=4,
                   help="davon beschriftet")
    p.add_argument("--peak-span-db", type=float, default=40.0,
                   help="Marker nur fuer Peaks, die weniger als so viel dB "
                        "unter dem eigenen Maximum der Spur liegen")
    p.add_argument("--peak-floor-db", type=float, default=None,
                   help="absolute Schwelle statt --peak-span-db")
    p.add_argument("--ymin", type=float, default=-110.0)
    p.add_argument("--relative", action="store_true",
                   help="jede Spur auf ihr eigenes Maximum normieren "
                        "(versteckt genau den gesuchten Unterschied)")
    p.add_argument("--window", default="kaiser",
                   choices=["hann", "blackmanharris", "kaiser"])
    p.add_argument("--kaiser-beta", type=float, default=12.0)
    p.add_argument("--trim", type=float, default=0.01)
    p.add_argument("--tau-aux-ns", type=float, default=None,
                   help="ueberschreibt tau_aux aus den Metadaten")
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

    files = series_files(a.raw_dir, a.series)
    if not files:
        sys.exit("keine Dateien fuer series_id %r in %s" % (a.series, a.raw_dir))
    print("%d Laeufe: %s\n" % (len(files), a.series))

    runs = []
    for run, path in files:
        r = reflectogram(path, a)
        r["run"] = run
        r["voa"] = voa_label(r["meta"])
        runs.append(r)
        print("  Lauf %2d  VOA %4s  tau_aux %.4f ns  %9.0f Fringes  "
              "Aux-Kontrast min %3.0f %%  P_meas %.4f mW"
              % (run, str(r["voa"]), r["tau_aux"] * 1e9, r["fringes"],
                 r["amp_min"] * 100, r["p_meas"]))

    # gemeinsame Referenz -- der staerkste Peak der ganzen Serie
    ref = max(runs, key=lambda r: r["R"].max())
    for r in runs:
        base = r["R"].max() if a.relative else ref["R"].max()
        r["db"] = 20 * np.log10(r["R"] / base + 1e-15)

    z_end = min(r["z"][-1] for r in runs)
    if a.zoom:
        windows = []
        for s in a.zoom:
            lo, hi = (float(v) * 1e-3 for v in s.split(":"))
            windows.append((lo, hi))
    else:
        # Die Zoomleiter haengt an dem, was die Serie zeigt: dem staerksten
        # Peak ueberhaupt. Ein fester Ausschnitt wie "die ersten 20 mm"
        # traefe ihn bei dieser Serie um 1.5 m daneben.
        zp = ref["z"][int(np.argmax(ref["R"]))]
        windows = [(0.0, z_end),
                   (max(zp - 0.100, 0.0), min(zp + 0.100, z_end)),
                   (max(zp - 0.005, 0.0), min(zp + 0.005, z_end)),
                   (0.0, 0.010)]
        print("\nstaerkster Peak der Serie bei %.3f mm (Lauf %d) "
              "-- Zoomleiter richtet sich danach"
              % (zp * 1e3, ref["run"]))

    for win in windows:
        for r in runs:
            m = (r["z"] >= win[0]) & (r["z"] <= win[1])
            r["nf"] = pra.noise_floor(r["db"][m])
        tag = "%.0f-%.0fmm" % (win[0] * 1e3, win[1] * 1e3)
        stack_plot(runs, ref, win, a,
                   os.path.join(out_dir, "%s_stack_%s.png" % (a.series, tag)))
    overlay_plot(runs, windows[min(1, len(windows) - 1)], a,
                 os.path.join(out_dir, "%s_overlay.png" % a.series))

    # Kennzahlen je Lauf, ueber die ganze Achse
    zfull = z_end
    csv = os.path.join(out_dir, "%s_voa_summary.csv" % a.series)
    with open(csv, "w", newline="") as f:
        f.write("run,voa_setting,tau_aux_ns,fringes,aux_contrast_min,"
                "P_meas_mW,P_aux_mW,max_dB,z_max_mm,noise_floor_dB,"
                "dynamic_range_dB\n")
        print("\n  Lauf  VOA   max dB    bei mm    Boden dB   Dynamik dB")
        for r in runs:
            keep = r["z"] <= zfull
            nf = pra.noise_floor(r["db"][keep])
            j = int(np.argmax(r["db"][keep]))
            f.write("%d,%s,%.4f,%.0f,%.4f,%.6f,%.6f,%.2f,%.4f,%.2f,%.2f\n"
                    % (r["run"], r["voa"], r["tau_aux"] * 1e9, r["fringes"],
                       r["amp_min"], r["p_meas"], r["p_aux"],
                       r["db"][keep].max(), r["z"][j] * 1e3, nf,
                       r["db"][keep].max() - nf))
            print("  %4d  %3s  %7.2f  %8.3f  %10.2f  %11.2f"
                  % (r["run"], str(r["voa"]), r["db"][keep].max(),
                     r["z"][j] * 1e3, nf, r["db"][keep].max() - nf))
    print("\nwrote: %s" % csv)


if __name__ == "__main__":
    main()
