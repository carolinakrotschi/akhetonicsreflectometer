#!/usr/bin/env python3
"""Fibre length from OFDR runs with and without the fibre -- reflectogram,
peak positions, length per run, comparability check and a zoom-ladder plot.

Built 2026-10-01 for the first Mk2-LOG measurement (2 m fibre, open end).

Method
    Every raw scan goes through process_reflectogram_aux (scan_io loading,
    sweep crop, tau_aux from the sweep, balanced subtraction, aux resampling,
    Kaiser FFT), CSV + sidecar into --out-dir. Peaks are located by parabolic
    interpolation on the dB curve, so a position is not quantised to the
    16.6 um bin.

    length = z(fibre end) - z(connector)

    computed two ways:
      * against the connector reflex IN THE SAME RUN (it survives at about
        -20 dB with the fibre attached). This cancels the small common scale
        drift between runs (the aux delay wanders by ~1e-4 over half an hour:
        every peak moves in proportion to z), so it is the better number;
      * against the mean connector position of the runs without fibre.

    --n-end K looks for the K strongest sub-peaks inside the end window and
    reports a length for each (the 2026-10-01 fibre end is three peaks
    ~1.6 mm apart).

Scale: positions use n_g = 1.468 (process_reflectogram_aux.NG); --ng
rescales the lengths (L is proportional to 1/n_g). The aux delay is scaled
from the commanded sweep endpoints, which the EXFO hits to ~20 pm in 50 nm,
i.e. a common +-0.04 % on every length (+-0.8 mm on 2 m) that no repeat
reduces.

Example
    python tools/measure_fiber_length.py \\
        --without "raw_data/sweep_ofdr_mk2_nofiber_run*.npz" \\
        --with    "raw_data/sweep_ofdr_mk2_2m_run*.npz" \\
        --conn-window 2700 2720 --end-window 4800 4815 --n-end 3 \\
        --out-dir results/2026-10-01/mk2_reference --tag mk2_2m
"""

import argparse
import glob
import os
import sys

import numpy as np
from scipy.signal import find_peaks

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import process_reflectogram_aux as P  # noqa: E402
from scan_io import check_comparable, provenance_for  # noqa: E402


def reflectogram_csv(path, out_dir):
    """Raw scan -> its _reflectogram.csv in out_dir (computed if missing)."""
    if path.endswith("_reflectogram.csv"):
        return path
    stem = os.path.splitext(os.path.basename(path))[0]
    prefix = os.path.join(out_dir, stem)
    csv = prefix + "_reflectogram.csv"
    if not os.path.exists(csv):
        ns = P.build_argparser().parse_args([path, "--out", prefix])
        P.process(ns)
    return csv


def load_csv(csv):
    d = np.genfromtxt(csv, delimiter=",", names=True)
    return d["distance_m"] * 1e3, d["amplitude_dB"]          # mm, dB


def peaks_in(z, db, lo, hi, n=1, min_sep_mm=0.5):
    """The n strongest peaks in [lo, hi] mm, parabolic-interpolated,
    sorted by position: [(z_mm, dB), ...]."""
    m = np.flatnonzero((z > lo) & (z < hi))
    dz = z[1] - z[0]
    k, _ = find_peaks(db[m], distance=max(1, int(min_sep_mm / dz)))
    k = k[np.argsort(-db[m][k])][:n]
    out = []
    for j in sorted(k):
        i = m[j]
        a, b, c = db[i - 1], db[i], db[i + 1]
        den = a - 2 * b + c
        d = 0.5 * (a - c) / den if den else 0.0
        out.append((z[i] + d * dz, b - 0.25 * (a - c) * d))
    return out


def ladder(x_end, xmin, xmax, w_last):
    """Four windows: whole axis, then geometric zoom onto x_end."""
    span = xmax - xmin
    widths = span * (w_last / span) ** (np.arange(4) / 3.0)
    wins = [(xmin, xmax)]
    for w in widths[1:]:
        lo = max(xmin, min(x_end - w / 2, xmax - w))
        wins.append((lo, lo + w))
    return wins


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--without", required=True, help="glob: runs without fibre")
    ap.add_argument("--with", dest="with_", required=True, help="glob: runs with fibre")
    ap.add_argument("--conn-window", nargs=2, type=float, required=True,
                    metavar="MM", help="search window of the connector reflex")
    ap.add_argument("--end-window", nargs=2, type=float, required=True,
                    metavar="MM", help="search window of the fibre end")
    ap.add_argument("--n-end", type=int, default=1,
                    help="sub-peaks to resolve inside the end window")
    ap.add_argument("--ng", type=float, default=P.NG,
                    help="group index for the reported length (default %g)" % P.NG)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tag", default="fiber_length")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    scale = P.NG / a.ng

    runs = []
    for kind, pat in (("without", a.without), ("with", a.with_)):
        files = sorted(glob.glob(pat))
        if not files:
            sys.exit("no files for %r" % pat)
        for f in files:
            csv = reflectogram_csv(f, a.out_dir)
            z, db = load_csv(csv)
            r = dict(kind=kind, name=os.path.basename(f), csv=csv, z=z, db=db,
                     prov=provenance_for(csv))
            r["conn"] = peaks_in(z, db, *a.conn_window)[0]
            r["end"] = (peaks_in(z, db, *a.end_window, n=a.n_end)
                        if kind == "with" else [])
            runs.append(r)
    check_comparable([(r["name"], r["prov"]) for r in runs])

    ref = [r for r in runs if r["kind"] == "without"]
    fib = [r for r in runs if r["kind"] == "with"]
    zc_ref = float(np.mean([r["conn"][0] for r in ref]))

    lines = ["run,kind,connector_mm,connector_dB" + "".join(
        ",end%s_mm,end%s_dB,L%s_own_m,L%s_ref_m" % ((k,) * 4)
        for k in "ABCDEFG"[:a.n_end])]
    print("\n| run | connector (mm) | dB |" + "".join(
        " end %s (mm) | dB | L %s (m) |" % (k, k) for k in "ABCDEFG"[:a.n_end]))
    print("|---|---|---|" + "---|---|---|" * a.n_end)
    for r in runs:
        row = "| %s | %.3f | %.1f |" % (r["name"][:40], r["conn"][0], r["conn"][1])
        csvrow = "%s,%s,%.4f,%.2f" % (r["name"], r["kind"], r["conn"][0], r["conn"][1])
        for ze, de in r["end"]:
            l_own = (ze - r["conn"][0]) * 1e-3 * scale
            l_ref = (ze - zc_ref) * 1e-3 * scale
            row += " %.3f | %.1f | %.5f |" % (ze, de, l_own)
            csvrow += ",%.4f,%.2f,%.6f,%.6f" % (ze, de, l_own, l_ref)
        if not r["end"]:
            row += " - | | - |" * a.n_end
        print(row)
        lines.append(csvrow)
    csv_out = os.path.join(a.out_dir, "%s_lengths.csv" % a.tag)
    with open(csv_out, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\nn_g = %g  (lengths scale with 1/n_g)" % a.ng)
    print("connector without fibre: %.4f mm, spread %.3f mm over %d runs"
          % (zc_ref, np.ptp([r["conn"][0] for r in ref]), len(ref)))
    summary = []
    for k in range(a.n_end):
        own = np.array([(r["end"][k][0] - r["conn"][0]) for r in fib]) * 1e-3 * scale
        rf = np.array([(r["end"][k][0] - zc_ref) for r in fib]) * 1e-3 * scale
        sd = lambda x: x.std(ddof=1) * 1e3 if x.size > 1 else float("nan")  # noqa: E731
        summary.append((own.mean(), sd(own), rf.mean(), sd(rf)))
        print("end %s: L = %.5f m (sd %.3f mm) vs own connector, %.5f m "
              "(sd %.3f mm) vs reference runs"
              % ("ABCDEFG"[k], own.mean(), sd(own), rf.mean(), sd(rf)))
    print("systematic: +-0.04 %% scale from the sweep endpoints = +-%.1f mm"
          % (0.0004 * summary[0][0] * 1e3))
    print("wrote: %s" % csv_out)

    # ------------------------------------------------------------- plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # x = distance from the connector (reference runs), so the fibre end
    # reads as the fibre length; the top axis gives the absolute z.
    x0 = zc_ref
    ends = [float(np.mean([r["end"][j][0] for r in fib])) - x0
            for j in range(a.n_end)]
    xmin = -x0
    xmax = min(r["z"][-1] for r in runs) - x0
    wins = ladder(ends[0], xmin, xmax, max(a.end_window[1] - a.end_window[0], 5.0))
    fig, axs = plt.subplots(4, 1, figsize=(12, 13.5))
    cmap = plt.get_cmap("viridis")
    for k, (lo, hi) in enumerate(wins):
        ax = axs[k]
        for r in ref:
            m = (r["z"] - x0 >= lo) & (r["z"] - x0 <= hi)
            ax.plot(r["z"][m] - x0, r["db"][m], color="0.55", lw=0.6, ls="--",
                    label="without fibre" if r is ref[0] else None)
        for i, r in enumerate(fib):
            m = (r["z"] - x0 >= lo) & (r["z"] - x0 <= hi)
            ax.plot(r["z"][m] - x0, r["db"][m], color=cmap(i / max(1, len(fib) - 1)),
                    lw=0.6, alpha=0.85, label=r["name"][:32] if k == 0 else None)
        if k < 2:
            # The floor steps down where the glass ends: Rayleigh backscatter
            # of the fibre (verified 2026-10-01 -- the step moves with the end).
            ax.axvspan(0, ends[0], color="gold", alpha=0.12, lw=0,
                       label="fibre backscatter (with fibre)" if k == 0 else None)
        if a.n_end > 1 and hi - lo > 20 * (ends[-1] - ends[0]):
            # too wide to tell the sub-peaks apart: one marker for the cluster
            marks = [(ends[0], "fibre end (%s-%s)" % ("A", "ABCDEFG"[a.n_end - 1]))]
        else:
            marks = [(x, "end %s" % "ABCDEFG"[j]) for j, x in enumerate(ends)]
        for x, lab in [(0.0, "connector")] + marks:
            if lo <= x <= hi:
                ax.axvline(x, color="crimson", lw=0.7, ls=":")
                ax.annotate("%s\n%.2f mm" % (lab, x), (x, 2), fontsize=7,
                            ha="center", va="bottom", color="crimson")
        ax.set_xlim(lo, hi)
        ax.set_ylim(-90, 14)
        ax.set_ylabel("dB rel. strongest peak")
        ax.grid(alpha=0.3)
        top = ax.secondary_xaxis("top", functions=(lambda x: x + x0,
                                                   lambda x: x - x0))
        top.tick_params(labelsize=7)
        if k == 0:
            top.set_xlabel("absolute z (mm)", fontsize=8)
        ax.set_title("%.0f to %.0f mm from connector  (%.0f mm wide)"
                     % (lo, hi, hi - lo), fontsize=9)
    axs[-1].set_xlabel("distance from connector (mm, n_g = %.3f); connector "
                       "= %.2f mm absolute" % (P.NG, x0))
    if scale != 1.0:
        print("note: plot axis uses n_g %g, lengths in the table use %g"
              % (P.NG, a.ng))
    axs[0].legend(fontsize=7, ncol=3, loc="lower right")
    L = summary[min(1, a.n_end - 1)]
    fig.suptitle("%s: L = %.4f m (sd %.3f mm, n_g %g), connector %.2f mm"
                 % (a.tag, L[0], L[1], a.ng, zc_ref), fontsize=11)
    fig.tight_layout()
    png = os.path.join(a.out_dir, "%s_zoom4.png" % a.tag)
    fig.savefig(png, dpi=140)
    print("wrote: %s" % png)


if __name__ == "__main__":
    main()
