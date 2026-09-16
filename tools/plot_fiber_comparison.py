#!/usr/bin/env python3
"""
Fiber-comparison plots used in the 2026-08-31 progress report (KW36).

Reads the already-processed per-scan reflectogram CSVs under results/ (each
has columns distance_m, amplitude_dB, normalized to 0 dB at that scan's own
maximum) and produces:

  - results/<date>/fiber_comparison_relative_to_connector.png
    every scan's dominant (0 dB) connector peak shifted to x=0, calibrated
    tau_aux axis.
  - results/<date>/fiber_comparison_extended_with_2026-08-31.png
    same data, absolute (not connector-relative) x-axis.
  - results/<date>/fiber_comparison_aux_corrected.png
    same scans, but read from results/<date>/aux_corrected/ (dL forced to
    each scan's filename-indicated nominal aux length instead of the
    fringe-counting-calibrated tau_aux) -- see logs/2026-08-31.md point 4.
  - results/<date>/<date>_vs_2026-08-20_nofiber_comparison.png
    the five 2026-08-20 no-fiber connector-swap controls vs. today's two
    repeat no-fiber scans, to check whether today's connector position is a
    new state or matches where the bench was left last session.

"First peak" = each curve's own argmax (its connector reflection).
"Important peak" = the strongest peak beyond a guard band around the first
peak (found with scipy.signal.find_peaks on the dB trace, floor -50 dB),
i.e. the fiber/aux secondary reflection -- see HANDOVER.md and
logs/2026-08-31.md for the double-bounce interpretation.

Usage: python tools/plot_fiber_comparison.py
"""

import sys
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from process_reflectogram_aux import noise_floor   # noqa: E402
DATE = "2026-08-31"
OUT = ROOT / "results" / DATE

PEAK_FLOOR_DB = -50.0
GUARD_MM = 50.0       # exclude the first peak's own sidelobes from the search
SEARCH_MAX_MM = 4500.0
MIN_PEAK_SPACING_M = 200e-6  # matches process_reflectogram_aux.py's find_peaks distance

# category -> matplotlib color sequence, one entry per scan in that category,
# in the order scans are listed below (gray=no fiber, red=1m, blue=3m, purple=combo)
CATEGORY_COLORS = {
    "nofiber": ["#bdbdbd", "#757575", "#4d4d4d", "#262626"],
    "1m": ["#f4a582", "#e8836a", "#d6604d", "#b2182b", "#67001f"],
    "3m": ["#92c5de", "#4393c3", "#2166ac", "#053061"],
    "combo": ["#762a83"],
}

# (label, category, path to calibrated reflectogram CSV, path to
#  aux-corrected reflectogram CSV) -- both relative to ROOT.
SCANS = [
    ("No fiber (control)", "nofiber",
     "results/2026-08-20/2026-08-20_nofiberattheend_trimmed1530_reflectogram_ch1ch3meas_wide.csv",
     "results/2026-08-31/aux_corrected/nofiber_control_reflectogram.csv"),
    ("No fiber, 2026-08-20 last scan (connector #2, settle-margin test)", "nofiber",
     "results/2026-08-20/2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleanedandnewconnectorneueeinschwingzeit_full_reflectogram_ch1ch3meas_wide_reflectogram.csv",
     "results/2026-08-31/aux_corrected/nofiber_0820last_reflectogram.csv"),
    ("No fiber, 2026-08-31 (aux stretched to ~4.49 m)", "nofiber",
     "results/2026-08-31/2026-08-31-ohnefiberaux5m_reflectogram_reflectogram.csv",
     "results/2026-08-31/aux_corrected/nofiber_aux5m_reflectogram.csv"),
    ("No fiber, 2026-08-31 (aux stretched to ~5.39 m)", "nofiber",
     "results/2026-08-31/2026-08-31-ohnefiberaux6m_reflectogram_reflectogram.csv",
     "results/2026-08-31/aux_corrected/nofiber_aux6m_reflectogram.csv"),
    ("1 m fiber #1", "1m",
     "results/2026-08-20/2026-08-19_deltaLwirklich1m_trimmed1530_reflectogram_ch1ch3meas_wide.csv",
     "results/2026-08-31/aux_corrected/1mfiber1_reflectogram.csv"),
    ("1 m fiber #2", "1m",
     "results/2026-08-20/2026-08-19_neuefaserauch1mlang_trimmed1530_reflectogram_ch1ch3meas_wide.csv",
     "results/2026-08-31/aux_corrected/1mfiber2_reflectogram.csv"),
    ("1 m fiber #3", "1m",
     "results/2026-08-20/2026-08-19_neueneuefaserauch1mlang_trimmed1530_reflectogram_ch1ch3meas_wide.csv",
     "results/2026-08-31/aux_corrected/1mfiber3_reflectogram.csv"),
    ("1 m fiber #4", "1m",
     "results/2026-08-20/2026-08-19_neueneueneuefaserauch1mlang_trimmed1530_reflectogram_ch1ch3meas_wide.csv",
     "results/2026-08-31/aux_corrected/1mfiber4_reflectogram.csv"),
    ("1 m fiber, 2026-08-31 (aux 4.49 m)", "1m",
     "results/2026-08-31/2026-08-31-1mfiberaux5m_reflectogram_reflectogram.csv",
     "results/2026-08-31/aux_corrected/1mfiber_aux5m_reflectogram.csv"),
    ("3m fiber", "3m",
     "results/2026-08-20/2026-08-20_3mfiber_trimmed1530_reflectogram_ch1ch3meas_wide.csv",
     "results/2026-08-31/aux_corrected/3mfiber_reflectogram.csv"),
    ("3m fiber, 2026-08-31 (aux 4.49 m)", "3m",
     "results/2026-08-31/2026-08-31-3mfiberaux5m_reflectogram_reflectogram.csv",
     "results/2026-08-31/aux_corrected/3mfiber_aux5m_reflectogram.csv"),
    ("3m fiber, 2026-08-31 (aux 3.58 m)", "3m",
     "results/2026-08-31/2026-08-31-3mfiberaux4m_reflectogram_reflectogram.csv",
     "results/2026-08-31/aux_corrected/3mfiber_aux4m_reflectogram.csv"),
    ("3m + 1m fiber in series, 2026-08-31 (aux 3.58 m)", "combo",
     "results/2026-08-31/2026-08-31-3mplus1mfiberaux4m_reflectogram_reflectogram.csv",
     "results/2026-08-31/aux_corrected/3mplus1mfiber_aux4m_reflectogram.csv"),
]

# 2026-08-20 connector-swap controls vs. today's two repeats, for the
# separate connector-family comparison plot.
NOFIBER_FAMILY_SCANS = [
    ("2026-08-20 no-fiber, Thorlabs connector",
     "results/2026-08-20/2026-08-20_nofiberattheend_trimmed1530_reflectogram_ch1ch3meas_wide.csv",
     "#4c72b0", 0.7, 1),
    ("2026-08-20 no-fiber, connector #1",
     "results/2026-08-20/2026-08-20_nofiberattheendandnotthorlabsconnector_trimmed1530_reflectogram_ch1ch3meas_wide_reflectogram.csv",
     "#dd8452", 0.7, 1),
    ("2026-08-20 no-fiber, connector #1 (cleaned)",
     "results/2026-08-20/2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleaned_trimmed1530_reflectogram_ch1ch3meas_wide_reflectogram.csv",
     "#55a868", 0.7, 1),
    ("2026-08-20 no-fiber, connector #2",
     "results/2026-08-20/2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleanedandnewconnector_trimmed1533_reflectogram_ch1ch3meas_wide_reflectogram.csv",
     "#8c4b4b", 0.7, 1),
    ("2026-08-20 no-fiber, connector #2 (settle-margin test)",
     "results/2026-08-20/2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleanedandnewconnectorneueeinschwingzeit_full_reflectogram_ch1ch3meas_wide_reflectogram.csv",
     "#8172b2", 0.7, 1),
    ("2026-08-31 no-fiber (today, 1st scan)",
     "results/2026-08-31/2026-08-31-ohnefiberallesnormalerstemessungdestages_reflectogram_reflectogram.csv",
     "black", 1.0, 2),
    ("2026-08-31 no-fiber (today, 2nd scan)",
     "results/2026-08-31/2026-08-31-ohnefiberallesnormalzweitemessungdestages_reflectogram_reflectogram.csv",
     "magenta", 1.0, 2),
]

REFLECTOR_587_MM = 587.82   # Thorlabs/connector-#2 family reference (own peak, no-fiber control)
REFLECTOR_518_MM = 518.0    # second fixed internal reflection, see HANDOVER.md


def load_csv(path):
    d = np.genfromtxt(ROOT / path, delimiter=",", skip_header=1)
    return d[:, 0] * 1000.0, d[:, 1]  # distance in mm, amplitude in dB


def find_first_and_important_peak(z_mm, db):
    """Own connector peak (argmax) + strongest peak beyond the guard band."""
    i0 = int(np.argmax(db))
    first_mm = z_mm[i0]
    rel_mm = z_mm - first_mm

    dz_m = (z_mm[1] - z_mm[0]) / 1000.0
    distance = max(3, int(MIN_PEAK_SPACING_M / dz_m))
    pk, _ = find_peaks(db, height=PEAK_FLOOR_DB, distance=distance)

    candidates = [
        (rel_mm[p], db[p]) for p in pk
        if GUARD_MM < rel_mm[p] < SEARCH_MAX_MM
    ]
    if not candidates:
        return first_mm, None, None
    important_rel_mm, important_db = max(candidates, key=lambda c: c[1])
    return first_mm, important_rel_mm, important_db


def legend_label(label, category, first_mm, important_rel_mm):
    base = f"{label} (first peak {first_mm:.2f} mm"
    # no-fiber scans have no real secondary reflector -- whatever find_peaks
    # turns up beyond the guard band there is just noise, not worth showing
    if important_rel_mm is None or category == "nofiber":
        return base + ")"
    return base + f"; important peak {important_rel_mm:+.1f} mm relative)"


def plot_comparison(csv_key, mode, out_path, title, subtitle, xlabel):
    """mode: 'relative' (x-axis re-referenced to each scan's first peak) or
    'absolute' (raw distance)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(20, 8))
    color_counts = {cat: 0 for cat in CATEGORY_COLORS}
    floors = []

    for label, category, csv_calibrated, csv_aux_corrected in SCANS:
        csv_path = csv_calibrated if csv_key == "csv" else csv_aux_corrected
        z_mm, db = load_csv(csv_path)
        first_mm, important_rel_mm, _ = find_first_and_important_peak(z_mm, db)

        x = (z_mm - first_mm) if mode == "relative" else z_mm
        color = CATEGORY_COLORS[category][color_counts[category]]
        color_counts[category] += 1
        lw = 2.2 if category == "combo" else 0.8

        nf = noise_floor(db)
        floors.append(nf)
        ax.plot(x, db, color=color, lw=lw, alpha=0.85,
                label=legend_label(label, category, first_mm, important_rel_mm)
                + f"  | RMS floor {nf:.1f} dB")
        ax.axhline(nf, color=color, ls="--", lw=0.9, alpha=0.5)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Amplitude (dB relative to each scan's own maximum)")
    ax.set_title(f"{title}\n{subtitle}" if subtitle else title, fontsize=11)
    ax.legend(fontsize=8, title="Scan (fiber length) -- gray=no fiber, red=1m, blue=3m, purple=3m+1m combo")
    ax.set_ylim(min([-60.0] + [f - 5.0 for f in floors]), 5)
    if mode == "relative":
        ax.set_xlim(-50, 4500)
    else:
        ax.set_xlim(0, 4500)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def plot_nofiber_family(out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    ranges = [(0, 750), (0, 3000)]
    titles = ["0-750 mm (587mm / 518mm family)", "0-3000 mm (wider view)"]

    for ax, (xmin, xmax), sub_title in zip(axes, ranges, titles):
        for label, csv_path, color, alpha, lw in NOFIBER_FAMILY_SCANS:
            z_mm, db = load_csv(csv_path)
            nf = noise_floor(db)
            ax.plot(z_mm, db, color=color, alpha=alpha, lw=lw,
                    label=f"{label}  | RMS floor {nf:.1f} dB")
            ax.axhline(nf, color=color, ls="--", lw=0.9, alpha=0.5)
        ax.axvline(REFLECTOR_587_MM, color="gray", ls="--", lw=1)
        ax.axvline(REFLECTOR_518_MM, color="gray", ls=":", lw=1)
        ax.set_xlim(xmin, xmax)
        ax.set_title(sub_title)
        ax.set_xlabel("Distance (mm, reflection convention)")
        ax.set_ylabel("Amplitude (dB rel. own scan's max)")
        ax.legend(fontsize=8)

    fig.suptitle("No-fiber control: today vs. 2026-08-20 "
                  "(dashed=587mm family, dotted=518mm family)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    plot_comparison(
        csv_key="csv", mode="relative",
        out_path=OUT / "fiber_comparison_relative_to_connector.png",
        title="Fibre lengths on the fringe-counting-calibrated tau_aux axis "
              "(NO aux correction)",
        subtitle="",
        xlabel="Distance relative to each scan's own connector peak (mm)",
    )

    plot_comparison(
        csv_key="csv", mode="absolute",
        out_path=OUT / f"fiber_comparison_extended_with_{DATE}.png",
        title="Fiber comparison -- all scans referenced (aux = Ch2/Ch4, measurement = Ch1/Ch3)",
        subtitle="587mm connector family held through 2026-08-20's main sequence; the "
                  "connector's position later shifted to ~495-496mm (last 2026-08-20 scan "
                  f"and all of {DATE}, unrelated to aux length) -- see logs/{DATE}.md",
        xlabel="Distance (mm, one-way / reflection convention)",
    )

    plot_comparison(
        csv_key="aux_csv", mode="relative",
        out_path=OUT / "fiber_comparison_aux_corrected.png",
        title="Fibre lengths with dL forced to the nominal aux length "
              "(aux-corrected)",
        subtitle="",
        xlabel="Distance relative to each scan's own connector peak (mm)",
    )

    plot_nofiber_family(OUT / f"{DATE}_vs_2026-08-20_nofiber_comparison.png")


if __name__ == "__main__":
    main()
