#!/usr/bin/env python3
"""
Compare the 6-32cm "dirt band" across several scans on one plot.

Each curve is normalized to its own peak (dB rel. own max), so this shows
whether the band's brightness changes with what's connected at the test
port -- it does NOT show absolute power differences between scans (see
raw_data/README.md for that check).

Usage: python tools/plot_band_comparison.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import diagnose_artifacts as da

RAW = Path(__file__).resolve().parent.parent / "raw_data"
OUT = Path(__file__).resolve().parent.parent / "results" / "2026-08-18"

FILES = [
    (RAW / "2026-08-18_open_end.json", "open end (no terminator)"),
    (RAW / "2026-08-18_no_fiber.json", "no test fiber"),
    (RAW / "2026-08-18_terminator1.json", "terminator 1"),
    (RAW / "2026-08-18_terminator2.json", "terminator 2"),
    (RAW / "2026-08-18_10db_coupler.json", "10dB coupler"),
]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(12, 5.5))
    print(f"{'condition':25s} {'median dB':>10s} {'max dB':>10s}")
    for path, label in FILES:
        z, R = da.spectrum_diag(*da.preprocess(*da.load(str(path))))
        db = 20 * np.log10(R / R.max() + 1e-15)
        m = (z >= 0.06) & (z <= 0.32)
        ax.plot(z[m] * 100, db[m], lw=0.6, label=label, alpha=0.8)
        print(f"{label:25s} {np.median(db[m]):10.1f} {db[m].max():10.1f}")

    ax.set_xlabel("apparent distance (cm)")
    ax.set_ylabel("dB rel. own max")
    ax.set_title("6-32 cm dirt band: comparison across scans")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()

    OUT.mkdir(parents=True, exist_ok=True)
    out_path = OUT / "band_6_32cm_comparison.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
