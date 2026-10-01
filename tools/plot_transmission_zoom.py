#!/usr/bin/env python3
"""
Zoom ladder for a transmission scan, with an instrument reference overlaid.

The reference is a reflectogram taken with the SAME sweep settings (span,
speed, rate); its strongest peak is shifted onto the main peak of the scan,
so whatever sits at the same offset in both curves is the instrument skirt,
not the DUT.

    python tools/plot_transmission_zoom.py scan_reflectogram.csv \
        --ref ref_reflectogram.csv --out results/<date>/<name>_zoom4.png
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    z, a = np.loadtxt(path, delimiter=",", skiprows=1, unpack=True)
    return z * 1e3, a


def main():
    p = argparse.ArgumentParser()
    p.add_argument("scan")
    p.add_argument("--ref", required=True)
    p.add_argument("--ref-label", default="instrument reference (same sweep)")
    p.add_argument("--mark", type=float, nargs="*", default=[],
                   help="offsets in mm from the main peak to mark")
    p.add_argument("--title", default="")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    z, s = load(a.scan)
    zr, r = load(a.ref)
    z0 = z[np.argmax(s)]
    zr_shift = zr - zr[np.argmax(r)] + z0

    lo = z[1] - z[0]
    panels = [(z[0], z[-1]), (z0 - 500, z0 + 500), (z0 - 5, z0 + 5),
              (z0 - 1, z0 + 2.5)]
    fig, ax = plt.subplots(len(panels), 1, figsize=(12, 13))
    for k, (x0, x1) in enumerate(panels):
        m = (z >= x0) & (z <= x1)
        ax[k].plot(z[m], s[m], lw=0.6 if k == 0 else 0.9, label="transmission")
        if k >= 2:
            mr = (zr_shift >= x0) & (zr_shift <= x1)
            ax[k].plot(zr_shift[mr], r[mr], lw=0.8, alpha=0.8,
                       label=a.ref_label)
            for d in a.mark:
                ax[k].axvline(z0 + d, color="C3", ls=":", lw=1)
        ax[k].set_xlim(x0, x1)
        ax[k].set_ylim(-100 if k < 2 else -60, 3)
        ax[k].set_ylabel("dB rel. main peak")
        ax[k].grid(alpha=0.3)
        ax[k].set_title(f"{x0:.1f} - {x1:.1f} mm", fontsize=9)
        if k >= 2:
            ax[k].legend(fontsize=8, loc="upper right")
    ax[-1].set_xlabel(f"z (mm, absolute, reflection convention, bin {lo*1e3:.2f} um)")
    fig.suptitle(a.title or a.scan, fontsize=10)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote:", a.out, f"(main peak {z0:.4f} mm)")


if __name__ == "__main__":
    main()
