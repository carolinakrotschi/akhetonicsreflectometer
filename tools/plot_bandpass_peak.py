#!/usr/bin/env python3
"""
Find and plot the transmission peak of a physical bandpass filter placed in
the beam path, from raw EXFO/CoreDAQ scans.

Unlike the OFDR reflectogram scripts, this does not run the interferometric
pipeline at all -- when a bandpass filter blocks the sweep almost
everywhere, ALL FOUR channels (aux pair and measurement pair alike) read
near-zero power except inside the passband, where they all show power
simultaneously (each scaled by its usual channel split ratio). So the peak
is found directly on raw power vs. wavelength, on the strongest channel
(Ch4, or Ch2 as a cross-check), lightly smoothed only to denoise -- no
interferometric processing needed.

Usage:
    python tools/plot_bandpass_peak.py raw_data/2026-09-01-bandpassfilter1550nm_*.json
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d


def load_wl_channel(path, channel="Ch4 [mW]"):
    with open(path) as f:
        doc = json.load(f)
    e = doc["data"][0]
    wl = np.asarray(e["Wavelength [nm]"], float)
    x = np.asarray(e[channel], float)
    return wl, x


def find_peak(wl, x, smooth_pts=100):
    sm = uniform_filter1d(x, smooth_pts)
    i = int(np.argmax(sm))
    half = sm[i] / 2
    l = r = i
    while l > 0 and sm[l] > half:
        l -= 1
    while r < len(sm) - 1 and sm[r] > half:
        r += 1
    fwhm_pm = (wl[r] - wl[l]) * 1000
    return wl[i], sm, fwhm_pm


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scans", nargs="+")
    p.add_argument("--channel", default="Ch4 [mW]")
    p.add_argument("--center-nm", type=float, default=None,
                   help="nominal center wavelength for the pm x-axis "
                        "(default: mean of the found peaks)")
    p.add_argument("--out", default=None)
    args = p.parse_args()

    results = []
    for scan in args.scans:
        wl, x = load_wl_channel(scan, args.channel)
        peak_wl, sm, fwhm_pm = find_peak(wl, x)
        results.append((scan, wl, sm, peak_wl, fwhm_pm))
        print(f"{scan}: peak at {peak_wl:.4f} nm, FWHM {fwhm_pm:.1f} pm, "
              f"peak power {sm.max()*1e3:.2f} uW")

    peaks = np.array([r[3] for r in results])
    print(f"\n{len(results)} scan(s): peak {peaks.mean():.4f} +/- "
          f"{peaks.std():.4f} nm (range {peaks.min():.4f}-{peaks.max():.4f} nm)")

    center = args.center_nm if args.center_nm is not None else peaks.mean()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 5))
    for scan, wl, sm, peak_wl, fwhm_pm in results:
        ax.plot((wl - center) * 1000, sm * 1e3, lw=0.8,
                label=f"{Path(scan).stem}  (peak {peak_wl:.3f} nm, "
                      f"FWHM {fwhm_pm:.0f} pm)")
    ax.set_xlabel(f"Wavelength offset from {center:.3f} nm (pm)")
    ax.set_ylabel(f"{args.channel} (uW, smoothed)")
    ax.set_title("Bandpass filter transmission peak")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    max_fwhm_pm = max(r[4] for r in results)
    ax.set_xlim(-4 * max_fwhm_pm, 4 * max_fwhm_pm)
    fig.tight_layout()

    prefix = args.out or "bandpass_peak"
    out_path = f"{prefix}.png"
    fig.savefig(out_path, dpi=150)
    print(f"\nwrote: {out_path}")


if __name__ == "__main__":
    main()
