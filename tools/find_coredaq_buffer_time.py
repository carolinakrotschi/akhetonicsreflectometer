#!/usr/bin/env python3
"""
Find the coreDAQ-LIN "buffer time" -- the delay between coreDAQ recording
start and the laser sweep actually beginning to ramp -- from a bandpass
filter response-sweep recording (LinaWindow / Interface/OBR.py).

Why this works: point a narrow bandpass filter at the sweep's *start*
wavelength (Sweep start = filter's passband). If coreDAQ starts recording
before the sweep physically begins, the laser sits still at the start
wavelength for a while, so Ch2/Ch4 (the strong pair) stay pinned near
their transmission maximum for that entire dead period -- then drop
sharply within milliseconds once the sweep starts moving the wavelength
out of the (typically sub-nm-wide) passband. That drop's location in
*recording time* (not the GUI's naive linear-model wavelength axis, which
assumes no dead time and is exactly what this dead time invalidates) is
the buffer time.

Usage:
    python tools/find_coredaq_buffer_time.py raw_data/2026-09-01-bandpassfilter1527nm.json --rate 100000
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scan")
    p.add_argument("--rate", type=float, required=True,
                   help="coreDAQ sample rate in Hz")
    p.add_argument("--smooth", type=int, default=500)
    p.add_argument("--out", default=None)
    args = p.parse_args()

    with open(args.scan) as f:
        doc = json.load(f)
    e = doc["data"][0]
    ch2 = np.asarray(e["Ch2 [mW]"], float)
    ch4 = np.asarray(e["Ch4 [mW]"], float)
    n = len(ch2)
    t = np.arange(n) / args.rate

    env = uniform_filter1d(ch2 + ch4, args.smooth)
    i_peak = int(np.argmax(env))
    plateau = env[:i_peak].mean() if i_peak > args.smooth else env.max()

    thresholds = {}
    for frac in (0.9, 0.5, 0.1):
        thresh = plateau * frac
        above = env > thresh
        idxs = np.nonzero(above)[0]
        i_last = idxs[idxs > i_peak].max() if (idxs > i_peak).any() else idxs.max()
        if env[i_last + 1] == env[i_last]:
            t_cross = t[i_last]
        else:
            t0, t1 = t[i_last], t[i_last + 1]
            v0, v1 = env[i_last], env[i_last + 1]
            t_cross = t0 + (thresh - v0) * (t1 - t0) / (v1 - v0)
        thresholds[frac] = t_cross

    print(f"{args.scan}: plateau level {plateau:.5f} mW (Ch2+Ch4, smoothed)")
    for frac, tc in thresholds.items():
        print(f"  {frac*100:.0f}% falling-edge crossing: {tc:.5f} s")
    print(f"\nbuffer time (onset of sweep motion) ~= {thresholds[0.9]:.3f} s "
          f"(90% edge) .. {thresholds[0.5]:.3f} s (50% edge)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(t, env, lw=0.8)
    for frac, tc in thresholds.items():
        ax.axvline(tc, color="gray", ls="--", lw=0.8)
        ax.annotate(f"{frac*100:.0f}%: {tc:.3f}s", (tc, plateau * frac),
                    textcoords="offset points", xytext=(5, 5), fontsize=8)
    ax.set_xlabel("recording time (s)")
    ax.set_ylabel("Ch2+Ch4 (mW, smoothed)")
    ax.set_title(f"{args.scan} -- coreDAQ buffer-time estimate")
    ax.grid(alpha=0.3)
    fig.tight_layout()

    out_path = args.out or (Path(args.scan).stem + "_buffertime.png")
    fig.savefig(out_path, dpi=150)
    print(f"\nwrote: {out_path}")


if __name__ == "__main__":
    main()
