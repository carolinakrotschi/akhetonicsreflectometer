#!/usr/bin/env python3
"""
Just the reflectogram plot + CSV -- no peak list, no truth comparison.

Runs the same core pipeline as process_reflectogram_aux.py (load, balanced
subtraction, aux-phase resampling, window, FFT -- imported from there, not
duplicated) and writes the reflectogram CSV plus one PNG: distance vs.
amplitude, same axes/title style as that script's plot, with the aux arm
length difference (dL_aux) annotated top-right.

Usage:
    python tools/plot_reflectogram_only.py raw_data/scan.json --tau-aux-ns 17.6058
    python tools/plot_reflectogram_only.py raw_data/scan.json --dl 4 --zmax 2.5
"""

import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal.windows import kaiser, hann, blackmanharris

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import process_reflectogram_aux as pra


def build_argparser():
    p = pra.build_argparser()
    p.description = __doc__
    return p


def make_plot(a):
    if a.tau_aux_ns is not None:
        tau_aux = a.tau_aux_ns * 1e-9
    elif a.dl is not None:
        tau_aux = pra.NG * a.dl / pra.C
    else:
        sys.exit("specify --tau-aux-ns or --dl")

    dl_aux = pra.C * tau_aux / pra.NG

    ch1_, ch2_, ch3_, ch4_, meta = pra.load(a.scan)
    chans = {1: ch1_, 2: ch2_, 3: ch3_, 4: ch4_}
    ch1, ch2 = chans[a.aux_a], chans[a.aux_b]
    ch3, ch4 = chans[a.meas_a], chans[a.meas_b]

    if a.single:
        aux_ch = ch2 if np.median(ch2) > np.median(ch1) else ch1
        meas_ch = ch4 if np.median(ch4) > np.median(ch3) else ch3
        seg_a = max(len(aux_ch) // 32, 256)
        seg_m = max(len(meas_ch) // 32, 256)
        aux = aux_ch - uniform_filter1d(aux_ch, seg_a)
        meas = meas_ch - uniform_filter1d(meas_ch, seg_m)
    else:
        aux, _ = pra.balanced(ch1, ch2)
        meas, _ = pra.balanced(ch3, ch4)

    y, dnu, span_nu, _ = pra.resample_on_aux(meas, aux, tau_aux, a.trim)
    m = len(y)

    t = np.linspace(-1, 1, m)
    y = y - np.polyval(np.polyfit(t, y, 5), t)

    dz_bin = pra.C / (2 * pra.NG * span_nu)
    z_nyq = pra.C / (4 * pra.NG * dnu)

    win = {"hann": hann(m), "blackmanharris": blackmanharris(m),
           "kaiser": kaiser(m, a.kaiser_beta)}[a.window]
    R = np.abs(np.fft.rfft(y * win))
    z = np.arange(len(R)) * pra.C / (2 * pra.NG * dnu * m)
    db = 20 * np.log10(R / R.max() + 1e-15)

    zmax = a.zmax if a.zmax else z[-1]
    keep = z <= zmax

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nf = pra.noise_floor(db[keep])
    print(f"RMS noise floor {nf:.1f} dB   -> dynamic range {-nf:.1f} dB")
    fig, ax = plt.subplots(figsize=(11, 5))
    ax.plot(z[keep], db[keep], lw=0.6)
    ax.axhline(nf, color="darkviolet", ls="--", lw=1.1, zorder=4)
    ax.annotate(f"RMS noise floor {nf:.1f} dB  (dynamic range {-nf:.1f} dB)",
                (z[keep][-1], nf), fontsize=8.5, ha="right", va="bottom",
                color="darkviolet",
                bbox=dict(fc="white", ec="none", alpha=.75, pad=1.0))
    ax.set_xlabel("Distance (m, one-way / reflection convention)")
    ax.set_ylabel("Amplitude (dB rel. maximum)")
    ax.set_title(f"{a.scan} | aux-referenced, {a.window}, "
                 f"dz {dz_bin*1e6:.1f} um, Nyquist {z_nyq:.2f} m")
    ax.grid(alpha=0.3)
    ax.set_ylim(max(-110, db[keep].min() - 5), 5)
    ax.text(0.99, 0.97, f"dL$_{{aux}}$ = {dl_aux:.4f} m",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10, bbox=dict(boxstyle="round", facecolor="white",
                                    edgecolor="gray", alpha=0.85))
    fig.tight_layout()

    prefix = a.out or a.scan.rsplit(".", 1)[0]

    csv_path = f"{prefix}_reflectogram.csv"
    np.savetxt(csv_path, np.column_stack([z[keep], db[keep]]),
               delimiter=",", header="distance_m,amplitude_dB", comments="")
    print(f"wrote: {csv_path}")

    png_path = f"{prefix}_reflectogram.png"
    fig.savefig(png_path, dpi=150)
    print(f"wrote: {png_path}")


def main():
    a = build_argparser().parse_args()
    make_plot(a)


if __name__ == "__main__":
    main()
