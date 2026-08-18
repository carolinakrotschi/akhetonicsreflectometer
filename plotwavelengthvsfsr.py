#!/usr/bin/env python3
"""Wavelength vs. FSR: find fringe maxima in the raw signal (scipy).

The raw signal is not a clean single tone -- it sits (see CHANGES.md) close
to the Nyquist limit and several reflections overlap, so the sign flips
almost every sample. Searching for peaks directly on that gives nothing but
noise.

So: first find the dominant fringe frequency via FFT and isolate it with a
Gaussian window in frequency space (a smooth bandpass, no filter ringing),
THEN use scipy.signal.find_peaks to find the maxima of the cleaned single
tone. The wavelength spacing between consecutive maxima is the local FSR
(Free Spectral Range).

For an interferometer with a fixed path difference OPD:
    FSR(lambda) = lambda^2 / (n_g * OPD)
Over the ~50 nm bandwidth here, FSR only changes by ~7 % because of this --
the curve looks practically linear. That's exactly what gets plotted,
including a linear least-squares fit.

    python3 plotwavelengthvsfsr.py                  # all files in "test data"
    python3 plotwavelengthvsfsr.py file.json         # a single file
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from fft_reflectometer_fixed import (
    DATA_DIR,
    N_GROUP_DEFAULT,
    load_sweep,
    clean,
    balanced,
    detrend,
)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#9a9993"
S1 = "#2a78d6"     # filtered fringe
S2 = "#eb6834"     # peaks / fit
RAW = "#c7c6c1"    # raw signal, background context only


# --------------------------------------------------------- isolate fringe
def isolate_fringe(sig: np.ndarray, d_lam: float, bw_frac: float = 0.4):
    """Find the dominant frequency via FFT, then isolate it with a Gaussian
    window in frequency space -> a smooth single tone without the
    overlapping secondary reflections / noise beyond Nyquist."""
    n = sig.size
    win = np.hanning(n)
    f = np.fft.rfftfreq(n, d=d_lam)                       # cycles/nm

    S_win = np.fft.rfft(sig * win)
    skip = max(2, int(0.02 * f.size))                     # exclude DC/envelope leftovers
    f0 = float(f[skip + np.argmax(np.abs(S_win[skip:]))])

    S_raw = np.fft.rfft(sig)
    mask = np.exp(-0.5 * ((f - f0) / (bw_frac * f0)) ** 2)
    filtered = np.fft.irfft(S_raw * mask, n=n)
    return filtered, f0


def find_fringe_peaks(filtered: np.ndarray, d_lam: float, f0: float):
    period_samples = 1.0 / (f0 * d_lam)
    dist = max(1, int(0.5 * period_samples))
    prom = 0.4 * np.std(filtered)
    idx, _ = find_peaks(filtered, distance=dist, prominence=prom)
    return idx


def fsr_from_peaks(wl: np.ndarray, idx: np.ndarray):
    peak_wl = wl[idx]
    mid_wl = 0.5 * (peak_wl[:-1] + peak_wl[1:])
    fsr_pm = np.diff(peak_wl) * 1e3
    return mid_wl, fsr_pm


# ----------------------------------------------------------------- plot
def plot(path: Path, wl, raw_sig, filtered, idx, mid_wl, fsr_pm, f0, ch: str, fit):
    fig = plt.figure(figsize=(12, 8), facecolor=SURFACE)
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.1], hspace=0.38, wspace=0.28)

    period_nm = 1.0 / f0
    span_zoom = 25 * period_nm                          # ~25 fringes per zoom window
    n = wl.size
    centers_idx = [int(n * 0.08), int(n * 0.5), int(n * 0.92)]

    for col, ci in enumerate(centers_idx):
        ax = fig.add_subplot(gs[0, col])
        ax.set_facecolor(SURFACE)
        c = wl[ci]
        m = (wl >= c - span_zoom / 2) & (wl <= c + span_zoom / 2)
        ax.plot(wl[m], raw_sig[m], color=RAW, lw=0.8, label="raw signal")
        ax.plot(wl[m], filtered[m], color=S1, lw=1.4, label="isolated fringe")
        sel = idx[m[idx]] if idx.size else idx
        ax.plot(wl[sel], filtered[sel], "o", color=S2, ms=5, zorder=5, label="peak")
        ax.set_title(f"~{c:.1f} nm", color=INK, fontsize=10, loc="left")
        ax.grid(True, color=MUTED, alpha=0.25, lw=0.6)
        for s in ax.spines.values():
            s.set_color(MUTED)
        ax.tick_params(colors=INK2, labelsize=8)
        if col == 0:
            ax.set_ylabel("Signal", color=INK2)
            ax.legend(fontsize=7, frameon=False, loc="upper right")
        ax.set_xlabel("Wavelength [nm]", color=INK2, fontsize=9)

    ax2 = fig.add_subplot(gs[1, :])
    ax2.set_facecolor(SURFACE)
    ax2.plot(mid_wl, fsr_pm, "o", color=S1, ms=4.5, alpha=0.85, label="measured FSR")
    if fit is not None:
        slope, intercept, r2 = fit
        xx = np.array([mid_wl.min(), mid_wl.max()])
        ax2.plot(xx, slope * xx + intercept, color=S2, lw=2, ls="--",
                 label=f"Fit: FSR = {slope:.5f}*lambda + {intercept:.3f}  (R^2={r2:.3f})")
    ax2.set_xlabel("Wavelength [nm]", color=INK2)
    ax2.set_ylabel("FSR [pm]", color=INK2)
    ax2.set_title(f"Wavelength vs. FSR — {ch}", color=INK, fontsize=12, loc="left")
    ax2.grid(True, color=MUTED, alpha=0.25, lw=0.6)
    for s in ax2.spines.values():
        s.set_color(MUTED)
    ax2.tick_params(colors=INK2)
    ax2.legend(fontsize=9, frameon=False, loc="upper left")

    fig.savefig(path, dpi=170, facecolor=SURFACE)
    plt.close(fig)


# ----------------------------------------------------------------- main
def process(path: Path, channel: str | None, n_g: float, out: Path):
    wl, chans = load_sweep(path)
    wl, chans = clean(wl, chans)
    name, y = balanced(chans, channel)
    raw_sig = detrend(y)
    d_lam = float(np.median(np.diff(wl)))

    filtered, f0 = isolate_fringe(raw_sig, d_lam)
    idx = find_fringe_peaks(filtered, d_lam, f0)

    print(f"File      : {path.name}")
    print(f"  Channel          : {name}")
    print(f"  Samples          : {wl.size}  (step {d_lam*1e3:.4f} pm)")
    print(f"  Dominant fringe  : {1/f0*1e3:.4f} pm  ({f0:.4f} 1/nm)")

    if idx.size < 3:
        print("  -> too few peaks found, skipping.\n")
        return

    mid_wl, fsr_pm = fsr_from_peaks(wl, idx)
    slope, intercept = np.polyfit(mid_wl, fsr_pm, 1)
    pred = slope * mid_wl + intercept
    ss_res = np.sum((fsr_pm - pred) ** 2)
    ss_tot = np.sum((fsr_pm - fsr_pm.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    lam0 = float(mid_wl.mean())
    fsr0_nm = (slope * lam0 + intercept) * 1e-3
    opd_est_m = (lam0 ** 2) / (n_g * fsr0_nm) * 1e-9 if fsr0_nm > 0 else float("nan")

    print(f"  Peaks found      : {idx.size}")
    print(f"  FSR (mean)       : {np.mean(fsr_pm):.4f} pm  (std {np.std(fsr_pm):.4f} pm)")
    print(f"  Fit FSR(lambda)  : {slope:.5f}*lambda + {intercept:.3f}   (R^2 = {r2:.3f})")
    print(f"  -> OPD at lambda0={lam0:.1f} nm : {opd_est_m*100:.2f} cm"
          f"  (reflector distance {opd_est_m*50:.2f} cm)")

    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{path.stem}_wl_vs_fsr.png"
    plot(p, wl, raw_sig, filtered, idx, mid_wl, fsr_pm, f0, name, (slope, intercept, r2))
    print(f"  -> {p}\n")


def main():
    ap = argparse.ArgumentParser(description="Find peaks in the raw signal and plot wavelength vs. FSR.")
    ap.add_argument("filename", nargs="?")
    ap.add_argument("--channel", help="force a single channel instead of balanced detection")
    ap.add_argument("--group-index", type=float, default=N_GROUP_DEFAULT)
    ap.add_argument("--output-dir", default="results")
    a = ap.parse_args()

    paths = [DATA_DIR / a.filename] if a.filename else sorted(DATA_DIR.glob("*.json"))
    out = Path(__file__).resolve().parent / a.output_dir
    for p in paths:
        process(p, a.channel, a.group_index, out)


if __name__ == "__main__":
    main()
