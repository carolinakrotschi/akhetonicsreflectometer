#!/usr/bin/env python3
"""OFDR evaluation of a wavelength sweep (EXFO JSON).

Corrected version. See CHANGES.md for the main differences from the old
version.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).resolve().parent / "test data"
N_GROUP_DEFAULT = 1.468          # SMF-28 @1550 nm, group index
C0 = 299_792_458.0               # m/s

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#9a9993"
S1 = "#2a78d6"                   # Ch1 / signal
S2 = "#eb6834"                   # Ch2 / warning marker


# ------------------------------------------------------------------ loading
def load_sweep(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    with path.open("r", encoding="utf-8") as fh:
        raw = json.load(fh)
    block = raw["data"][0]
    arrays = {
        k: np.asarray(v, dtype=float)
        for k, v in block.items()
        if isinstance(v, list) and v and isinstance(v[0], (int, float))
    }
    wl_key = next(k for k in arrays if k.lower().startswith("wavelength"))
    wl = arrays.pop(wl_key)
    return wl, arrays


def clean(wl: np.ndarray, chans: dict[str, np.ndarray]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Drop non-monotonic outliers (the EXFO appends a jump back to the
    sweep's start wavelength at the end -- that wrecks mean(diff))."""
    keep = np.ones(wl.size, dtype=bool)
    step = np.median(np.diff(wl))
    bad = np.where(np.abs(np.diff(wl) - step) > 50 * abs(step))[0]
    for i in bad:
        keep[i + 1] = False
    if (~keep).any():
        print(f"  [clean] dropped {int((~keep).sum())} outlier(s) on the wavelength axis")
    return wl[keep], {k: v[keep] for k, v in chans.items()}


# --------------------------------------------------------------- signal path
def balanced(chans: dict[str, np.ndarray], channel: str | None) -> tuple[str, np.ndarray]:
    """Two complementary coupler ports -> (a-b)/(a+b).

    That removes the laser's power envelope exactly instead of approximately
    and doubles the fringe contrast.
    """
    if channel:
        return channel, chans[channel]
    ports = [k for k in chans if k.lower().startswith("ch")]
    if len(ports) >= 2:
        a, b = chans[ports[0]], chans[ports[1]]
        if np.corrcoef(a, b)[0, 1] < -0.3:          # actually anti-phase?
            return f"({ports[0]}-{ports[1]})/({ports[0]}+{ports[1]})", (a - b) / (a + b)
    return ports[0], chans[ports[0]]


def detrend(y: np.ndarray, width: int = 1001) -> np.ndarray:
    w = min(width | 1, y.size // 4 | 1)
    ker = np.ones(w) / w
    env = np.convolve(np.pad(y, (w // 2, w // 2), mode="edge"), ker, mode="valid")[: y.size]
    return y - env


def reflectogram(wl: np.ndarray, y: np.ndarray, n_g: float, zero_pad: int = 4):
    """FFT against 1/lambda -> the peak frequency is directly n_g * OPD.

    The device's sampling grid is uniform in lambda; the reported
    wavelengths are rounded to 1 pm. For the phase, the ideal grid is a
    better assumption than the rounded values.
    """
    lam = np.linspace(wl[0], wl[-1], wl.size)        # nm, ideally uniform
    k = 1.0 / lam                                    # 1/nm
    order = np.argsort(k)
    ku = np.linspace(k[order][0], k[order][-1], k.size)
    yu = np.interp(ku, k[order], y[order])

    sig = detrend(yu)
    n = sig.size
    spec = np.abs(np.fft.rfft(sig * np.hanning(n), n=zero_pad * n))
    f = np.fft.rfftfreq(zero_pad * n, d=ku[1] - ku[0])   # unit: nm

    opd_m = f / n_g * 1e-9                            # optical path difference [m]
    return lam, opd_m, spec


def limits(lam: np.ndarray, n_g: float) -> dict[str, float]:
    d_lam = (lam[-1] - lam[0]) / (lam.size - 1)       # nm
    lam0 = float(lam.mean())
    span = float(lam[-1] - lam[0])
    B_hz = C0 * span * 1e-9 / (lam0 * 1e-9) ** 2
    return {
        "d_lam_pm": d_lam * 1e3,
        "lam0": lam0,
        "span_nm": span,
        "B_THz": B_hz / 1e12,
        "resolution_um": C0 / (2 * n_g * B_hz) * 1e6,
        "opd_max_m": lam0 ** 2 / (2 * n_g * d_lam) * 1e-9,
    }


# ----------------------------------------------------------------- plot
def plot(path: Path, lam, opd_m, spec, lim, n_g, round_trip: bool, ch: str):
    scale = 0.5 if round_trip else 1.0
    x = opd_m * scale * 100                            # cm
    xmax = lim["opd_max_m"] * scale * 100
    label = "Reflector distance [cm]" if round_trip else "Arm length difference [cm]"

    fig, ax = plt.subplots(figsize=(11, 5.2), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    m = x <= xmax * 1.02
    ax.plot(x[m], spec[m], color=S1, lw=1.0)
    ax.axvline(xmax, color=S2, lw=2, ls="--")
    ax.text(xmax, ax.get_ylim()[1] * 0.97, f"  Nyquist {xmax:.1f} cm",
            color=S2, va="top", ha="left", fontsize=10, fontweight="bold")
    ax.set_xlim(0, xmax * 1.02)
    ax.set_xlabel(label, color=INK2)
    ax.set_ylabel("Amplitude", color=INK2)
    ax.set_title(f"Reflectogram — {ch}   ·   resolution {lim['resolution_um']:.0f} um",
                 color=INK, fontsize=12, loc="left")
    ax.grid(True, color=MUTED, alpha=0.25, lw=0.6)
    for s in ax.spines.values():
        s.set_color(MUTED)
    ax.tick_params(colors=INK2)
    fig.tight_layout()
    fig.savefig(path, dpi=170, facecolor=SURFACE)
    plt.close(fig)


# ----------------------------------------------------------------- main
def process(path: Path, channel: str | None, n_g: float, out: Path, round_trip: bool):
    wl, chans = load_sweep(path)
    wl, chans = clean(wl, chans)
    name, y = balanced(chans, channel)
    lam, opd_m, spec = reflectogram(wl, y, n_g)
    lim = limits(lam, n_g)
    scale = 0.5 if round_trip else 1.0

    print(f"File      : {path.name}")
    print(f"  Channel          : {name}")
    print(f"  Samples          : {lam.size}")
    print(f"  lambda           : {lam[0]:.3f} .. {lam[-1]:.3f} nm  (span {lim['span_nm']:.3f} nm)")
    print(f"  Step size        : {lim['d_lam_pm']:.4f} pm")
    print(f"  Bandwidth        : {lim['B_THz']:.3f} THz")
    print(f"  Resolution dz    : {lim['resolution_um']:.1f} um")
    print(f"  NYQUIST LIMIT    : OPD {lim['opd_max_m']*100:.2f} cm"
          f"  =  reflector distance {lim['opd_max_m']*50:.2f} cm")

    band = opd_m * scale <= lim["opd_max_m"] * scale
    valid = band & (opd_m * scale > 0.002)
    i = int(np.argmax(np.where(valid, spec, 0)))
    d_peak = opd_m[i] * scale
    print(f"  Peak             : {d_peak*100:.3f} cm"
          f"   ({'reflector distance' if round_trip else 'arm length difference'})")

    # energy distribution -> aliasing warning
    e_tot = (spec[valid] ** 2).sum()
    top = valid & (opd_m * scale > 0.8 * lim["opd_max_m"] * scale)
    frac = (spec[top] ** 2).sum() / e_tot
    if frac > 0.25 or d_peak > 0.8 * lim["opd_max_m"] * scale:
        print(f"\n  !!! WARNING: {frac*100:.0f} % of the energy is in the top 20 %"
              f" before the Nyquist limit.\n"
              f"      That is the signature of aliasing. The displayed length is"
              f" probably WRONG.\n"
              f"      Remedy: reduce the wavelength step delta-lambda"
              f" (halve it -> doubles the range).\n"
              f"      Cross-check: record a sweep with half the delta-lambda."
              f" If the peak moves, it was folded.\n")

    out.mkdir(parents=True, exist_ok=True)
    p = out / f"{path.stem}_reflectogram.png"
    plot(p, lam, opd_m, spec, lim, n_g, round_trip, name)
    np.savez_compressed(out / f"{path.stem}_reflectogram.npz",
                        opd_m=opd_m, amplitude=spec, **lim)
    print(f"  -> {p}\n")


def main():
    ap = argparse.ArgumentParser(description="OFDR evaluation of an EXFO wavelength sweep.")
    ap.add_argument("filename", nargs="?")
    ap.add_argument("--channel", help="force a single channel instead of balanced detection")
    ap.add_argument("--group-index", type=float, default=N_GROUP_DEFAULT)
    ap.add_argument("--single-pass", action="store_true",
                    help="arms are single-pass (no circulator/no reflection) "
                         "-> x-axis is the arm length difference instead of reflector distance")
    ap.add_argument("--output-dir", default="results")
    a = ap.parse_args()

    paths = [DATA_DIR / a.filename] if a.filename else sorted(DATA_DIR.glob("*.json"))
    out = Path(__file__).resolve().parent / a.output_dir
    for p in paths:
        process(p, a.channel, a.group_index, out, round_trip=not a.single_pass)


if __name__ == "__main__":
    main()
