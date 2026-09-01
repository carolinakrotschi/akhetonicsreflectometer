"""
One overview figure of everything the Lina wavelength-axis investigation
measured. Reads only saved files (no hardware), so it can be re-run any time
to regenerate the summary as more data is collected.

  python lina/scripts/lina_cal_summary.py [--data lina_cal_data] [--out <png>]

Panels
------
  1  TOF1550 passband measured in CW              (the marker itself)
  2  raw capture in time, with the sweep window   (what the buffer contains)
  3  lambda(t) fit through the filter markers     (the calibration)
  4  fit residuals per polynomial degree          (why degree 2)
  5  axis error: naive linear vs calibrated       (what it is worth)
  6  marker position across repeated sweeps       (stability in time)
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

from lina.analysis.lina_wl_cal import (
    fit_time_to_wavelength, wavelength_axis, sweep_window,
)
from lina_sweep_test import (
    load_raw, find_buffer_peaks, select_forward_peak,
)

MARKER_NM = 1550.0        # the filter position used for the single-sweep tests


def newest(pattern):
    hits = sorted(glob.glob(pattern), key=os.path.getmtime)
    return hits[-1] if hits else None


def panel_passband(ax, data_dir):
    path = newest(os.path.join(data_dir, "filter_cw_*.csv"))
    if not path:
        ax.set_title("no filter-cw data")
        return
    d = np.genfromtxt(path, delimiter=",", names=True)
    wl = d["wavelength_nm"]
    cols = [c for c in d.dtype.names if c.endswith("_mw")]
    for c in cols:
        p = d[c]
        if p.max() < 1e-3:        # skip the near-dark channels
            continue
        ax.plot(wl, p * 1000, marker=".", ms=3, lw=0.9, label=c.replace("_mw", ""))
    ax.axvline(MARKER_NM, color="k", ls=":", lw=1.0, label=f"set {MARKER_NM:.0f} nm")
    ax.set_xlabel("laser CW wavelength (nm)")
    ax.set_ylabel("power (µW)")
    ax.set_title("1 · TOF1550 passband in CW\n(the wavelength marker, FWHM ≈ 0.18 nm)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)


def panel_raw(ax, data_dir, cal):
    path = newest(os.path.join(data_dir, "sweep_t100k_*.npz")) or \
           newest(os.path.join(data_dir, "sweep_*.npz"))
    if not path:
        ax.set_title("no raw capture")
        return None
    res = load_raw(path)
    rate, n = res["rate_hz"], res["n_pts"]
    ch = max(res["traces_mw"], key=lambda c: np.nanmax(res["traces_mw"][c]))
    tr = res["traces_mw"][ch]
    t = np.arange(n) / rate
    ax.plot(t, tr * 1000, lw=0.8, color="tab:blue", label=f"Ch{ch}")
    peaks, _, _ = find_buffer_peaks(tr, rate_hz=rate)
    fwd = select_forward_peak(peaks)
    for p in peaks:
        is_fwd = p is fwd
        ax.axvline(p["index"] / rate, color="tab:red" if is_fwd else "tab:gray",
                   ls="--", lw=0.9)
        ax.annotate(f"{'sweep crossing' if is_fwd else 'return slew'}\n"
                    f"{p['index'] / rate:.2f} s, {p['width_s'] * 1000:.0f} ms",
                    (p["index"] / rate, tr.max() * 1000),
                    fontsize=6.5, ha="left" if is_fwd else "right", va="top")
    if cal is not None:
        lo, hi = sweep_window(cal, n, rate)
        ax.axvspan(lo / rate, hi / rate, color="tab:green", alpha=0.12,
                   label=f"measured sweep window\n{lo / rate:.2f}–{hi / rate:.2f} s")
        ax.axvspan(hi / rate, n / rate, color="tab:red", alpha=0.10,
                   label=f"post-sweep tail (cropped)\n{(n / rate - hi / rate):.2f} s")
    ax.set_xlabel("time in capture (s)")
    ax.set_ylabel("power (µW)")
    ax.set_title(f"2 · raw capture, {res['rate_hz'] // 1000} kHz, "
                 f"{n / 1e6:.1f} M samples\n"
                 f"{res['actual_start_nm']:.0f}→{res['actual_stop_nm']:.0f} nm @ "
                 f"{res['actual_speed_nm_s']:.0f} nm/s")
    ax.legend(fontsize=6.5, loc="upper left")
    ax.grid(alpha=0.3)
    return res


def panel_fit(ax, cal):
    pts = cal["points"]
    t = np.array([p["time_s"] for p in pts])
    wl = np.array([p["true_nm"] for p in pts])
    meta = cal["sweep"]
    tt = np.linspace(0, t.max() * 1.05, 300)
    ax.plot(t, wl, "o", ms=7, color="k", zorder=5,
            label=f"{len(pts)} filter markers")
    ax.plot(tt, meta["actual_start_nm"] + meta["actual_speed_nm_s"] * tt,
            ls="--", lw=1.2, color="tab:red",
            label=f"assumed: {meta['actual_speed_nm_s']:.2f} nm/s commanded")
    f2 = fit_time_to_wavelength(t, wl, degree=2)
    ax.plot(tt, np.polyval(f2["coeffs"], tt), lw=1.2, color="tab:blue",
            label=f"measured λ(t), degree 2\n({f2['residual_pm_rms']:.1f} pm rms)")
    ax.set_xlabel("time from trigger edge (s)")
    ax.set_ylabel("wavelength (nm)")
    ax.set_title("3 · the calibration\nλ(t), in seconds → rate-independent")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)


def panel_residuals(ax, cal):
    pts = cal["points"]
    t = np.array([p["time_s"] for p in pts])
    wl = np.array([p["true_nm"] for p in pts])
    for d, col in zip((1, 2, 3), ("tab:red", "tab:blue", "tab:green")):
        if d > len(pts) - 1:
            break
        f = fit_time_to_wavelength(t, wl, degree=d)
        ax.plot(wl, (wl - np.polyval(f["coeffs"], t)) * 1000, "o-", ms=5,
                lw=1.0, color=col,
                label=f"degree {d}: {f['residual_pm_rms']:.1f} pm rms")
    ax.axhline(0, color="k", lw=0.8)
    ax.axhspan(-32, 32, color="tab:green", alpha=0.10,
               label="filter's own setting error (±32 pm)")
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("residual (pm)")
    ax.set_title("4 · why degree 2\ndegree 1 leaves a systematic arch")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)


def panel_axis_error(ax, res, cal):
    """What each axis choice does to a marker of known wavelength."""
    if res is None or cal is None:
        ax.set_title("no capture to compare")
        return
    rate, n = res["rate_hz"], res["n_pts"]
    ch = max(res["traces_mw"], key=lambda c: np.nanmax(res["traces_mw"][c]))
    tr = res["traces_mw"][ch]
    peaks, _, _ = find_buffer_peaks(tr, rate_hz=rate)
    fwd = select_forward_peak(peaks)
    lo, hi = fwd["start"], fwd["stop"] + 1
    w = np.clip(tr[lo:hi] - np.median(tr), 0, None)
    apparent = lambda axis: float(np.sum(axis[lo:hi] * w) / np.sum(w))
    a0, a1 = res["actual_start_nm"], res["actual_stop_nm"]

    labels, errors = [], []
    labels.append("linear axis\n(old default)")
    errors.append(apparent(np.linspace(a0, a1, n)) - MARKER_NM)
    for bt in (1.0, 2.0, 3.0):
        crop = int(round(bt * rate))
        axis = np.concatenate([np.full(crop, np.nan),
                               np.linspace(a0, a1, n - crop)])
        labels.append(f"+ 'Buffer time'\n{bt:.0f} s (removed)")
        errors.append(apparent(axis) - MARKER_NM)
    labels.append("measured λ(t)\n(new)")
    errors.append(apparent(wavelength_axis(cal, n, rate)) - MARKER_NM)

    colors = ["tab:red"] * (len(errors) - 1) + ["tab:green"]
    bars = ax.bar(range(len(errors)), errors, color=colors, alpha=0.85)
    for b, e in zip(bars, errors):
        ax.annotate(f"{e * 1000:+.0f} pm" if abs(e) < 0.1 else f"{e:+.2f} nm",
                    (b.get_x() + b.get_width() / 2, e), ha="center",
                    va="top" if e < 0 else "bottom", fontsize=7.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=6.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_ylabel(f"error at the {MARKER_NM:.0f} nm marker (nm)")
    ax.set_title("5 · what it is worth\nsame data, different wavelength axis")
    ax.grid(alpha=0.3, axis="y")


def panel_stability(ax, data_dir, cal):
    paths = sorted(glob.glob(os.path.join(data_dir, "repeat_*.csv")),
                   key=os.path.getmtime)
    if not paths:
        ax.set_title("6 · no repeat data yet")
        return
    speed = cal["sweep"]["actual_speed_nm_s"] if cal else 5.0
    for path in paths:
        d = np.genfromtxt(path, delimiter=",", names=True)
        el = np.atleast_1d(d["elapsed_s"]) / 60.0
        ts = np.atleast_1d(d["time_s"])
        dev = (ts - ts.mean()) * speed * 1000
        tag = os.path.basename(path).split("_")[1]
        ax.plot(el, dev, "o-", ms=5,
                label=f"{tag}: {len(ts)} runs over {el.max():.0f} min, "
                      f"{dev.std(ddof=1) if dev.size > 1 else 0:.1f} pm rms")
    ax.axhline(0, color="k", lw=0.8)
    ax.axhspan(-5, 5, color="tab:green", alpha=0.10,
               label="±5 pm (the fit's own residual)")
    ax.set_xlabel("elapsed time (min)")
    ax.set_ylabel("marker deviation (pm)")
    ax.set_title("6 · stability in time\nis the calibration still valid later?")
    ax.legend(fontsize=6.5)
    ax.grid(alpha=0.3)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="lina_cal_data")
    ap.add_argument("--cal", default=None,
                    help="calibration JSON (default: newest in --data)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    cal_path = a.cal or newest(os.path.join(a.data, "lina_wl_cal*.json"))
    cal = json.load(open(cal_path)) if cal_path else None
    if cal is None:
        sys.exit(f"no calibration JSON found in {a.data}")
    print(f"[in] calibration: {cal_path}")

    fig, axes = plt.subplots(2, 3, figsize=(17.5, 9.5))
    panel_passband(axes[0, 0], a.data)
    res = panel_raw(axes[0, 1], a.data, cal)
    panel_fit(axes[0, 2], cal)
    panel_residuals(axes[1, 0], cal)
    panel_axis_error(axes[1, 1], res, cal)
    panel_stability(axes[1, 2], a.data, cal)

    d = cal.get("derived", {})
    fig.suptitle(
        f"Lina wavelength-axis calibration — exfo-1 + coredaq-1 + TOF1550   |   "
        f"sweep occupies t = {d.get('buffer_time_s', 0):.3f}–"
        f"{d.get('sweep_end_s', 0):.3f} s of the armed buffer   |   "
        f"effective speed {d.get('effective_speed_nm_s', float('nan')):.3f} nm/s "
        f"vs {cal['sweep']['actual_speed_nm_s']:.3f} commanded",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = a.out or os.path.join(a.data, "SUMMARY_lina_wl_cal.png")
    fig.savefig(out, dpi=130)
    print(f"[out] {out}")


if __name__ == "__main__":
    main()
