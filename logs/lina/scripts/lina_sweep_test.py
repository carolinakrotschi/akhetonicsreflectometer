"""
Headless Lina (EXFO T200S + LINEAR coreDAQ) sweep test & wavelength-axis
calibration, extracted from LabGUI/dashboards/instruments/lina_window.py
(``_CoreDAQSweepWorker``, stepped=False path) so the trigger-synchronisation
problem can be investigated without the GUI.

WHY THIS EXISTS
---------------
In LINEAR mode the coreDAQ capture is started by ONE rising edge on the EXFO
"Trig out" BNC (OUTPut:SYNChronization:STATe ON goes high while the laser is
sweeping) and then free-runs at its own sample rate, so every sample carries a
TIME, not a wavelength. The GUI assumed the buffer maps linearly onto
[start, stop] nm; measured on exfo-1 + coredaq-1 that is wrong by several nm.

What the capture actually contains (measured, 1520->1570 nm @ 5 nm/s, 12 s
armed buffer):

  * the sweep, occupying t = 0.004 .. 9.889 s — it starts essentially AT the
    trigger edge: no laser settling is captured, because the EXFO holds its
    sync output until it has settled (pre-parking and settling 5 s before
    arming moves the marker by <1 pm)
  * a 2.1 s post-sweep tail, including the laser slewing back to its parked
    wavelength at t ~ 10.12 s, re-crossing the same wavelengths BACKWARDS

So the dead time is at the END, and the GUI's "Buffer time" front-crop
corrects the wrong end — it deletes real sweep data and makes the error worse
(1550 nm marker apparent at 1544.6 nm with no crop, 1539.5 nm with 2 s).

This script measures the mapping instead, using the Thorlabs TOF1550 tunable
filter in the path as a wavelength marker: it only transmits near its own
centre wavelength (0.21 nm passband), so its transmission peak marks a KNOWN
wavelength at a MEASURED time. Parking it at several wavelengths and fitting
gives lambda(t) — in SECONDS from the trigger edge, so the result is
rate-independent: a fit measured at 5 kHz recovered a 1550.000 nm marker to
-3 pm in a 100 kHz / 1.2 M-sample capture. Degree 2 fits to ~3 pm rms
(degree 1 leaves a systematic ~27 pm arch).

The filter is driven over its own serial port (lina/interface/TunableFilters.py),
so calibration runs unattended. The fit format and its application live in
lina/analysis/lina_wl_cal.py, shared with the GUI window.

SUBCOMMANDS
-----------
  sweep        One triggered capture: RAW, uncropped buffer + diagnostics.
               With --cal, also reports where the marker lands on the
               calibrated axis vs the naive linear one.
  filter-cw    Step the laser in CW across the filter passband -> the filter's
               TRUE peak wavelength (ground truth, independent of the filter's
               own <32 pm setting error).
  calibrate    For each filter position: optional CW ground truth + a sweep,
               locating the marker in time. Fits lambda(t), reports the
               effective sweep rate and dead time, and publishes the result
               where the Lina GUI window looks it up.
  repeat       N identical sweeps with the filter parked — run-to-run jitter
               and drift of the mapping (measured: 0.8 pm rms over 1.7 min).
  refit        Re-fit a saved calibration at another degree (no hardware).
  plot-cal     Plot a saved calibration: fit, residuals, naive-axis error.
  apply        Re-plot a saved raw capture on a calibrated axis.

EXAMPLES
--------
  python lina/scripts/lina_sweep_test.py filter-cw --filter-nm 1550
  python lina/scripts/lina_sweep_test.py --channels 2 calibrate \
      --start 1520 --stop 1570 --speed 5 --rate 5000 --trust-filter \
      --filter-positions 1530,1535,1540,1545,1550,1555,1560,1565
  python lina/scripts/lina_sweep_test.py --channels 2 sweep \
      --start 1520 --stop 1570 --speed 5 --rate 100000 --filter-nm 1550 \
      --cal Interface/config/calibrations/lina/lina_wl_cal_1520.0-1570.0nm_5.00nms.json
  python lina/scripts/lina_sweep_test.py --channels 2 repeat \
      --filter-nm 1550 --runs 6 --interval-s 180

NOTE on the sample rate: keep calibration runs at a few kHz. The calibration
is in seconds, so it applies at 100 kHz anyway, and the smaller USB transfer
is far less likely to stall (see run_sweep's retry comment).

NOTE: close the LabGUI Lina/coreDAQ windows first — LabDevice allows only one
open instance per serial id, and the coreDAQ USB-CDC port is exclusive.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from datetime import datetime

# repo root (for Interface/, LabGUI/) and the lina package root
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import numpy as np

from Interface.Lasers import EXFO
from Interface.PowerMeters import CoreDAQ
from lina.interface.TunableFilters import TOF1550
# Fit / axis / crop logic lives in one place, shared with the Lina GUI window.
from lina.analysis.lina_wl_cal import (
    fit_time_to_wavelength, wavelength_axis, sweep_window,
    save_calibration, describe as describe_cal, CAL_DIR,
)


# ---------------------------------------------------------------------------
# Device helpers
# ---------------------------------------------------------------------------

def open_devices(exfo_name: str, coredaq_name: str, filter_name: str = None):
    laser = EXFO(exfo_name)
    laser.open()
    pm = CoreDAQ(coredaq_name)
    pm.open()
    print(f"[dev] {exfo_name} @ {laser.ip}:{laser.port} | "
          f"{coredaq_name} {pm.frontend()} @ {pm.get_sample_rate_hz()} Hz")
    tof = None
    if filter_name:
        tof = TOF1550(filter_name)
        tof.open()
        print(f"[dev] {filter_name} @ {tof.port}, center {tof.get_wavelength():.3f} nm")
    return laser, pm, tof


def close_devices(laser, pm, tof=None):
    for dev in (tof, pm, laser):
        try:
            if dev is not None:
                dev.close()
                dev._deregister()
        except Exception as e:
            print(f"[dev] close failed: {e}")


def park(laser, wl_nm: float, settle_s: float):
    """Move to a CW wavelength and let the laser settle before anything is armed."""
    laser.set_wavelength(wl_nm)
    time.sleep(settle_s)


def read_power_mw(pm, channels, n_avg: int = 4) -> dict:
    out = {}
    for ch in channels:
        vals = [pm.get_power(ch) for _ in range(n_avg)]
        out[ch] = float(np.mean(vals))
    return out


# ---------------------------------------------------------------------------
# Triggered sweep capture — the LINEAR/continuous path, raw and uncropped
# ---------------------------------------------------------------------------

def run_sweep(laser, pm, *, start_nm, stop_nm, speed_nm_s, channels,
              rate_hz, n_samples=None, pad_s=2.0, park_settle_s=5.0,
              restore_wl=True, attempts=3):
    """Arm the coreDAQ on the EXFO sweep trigger, sweep, and return the RAW buffer.

    Unlike the GUI worker this crops nothing and builds no wavelength axis —
    the whole point is to see where in the buffer the sweep actually lives.

    Returns a dict with the raw per-channel mW arrays and the sweep metadata.
    """
    initial_wl = laser.get_wavelength()
    laser.enable("ON")

    # 1. Program the sweep, then read back what the EXFO actually accepted
    #    (it clamps/rounds start, stop and speed).
    laser.set_sweep_parameters(start_wavelength=start_nm,
                               stop_wavelength=stop_nm,
                               speed_nm_s=speed_nm_s, cycles=1)
    actual_start = float(laser.get_sweep_start_wavelength()) * 1e9
    actual_stop = float(laser.get_sweep_stop_wavelength()) * 1e9
    actual_speed = float(laser.get_sweep_speed()) * 1e9
    sweep_s = abs(actual_stop - actual_start) / max(actual_speed, 1e-6)

    # 2. Buffer size: sweep duration + pad, so the leading dead time and the
    #    start of the return slew are both inside the capture and visible.
    pm.set_sample_rate_hz(int(rate_hz))
    rate = int(pm.get_sample_rate_hz())
    if n_samples is None:
        n_frames = int(math.ceil((sweep_s + pad_s) * rate))
    else:
        n_frames = int(n_samples)
    try:
        max_frames = int(pm.max_capture_frames(channels))
        if n_frames > max_frames:
            print(f"[cap] {n_frames} frames > device max {max_frames}; capping")
            n_frames = max_frames
    except Exception:
        pass

    print(f"[cap] {actual_start:.3f} -> {actual_stop:.3f} nm @ "
          f"{actual_speed:.3f} nm/s = {sweep_s:.2f} s sweep; arming "
          f"{n_frames} frames @ {rate} Hz = {n_frames / rate:.2f} s, ch={channels}")

    # 3.-6. Arm, sweep, read — retried as a whole. The USB-CDC read of the
    #       buffer intermittently trickles and dies part-way (the py_coreDAQ
    #       transport's own 30 s budget expires with only 70-90 % received),
    #       and a failed read consumes the arming, so the only recovery is to
    #       re-arm and re-sweep. Keep the transfer small (a lower --rate) to
    #       make this rare.
    caps, n_pts, blocked_s, acq_state, captured = {}, 0, float("nan"), None, None
    last_err = None
    for attempt in range(1, int(attempts) + 1):
        try:
            # Park at the start wavelength and let it settle BEFORE arming, so
            # the parking slew is not part of the measurement window at all.
            park(laser, actual_start, park_settle_s)

            # Sync out ON (goes high while sweeping), then arm on the rising edge.
            laser.trig_out(state="ON", mode="SWEEP",
                           trigger_in=False, trigger_out=True)
            # The driver may arm one frame fewer to keep the transfer size
            # off a 512-byte USB packet boundary (see arm_sweep_capture) —
            # use what it actually armed.
            n_frames = pm.arm_sweep_capture(
                n_frames, channels=channels, trigger_rising=True,
                use_trigger=True, stepped=False)

            # Sweep. Blocks until the EXFO reports it is no longer sweeping or
            # setting. No coreDAQ USB I/O in this window — a command mid-DMA
            # corrupts samples.
            t0 = time.monotonic()
            laser.start_continuous_sweep()
            blocked_s = time.monotonic() - t0
            print(f"[cap] start_continuous_sweep() blocked {blocked_s:.2f} s "
                  f"(sweep itself {sweep_s:.2f} s -> {blocked_s - sweep_s:+.2f} s "
                  f"of settle/return inside the trigger-high window)")

            # The capture free-runs for n_frames/rate seconds from the TRIGGER
            # EDGE, independent of how long the sweep takes. If the armed
            # window is longer than the sweep (it must be, to see the return
            # slew), the buffer is still filling when start_continuous_sweep()
            # returns — reading it then stalls the transfer part-way, because
            # the remaining frames do not exist yet. So wait the window out.
            capture_s = n_frames / rate
            remaining = capture_s - blocked_s + 0.3
            if remaining > 0:
                print(f"[cap] waiting {remaining:.2f} s for the armed window "
                      f"({capture_s:.2f} s) to finish filling…")
                time.sleep(remaining)

            acq_state = pm._acq_state()
            labels = {0: "idle", 1: "ARMED — never triggered!",
                      2: "acquiring", 4: "done"}
            print(f"[cap] acq_state={acq_state} "
                  f"({labels.get(acq_state, 'unknown/firmware-specific')})")

            # Never finalize_capture() on LINEAR — stop_capture() zeroes the
            # driver's armed-frame bookkeeping that the pre-v4.3 collect path
            # needs, turning the read into "no capture was armed".
            time.sleep(0.5)
            captured = pm.captured_frames()
            if captured == -1 and acq_state != 1:
                captured = n_frames  # FRAMES? unsupported on this firmware
            caps = (pm.read_sweep_capture(captured, channels=channels)
                    if captured and captured > 0 else {})
            n_pts = min((len(v) for v in caps.values()), default=0)
            if n_pts == 0:
                raise RuntimeError(
                    f"coreDAQ captured 0 frames [armed={n_frames} rate={rate} "
                    f"acq_state={acq_state} captured={captured}] — check the "
                    f"EXFO Trig Out BNC is wired to the coreDAQ trigger input.")
            break
        except Exception as e:
            last_err = e
            print(f"[cap] attempt {attempt}/{attempts} failed: "
                  f"{type(e).__name__}: {e}")
            try:
                pm.abort_capture()
            except Exception:
                pass
            if attempt == attempts:
                try:
                    laser.trig_out(state="OFF", mode="SWEEP", trigger_out=False)
                except Exception:
                    pass
                if restore_wl:
                    try:
                        laser._write(f"SOUR:WAV:CW {initial_wl}NM")
                    except Exception:
                        pass
                raise RuntimeError(
                    f"capture failed after {attempts} attempt(s): "
                    f"{type(last_err).__name__}: {last_err}\n"
                    f"       {n_frames} frames x {len(channels)} ch = "
                    f"~{n_frames * 2 * len(channels) / 1e6:.2f} MB per read. "
                    f"Lower --rate/--pad or use fewer channels to shrink it."
                ) from last_err
            # Wait out the EXFO before re-arming: after an aborted sweep it
            # stays busy for a few seconds and rejects the trigger-output
            # reconfiguration with Metrino.T200S.BusyException, which would
            # burn a retry on the laser rather than on the failure we are
            # actually retrying (observed: 1 s was not enough).
            time.sleep(5.0)

    try:
        laser.trig_out(state="OFF", mode="SWEEP", trigger_out=False)
    except Exception:
        pass
    if restore_wl:
        try:
            laser._write(f"SOUR:WAV:CW {initial_wl}NM")
        except Exception:
            pass

    print(f"[cap] captured {n_pts} frames/ch ({n_pts / rate:.2f} s)")
    return {
        "traces_mw": {int(ch): np.asarray(caps[ch], dtype=float)[:n_pts]
                      for ch in caps},
        "n_pts": n_pts,
        "rate_hz": rate,
        "actual_start_nm": actual_start,
        "actual_stop_nm": actual_stop,
        "actual_speed_nm_s": actual_speed,
        "sweep_s": sweep_s,
        "blocked_s": blocked_s,
        "requested_start_nm": float(start_nm),
        "requested_stop_nm": float(stop_nm),
        "requested_speed_nm_s": float(speed_nm_s),
        "acq_state": acq_state,
    }


# ---------------------------------------------------------------------------
# CW step scan — the filter's true passband peak
# ---------------------------------------------------------------------------

def run_filter_cw(laser, pm, *, center_nm, span_nm, step_nm, channels,
                  settle_s=0.6, n_avg=4):
    """Step the laser across the filter passband in CW and record power.

    This is the ground truth: the x axis is the laser's own CW wavelength
    setting, so the resulting peak is a real wavelength — no dependence on any
    sweep timing, and independent of the filter's own setting error.
    """
    laser.enable("ON")
    wls = np.arange(center_nm - span_nm / 2, center_nm + span_nm / 2 + step_nm / 2,
                    step_nm)
    power = {ch: np.zeros(wls.size) for ch in channels}
    for i, wl in enumerate(wls):
        laser.set_wavelength(float(wl))
        time.sleep(settle_s)
        vals = read_power_mw(pm, channels, n_avg)
        for ch in channels:
            power[ch][i] = vals[ch]
        print(f"\r[cw] {i + 1}/{wls.size}  {wl:.3f} nm  " +
              "  ".join(f"Ch{ch} {vals[ch]:.4g}" for ch in channels) + "      ",
              end="", flush=True)
    print()
    peaks = {}
    for ch in channels:
        pk, fwhm = peak_wavelength(wls, power[ch])
        contrast = float(power[ch].max()) / max(float(np.median(power[ch])), 1e-15)
        peaks[ch] = {"peak_nm": pk, "fwhm_nm": fwhm, "contrast": contrast,
                     "max_mw": float(power[ch].max()),
                     "median_mw": float(np.median(power[ch]))}
        print(f"[cw] Ch{ch}: peak {pk:.4f} nm, FWHM {fwhm:.3f} nm, "
              f"max {power[ch].max():.4g} mW, baseline {np.median(power[ch]):.4g} mW, "
              f"contrast {contrast:.1f}x")
    # Pick by absolute peak power, not contrast: an unconnected channel can
    # show huge contrast on nanowatt-level leakage.
    best = max(channels, key=lambda c: peaks[c]["max_mw"])
    return {"wavelengths_nm": wls, "power_mw": power, "peaks": peaks,
            "best_channel": best,
            "peak_nm": peaks[best]["peak_nm"], "fwhm_nm": peaks[best]["fwhm_nm"]}


def peak_wavelength(x, y):
    """Centroid of the above-half-max region — robust for a filter passband."""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    base = float(np.median(y))
    amp = float(y.max()) - base
    if amp <= 0:
        return float(x[int(np.argmax(y))]), float("nan")
    mask = y >= base + amp / 2
    peak = float(np.sum(x[mask] * (y[mask] - base)) / np.sum(y[mask] - base))
    fwhm = float(x[mask].max() - x[mask].min())
    return peak, fwhm


# ---------------------------------------------------------------------------
# Peak location inside a raw sweep buffer
# ---------------------------------------------------------------------------

def find_buffer_peaks(trace, *, rate_hz, min_separation_s=0.05, rel_height=0.3):
    """Locate filter transmission peaks in a raw capture, in sample index.

    Above-threshold runs closer together than ``min_separation_s`` are MERGED
    into one peak: the crossing is not monotonic (noise, and TIA range
    switching near full scale, dip below the threshold part-way up), and
    treating those dips as peak boundaries otherwise splits one crossing into
    a sliver plus a missed maximum.

    The slow forward sweep crossing and the fast post-sweep return slew appear
    as separate peaks — the return is 10x narrower at the same height, so
    ``area_mw_s`` (not height) is what distinguishes them; see
    :func:`select_forward_peak`.
    """
    trace = np.asarray(trace, dtype=float)
    base = float(np.median(trace))
    amp = float(np.nanmax(trace)) - base
    if amp <= 0:
        return [], base, amp
    thresh = base + rel_height * amp
    idx = np.flatnonzero(trace >= thresh)
    if idx.size == 0:
        return [], base, amp
    gap = max(int(min_separation_s * rate_hz), 1)
    groups = np.split(idx, np.flatnonzero(np.diff(idx) > gap) + 1)
    peaks = []
    for g in groups:
        lo, hi = int(g[0]), int(g[-1])
        seg = trace[lo:hi + 1] - base
        w = np.clip(seg, 0.0, None)
        if w.sum() <= 0:
            continue
        centroid = float(np.sum(np.arange(lo, hi + 1) * w) / w.sum())
        peaks.append({"index": centroid,
                      "peak_index": int(lo + int(np.argmax(seg))),
                      "start": lo, "stop": hi,
                      "width_s": (hi - lo + 1) / rate_hz,
                      "amplitude_mw": float(seg.max()),
                      "area_mw_s": float(w.sum() / rate_hz)})
    return peaks, base, amp


def select_forward_peak(peaks):
    """The forward-sweep crossing: the peak with the largest integrated area.

    The laser crosses the passband slowly while sweeping (at 5 nm/s a 0.18 nm
    passband takes ~36 ms) and very fast on the way back (a few ms), so the
    forward crossing carries far more area even when both reach a similar
    height.
    """
    return max(peaks, key=lambda p: p["area_mw_s"]) if peaks else None


# ---------------------------------------------------------------------------
# Calibration reporting (the fit itself lives in lina/analysis/lina_wl_cal.py)
# ---------------------------------------------------------------------------

def report_fit(fit, meta):
    """Translate the fit into the numbers the GUI needs: dead time and rate."""
    coeffs = fit["coeffs"]
    rate = meta["rate_hz"]
    eff_speed = coeffs[-2] if len(coeffs) >= 2 else float("nan")
    # Time at which the fit reaches the actual sweep start / stop wavelength.
    real = lambda rs: [float(r.real) for r in np.atleast_1d(rs)
                       if abs(np.imag(r)) < 1e-6 and r.real > -1.0]
    offset = lambda target: np.array(coeffs, dtype=float) - np.array(
        [0.0] * (len(coeffs) - 1) + [target])
    t_start = min(real(np.roots(offset(meta["actual_start_nm"]))), default=float("nan"))
    t_stop = min(real(np.roots(offset(meta["actual_stop_nm"]))), default=float("nan"))
    print("\n=== calibration result ===")
    print(f"  fit degree           : {fit['degree']}  ({fit['units']})")
    print(f"  residual             : {fit['residual_pm_rms']:.1f} pm rms, "
          f"{fit['residual_pm_max']:.1f} pm max")
    print(f"  effective sweep rate : {eff_speed:.4f} nm/s  "
          f"(commanded {meta['actual_speed_nm_s']:.4f} nm/s, "
          f"{100 * eff_speed / max(meta['actual_speed_nm_s'], 1e-9):.1f} %)")
    print(f"  = {eff_speed / rate * 1000:.4f} pm per sample at {rate} Hz")
    if not math.isnan(t_start):
        print(f"  trigger dead time    : {t_start:.4f} s before the sweep "
              f"reaches {meta['actual_start_nm']:.3f} nm   <-- GUI 'Buffer time'")
    if not math.isnan(t_stop):
        print(f"  sweep ends at        : {t_stop:.4f} s; everything after that "
              f"is the return slew")
    return {"effective_speed_nm_s": eff_speed,
            "buffer_time_s": None if math.isnan(t_start) else t_start,
            "sweep_end_s": None if math.isnan(t_stop) else t_stop}


# ---------------------------------------------------------------------------
# Saving / plotting
# ---------------------------------------------------------------------------

def stamp():
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def save_raw(res, out_dir, tag="", write_csv=False):
    """Save a raw capture. NPZ always; CSV only on request.

    A 1.2 M-sample capture is ~30 kB as compressed NPZ but ~25 MB as CSV text,
    so the CSV is opt-in (``--csv``) — it is only worth writing when the data
    has to leave Python.
    """
    os.makedirs(out_dir, exist_ok=True)
    name = f"sweep_{tag + '_' if tag else ''}{stamp()}"
    npz = os.path.join(out_dir, name + ".npz")
    meta = {k: v for k, v in res.items() if k != "traces_mw"}
    np.savez_compressed(npz, meta=json.dumps(meta),
                        **{f"ch{ch}": tr for ch, tr in res["traces_mw"].items()})
    print(f"[out] {npz}")
    if write_csv:
        csv = os.path.join(out_dir, name + ".csv")
        chs = sorted(res["traces_mw"])
        cols = [np.arange(res["n_pts"]) / res["rate_hz"]] + \
               [res["traces_mw"][ch] for ch in chs]
        np.savetxt(csv, np.column_stack(cols), delimiter=",",
                   header="time_s," + ",".join(f"ch{ch}_mw" for ch in chs),
                   comments="")
        print(f"[out] {csv}")
    return npz


def load_raw(npz_path):
    d = np.load(npz_path, allow_pickle=False)
    meta = json.loads(str(d["meta"]))
    traces = {int(k[2:]): d[k] for k in d.files if k.startswith("ch")}
    meta["traces_mw"] = traces
    return meta


def plot_filter_cw(res, out_dir, tag="", filter_nm=None):
    import matplotlib.pyplot as plt
    wls = res["wavelengths_nm"]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for ch, p in sorted(res["power_mw"].items()):
        ax.plot(wls, p, marker=".", ms=3, lw=0.9, label=f"Ch{ch}")
        ax.axvline(res["peaks"][ch]["peak_nm"], lw=0.7, ls="--", alpha=0.5)
    if filter_nm is not None:
        ax.axvline(filter_nm, color="k", lw=1.0, ls=":",
                   label=f"filter set {filter_nm:.3f} nm")
    ax.set_xlabel("laser CW wavelength (nm)")
    ax.set_ylabel("power (mW)")
    ax.set_title("TOF1550 passband, measured in CW (ground truth)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, f"filter_cw_{tag + '_' if tag else ''}{stamp()}.png")
    fig.savefig(png, dpi=130)
    print(f"[out] {png}")
    try:
        plt.show()
    except Exception:
        pass


def plot_calibration(cal, out_dir, tag="", repeat_csv=None, max_degree=3):
    """Visualise the calibration: fit, residuals per degree, and the error the
    naive linear axis makes — the numbers that justify the fit degree."""
    import matplotlib.pyplot as plt

    pts = cal["points"]
    t = np.array([p["time_s"] for p in pts])
    wl = np.array([p["true_nm"] for p in pts])
    meta = cal["sweep"]
    start, stop = meta["actual_start_nm"], meta["actual_stop_nm"]
    speed = meta["actual_speed_nm_s"]

    n = 2 if repeat_csv is None else 3
    fig, axes = plt.subplots(2, n, figsize=(5.5 * n, 8))

    # (a) the fit itself
    ax = axes[0, 0]
    tt = np.linspace(0, max(t.max() * 1.05, meta["sweep_s"]), 400)
    ax.plot(t, wl, "o", ms=7, color="k", label="filter markers (measured)", zorder=5)
    ax.plot(tt, start + speed * tt, lw=1.2, ls="--", color="tab:red",
            label=f"naive linear axis ({speed:.3f} nm/s commanded)")
    for d in range(1, min(max_degree, len(pts) - 1) + 1):
        f = fit_time_to_wavelength(t, wl, degree=d)
        ax.plot(tt, np.polyval(f["coeffs"], tt), lw=1.0,
                label=f"degree {d} fit ({f['residual_pm_rms']:.1f} pm rms)")
    ax.set_xlabel("time from trigger edge (s)")
    ax.set_ylabel("wavelength (nm)")
    ax.set_title("sample time -> wavelength")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # (b) residuals per degree — shows WHY degree 2 is needed (the arch)
    ax = axes[1, 0]
    for d in range(1, min(max_degree, len(pts) - 1) + 1):
        f = fit_time_to_wavelength(t, wl, degree=d)
        resid = (wl - np.polyval(f["coeffs"], t)) * 1000
        ax.plot(wl, resid, "o-", ms=5, lw=1.0,
                label=f"degree {d}: {f['residual_pm_rms']:.1f} pm rms")
    ax.axhline(0, color="k", lw=0.8)
    ax.axhspan(-32, 32, color="tab:green", alpha=0.10,
               label="TOF1550 setting error (±32 pm spec)")
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("fit residual (pm)")
    ax.set_title("residuals — degree 1 leaves a systematic arch")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)

    # (c) what the uncorrected axis costs
    ax = axes[0, 1]
    naive_err = (start + speed * t - wl) * 1000
    ax.plot(wl, naive_err, "o-", ms=5, color="tab:red")
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("naive-axis error (pm)")
    ax.set_title(f"error of the uncorrected linear axis\n"
                 f"(up to {np.abs(naive_err).max():.0f} pm)")
    ax.grid(alpha=0.3)

    # (d) residual vs degree
    ax = axes[1, 1]
    degs, rms = [], []
    for d in range(1, min(max_degree + 1, len(pts) - 1) + 1):
        f = fit_time_to_wavelength(t, wl, degree=d)
        degs.append(d)
        rms.append(f["residual_pm_rms"])
    ax.bar([str(d) for d in degs], rms, color="tab:blue", alpha=0.8)
    for d, r in zip(range(len(degs)), rms):
        ax.annotate(f"{r:.1f}", (d, r), ha="center", va="bottom", fontsize=8)
    ax.axhline(32, color="tab:green", ls="--", lw=1.0,
               label="filter setting error (32 pm)")
    ax.set_xlabel("polynomial degree")
    ax.set_ylabel("residual (pm rms)")
    ax.set_yscale("log")
    ax.set_title("how much fit complexity is justified")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, axis="y")

    # (e/f) repeatability, if a repeat CSV was given
    if repeat_csv is not None:
        rep = np.genfromtxt(repeat_csv, delimiter=",", names=True)
        el = np.atleast_1d(rep["elapsed_s"]) / 60.0
        ts = np.atleast_1d(rep["time_s"])
        pm_err = (ts - ts.mean()) * speed * 1000
        ax = axes[0, 2]
        ax.plot(el, pm_err, "o-", ms=6)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel("elapsed time (min)")
        ax.set_ylabel("peak position deviation (pm)")
        ax.set_title(f"repeatability of the marker\n"
                     f"({pm_err.std(ddof=1) if pm_err.size > 1 else 0:.1f} pm rms, "
                     f"{pm_err.max() - pm_err.min():.1f} pm p-p)")
        ax.grid(alpha=0.3)

        ax = axes[1, 2]
        ax.plot(el, np.atleast_1d(rep["amplitude_mw"]) * 1000, "o-", ms=6,
                label="peak height (µW)")
        ax.plot(el, np.atleast_1d(rep["width_s"]) * 1000, "s-", ms=6,
                label="peak width (ms)")
        ax.set_xlabel("elapsed time (min)")
        ax.set_title("marker shape stability")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, f"calibration_{tag + '_' if tag else ''}{stamp()}.png")
    fig.savefig(png, dpi=130)
    print(f"[out] {png}")
    try:
        plt.show()
    except Exception:
        pass
    return png


def plot_raw(res, out_dir, tag="", cal=None):
    import matplotlib
    matplotlib.use("Agg") if os.environ.get("LINA_NO_DISPLAY") else None
    import matplotlib.pyplot as plt

    rate = res["rate_hz"]
    n = res["n_pts"]
    t = np.arange(n) / rate
    fig, axes = plt.subplots(2 if cal else 1, 1, figsize=(11, 7 if cal else 4.5))
    axes = np.atleast_1d(axes)

    ax = axes[0]
    for ch, tr in sorted(res["traces_mw"].items()):
        ax.plot(t, tr, lw=0.8, label=f"Ch{ch}")
        peaks, base, _ = find_buffer_peaks(tr, rate_hz=rate)
        for k, p in enumerate(peaks):
            ax.axvline(p["index"] / rate, ls="--", lw=0.8, color="r" if k == 0 else "gray")
            ax.annotate(f"{'fwd' if k == 0 else 'return'} @ {p['index'] / rate:.3f} s",
                        (p["index"] / rate, tr.max()), fontsize=7,
                        rotation=90, va="top")
    ax.axvspan(0, min(res["sweep_s"], t[-1]), color="tab:green", alpha=0.08,
               label="nominal sweep duration (from t=0)")
    ax.set_xlabel("time in capture (s)")
    ax.set_ylabel("power (mW)")
    ax.set_title(f"RAW capture — {res['actual_start_nm']:.2f} -> "
                 f"{res['actual_stop_nm']:.2f} nm @ "
                 f"{res['actual_speed_nm_s']:.2f} nm/s, {rate} Hz")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    if cal is not None:
        wl = wavelength_axis(cal, n, rate)
        ax2 = axes[1]
        for ch, tr in sorted(res["traces_mw"].items()):
            ax2.plot(wl, tr, lw=0.8, label=f"Ch{ch}")
        ax2.axvline(res["actual_start_nm"], color="k", ls=":", lw=0.8)
        ax2.axvline(res["actual_stop_nm"], color="k", ls=":", lw=0.8)
        ax2.set_xlabel("calibrated wavelength (nm)")
        ax2.set_ylabel("power (mW)")
        ax2.set_title("calibrated axis (dotted = actual sweep range)")
        ax2.legend(fontsize=8)
        ax2.grid(alpha=0.3)

    fig.tight_layout()
    os.makedirs(out_dir, exist_ok=True)
    png = os.path.join(out_dir, f"sweep_{tag + '_' if tag else ''}{stamp()}.png")
    fig.savefig(png, dpi=130)
    print(f"[out] {png}")
    try:
        plt.show()
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Subcommands
# ---------------------------------------------------------------------------

def cmd_sweep(a):
    laser = pm = tof = None
    try:
        laser, pm, tof = open_devices(a.exfo, a.coredaq,
                                      a.filter if a.filter_nm is not None else None)
        if a.filter_nm is not None:
            print(f"[tof] filter set to {tof.set_wavelength(a.filter_nm):.3f} nm")
        res = run_sweep(laser, pm, start_nm=a.start, stop_nm=a.stop,
                        speed_nm_s=a.speed, channels=a.channels,
                        rate_hz=a.rate, n_samples=a.samples, pad_s=a.pad,
                        park_settle_s=a.park_settle, restore_wl=not a.no_restore)
    finally:
        close_devices(laser, pm, tof)

    for ch, tr in sorted(res["traces_mw"].items()):
        peaks, base, amp = find_buffer_peaks(tr, rate_hz=res["rate_hz"])
        print(f"[peak] Ch{ch}: baseline {base:.6g} mW, amplitude {amp:.6g} mW, "
              f"{len(peaks)} peak(s)")
        for k, p in enumerate(peaks):
            print(f"        #{k}: index {p['index']:.1f} "
                  f"({p['index'] / res['rate_hz']:.4f} s), "
                  f"width {p['width_s'] * 1000:.1f} ms, "
                  f"amp {p['amplitude_mw']:.6g} mW")
    cal = json.load(open(a.cal)) if a.cal else None
    if cal is not None:
        # Independent check of the calibration: the filter sits at a known
        # wavelength, so the calibrated axis must put the peak there — even
        # when this capture's sample rate differs from the calibration's,
        # because the fit is in seconds, not samples.
        wl = wavelength_axis(cal, res["n_pts"], res["rate_hz"])
        print(f"\n[cal] applying {a.cal} (fit degree {cal['fit']['degree']}, "
              f"measured at {cal['sweep']['rate_hz']} Hz)")
        for ch, tr in sorted(res["traces_mw"].items()):
            peaks, _, _ = find_buffer_peaks(tr, rate_hz=res["rate_hz"])
            fwd = select_forward_peak(peaks)
            if fwd is None:
                continue
            lo, hi = fwd["start"], fwd["stop"] + 1
            w = np.clip(tr[lo:hi] - np.median(tr), 0, None)
            peak_nm = float(np.sum(wl[lo:hi] * w) / np.sum(w))
            naive = np.linspace(res["actual_start_nm"], res["actual_stop_nm"],
                                res["n_pts"])
            naive_nm = float(np.sum(naive[lo:hi] * w) / np.sum(w))
            print(f"[cal] Ch{ch}: peak at {peak_nm:.4f} nm calibrated"
                  + (f" ({(peak_nm - a.filter_nm) * 1000:+.0f} pm vs filter "
                     f"{a.filter_nm:.3f} nm)" if a.filter_nm is not None else "")
                  + f"  |  naive linear axis would say {naive_nm:.4f} nm"
                  + (f" ({(naive_nm - a.filter_nm) * 1000:+.0f} pm)"
                     if a.filter_nm is not None else ""))
    save_raw(res, a.out, a.tag, write_csv=a.csv)
    if not a.no_plot:
        plot_raw(res, a.out, a.tag, cal)


def cmd_filter_cw(a):
    laser = pm = tof = None
    try:
        laser, pm, tof = open_devices(a.exfo, a.coredaq,
                                      a.filter if a.filter_nm is not None else None)
        center = a.center
        if a.filter_nm is not None:
            rb = tof.set_wavelength(a.filter_nm)
            print(f"[tof] filter set to {a.filter_nm:.3f} nm (readback {rb:.3f})")
            if center is None:
                center = rb
        if center is None:
            sys.exit("give --center (CW scan centre) or --filter-nm")
        res = run_filter_cw(laser, pm, center_nm=center, span_nm=a.span,
                            step_nm=a.step, channels=a.channels,
                            settle_s=a.settle, n_avg=a.avg)
        if not a.no_restore:
            laser.set_wavelength(center)
        print(f"\n[cw] best channel Ch{res['best_channel']}; "
              f"true peak {res['peak_nm']:.4f} nm"
              + (f", filter was set to {a.filter_nm:.3f} nm -> offset "
                 f"{(res['peak_nm'] - a.filter_nm) * 1000:+.0f} pm"
                 if a.filter_nm is not None else ""))
    finally:
        close_devices(laser, pm, tof)

    os.makedirs(a.out, exist_ok=True)
    csv = os.path.join(a.out, f"filter_cw_{a.tag + '_' if a.tag else ''}{stamp()}.csv")
    chs = sorted(res["power_mw"])
    np.savetxt(csv, np.column_stack([res["wavelengths_nm"]] +
                                     [res["power_mw"][ch] for ch in chs]),
               delimiter=",",
               header="wavelength_nm," + ",".join(f"ch{ch}_mw" for ch in chs),
               comments="")
    print(f"[out] {csv}")
    if not a.no_plot:
        plot_filter_cw(res, a.out, a.tag, a.filter_nm)


def cmd_calibrate(a):
    positions = [float(p) for p in a.filter_positions.split(",") if p.strip()]
    if len(positions) < 2:
        sys.exit("calibrate needs at least 2 filter positions")

    laser = pm = tof = None
    points = []
    meta = None
    try:
        laser, pm, tof = open_devices(a.exfo, a.coredaq, a.filter)
        for pos in positions:
            print(f"\n================ filter {pos:.3f} nm ================")
            rb = tof.set_wavelength(pos)
            print(f"[tof] set {pos:.3f} nm, readback {rb:.3f} nm")

            # Stage A — ground truth: where is the passband REALLY? The filter's
            # own setting error (<32 pm spec) never enters the calibration.
            true_nm = pos
            cw = None
            if not a.trust_filter:
                cw = run_filter_cw(laser, pm, center_nm=rb, span_nm=a.cw_span,
                                   step_nm=a.cw_step, channels=a.channels,
                                   settle_s=a.settle, n_avg=a.avg)
                true_nm = cw["peak_nm"]
                print(f"[cal] filter set {pos:.3f} nm -> true peak "
                      f"{true_nm:.4f} nm ({(true_nm - pos) * 1000:+.0f} pm)")
            if not (min(a.start, a.stop) <= true_nm <= max(a.start, a.stop)):
                print(f"[warn] {true_nm:.3f} nm is outside the sweep range — "
                      f"cannot be located in the buffer, skipping")
                continue

            # Stage B — where does that wavelength land in the capture?
            res = run_sweep(laser, pm, start_nm=a.start, stop_nm=a.stop,
                            speed_nm_s=a.speed, channels=a.channels,
                            rate_hz=a.rate, n_samples=a.samples, pad_s=a.pad,
                            park_settle_s=a.park_settle, restore_wl=False)
            meta = {k: v for k, v in res.items() if k != "traces_mw"}
            trace = res["traces_mw"][a.channels[0]]
            peaks, base, amp = find_buffer_peaks(trace, rate_hz=res["rate_hz"])
            if not peaks:
                print(f"[warn] no filter peak in the capture "
                      f"(baseline {base:.4g} mW, amplitude {amp:.4g} mW) — skipping")
                save_raw(res, a.out, f"cal_{pos:.0f}_nopeak")
                continue
            fwd = select_forward_peak(peaks)
            print(f"[cal] true {true_nm:.4f} nm at index {fwd['index']:.1f} "
                  f"({fwd['index'] / res['rate_hz']:.4f} s), width "
                  f"{fwd['width_s'] * 1000:.1f} ms"
                  + (f"; {len(peaks) - 1} other peak(s) = return slew inside the "
                     f"capture window" if len(peaks) > 1 else ""))
            points.append({"filter_set_nm": pos, "filter_readback_nm": rb,
                           "true_nm": true_nm,
                           "fwhm_nm": (cw["fwhm_nm"] if cw else None),
                           "index": fwd["index"],
                           "time_s": fwd["index"] / res["rate_hz"],
                           "n_peaks": len(peaks), "timestamp": stamp()})
            save_raw(res, a.out, f"cal_{pos:.0f}")
    except KeyboardInterrupt:
        print("\n[cal] aborted by user")
    finally:
        close_devices(laser, pm, tof)

    if len(points) < 2:
        sys.exit("not enough usable calibration points")

    fit = fit_time_to_wavelength([p["time_s"] for p in points],
                                 [p["true_nm"] for p in points],
                                 degree=a.degree)
    derived = report_fit(fit, meta)
    cal = {"created": stamp(), "exfo": a.exfo, "coredaq": a.coredaq,
           "channel": a.channels[0], "sweep": meta, "points": points,
           "fit": fit, "derived": derived}
    os.makedirs(a.out, exist_ok=True)
    path = os.path.join(a.out, a.cal_name)
    with open(path, "w") as fp:
        json.dump(cal, fp, indent=2)
    print(f"[out] {path}")
    # Publish into the shared store the GUI looks up by sweep configuration,
    # unless asked not to (e.g. an exploratory run that should not become the
    # calibration the Lina window starts using).
    if not a.no_publish:
        published = save_calibration(cal)
        print(f"[out] published for the GUI: {published}")
        print(f"      {describe_cal(cal)}")
    if not a.no_plot:
        plot_calibration(cal, a.out, a.tag)


def cmd_repeat(a):
    """Same sweep N times with the filter parked — is the mapping stable in time?

    The filter peak sits at a FIXED wavelength, so its sample index must come
    out the same every run. Any scatter is trigger/sweep-start jitter, and any
    monotonic drift means a calibration measured once will go stale.
    """
    laser = pm = tof = None
    rows = []
    try:
        laser, pm, tof = open_devices(a.exfo, a.coredaq, a.filter)
        rb = tof.set_wavelength(a.filter_nm)
        print(f"[tof] filter parked at {a.filter_nm:.3f} nm (readback {rb:.3f})")
        t_start = time.time()
        for k in range(a.runs):
            print(f"\n---------------- run {k + 1}/{a.runs} ----------------")
            res = run_sweep(laser, pm, start_nm=a.start, stop_nm=a.stop,
                            speed_nm_s=a.speed, channels=a.channels,
                            rate_hz=a.rate, n_samples=a.samples, pad_s=a.pad,
                            park_settle_s=a.park_settle, restore_wl=False)
            trace = res["traces_mw"][a.channels[0]]
            peaks, base, amp = find_buffer_peaks(trace, rate_hz=res["rate_hz"])
            if not peaks:
                print(f"[warn] run {k + 1}: no peak found "
                      f"(baseline {base:.4g} mW, amplitude {amp:.4g} mW)")
                continue
            fwd = select_forward_peak(peaks)
            idx = fwd["index"]
            rows.append({"run": k + 1, "elapsed_s": time.time() - t_start,
                         "index": idx, "time_s": idx / res["rate_hz"],
                         "width_s": fwd["width_s"],
                         "amplitude_mw": fwd["amplitude_mw"],
                         "n_peaks": len(peaks),
                         "blocked_s": res["blocked_s"]})
            print(f"[rep] run {k + 1}: peak index {idx:.1f} "
                  f"({idx / res['rate_hz']:.4f} s), amp "
                  f"{fwd['amplitude_mw']:.4g} mW, width "
                  f"{fwd['width_s'] * 1000:.1f} ms, {len(peaks)} peak(s)")
            if a.save_raw:
                save_raw(res, a.out, f"rep{k + 1}", write_csv=a.csv)
            if a.interval_s and k < a.runs - 1:
                print(f"[rep] waiting {a.interval_s:.0f} s before the next run…")
                time.sleep(a.interval_s)
    except KeyboardInterrupt:
        print("\n[rep] aborted by user")
    finally:
        close_devices(laser, pm, tof)

    if not rows:
        sys.exit("no usable runs")
    idx = np.array([r["index"] for r in rows])
    rate = a.rate
    print("\n=== repeatability ===")
    print(f"  runs                 : {len(rows)} over "
          f"{rows[-1]['elapsed_s'] / 60:.1f} min")
    print(f"  peak index mean      : {idx.mean():.1f} ({idx.mean() / rate:.4f} s)")
    print(f"  peak index spread    : {idx.std(ddof=1) if idx.size > 1 else 0:.1f} "
          f"samples rms, {idx.max() - idx.min():.1f} peak-to-peak")
    # Convert index scatter into the wavelength error it causes.
    nm_per_sample = a.speed / rate
    print(f"  -> wavelength error  : {idx.std(ddof=1) * nm_per_sample * 1000 if idx.size > 1 else 0:.1f} pm rms, "
          f"{(idx.max() - idx.min()) * nm_per_sample * 1000:.1f} pm peak-to-peak "
          f"(at {a.speed} nm/s, {rate} Hz -> {nm_per_sample * 1000:.2f} pm/sample)")
    if idx.size >= 3:
        el = np.array([r["elapsed_s"] for r in rows])
        slope = np.polyfit(el, idx, 1)[0]
        print(f"  drift                : {slope * nm_per_sample * 1000 * 60:+.1f} pm/min "
              f"({slope:+.3f} samples/s of elapsed time)")
    os.makedirs(a.out, exist_ok=True)
    csv = os.path.join(a.out, f"repeat_{a.tag + '_' if a.tag else ''}{stamp()}.csv")
    keys = list(rows[0])
    np.savetxt(csv, np.array([[r[k] for k in keys] for r in rows]),
               delimiter=",", header=",".join(keys), comments="")
    print(f"[out] {csv}")


def cmd_refit(a):
    """Re-fit an existing calibration at a different polynomial degree.

    Choosing the degree needs no new hardware runs — the (time, wavelength)
    points are already in the JSON.
    """
    cal = json.load(open(a.cal))
    pts = cal["points"]
    t = [p["time_s"] for p in pts]
    wl = [p["true_nm"] for p in pts]
    for d in range(1, min(a.max_degree, len(pts) - 1) + 1):
        f = fit_time_to_wavelength(t, wl, degree=d)
        print(f"  degree {d}: {f['residual_pm_rms']:6.1f} pm rms, "
              f"{f['residual_pm_max']:6.1f} pm max")
    cal["fit"] = fit_time_to_wavelength(t, wl, degree=a.degree)
    cal["derived"] = report_fit(cal["fit"], cal["sweep"])
    cal["refit"] = stamp()
    out = a.out_cal or a.cal
    with open(out, "w") as fp:
        json.dump(cal, fp, indent=2)
    print(f"[out] {out}")


def cmd_apply(a):
    res = load_raw(a.raw)
    cal = json.load(open(a.cal))
    plot_raw(res, a.out, a.tag, cal)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exfo", default="exfo-1")
    p.add_argument("--coredaq", default="coredaq-1")
    p.add_argument("--filter", default="tof1550-1",
                   help="TOF1550 device name in TunableFilter_config.json")
    p.add_argument("--channels", default="1",
                   help="comma-separated 1-based coreDAQ channels (default 1)")
    p.add_argument("--out", default="lina_cal_data", help="output directory")
    p.add_argument("--tag", default="", help="filename tag")
    p.add_argument("--no-restore", action="store_true",
                   help="leave the laser at the sweep start instead of restoring")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_sweep_args(sp):
        sp.add_argument("--start", type=float, default=1520.0)
        sp.add_argument("--stop", type=float, default=1570.0)
        sp.add_argument("--speed", type=float, default=5.0, help="nm/s")
        sp.add_argument("--rate", type=int, default=100_000,
                        help="coreDAQ sample rate (Hz, default 100k)")
        sp.add_argument("--samples", type=int, default=None,
                        help="frames/channel (default: (sweep+pad) x rate)")
        sp.add_argument("--pad", type=float, default=2.0,
                        help="extra seconds of buffer beyond the sweep (default 2)")
        sp.add_argument("--park-settle", type=float, default=5.0,
                        help="seconds to settle at the start wavelength before arming")

    sp = sub.add_parser("sweep", help="one raw triggered capture")
    add_sweep_args(sp)
    sp.add_argument("--filter-nm", type=float, default=None,
                    help="park the TOF1550 here first (omit = leave it alone)")
    sp.add_argument("--cal", default=None, help="calibration JSON to overlay")
    sp.add_argument("--csv", action="store_true",
                    help="also write the capture as CSV text (~25 MB per "
                         "1.2 M-sample capture; NPZ is always written)")
    sp.add_argument("--no-plot", action="store_true")
    sp.set_defaults(func=cmd_sweep)

    sp = sub.add_parser("filter-cw", help="CW step scan of the filter passband")
    sp.add_argument("--filter-nm", type=float, default=None,
                    help="set the TOF1550 here first, and centre the scan on it")
    sp.add_argument("--center", type=float, default=None,
                    help="scan centre (default: the filter readback)")
    sp.add_argument("--span", type=float, default=2.0)
    sp.add_argument("--step", type=float, default=0.025)
    sp.add_argument("--settle", type=float, default=0.6)
    sp.add_argument("--avg", type=int, default=4)
    sp.add_argument("--no-plot", action="store_true")
    sp.set_defaults(func=cmd_filter_cw)

    sp = sub.add_parser("repeat", help="N identical sweeps — stability in time")
    add_sweep_args(sp)
    sp.add_argument("--filter-nm", type=float, required=True,
                    help="wavelength to park the filter at (the fixed marker)")
    sp.add_argument("--runs", type=int, default=5)
    sp.add_argument("--interval-s", type=float, default=0.0,
                    help="wait between runs, to probe longer-term drift")
    sp.add_argument("--save-raw", action="store_true")
    sp.add_argument("--csv", action="store_true",
                    help="write CSV text alongside the NPZ for saved captures")
    sp.set_defaults(func=cmd_repeat)

    sp = sub.add_parser("calibrate", help="multi-position filter calibration")
    add_sweep_args(sp)
    sp.add_argument("--filter-positions", required=True,
                    help="comma-separated filter wavelengths, e.g. 1535,1545,1555")
    sp.add_argument("--trust-filter", action="store_true",
                    help="use the filter's own setting as truth (spec <32 pm) and "
                         "skip the CW ground-truth scan at each point — much faster")
    sp.add_argument("--cw-span", type=float, default=2.0,
                    help="CW scan span around each filter position (nm)")
    sp.add_argument("--cw-step", type=float, default=0.025)
    sp.add_argument("--settle", type=float, default=0.6)
    sp.add_argument("--avg", type=int, default=4)
    sp.add_argument("--degree", type=int, default=2,
                    help="time->nm polynomial degree (2 fits the measured "
                         "sweep non-uniformity to ~5 pm; 1 leaves ~30 pm)")
    sp.add_argument("--cal-name", default="lina_wl_cal.json")
    sp.add_argument("--no-publish", action="store_true",
                    help=f"don't copy the result into {CAL_DIR} (where the "
                         f"Lina GUI window looks it up)")
    sp.add_argument("--no-plot", action="store_true")
    sp.set_defaults(func=cmd_calibrate)

    sp = sub.add_parser("plot-cal", help="plot a saved calibration, no hardware")
    sp.add_argument("--cal", required=True)
    sp.add_argument("--repeat-csv", default=None,
                    help="a repeat-mode CSV, to add the stability panels")
    sp.add_argument("--max-degree", type=int, default=3)
    sp.set_defaults(func=lambda a: plot_calibration(
        json.load(open(a.cal)), a.out, a.tag, a.repeat_csv, a.max_degree))

    sp = sub.add_parser("refit", help="re-fit a saved calibration, no hardware")
    sp.add_argument("--cal", required=True)
    sp.add_argument("--degree", type=int, default=2)
    sp.add_argument("--max-degree", type=int, default=4,
                    help="also report residuals up to this degree")
    sp.add_argument("--out-cal", default=None, help="write here (default: in place)")
    sp.set_defaults(func=cmd_refit)

    sp = sub.add_parser("apply", help="re-plot a saved raw sweep with a calibration")
    sp.add_argument("--raw", required=True)
    sp.add_argument("--cal", required=True)
    sp.add_argument("--no-plot", action="store_true")
    sp.set_defaults(func=cmd_apply)

    a = p.parse_args()
    a.channels = [int(c) for c in str(a.channels).split(",") if c.strip()]
    a.func(a)


if __name__ == "__main__":
    main()
