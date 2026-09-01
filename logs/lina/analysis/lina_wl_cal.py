"""
Lina sweep wavelength-axis calibration — shared by the measurement script
(``scripts/lina_sweep_test.py``) and the GUI (``LabGUI/dashboards/instruments/
lina_window.py``).

WHAT IS CALIBRATED
------------------
A LINEAR-variant coreDAQ capture is started by ONE trigger edge from the EXFO
and then free-runs at its own sample rate, so every sample is tied to a TIME,
not to a wavelength. Mapping those samples onto wavelengths by assuming the
buffer spans exactly [start, stop] nm is wrong twice over:

  * the armed buffer is longer than the sweep (it has to be, or the end is
    lost), so the sweep occupies only the first ~82 % of it — measured: the
    sweep runs t = 0.004 .. 9.889 s of a 12.00 s buffer, and the laser's
    return slew to its parked wavelength shows up at ~10.12 s;
  * the EXFO does not sweep at exactly the commanded speed (measured ~+1..2 %)
    and is slightly non-uniform in time, so even a correctly cropped linear
    axis drifts by ~0.5 nm over a 50 nm sweep.

So the calibration is a polynomial  lambda(t)  in SECONDS from the trigger
edge, measured by parking a narrow tunable filter (Thorlabs TOF1550, 0.21 nm
passband) at several known wavelengths and recording at what time its
transmission peak appears in the capture. Being in seconds rather than sample
index, it applies unchanged at any sample rate: a fit measured at 5 kHz
reproduced a 1550.000 nm marker to -3 pm in a 100 kHz capture.

Measured accuracy: degree 2 -> ~3 pm rms (degree 1 leaves a systematic ~27 pm
arch; degree 3 gains nothing).

SCOPE — IMPORTANT
-----------------
The coefficients are absolute: they encode one specific sweep configuration
(start wavelength, stop, speed). They are NOT valid for a different range or
speed, so calibrations are stored one file per configuration and looked up by
matching those parameters. ``find_calibration`` returns None rather than
guessing when no file matches; run the calibration again for that
configuration (see the script's ``calibrate`` subcommand).
"""

from __future__ import annotations

import glob
import json
import os
from datetime import datetime

import numpy as np

# One file per sweep configuration, next to the other instrument calibrations.
# lina/calibrations/ — one file per sweep configuration, beside the scripts
# and data that produced them.
CAL_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'calibrations')

# Tolerances for deciding that a stored calibration describes the sweep about
# to be run. Tight on purpose — a calibration from a different configuration is
# worse than no calibration, because it is silently wrong.
TOL_NM = 0.05
TOL_SPEED_REL = 0.02


def cal_filename(start_nm: float, stop_nm: float, speed_nm_s: float) -> str:
    return (f"lina_wl_cal_{start_nm:.1f}-{stop_nm:.1f}nm"
            f"_{speed_nm_s:.2f}nms.json")


def fit_time_to_wavelength(times_s, wavelengths_nm, degree: int = 2) -> dict:
    """Fit capture time (s from the trigger edge) -> wavelength (nm)."""
    t = np.asarray(times_s, dtype=float)
    wl = np.asarray(wavelengths_nm, dtype=float)
    degree = int(min(degree, t.size - 1))
    coeffs = np.polyfit(t, wl, degree)
    resid = wl - np.polyval(coeffs, t)
    return {"coeffs": [float(c) for c in coeffs], "degree": degree,
            "units": "nm vs seconds from trigger edge",
            "n_points": int(t.size),
            "residual_pm_rms": float(np.sqrt(np.mean(resid ** 2)) * 1000),
            "residual_pm_max": float(np.max(np.abs(resid)) * 1000)}


def save_calibration(cal: dict, cal_dir: str = None) -> str:
    """Write a calibration under its configuration-derived name."""
    cal_dir = cal_dir or CAL_DIR
    os.makedirs(cal_dir, exist_ok=True)
    sweep = cal["sweep"]
    path = os.path.join(cal_dir, cal_filename(sweep["requested_start_nm"],
                                              sweep["requested_stop_nm"],
                                              sweep["requested_speed_nm_s"]))
    cal.setdefault("saved", datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    with open(path, "w") as fp:
        json.dump(cal, fp, indent=2)
    return path


def load_calibration(path: str) -> dict:
    with open(path) as fp:
        return json.load(fp)


def find_calibration(start_nm: float, stop_nm: float, speed_nm_s: float,
                     cal_dir: str = None):
    """Newest stored calibration matching this sweep configuration, or None.

    Matching is on the REQUESTED sweep parameters (what the user asked for),
    within :data:`TOL_NM` / :data:`TOL_SPEED_REL`.
    """
    cal_dir = cal_dir or CAL_DIR
    best, best_time = None, None
    for path in glob.glob(os.path.join(cal_dir, "lina_wl_cal_*.json")):
        try:
            cal = load_calibration(path)
            sweep = cal["sweep"]
            if (abs(float(sweep["requested_start_nm"]) - start_nm) <= TOL_NM
                    and abs(float(sweep["requested_stop_nm"]) - stop_nm) <= TOL_NM
                    and abs(float(sweep["requested_speed_nm_s"]) - speed_nm_s)
                    <= TOL_SPEED_REL * max(abs(speed_nm_s), 1e-9)):
                mtime = os.path.getmtime(path)
                if best_time is None or mtime > best_time:
                    best, best_time = cal, mtime
                    cal["_path"] = path
        except Exception:
            continue
    return best


def wavelength_axis(cal: dict, n_pts: int, rate_hz: float) -> np.ndarray:
    """Wavelength for every sample of a capture taken at ``rate_hz``.

    Rate-independent: the fit is in seconds, so the calibration may have been
    measured at a different sample rate than the capture it is applied to.
    """
    t = np.arange(int(n_pts)) / float(rate_hz)
    return np.polyval(cal["fit"]["coeffs"], t)


def sweep_window(cal: dict, n_pts: int, rate_hz: float):
    """``(lo, hi)`` sample indices of the part of the buffer that is the sweep.

    Everything before ``lo`` is trigger dead time and everything from ``hi`` on
    is the post-sweep tail — the laser slewing back to its parked wavelength,
    which crosses the same wavelengths backwards and must not be plotted as
    sweep data.
    """
    derived = cal.get("derived", {})
    t0 = derived.get("buffer_time_s") or 0.0
    t1 = derived.get("sweep_end_s")
    lo = max(int(round(max(t0, 0.0) * rate_hz)), 0)
    hi = int(round(t1 * rate_hz)) if t1 else int(n_pts)
    return lo, min(max(hi, lo + 1), int(n_pts))


def apply_calibration(cal: dict, traces: dict, rate_hz: float,
                      crop_to_sweep: bool = True):
    """Map a raw capture onto wavelengths.

    Args:
        cal: calibration dict (from :func:`find_calibration` / :func:`load_calibration`).
        traces: ``{channel: 1-D array}`` raw capture, all the same length.
        rate_hz: sample rate the capture was taken at.
        crop_to_sweep: drop the trigger dead time and the post-sweep return
            slew (see :func:`sweep_window`).

    Returns:
        ``(wavelengths_nm, {channel: array})`` — both cropped consistently.
    """
    n_pts = min((len(v) for v in traces.values()), default=0)
    if n_pts == 0:
        return np.empty(0), {ch: np.empty(0) for ch in traces}
    wl = wavelength_axis(cal, n_pts, rate_hz)
    out = {ch: np.asarray(v, dtype=float)[:n_pts] for ch, v in traces.items()}
    if crop_to_sweep:
        lo, hi = sweep_window(cal, n_pts, rate_hz)
        wl = wl[lo:hi]
        out = {ch: v[lo:hi] for ch, v in out.items()}
    return wl, out


def describe(cal: dict) -> str:
    """One-line summary for a status bar / plot label."""
    fit = cal.get("fit", {})
    sweep = cal.get("sweep", {})
    return (f"λ(t) deg {fit.get('degree', '?')}, "
            f"{fit.get('residual_pm_rms', float('nan')):.0f} pm rms, "
            f"{sweep.get('requested_start_nm', '?')}–"
            f"{sweep.get('requested_stop_nm', '?')} nm @ "
            f"{sweep.get('requested_speed_nm_s', '?')} nm/s")
