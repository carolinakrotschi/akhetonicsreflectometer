"""
Lina sweep wavelength-axis calibration — shared by the measurement script
(``scripts/lina_sweep_test.py``) and the GUI (``LabGUI/dashboards/instruments/
lina_window.py``).

WHAT IS CALIBRATED
------------------
A continuous coreDAQ capture (LINEAR or LOG frontend) is started by ONE
trigger edge from the EXFO and then free-runs at its own sample rate, so every sample is tied to a TIME,
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


# ---------------------------------------------------------------------------
# Aux-referenced axis — for sweep configurations that have no lambda(t) file
# ---------------------------------------------------------------------------
#
# A lambda(t) calibration describes ONE configuration and is measured against a
# tunable filter, so a sweep outside the filter's range (TOF1550: 1527-1567 nm)
# can never get one. Those sweeps used to fall back to the linear axis, which
# is wrong twice over (see the module docstring).
#
# The aux MZI answers both objections without a calibration run, because it
# MEASURES the sweep instead of predicting it:
#
#   * its unwrapped fringe phase is proportional to optical frequency, so the
#     phase says how far through the sweep — in FREQUENCY — each sample is;
#   * its fringes exist only while the laser is moving, so where they stop is
#     where the sweep stops. That is the crop, measured rather than assumed.
#
# The absolute scale comes from the commanded endpoints. That split is
# deliberate. An earlier design took the scale from a STORED aux delay, which
# would have been wrong: measured 2026-09-15, captures on this bench imply
# tau_aux 20.64 ns against the 20.39 ns measured on 2026-09-01 — the aux arm
# had been re-patched and nothing noticed. The endpoints are an instrument
# property and cannot go stale that way, so tau_aux is computed from every
# sweep and reported as a CHECK, never used to set the axis.
#
# Validated on two fibres swept over 1520-1570 nm and over 1505-1625 nm (2.4x
# wider): the OFDR peaks agree between the two windows to under 0.04 mm on
# 1-5 m, i.e. 0.001-0.003 %.

C_VAC = 299_792_458.0        # m/s

# The reference/aux MZI pair on this setup — same convention as the OFDR
# evaluation in lina_window.py and lina/scripts/lina_ofdr_test.py.
AUX_CHANNELS = (2, 4)

# Refusal thresholds. An aux axis that cannot be trusted must fall back to the
# linear axis and say why, rather than quietly produce a wrong one.
AUX_MIN_FRINGES = 1000.0      # fewer than this is noise, not a sweep
AUX_MAX_NONMONOTONIC = 0.01   # 1 %, matching the OFDR resampling guard
AUX_MAX_NYQUIST_FRAC = 0.45   # above this the measured line is unusable

# A cleanly aliased fringe CANNOT be detected from the capture alone: folding a
# pure tone gives another pure tone, below Nyquist, with a smooth monotonic
# phase -- every intrinsic test passes. What does catch it is the aux delay the
# sweep implies: aliasing throws it far off the bench's known value. So a gross
# mismatch against the stored reference is a refusal, and a small one is only a
# warning (the bench really can be re-patched, and the axis stays usable).
AUX_TAU_DEVIATION_WARN = 0.002
AUX_TAU_DEVIATION_REFUSE = 0.10

# Sweep-window detection: the fringe envelope, smoothed over this many samples,
# must stay above this fraction of its median.
#
# Both numbers are measured, not guessed. A boxcar smear puts the detected edge
# about half its own width past the true one, and samples of dead tail inside
# the window are NOT harmless: the Hilbert phase of high-passed noise keeps
# advancing, so they inject fringes that never happened. At 20001 that was
# +5628 samples of tail and a 0.69 % error in the implied aux delay -- enough
# to trip the staleness check on every sweep. At 2001 it is +126 samples and
# 0.004 %.
#
# The threshold stays at 0.25 rather than the 0.5 that would centre the edge
# exactly: over a 120 nm sweep the real fringe envelope varies by ~50 % end to
# end, and at 0.5 the mask breaks up INSIDE the sweep, which costs far more
# than a slightly late edge (measured on a 1505-1625 nm capture: 19.89 ns
# against 20.65 ns for every other setting).
AUX_ENVELOPE_WIN = 2001
AUX_ENVELOPE_FRAC = 0.25

# Stored aux delay — a reference for the staleness check only.
AUX_REF_FILE = os.path.join(CAL_DIR, 'aux_mzi.json')


class AuxAxisUnavailable(RuntimeError):
    """The aux fringe cannot carry a wavelength axis.

    Carries the reason, so a caller that falls back to the linear axis can say
    WHY instead of silently degrading.
    """


def _moving_mean(x, win: int):
    """Centred moving average in O(n).

    np.convolve would be O(n*win); at 2.5 M samples and a 20 k window that is
    5e10 operations and unusable inside a GUI worker.
    """
    x = np.asarray(x, dtype=float)
    win = int(max(1, win))
    c = np.cumsum(np.concatenate(([0.0], x)))
    half = win // 2
    i = np.arange(x.size)
    lo = np.clip(i - half, 0, x.size)
    hi = np.clip(i + half + 1, 0, x.size)
    return (c[hi] - c[lo]) / (hi - lo)


def aux_sweep_window(aux_ac, dc=None, win: int = AUX_ENVELOPE_WIN,
                     frac: float = AUX_ENVELOPE_FRAC):
    """Longest stretch of the capture where the aux is actually fringing.

    Fringes exist only while the laser moves, so this is the sweep. Measured
    on a 1505-1625 nm capture: the envelope drops by a factor ~500 within a
    fraction of a nm at t = 24.07 s of a 25.20 s buffer — the laser had
    finished and the buffer was still recording.

    ``dc`` (the raw, un-high-passed aux level) only decides what to CALL the
    tail: an abrupt loss of fringes with the light still there is the end of
    the sweep, whereas losing both is the laser running out of power at a band
    edge. The window is the same either way; the distinction goes in the note,
    because the second case means the sweep was cut short and the commanded
    endpoints no longer describe the cropped region.

    Returns:
        ``(lo, hi, note)`` — half-open sample range, plus that diagnosis.
    """
    env = _moving_mean(np.abs(aux_ac - aux_ac.mean()), win)
    med = float(np.median(env))
    if med <= 0:
        raise AuxAxisUnavailable("aux interferogram is flat - no fringes at all")
    mask = env > frac * med
    edges = np.flatnonzero(np.diff(np.concatenate(([0], mask.view(np.int8), [0]))))
    if edges.size < 2:
        raise AuxAxisUnavailable("no contiguous fringing region in the capture")
    runs = edges.reshape(-1, 2)
    lo, hi = runs[np.argmax(runs[:, 1] - runs[:, 0])]
    lo, hi = int(lo), int(hi)

    note = ""
    if dc is not None and hi < aux_ac.size - win:
        dc_sweep = float(np.mean(dc[lo:hi]))
        dc_tail = float(np.mean(dc[hi:]))
        if dc_sweep > 0 and dc_tail > 0.5 * dc_sweep:
            note = "tail is a parked laser (fringes stop, light stays)"
        else:
            note = ("tail loses light as well as fringes - the sweep may have "
                    "been cut short at the band edge, so the axis scale is "
                    "suspect")
    return lo, hi, note


def load_aux_reference(path: str = None):
    """Stored aux delay, used only for the staleness check. None when absent."""
    try:
        with open(path or AUX_REF_FILE, 'r') as fh:
            return json.load(fh)
    except Exception:
        return None


def aux_wavelength_axis(traces: dict, start_nm: float, stop_nm: float,
                        rate_hz: float, aux_channels=AUX_CHANNELS,
                        tau_ref_s: float = None, log=print):
    """Wavelength for every sample, from the aux fringe phase.

    For sweep configurations that have no lambda(t) calibration. The aux
    supplies the SHAPE of the axis and the crop; the commanded endpoints
    supply the SCALE.

    Args:
        traces: ``{channel: 1-D array}`` raw capture, aux pair included.
        start_nm / stop_nm: the range the laser was actually set to sweep (read
            back from the instrument, not the user's request).
        rate_hz: sample rate of the capture.
        tau_ref_s: aux delay to check against, for callers that have one.
            Defaults to the stored reference. This is the ONLY thing that
            catches an aliased fringe — see AUX_TAU_DEVIATION_REFUSE — so
            without it a sweep too fast for the sample rate is accepted and
            silently wrong.

    Returns:
        ``(wavelengths_nm, {channel: array}, diag)`` — traces cropped to the
        sweep, consistently with the axis.

    Raises:
        AuxAxisUnavailable: the fringe cannot carry an axis. Fall back to the
            linear axis and report ``str(exc)``.
    """
    from lina.analysis.lina_ofdr import (combine_pair, highpass, aux_phase,
                                         aux_quality)

    missing = [ch for ch in aux_channels if ch not in traces]
    if missing:
        raise AuxAxisUnavailable(
            "aux channels %s not in the capture (missing %s) - all four "
            "channels must be captured in one sweep"
            % (list(aux_channels), missing))
    n_pts = min((len(v) for v in traces.values()), default=0)
    if n_pts < AUX_MIN_FRINGES:
        raise AuxAxisUnavailable("capture is only %d samples" % n_pts)
    cut = {ch: np.asarray(v, dtype=float)[:n_pts] for ch, v in traces.items()}

    def quiet(*_a, **_k):
        pass

    aux_ac = highpass(
        combine_pair(cut, list(aux_channels), "auto", "aux", log=quiet),
        rate_hz, 500.0)
    lo, hi, tail_note = aux_sweep_window(aux_ac, np.abs(cut[aux_channels[0]]))
    win_ac = aux_ac[lo:hi]
    if win_ac.size < AUX_MIN_FRINGES:
        raise AuxAxisUnavailable(
            "the fringing region is only %d samples long" % win_ac.size)

    f_aux, tone = aux_quality(win_ac, rate_hz)
    if f_aux > AUX_MAX_NYQUIST_FRAC * rate_hz:
        raise AuxAxisUnavailable(
            "aux fringe %.1f kHz is %.0f %% of the sample rate - aliased; "
            "sweep slower or sample faster"
            % (f_aux / 1e3, 100 * f_aux / rate_hz))

    phase = aux_phase(win_ac)
    d = np.diff(phase)
    rising = float(np.median(d)) >= 0
    bad = float(np.mean(d <= 0)) if rising else float(np.mean(d >= 0))
    if bad > AUX_MAX_NONMONOTONIC:
        raise AuxAxisUnavailable(
            "%.1f %% of the aux phase reverses - the fringe is not resolved"
            % (100 * bad))
    # Isolated noise reversals survive the 1 % guard above but would make the
    # wavelength axis non-monotonic, which breaks plotting and any downstream
    # interpolation. Same cumulative fix as resample_on_aux().
    phase = (np.maximum.accumulate(phase) if rising
             else np.minimum.accumulate(phase))
    span = float(phase[-1] - phase[0])
    fringes = abs(span) / (2 * np.pi)
    if fringes < AUX_MIN_FRINGES or span == 0.0:
        raise AuxAxisUnavailable(
            "only %.0f aux fringes in the sweep" % fringes)

    # Shape from the phase, scale from the endpoints.
    nu_start = C_VAC / (float(start_nm) * 1e-9)
    nu_stop = C_VAC / (float(stop_nm) * 1e-9)
    nu = nu_start + (phase - phase[0]) / span * (nu_stop - nu_start)
    wl_nm = C_VAC / nu * 1e9

    # Every sweep measures the aux delay in passing. It is NOT used above — it
    # is the check that this is still the bench the stored value describes.
    tau_implied = fringes / abs(nu_stop - nu_start)
    diag = {"fringes": fringes, "tau_aux_implied_s": tau_implied,
            "aux_fringe_hz": f_aux, "aux_tone_fraction": tone,
            "nonmonotonic_fraction": bad, "window": (lo, hi),
            "n_points": int(hi - lo), "tail_note": tail_note,
            "start_nm": float(start_nm), "stop_nm": float(stop_nm)}
    tau_ref = tau_ref_s
    if tau_ref is None:
        tau_ref = (load_aux_reference() or {}).get("tau_aux_s")
    if tau_ref:
        dev = tau_implied / float(tau_ref) - 1.0
        diag["tau_aux_reference_s"] = float(tau_ref)
        diag["tau_aux_deviation"] = dev
        if abs(dev) > AUX_TAU_DEVIATION_REFUSE:
            raise AuxAxisUnavailable(
                "the sweep implies an aux delay of %.3f ns against the bench's "
                "%.3f ns (%+.0f %%) - too far off to be a re-patched arm; the "
                "fringe is probably aliased"
                % (tau_implied * 1e9, float(tau_ref) * 1e9, 100 * dev))
        if abs(dev) > AUX_TAU_DEVIATION_WARN:
            log("aux delay %.4f ns against the stored %.4f ns (%+.2f %%) - the "
                "aux arm or the laser's range has changed since that "
                "measurement" % (tau_implied * 1e9, float(tau_ref) * 1e9,
                                 100 * dev))
    if tail_note:
        log("sweep window %d..%d of %d samples - %s"
            % (lo, hi, n_pts, tail_note))
    return wl_nm, {ch: v[lo:hi] for ch, v in cut.items()}, diag


def describe_aux(diag: dict) -> str:
    """One-line summary for a status bar / plot label."""
    s = ("aux axis, %.0f fringes, tau_aux %.3f ns"
         % (diag["fringes"], diag["tau_aux_implied_s"] * 1e9))
    if "tau_aux_deviation" in diag:
        s += " (%+.2f %% vs stored)" % (100 * diag["tau_aux_deviation"])
    return s
