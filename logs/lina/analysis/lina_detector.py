"""
Which detector produced a Lina capture, and how much of it sat on a clamp.

Shared by the GUI worker (``LabGUI/dashboards/instruments/lina_window.py``) and
the scripts (``lina/scripts/lina_sweep_test.py:run_sweep``, and through it
``lina_chip_measure.py`` / ``lina_voa_series.py``), so a capture records the
same ``detector`` block no matter which path recorded it. Before this the
scripts recorded nothing, and a file could not say whether it came from the
Mk1 LINEAR or a Mk2 LOG unit.

Pure: reads only attributes the driver already holds, no device I/O, so it is
safe to call after the capture lock is released.

A LOG frontend pins out-of-range light to a floor and a ceiling (~3 mW)
WITHOUT raising over-range, so the fraction of samples on either clamp is the
only sign of it. A clamped fringe is distorted, and in OFDR that distortion
shows up as ghost peaks at multiples of the real delay.
"""

from __future__ import annotations

import numpy as np


def clamp_fractions(caps: dict, floor_mw: float, ceiling_mw: float):
    """Fraction of samples at the floor / at the ceiling, per channel.

    Returns ``({"Ch1": f, ...}, {"Ch1": f, ...})``.
    """
    at_floor, at_ceil = {}, {}
    for ch, arr in caps.items():
        a = np.asarray(arr, dtype=float)
        if a.size:
            at_floor[f"Ch{ch}"] = float(np.mean(a <= 1.001 * floor_mw))
            at_ceil[f"Ch{ch}"] = float(np.mean(a >= 0.999 * ceiling_mw))
    return at_floor, at_ceil


def build_detector(pm, caps: dict, *, frontend, generation, sample_rate_hz,
                   capture_mode: str) -> dict:
    """The ``detector`` header block for one capture.

    ``caps`` must be the RAW capture -- before any responsivity factor moves
    the values off the clamp.
    """
    det = {"device": getattr(pm, "id", None),
           "frontend": frontend, "generation": generation,
           "serial": getattr(pm, "_want_serial", "") or None,
           "sample_rate_hz": sample_rate_hz,
           "capture_mode": capture_mode}
    if str(frontend).upper() == "LOG" and caps:
        lo = float(getattr(pm, "_log_min_mw", 0.0) or 0.0)
        hi = float(getattr(pm, "_log_max_mw", float("inf")) or float("inf"))
        at_floor, at_ceil = clamp_fractions(caps, lo, hi)
        det.update({"log_floor_mw": lo, "log_ceiling_mw": hi,
                    "fraction_at_floor": at_floor,
                    "fraction_at_ceiling": at_ceil})
    return det


def clamp_summary(det: dict) -> list:
    """Human-readable list of the channels that touched a clamp ([] if none)."""
    out = [f"{k} {100 * v:.2g} % at ceiling"
           for k, v in (det.get("fraction_at_ceiling") or {}).items() if v > 0]
    out += [f"{k} {100 * v:.2g} % at floor"
            for k, v in (det.get("fraction_at_floor") or {}).items() if v > 0]
    return out
