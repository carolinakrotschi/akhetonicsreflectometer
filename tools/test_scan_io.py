#!/usr/bin/env python3
"""Checks for scan_io.py -- a Mk2 LOG capture must come back on the same
footing as the Mk1 scans.

Synthetic, no files needed, ~10 s:

    python tools/test_scan_io.py        (or: python -m pytest tools/test_scan_io.py)

What is pinned:
  * the GUI's R(lambda) correction is divided back out exactly
  * [uW] / [W] columns come back in mW
  * a 1 MHz capture is resampled to 100 kHz without moving a reflector
  * LOG clamp fractions are measured, and a ceiling clamp warns
  * the crop is applied to a Lina recorder NPZ and NOT to an already-cropped
    file that still carries "crop" in its header
  * a file without a detector header is Mk1 LINEAR
  * check_comparable flags mixed detectors and rates
"""

import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import process_reflectogram_aux as pra  # noqa: E402
from scan_io import check_comparable, load_scan  # noqa: E402

QUIET = dict(log=lambda *_: None)


def _fringes(rate, dur=2.0, f_aux=13e3, f_meas=4e3, seed=0):
    """Four channels: aux pair (complementary) and measurement pair, in mW."""
    rng = np.random.default_rng(seed)
    t = np.arange(int(rate * dur)) / rate
    ph_aux = 2 * np.pi * f_aux * t * (1 + 0.01 * np.sin(2 * np.pi * 0.7 * t))
    ph_meas = ph_aux * (f_meas / f_aux)
    c1 = 0.02 * (1 + 0.6 * np.cos(ph_meas))
    c3 = 0.02 * (1 - 0.6 * np.cos(ph_meas))
    c2 = 0.35 * (1 + 0.8 * np.cos(ph_aux))
    c4 = 0.35 * (1 - 0.8 * np.cos(ph_aux))
    ch = {1: c1, 2: c2, 3: c3, 4: c4}
    return {n: a + 1e-5 * rng.standard_normal(a.size) for n, a in ch.items()}, t


def _write_gui_json(path, ch, wl, header, unit="mW"):
    f = {"mW": 1.0, "uW": 1e3, "W": 1e-3}[unit]
    e = {"Device Description": {}, "Wavelength [nm]": list(map(float, wl))}
    for n, a in ch.items():
        e["Ch%d [%s]" % (n, unit)] = list(map(float, a * f))
    with open(path, "w") as fh:
        json.dump({"header": header, "data": [e]}, fh)


def _log_detector(rate, rc=None):
    det = {"device": "coredaq-1", "frontend": "LOG", "generation": "mk2",
           "serial": "SN0001", "sample_rate_hz": rate,
           "capture_mode": "continuous", "log_floor_mw": 1e-7,
           "log_ceiling_mw": 3.0, "wl_cal_applied": rc is not None}
    if rc:
        det["responsivity_correction"] = rc
    return det


def test_responsivity_undone_and_units():
    ch, t = _fringes(100e3, dur=0.5)
    wl = 1505 + 120 * t / t[-1]
    gw = np.linspace(1505, 1625, 256)
    gr = 1.0 + 0.0004 * (gw - 1550)                 # a few % across the band
    r_set = float(np.interp(1550, gw, gr))
    fac = r_set / np.interp(wl, gw, gr)
    stored = {n: a * fac for n, a in ch.items()}    # what the GUI writes
    rc = {"wl_set_nm": 1550.0, "r_set_a_per_w": r_set,
          "grid_wl_nm": gw.tolist(), "grid_r_a_per_w": gr.tolist()}
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "gui.json")
        _write_gui_json(p, stored, wl, {"detector": _log_detector(100000, rc)},
                        unit="uW")
        s = load_scan(p, **QUIET)
        assert s["prov"]["responsivity"] == "undone"
        for n in ch:
            assert np.allclose(s["ch"][n], ch[n], rtol=1e-9, atol=0)
        k = load_scan(p, undo_responsivity=False, **QUIET)
        assert k["prov"]["responsivity"] == "applied (kept)"
        assert np.allclose(k["ch"][2], stored[2], rtol=1e-9, atol=0)


def test_resample_1mhz_keeps_reflector_position():
    """Same fringe pattern at 100 kHz and at 1 MHz -> same reflectogram peak."""
    ref, t1 = _fringes(100e3)
    fast, t2 = _fringes(1e6)
    with tempfile.TemporaryDirectory() as d:
        p1, p2 = os.path.join(d, "slow.json"), os.path.join(d, "fast.json")
        _write_gui_json(p1, ref, 1520 + 50 * t1 / t1[-1], {"rate_hz": 100000})
        _write_gui_json(p2, fast, 1520 + 50 * t2 / t2[-1],
                        {"detector": _log_detector(1000000)})
        s1, s2 = load_scan(p1, **QUIET), load_scan(p2, **QUIET)
        assert s2["prov"]["rate_hz"] == 100000.0
        assert s2["prov"]["resampled"] == "1e+06 Hz -> 100000 Hz"
        assert abs(len(s2["ch"][1]) - len(s1["ch"][1])) <= 1
        assert len(s2["wl_nm"]) == len(s2["ch"][1])

        def peak(s):
            c = s["ch"]
            aux, _ = pra.balanced(c[2], c[4])
            meas, _ = pra.balanced(c[1], c[3])
            y, dnu, _, _ = pra.resample_on_aux(meas, aux, 20e-9)
            R = np.abs(np.fft.rfft((y - y.mean()) * np.hanning(y.size)))
            z = np.arange(R.size) * pra.C / (2 * pra.NG * dnu * y.size)
            return z[int(np.argmax(R[5:])) + 5], z[1] - z[0]

        (za, dz), (zb, _) = peak(s1), peak(s2)
        assert abs(za - zb) <= dz, (za, zb, dz)
        n = load_scan(p2, target_rate_hz=None, **QUIET)
        assert n["prov"]["rate_hz"] == 1e6 and "resampled" not in n["prov"]
        assert any("NOT comparable" in w for w in n["prov"]["warnings"])


def test_clamp_detected():
    ch, t = _fringes(100e3, dur=0.2)
    ch[2] = np.minimum(ch[2] * 10, 3.0)             # aux driven into the ceiling
    ch[1] = np.maximum(ch[1] - 0.01, 1e-7)          # fringe minima on the floor
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "clamp.json")
        _write_gui_json(p, ch, 1520 + t, {"detector": _log_detector(100000)})
        s = load_scan(p, **QUIET)
        c = s["prov"]["clamp"]
        assert c["Ch2"]["ceiling"] > 0.1 and c["Ch4"]["ceiling"] == 0.0
        assert c["Ch1"]["floor"] > 0.01
        assert any("CEILING" in w and w.startswith("Ch2") for w in s["prov"]["warnings"])


def test_crop_only_where_not_yet_cropped():
    ch, _ = _fringes(100e3, dur=0.2)
    n = len(ch[1])
    lo, hi = 1000, n - 3000
    wl = np.linspace(1520, 1570, hi - lo)
    with tempfile.TemporaryDirectory() as d:
        # Lina recorder: full buffer + axis of the swept part + crop
        p = os.path.join(d, "rec.npz")
        meta = {"crop": [lo, hi], "rate_hz": 100000}
        np.savez_compressed(p, meta=json.dumps(meta), wl_nm=wl,
                            **{"ch%d" % k: v for k, v in ch.items()})
        s = load_scan(p, **QUIET)
        assert s["prov"]["cropped"] == [lo, hi]
        assert np.array_equal(s["ch"][3], ch[3][lo:hi])
        # json_to_npz of an already-cropped JSON: header still says crop
        p2 = os.path.join(d, "conv.npz")
        cut = {k: v[lo:hi] for k, v in ch.items()}
        np.savez_compressed(p2, meta=json.dumps(meta), wavelength_nm=wl,
                            **{"ch%d" % k: v for k, v in cut.items()})
        s2 = load_scan(p2, **QUIET)
        assert "cropped" not in s2["prov"]
        assert np.array_equal(s2["ch"][3], cut[3])


def test_raw_buffer_cut_before_return_slew():
    """A raw run_sweep buffer: 9.9 s forward sweep, then a 3x faster return
    slew (as on the 2026-10-01 Mk2 captures). The slew must be cut off."""
    rate, n_fwd, n_ret = 100e3, 99_000, 6_000
    rng = np.random.default_rng(1)
    f = np.concatenate((np.full(n_fwd, 13e3), np.full(n_ret, 38e3)))
    ph = 2 * np.pi * np.cumsum(f) / rate
    ch = {2: 0.03 * (1 + 0.8 * np.cos(ph)), 4: 0.03 * (1 - 0.8 * np.cos(ph)),
          1: 0.007 + 1e-4 * rng.standard_normal(ph.size),
          3: 0.007 + 1e-4 * rng.standard_normal(ph.size)}
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "raw.npz")
        meta = {"rate_hz": 100000, "sweep_s": 1.0, "actual_start_nm": 1520.0,
                "actual_stop_nm": 1570.0}
        np.savez_compressed(p, meta=json.dumps(meta),
                            **{"ch%d" % k: v for k, v in ch.items()})
        s = load_scan(p, **QUIET)
        lo, hi = s["prov"]["cropped"]
        assert lo == 0 and abs(hi - n_fwd) < 60, (lo, hi)
        assert "return slew" in s["prov"]["crop_source"]
        assert len(load_scan(p, crop_sweep=False, **QUIET)["ch"][2]) == ph.size


def test_json_to_npz_keeps_unit_and_detector():
    import json_to_npz
    ch, t = _fringes(100e3, dur=0.2)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "gui.json")
        _write_gui_json(p, ch, 1520 + t, {"detector": _log_detector(100000)},
                        unit="uW")
        out = json_to_npz.convert(p, verify=True)
        a, b = load_scan(p, **QUIET), load_scan(out, **QUIET)
        for n in ch:
            assert np.array_equal(a["ch"][n], b["ch"][n])
        assert b["prov"]["detector_label"] == "Mk2 LOG SN0001"


def test_legacy_is_mk1_and_comparability():
    ch, t = _fringes(100e3, dur=0.2)
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "old.json")
        _write_gui_json(p, ch, 1520 + t, {"rate_hz": 100000})
        old = load_scan(p, **QUIET)["prov"]
        assert old["detector_inferred"] and old["detector_label"].startswith("Mk1 LINEAR")
        p2 = os.path.join(d, "new.json")
        _write_gui_json(p2, ch, 1520 + t, {"detector": _log_detector(100000)})
        new = load_scan(p2, **QUIET)["prov"]
        w = check_comparable([("old", old), ("new", new)], **QUIET)
        assert any("mixed detectors" in x for x in w)
        assert not any("sample rates" in x for x in w)
        assert check_comparable([("a", old), ("b", None)], **QUIET) == []


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok  ", name)
    print("all passed")
