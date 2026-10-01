"""One loader for every raw scan, so old and new captures stay comparable.

Why this exists
---------------
Since 2026-09-17 the coreDAQs are Mk2 **LOG** units; everything before was the
Mk1 **LINEAR** demo unit. The file format did not change, but four things about
the data did, and each one would make a new scan quietly incomparable with the
old ones if every tool kept reading the file its own way:

    1. R(lambda) correction. The Lina GUI now divides each sample by the
       detector's wavelength-dependent responsivity (a few % across a
       1505-1625 nm band). The scripts and every Mk1 file do not. Here it is
       DIVIDED BACK OUT by default, from the grid the GUI stores in
       header.detector.responsivity_correction -> all scans are "raw".
    2. Sample rate. Mk2 HIGH can sample at 1 MHz; every old scan is 100 kHz.
       More samples lower the FFT noise floor by themselves (~10 dB for 10x),
       which would read as an improvement that is not one. Here a faster
       capture is low-pass filtered and resampled to REFERENCE_RATE_HZ by
       default.
    3. Clamping. A LOG input pins light above ~3 mW (and dark samples) to a
       fixed value WITHOUT an over-range flag. A clamped fringe makes ghost
       peaks at multiples of the real delay. Here the fraction on each clamp
       is measured on exactly the samples the analysis uses, and reported.
    4. Units. The GUI auto-scales its columns to [uW]/[mW]/[W]; every tool
       used to hard-code "[mW]". Here everything comes back in mW.

Files without a ``detector`` header are Mk1 LINEAR: the Mk2 units could not be
read at all before py-coredaq 2.4.0, which arrived together with that header.

What it does NOT change for an old file: nothing. Same samples, same crop,
same order -- checked bit-for-bit against the loaders it replaces.

Crop rule. A Lina recorder NPZ stores the whole buffer plus the wavelength
axis of the swept part (``wl_nm``) and its position (``meta["crop"]``); it is
cropped here. An NPZ made by json_to_npz.py, and every JSON, is already
cropped -- its header still carries ``crop`` from the recording, so the crop is
applied only when the channel is longer than the axis and the crop spans the
axis exactly. That keeps an already-cropped file from being cropped twice.

Usage
    from scan_io import load_scan
    s = load_scan("raw_data/x.json")
    s["ch"][1], s["wl_nm"], s["meta"], s["prov"]   # prov = provenance dict
    print(describe(s["prov"]))
"""

import json
import os
import re
from fractions import Fraction

import numpy as np

# Every scan up to 2026-09-17 was taken at 100 kHz. Analyses compare at this
# rate unless told otherwise.
REFERENCE_RATE_HZ = 100_000

# Clamp tolerances -- the same ones the recorder uses (lina_detector.py).
_FLOOR_TOL = 1.001
_CEIL_TOL = 0.999
# Warn above these fractions. Any sample on the ceiling is a distorted fringe;
# a few on the floor at fringe minima is normal for a dark channel.
CLAMP_WARN_CEILING = 0.0
CLAMP_WARN_FLOOR = 1e-3
# Without a LOG header there is no clamp level to test against, only the
# channel's own maximum -- and a 6-significant-figure JSON repeats its maximum
# a few hundred times on a 2.4 M-sample Mk1 scan (0.02 %) with no saturation
# at all. So the data-only check counts only a flat top above this.
PLATEAU_MIN = 1e-3

_UNIT_TO_MW = {"w": 1e3, "mw": 1.0, "uw": 1e-3, "µw": 1e-3, "μw": 1e-3,
               "nw": 1e-6, "pw": 1e-9}
_CH_KEY = re.compile(r"^Ch\s*(\d)\s*\[([^\]]+)\]$")


# ------------------------------------------------------------------ helpers
def _json_meta(d):
    if "meta" not in d.files:
        return {}
    try:
        return json.loads(str(d["meta"]))
    except Exception:
        return {}


def _channels_from_json_entry(e):
    """{n: array_in_mW} from 'ChN [unit]' columns, any power unit."""
    ch = {}
    for k, v in e.items():
        m = _CH_KEY.match(k)
        if not m or not isinstance(v, list):
            continue
        unit = m.group(2).strip().lower()
        if unit not in _UNIT_TO_MW:
            raise ValueError("unknown power unit in column %r" % k)
        a = np.asarray(v, float)
        f = _UNIT_TO_MW[unit]
        ch[int(m.group(1))] = a if f == 1.0 else a * f
    return ch


def detector_of(header):
    """(detector dict, inferred?) -- a file without one is Mk1 LINEAR."""
    det = header.get("detector") if isinstance(header, dict) else None
    if det:
        return dict(det), False
    return {"frontend": "LINEAR", "generation": "mk1"}, True


def detector_label(det, inferred=False):
    gen = str(det.get("generation") or "?").capitalize()
    fe = str(det.get("frontend") or "?").upper()
    sn = det.get("serial")
    s = "%s %s" % (gen, fe) + (" %s" % sn if sn else "")
    return s + (" (inferred)" if inferred else "")


def rate_of(header):
    """Native sample rate in Hz, or None when the file does not say."""
    for src in (header.get("rate_hz"),
                (header.get("detector") or {}).get("sample_rate_hz")):
        if src:
            return float(src)
    return None


def clamp_fractions(ch, det):
    """{ChN: {"floor": f, "ceiling": f}} on exactly these samples.

    LOG: against the stored floor/ceiling. Otherwise a data-only check -- the
    fraction sitting on the channel's own maximum, i.e. a flat-topped
    (saturated) plateau; a smooth signal touches its maximum once.
    """
    out = {}
    lo = det.get("log_floor_mw")
    hi = det.get("log_ceiling_mw")
    is_log = str(det.get("frontend", "")).upper() == "LOG" and lo is not None
    for n, a in ch.items():
        if not a.size:
            continue
        if is_log:
            out["Ch%d" % n] = {
                "floor": float(np.mean(a <= _FLOOR_TOL * lo)),
                "ceiling": float(np.mean(a >= _CEIL_TOL * hi)) if hi else 0.0}
        else:
            mx = float(a.max())
            k = int(np.count_nonzero(a >= mx - 1e-9 * abs(mx))) if mx > 0 else 0
            f = k / a.size
            out["Ch%d" % n] = {"floor": 0.0,          # >= 20: not 2 of 600
                               "ceiling": f if f > PLATEAU_MIN and k >= 20 else 0.0}
    return out


def _undo_responsivity(ch, wl, rc):
    """Divide the GUI's R(lambda_set)/R(lambda) factor back out."""
    gw = np.asarray(rc["grid_wl_nm"], float)
    gr = np.asarray(rc["grid_r_a_per_w"], float)
    r_l = np.interp(wl, gw, gr)
    with np.errstate(divide="ignore", invalid="ignore"):
        fac = np.where(r_l > 0, float(rc["r_set_a_per_w"]) / r_l, 1.0)
    return {n: a / fac for n, a in ch.items()}


def _band_edge_hz(x, rate, frac=0.99):
    """Frequency below which `frac` of the AC power of x lies."""
    seg = x[: min(x.size, 1 << 20)]
    seg = seg - np.convolve(seg, np.ones(201) / 201, mode="same")
    p = np.abs(np.fft.rfft(seg * np.hanning(seg.size))) ** 2
    f = np.fft.rfftfreq(seg.size, 1.0 / rate)
    p[f < 500] = 0
    c = np.cumsum(p)
    if c[-1] <= 0:
        return 0.0
    return float(f[np.searchsorted(c, frac * c[-1])])


def _resample(ch, wl, rate, target):
    """Anti-aliased resampling of every channel from `rate` to `target`."""
    from scipy.signal import resample_poly
    r = Fraction(int(round(target)), int(round(rate))).limit_denominator(1000)
    up, down = r.numerator, r.denominator
    out = {n: resample_poly(a, up, down) for n, a in ch.items()}
    m = min(len(a) for a in out.values())
    out = {n: a[:m] for n, a in out.items()}
    wl2 = None
    if wl is not None:
        t_old = np.arange(wl.size) / rate
        wl2 = np.interp(np.arange(m) / target, t_old, wl)
    return out, wl2


# --------------------------------------------------------------------- main
AUX_PAIR = (2, 4)            # aux MZI channels (HANDOVER: Ch2/Ch4)
# Sweep = where the aux fringe frequency stays within this fraction of its
# median. Measured 2026-10-01 on Mk2 buffers: 12.6-13.6 kHz during the sweep,
# 35-40 kHz on the return slew, ~0.4 kHz when parked (Mk1, 09-01).
SWEEP_FREQ_TOL = 0.35
SWEEP_SMOOTH = 2001          # samples; coarse window, edges refined after


def sweep_window(ch, rate, aux=AUX_PAIR):
    """(lo, hi, note): the part of a raw buffer where the laser sweeps forward.

    From the aux fringe FREQUENCY, not its envelope: on the Mk2 captures the
    return slew still fringes with half the envelope, so an envelope threshold
    keeps it. The frequency jumps by ~3x there and drops to ~0 when parked.
    """
    from scipy.ndimage import uniform_filter1d
    a = ch[aux[0]] - ch[aux[1]]
    a = a - uniform_filter1d(a, 401)
    X = np.fft.fft(a)
    n = a.size
    X[n // 2 + 1:] = 0
    X[1:n // 2] *= 2
    phi = np.unwrap(np.angle(np.fft.ifft(X)))
    f = np.abs(uniform_filter1d(np.gradient(phi), SWEEP_SMOOTH)) * rate / (2 * np.pi)
    f0 = float(np.median(f[n // 4: 3 * n // 4]))
    ok = np.abs(f - f0) < SWEEP_FREQ_TOL * f0
    edges = np.flatnonzero(np.diff(np.concatenate(([0], ok.view(np.int8), [0]))))
    runs = edges.reshape(-1, 2)
    lo, hi = (int(v) for v in runs[np.argmax(runs[:, 1] - runs[:, 0])])
    # Refine each edge on a 0.5 ms scale: the sweep ends where the fringe
    # frequency leaves +-50 % of its median -- a collapse when the laser stops
    # (Mk2 captures: it dips to ~5 kHz before the return slew) or a jump when
    # the slew starts directly. The tuning ripple inside the sweep stays
    # within ~+-25 %. The coarse 20 ms window
    # alone would put the edge ~10 ms off -- ~130 fringes, i.e. 0.1 % of the
    # aux delay when that is scaled from the commanded endpoints.
    fine = np.abs(uniform_filter1d(np.gradient(phi), 51)) * rate / (2 * np.pi)
    low = np.abs(fine - f0) > 0.5 * f0
    if hi < n:
        a0 = max(hi - SWEEP_SMOOTH, 0)
        k = np.flatnonzero(low[a0:hi + SWEEP_SMOOTH])
        hi = a0 + int(k[0]) if k.size else hi - SWEEP_SMOOTH // 2
    if lo > 0:
        b1 = min(lo + SWEEP_SMOOTH, n)
        k = np.flatnonzero(low[max(lo - SWEEP_SMOOTH, 0):b1])
        lo = max(lo - SWEEP_SMOOTH, 0) + int(k[-1]) + 1 if k.size else lo + SWEEP_SMOOTH // 2
    tail = f[hi:]
    what = ("tail: return slew" if tail.size and np.median(tail) > 1.5 * f0
            else "tail: parked" if tail.size else "no tail")
    return lo, hi, "%.3f-%.3f s, f_aux %.2f kHz, %s" % (
        lo / rate, hi / rate, f0 / 1e3, what)


def load_scan(path, *, undo_responsivity=True, target_rate_hz=REFERENCE_RATE_HZ,
              crop_sweep=True, log=print):
    """Read any raw scan (.json, .npz from json_to_npz, Lina recorder, sim).

    Args:
        undo_responsivity: divide a GUI R(lambda) correction back out
            (default). False keeps the samples as stored.
        target_rate_hz: resample a FASTER capture down to this rate (default
            100 kHz, the rate of every Mk1 scan). None keeps the native rate.
            A capture at or below the target, or of unknown rate, is never
            touched.
        crop_sweep: cut a raw run_sweep buffer (no axis, no crop) to the
            forward sweep, found from the aux fringe frequency (sweep_window).

    Returns dict:
        ch     {1..4: float64 array, mW}
        wl_nm  wavelength per sample or None
        meta   header + Device Description (+ truth_z/truth_db, step_us for sims)
        prov   provenance: what was done to get here, and what to watch for
    """
    prov = {"source": os.path.basename(path), "warnings": []}
    wl = None
    meta = {}

    if path.endswith(".npz"):
        d = np.load(path, allow_pickle=False)
        meta = _json_meta(d)
        ch = {int(k[2:]): np.asarray(d[k], float)
              for k in d.files if re.fullmatch(r"ch\d", k)}
        unit = str(d["unit"]) if "unit" in d.files else "mW"
        if unit.lower() != "mw":
            f = _UNIT_TO_MW[unit.lower()]
            ch = {n: a * f for n, a in ch.items()}
        for k in ("step_us", "lam0_nm", "speed_nms"):
            if k in d.files:
                meta[k] = d[k].item()
        if "truth_z" in d.files:
            meta["truth_z"] = d["truth_z"]
            meta["truth_db"] = d["truth_db"]
        if "wl_nm" in d.files:
            wl = np.asarray(d["wl_nm"], float)
            prov["format"] = "npz (Lina recorder)"
        elif "wavelength_nm" in d.files:
            wl = np.asarray(d["wavelength_nm"], float)
            prov["format"] = "npz (json_to_npz)"
        elif "lam_start_nm" in d.files and ch:
            wl = np.linspace(float(d["lam_start_nm"]), float(d["lam_stop_nm"]),
                             len(next(iter(ch.values()))))
            prov["format"] = "npz (json_to_npz, linear axis)"
        else:
            prov["format"] = "npz"
    else:
        with open(path) as f:
            doc = json.load(f)
        if len(doc["data"]) != 1:
            prov["warnings"].append("file has %d data entries, using entry 0"
                                    % len(doc["data"]))
        e = doc["data"][0]
        ch = _channels_from_json_entry(e)
        meta = dict(doc.get("header") or {})
        meta.update(e.get("Device Description") or {})
        if "Wavelength [nm]" in e:
            wl = np.asarray(e["Wavelength [nm]"], float)
        prov["format"] = "json"
    if not ch:
        raise ValueError("%s: no channel columns found" % path)

    # -- crop (Lina recorder NPZ only; see module docstring)
    n_ch = min(len(a) for a in ch.values())
    crop = meta.get("crop")
    if (crop and wl is not None and n_ch > wl.size
            and int(crop[1]) - int(crop[0]) == wl.size):
        lo, hi = int(crop[0]), int(crop[1])
        ch = {n: a[lo:hi] for n, a in ch.items()}
        prov["cropped"] = [lo, hi]
        prov["crop_source"] = "meta crop"
    elif (crop_sweep and wl is None and "sweep_s" in meta and "truth_z" not in meta
          and all(n in ch for n in AUX_PAIR)):
        # A raw lina_sweep_test / run_sweep buffer: the whole armed window,
        # no axis, no crop. Since the Mk2 change the EXFO's return slew lands
        # INSIDE that window (fast fringes, 35-40 kHz, after ~9.9 s of a 10 s
        # sweep) -- and a Hilbert phase cannot tell backward from forward, so
        # left in, it would be read as more sweep.
        lo, hi, note = sweep_window(ch, rate_of(meta) or REFERENCE_RATE_HZ)
        ch = {n: a[lo:hi] for n, a in ch.items()}
        prov["cropped"] = [lo, hi]
        prov["crop_source"] = "aux fringe frequency (%s)" % note

    # -- detector
    det, inferred = detector_of(meta)
    prov["detector"] = det
    prov["detector_inferred"] = inferred
    prov["detector_label"] = ("synthetic" if "truth_z" in meta
                              else detector_label(det, inferred))

    # -- R(lambda)
    rc = det.get("responsivity_correction")
    if det.get("wl_cal_applied"):
        if not undo_responsivity:
            prov["responsivity"] = "applied (kept)"
            prov["warnings"].append(
                "R(lambda) correction kept -- NOT comparable with raw scans "
                "(scripts, Mk1) in absolute power across the band")
        elif rc and wl is not None:
            ch = _undo_responsivity(ch, wl, rc)
            prov["responsivity"] = "undone"
        else:
            prov["responsivity"] = "applied (cannot undo)"
            prov["warnings"].append(
                "R(lambda) correction applied by the GUI and the file does "
                "not store the factor -- absolute power across the band is "
                "off by a few %% against raw scans")
    else:
        prov["responsivity"] = "none"

    # -- clamp, on the samples actually analysed (before any resampling)
    cf = clamp_fractions(ch, det)
    prov["clamp"] = cf
    is_log = str(det.get("frontend", "")).upper() == "LOG"
    for name, c in cf.items():
        if not is_log:
            if c["ceiling"] > 0:
                prov["warnings"].append(
                    "%s: flat top -- %.3g %% of samples on the channel "
                    "maximum (saturated?)" % (name, 100 * c["ceiling"]))
            continue
        if c["ceiling"] > CLAMP_WARN_CEILING:
            prov["warnings"].append(
                "%s: %.3g %% of samples on the CEILING clamp -- fringes "
                "distorted, expect ghost peaks at multiples of real delays; "
                "attenuate and re-measure" % (name, 100 * c["ceiling"]))
        if c["floor"] > CLAMP_WARN_FLOOR:
            prov["warnings"].append(
                "%s: %.3g %% of samples on the floor clamp (dark or "
                "near-dark channel)" % (name, 100 * c["floor"]))

    # -- rate
    rate = rate_of(meta)
    prov["rate_hz_native"] = rate
    prov["rate_hz"] = rate
    if target_rate_hz and rate and rate > target_rate_hz * 1.001:
        edge = max(_band_edge_hz(a, rate) for a in ch.values())
        ch, wl = _resample(ch, wl, rate, target_rate_hz)
        prov["rate_hz"] = float(target_rate_hz)
        prov["resampled"] = "%g Hz -> %g Hz" % (rate, target_rate_hz)
        if edge > 0.4 * target_rate_hz:
            prov["warnings"].append(
                "signal reaches %.0f kHz, above the %.0f kHz reference band "
                "-- resampled for comparability this content is lost; use "
                "--native-rate to see it (and then do not compare dB levels "
                "with 100 kHz scans)" % (edge / 1e3, target_rate_hz / 2e3))
    elif rate and rate > REFERENCE_RATE_HZ * 1.001:
        prov["warnings"].append(
            "native %g Hz kept -- noise floor in dB is NOT comparable with "
            "100 kHz scans (more samples alone lower it)" % rate)

    for w in prov["warnings"]:
        log("  WARNING [%s]: %s" % (prov["source"], w))
    return {"ch": ch, "wl_nm": wl, "meta": meta, "prov": prov}


def describe(prov):
    """One line: detector, rate, R(lambda), worst clamp."""
    worst = max(((c["ceiling"], c["floor"]) for c in prov.get("clamp", {}).values()),
                default=(0.0, 0.0))
    rate = prov.get("rate_hz")
    s = "%s | %s | R(lambda) %s | clamp ceil %.2g %% floor %.2g %%" % (
        prov.get("detector_label", "?"),
        ("%g kHz" % (rate / 1e3)) if rate else "rate ?",
        prov.get("responsivity", "?"), 100 * worst[0], 100 * worst[1])
    if prov.get("crop_source", "").startswith("aux"):
        s += " | sweep %s" % prov["crop_source"][len("aux fringe frequency ("):-1].split(",")[0]
    if prov.get("resampled"):
        s += " | resampled " + prov["resampled"]
    return s


def jsonable(prov):
    """prov with numpy scalars made plain, for a sidecar file."""
    return json.loads(json.dumps(prov, default=lambda o: o.item()
                                 if hasattr(o, "item") else str(o)))


# ---------------------------------------------------- comparability checks
def sidecar_path(prefix):
    return prefix + "_reflectogram.json"


def write_sidecar(prefix, prov, extra=None):
    out = dict(jsonable(prov))
    out.update(extra or {})
    with open(sidecar_path(prefix), "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    return sidecar_path(prefix)


def provenance_for(path):
    """Provenance of a reflectogram CSV (from its sidecar) or a raw scan.

    Returns None for a CSV written before sidecars existed -- all of those
    are Mk1 LINEAR at 100 kHz, which is what `check_comparable` assumes.
    """
    if path.endswith("_reflectogram.csv"):
        side = path[: -len(".csv")] + ".json"
        if os.path.exists(side):
            with open(side) as fh:
                return json.load(fh)
        return None
    if path.endswith(".csv"):
        return None
    # a raw scan: header only would be cheaper, but the clamp needs the data
    return jsonable(load_scan(path, log=lambda *_: None)["prov"])


_LEGACY = {"detector_label": "Mk1 LINEAR (inferred)", "rate_hz": 100000.0,
           "responsivity": "none", "clamp": {}, "warnings": []}


def check_comparable(items, log=print):
    """Print a provenance table for scans that are about to be compared and
    warn where they differ in a way that moves dB levels or peaks.

    items: list of (label, prov-or-None). Returns the list of warnings.
    """
    rows = [(lab, p or _LEGACY) for lab, p in items]
    log("\nComparability check (scan_io):")
    log("  %-34s %-26s %-9s %-14s %s" % ("scan", "detector", "rate",
                                        "R(lambda)", "clamp ceil/floor %"))
    for lab, p in rows:
        worst = max(((c["ceiling"], c["floor"])
                     for c in (p.get("clamp") or {}).values()), default=(0, 0))
        r = p.get("rate_hz")
        log("  %-34s %-26s %-9s %-14s %.2g / %.2g" % (
            str(lab)[:34], p.get("detector_label", "?")[:26],
            ("%g k" % (r / 1e3)) if r else "?", p.get("responsivity", "?"),
            100 * worst[0], 100 * worst[1]))
    warn = []
    gens = {p.get("detector_label", "?").replace(" (inferred)", "")
            .rsplit(" SN", 1)[0] for _, p in rows}
    if len(gens) > 1:
        warn.append("mixed detectors %s: peak POSITIONS compare; noise "
                    "floor / dynamic range are NOT yet shown to compare "
                    "(LOG noise is ~constant in dB, LINEAR noise constant in "
                    "mW) -- see logs/2026-10-01.md, Mk2 reference capture"
                    % sorted(gens))
    rates = {p.get("rate_hz") for _, p in rows if p.get("rate_hz")}
    if len(rates) > 1:
        warn.append("different sample rates %s Hz: dB levels not comparable"
                    % sorted(rates))
    resp = {p.get("responsivity") for _, p in rows}
    if resp - {"none", "undone"}:
        warn.append("some scans still carry an R(lambda) correction: "
                    "absolute powers across the band differ by a few %")
    for lab, p in rows:
        for c_name, c in (p.get("clamp") or {}).items():
            if c.get("ceiling", 0) > CLAMP_WARN_CEILING:
                warn.append("%s: %s on the ceiling clamp -- ghost peaks "
                            "possible" % (lab, c_name))
    for w in warn:
        log("  WARNING: " + w)
    if not warn:
        log("  ok: same detector class, rate and R(lambda) handling")
    return warn
