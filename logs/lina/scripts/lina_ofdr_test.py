"""
OFDR measurement + aux-referenced FFT for Lina (EXFO T200S + LINEAR coreDAQ).

Setup this assumes (override with --meas / --ref):

    Ch1, Ch3  measurement MZI — contains the test fibre whose length we want
    Ch2, Ch4  reference ("aux") MZI — the phase reference for calibration

Each MZI is read out on TWO channels. If they are the complementary outputs of
a 2x2 output coupler (balanced detection), their difference cancels the common
-mode intensity envelope and doubles the fringe amplitude. The script measures
the correlation between the two and only subtracts when they are genuinely
anti-correlated, reporting what it did (--balance never/always to force it).

WHAT THE AUX MZI IS FOR
-----------------------
The coreDAQ samples uniformly in TIME, but an FFT only gives a clean delay
(= distance) spectrum if the signal is uniform in OPTICAL FREQUENCY. The EXFO
is not: measured with a tunable filter as a marker, it sweeps 1-2 % faster
than commanded and slightly non-uniformly (see lina/analysis/lina_wl_cal.py). The
aux MZI has a FIXED delay, so its fringe phase is strictly proportional to
optical frequency:

    phi_aux(t) = 2*pi*tau_aux*nu(t)

Unwrapping that phase therefore measures nu(t) directly, whatever the laser
did. Resampling the measurement signal onto a uniform phi_aux grid makes it
uniform in nu, and the FFT then resolves delays properly. Without it, sweep
nonlinearity smears every peak.

DISTANCE AXIS
-------------
    dz     = c / (2 * n_g * dnu_total)          resolution (bin spacing)
    z_max  = dz * N/2                           Nyquist distance
             = c * f_s / (4 * n_g * dnu/dt)     equivalently, from the rate

The default is the REFLECTION convention (factor 2*n_g), because on this setup
the test fibre is traversed TWICE. That was established by measurement, not
assumed: nominally 1 m fibres gave a 2.079 m path contribution and a nominal
3 m fibre gave 6.051 m — a factor of two, consistently. So dividing by 2 makes
the axis read physical fibre length. (--convention transmission uses n_g if a
setup traverses the fibre once.)

READING FIBRE LENGTH DIRECTLY (--offset-m / --mirror)
The interferometer's OWN path imbalance adds to every fibre-dependent peak, so
even in the right convention a 1.040 m fibre sits at 0.574 + 1.040 = 1.614 m.
Measured across three states (two fibres and none):

    measurement pair:  peak = 0.574 m + L_fibre
    aux pair:          peak = 2.081 m - L_fibre     (moves the other way)

Those two offsets are just the peaks that remain with NO fibre connected. So
--offset-m 0.574 makes the measurement axis read fibre length, and
--offset-m 2.081 --mirror does the same for the aux pair. Verified: a 1 m
patch cord then peaks at 1.040 m and a 3 m one at 3.026 m, in both pairs
independently.

USAGE
-----
  # measure and process (filter must be OUT of the path)
  python lina/scripts/lina_ofdr_test.py measure --start 1520 --stop 1570 --speed 5

  # re-process a saved capture with different windowing, no hardware
  python lina/scripts/lina_ofdr_test.py process --raw lina_ofdr_data/ofdr_*.npz \
      --beta 12 --pad-factor 4

NOTE: 4 channels at 100 kHz for 10 s is a ~8 MB USB read; that transfer
intermittently stalls on this unit, so run_sweep retries the whole capture.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from scipy.signal import hilbert, find_peaks, butter, filtfilt

C_VAC = 299_792_458.0        # m/s
N_GROUP_SMF28 = 1.4682       # SMF-28 group index at 1550 nm


# ---------------------------------------------------------------------------
# Signal conditioning
# ---------------------------------------------------------------------------

def combine_pair(traces: dict, channels: list, mode: str, label: str):
    """Turn a channel pair into one AC interferogram.

    Balanced detection: the two outputs of a 2x2 coupler are complementary, so
    their difference removes the shared intensity envelope (and the laser's own
    power drift across the sweep) while doubling the fringe amplitude. Whether
    that is what is actually wired is decided from the data, not assumed.
    """
    arrs = [np.asarray(traces[ch], dtype=float) for ch in channels if ch in traces]
    if not arrs:
        raise SystemExit(f"{label}: none of channels {channels} present in the capture")
    if len(arrs) == 1:
        print(f"[sig] {label}: single channel {channels[0]}")
        return arrs[0] - arrs[0].mean()


    a, b = arrs[0], arrs[1]
    # Correlate the AC parts only; the DC envelope would dominate otherwise.
    ac = lambda v: v - np.mean(v)
    denom = np.std(ac(a)) * np.std(ac(b))
    corr = float(np.mean(ac(a) * ac(b)) / denom) if denom > 0 else 0.0
    balanced = (mode == "always") or (mode == "auto" and corr < -0.3)
    if balanced:
        sig = ac(a) - ac(b)
        how = "difference (balanced pair)"
    else:
        sig = ac(a) + ac(b) if corr > 0.3 else ac(a)
        how = ("sum (in-phase pair)" if corr > 0.3
               else f"channel {channels[0]} only (uncorrelated pair)")
    print(f"[sig] {label}: Ch{channels[0]}/Ch{channels[1]} correlation "
          f"{corr:+.3f} -> using {how}")
    return sig


def highpass(sig: np.ndarray, rate_hz: float, cutoff_hz: float = 500.0):
    """Strip the slow intensity envelope, keep the fringes.

    Balanced subtraction removes most of it, but not all: the two detectors
    have unequal responsivity, so a residual of the laser's power drift across
    the sweep survives. That residual is orders of magnitude larger than the
    fringe amplitude at DC, which (a) makes the aux spectrum look like it has
    no fringe at all and (b) biases the Hilbert phase, smearing every FFT peak.
    The fringes sit at kHz and the envelope at a few Hz, so they separate
    cleanly.
    """
    nyq = rate_hz / 2.0
    wn = min(max(cutoff_hz / nyq, 1e-6), 0.99)
    b, a_ = butter(4, wn, btype="highpass")
    return filtfilt(b, a_, sig)


# ---------------------------------------------------------------------------
# Aux-referenced resampling
# ---------------------------------------------------------------------------

def aux_phase(aux_sig: np.ndarray):
    """Unwrapped fringe phase of the aux interferogram.

    Uses the analytic signal (Hilbert transform). The aux fringe must be a
    clean single tone for this to be meaningful — the caller checks that its
    spectrum has one dominant peak.
    """
    analytic = hilbert(aux_sig - aux_sig.mean())
    return np.unwrap(np.angle(analytic))


def aux_quality(aux_sig: np.ndarray, rate_hz: float):
    """Dominant aux fringe frequency and how single-tone the aux really is."""
    spec = np.abs(np.fft.rfft(aux_sig - aux_sig.mean()))
    freqs = np.fft.rfftfreq(aux_sig.size, 1.0 / rate_hz)
    k = int(np.argmax(spec[1:]) + 1)
    total = float(np.sum(spec[1:] ** 2))
    # Energy within +-2 % of the dominant line, as a single-tone measure.
    band = (freqs > freqs[k] * 0.98) & (freqs < freqs[k] * 1.02)
    frac = float(np.sum(spec[band] ** 2) / total) if total > 0 else 0.0
    return float(freqs[k]), frac


def resample_on_aux(meas_sig: np.ndarray, phase: np.ndarray, n_out: int = None):
    """Interpolate the measurement onto a uniform aux-phase grid.

    Equal steps in aux phase are equal steps in optical frequency, which is
    what the FFT needs. Monotonicity of the unwrapped phase is required and
    checked — a phase that reverses means the aux fringe was not resolved
    (aliased) or the sweep direction changed mid-capture.
    """
    dphi = np.diff(phase)
    if np.median(dphi) < 0:            # descending sweep: flip to ascending
        phase, meas_sig = -phase[::-1], meas_sig[::-1]
        dphi = np.diff(phase)
    bad = float(np.mean(dphi <= 0))
    n_out = int(n_out or meas_sig.size)
    grid = np.linspace(phase[0], phase[-1], n_out)
    # np.interp needs an increasing x; enforce it cumulatively so isolated
    # non-monotonic samples (noise on the phase) don't corrupt the mapping.
    mono = np.maximum.accumulate(phase)
    out = np.interp(grid, mono, meas_sig)
    return out, bad, float(phase[-1] - phase[0])


# ---------------------------------------------------------------------------
# FFT to distance
# ---------------------------------------------------------------------------

def fft_to_distance(sig: np.ndarray, *, dnu_total_hz: float, n_group: float,
                    beta: float = 10.0, pad_factor: int = 1,
                    convention: str = "reflection"):
    """Windowed FFT of a frequency-uniform interferogram -> distance spectrum.

    Args:
        dnu_total_hz: total optical frequency span the samples cover. This
            alone sets the distance scale; the number of samples only sets how
            far the axis reaches.
        beta: Kaiser window beta. Higher = lower sidelobes, wider main lobe
            (10 ≈ -70 dB sidelobes; 12-14 for high dynamic range).
        pad_factor: zero-pad to this multiple — interpolates the peak shape,
            it does NOT improve resolution.
    """
    n = sig.size
    win = np.kaiser(n, beta)
    padded = int(n * max(pad_factor, 1))
    spec = np.fft.rfft((sig - sig.mean()) * win, n=padded)
    mag = np.abs(spec)

    # Delay per FFT bin: sampling in nu with step dnu = dnu_total/n gives a
    # delay axis of step 1/(padded*dnu).
    dnu = dnu_total_hz / n
    tau = np.arange(mag.size) / (padded * dnu)
    factor = (2.0 * n_group) if convention == "reflection" else n_group
    z = C_VAC * tau / factor

    dz = C_VAC / (factor * dnu_total_hz)          # bin spacing at pad_factor 1
    with np.errstate(divide="ignore"):
        db = 20 * np.log10(mag / mag.max())
    return {"z_m": z, "db": db, "mag": mag, "tau_s": tau,
            "dz_m": dz, "z_max_m": float(z[-1]), "dnu_total_hz": dnu_total_hz,
            "n_points": n, "padded": padded, "beta": beta,
            "convention": convention, "n_group": n_group}


def list_peaks(res, *, n_peaks=8, min_db=-45.0, exclude_below_m=0.02):
    """Strongest peaks, reported in both conventions so the axis label cannot
    mislead: an MZI imbalance of L fibre metres sits at L/2 on the reflection
    axis."""
    z, db = res["z_m"], res["db"]
    idx, _ = find_peaks(db, height=min_db,
                        distance=max(int(0.5 * res["padded"] / res["n_points"]), 1))
    idx = idx[np.abs(z[idx]) > exclude_below_m]
    idx = idx[np.argsort(db[idx])[::-1][:n_peaks]]
    idx = idx[np.argsort(z[idx])]
    rows = []
    for i in idx:
        refl = z[i] if res["convention"] == "reflection" else z[i] / 2
        rows.append({"z_axis_m": float(z[i]), "db": float(db[i]),
                     "reflection_m": float(refl),
                     "fibre_length_m": float(refl * 2),
                     "tau_ns": float(res["tau_s"][i] * 1e9)})
    return rows


# ---------------------------------------------------------------------------
# Optical frequency span
# ---------------------------------------------------------------------------

def dnu_from_wavelengths(start_nm: float, stop_nm: float) -> float:
    """Optical frequency span between two wavelengths (exact, not the c*dl/l^2
    approximation)."""
    return abs(C_VAC / (min(start_nm, stop_nm) * 1e-9)
               - C_VAC / (max(start_nm, stop_nm) * 1e-9))


def dnu_from_calibration(cal, n_pts: int, rate_hz: float):
    """Span actually covered by the capture, from the measured lambda(t).

    More trustworthy than the commanded sweep range: the calibration is what
    established that the EXFO does not sweep at the commanded speed. Returns
    None when no calibration matches.
    """
    from lina.analysis.lina_wl_cal import wavelength_axis, sweep_window
    lo, hi = sweep_window(cal, n_pts, rate_hz)
    wl = wavelength_axis(cal, n_pts, rate_hz)[lo:hi]
    return dnu_from_wavelengths(float(wl[0]), float(wl[-1])), (lo, hi)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_ofdr(res, meas_sig, aux_sig, rate_hz, phase, meta, out_png, peaks,
              xlim=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.4])

    # raw interferograms
    ax = fig.add_subplot(gs[0, 0])
    t = np.arange(meas_sig.size) / rate_hz
    ax.plot(t, meas_sig, lw=0.4, label="measurement MZI")
    ax.plot(t, aux_sig, lw=0.4, alpha=0.7, label="aux MZI")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("AC signal (mW)")
    ax.set_title("interferograms over the sweep")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # fringe zoom
    ax = fig.add_subplot(gs[0, 1])
    mid = meas_sig.size // 2
    span = min(int(0.004 * rate_hz), meas_sig.size // 4)
    sl = slice(mid, mid + span)
    ax.plot(t[sl] * 1000, meas_sig[sl], lw=0.8, label="measurement")
    ax.plot(t[sl] * 1000, aux_sig[sl], lw=0.8, alpha=0.8, label="aux")
    ax.set_xlabel("time (ms)")
    ax.set_title("fringe zoom — both must be resolved, not aliased")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    # aux phase nonlinearity = the sweep nonlinearity, measured
    ax = fig.add_subplot(gs[1, 0])
    lin = np.linspace(phase[0], phase[-1], phase.size)
    resid_frac = (phase - lin) / (phase[-1] - phase[0])
    ax.plot(t, resid_frac * 100, lw=0.8, color="tab:red")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("deviation (% of span)")
    ax.set_title(f"aux phase vs a linear ramp = the sweep's nonlinearity\n"
                 f"peak {np.abs(resid_frac).max() * 100:.2f} % of the frequency span")
    ax.grid(alpha=0.3)

    # local sweep rate from the aux phase
    ax = fig.add_subplot(gs[1, 1])
    k = max(int(0.01 * rate_hz), 10)
    inst = np.gradient(phase[::k], t[::k])
    inst = inst / np.median(inst) * 100
    ax.plot(t[::k], inst, lw=0.9, color="tab:blue")
    ax.axhline(100, color="k", lw=0.8, ls=":")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("instantaneous rate (% of median)")
    ax.set_title("sweep speed over time (from the aux fringe)")
    ax.grid(alpha=0.3)

    # the spectrum
    ax = fig.add_subplot(gs[2, :])
    ax.plot(res["z_m"], res["db"], lw=0.6)
    for p in peaks[:6]:
        ax.annotate(f"{p['z_axis_m']:.4f} m\n{p['db']:.1f} dB",
                    (p["z_axis_m"], p["db"]), fontsize=7,
                    xytext=(0, 12), textcoords="offset points", ha="center")
        ax.plot(p["z_axis_m"], p["db"], "v", ms=5, color="tab:red")
    ax.set_xlabel("test fibre length (m)" if res.get("axis_note")
                  else ("MZI path imbalance (m)"
                        if res["convention"] == "transmission"
                        else "distance (m, one-way / reflection convention)"))
    ax.set_ylabel("amplitude (dB rel. maximum)")
    ax.set_xlim(*(xlim if xlim else (0, res["z_max_m"])))
    if xlim:
        # Re-scale the dB reference to what is inside the window, so a peak
        # outside it cannot flatten everything shown here against the ceiling.
        vis = (res["z_m"] >= xlim[0]) & (res["z_m"] <= xlim[1])
        if vis.any():
            top = float(np.max(res["db"][vis]))
            ax.set_ylim(min(-90.0, top - 90), top + 8)
    ax.set_title(f"{meta.get('source', 'capture')} | aux-referenced, kaiser β="
                 f"{res['beta']:g}, dz {res['dz_m'] * 1e6:.1f} um, "
                 f"Nyquist {res['z_max_m']:.2f} m, n_g {res['n_group']:.4f}")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    print(f"[out] {out_png}")


def decimate_max(z, db, n_out=3000):
    """Thin a spectrum for display without losing peaks.

    Plotting every FFT bin draws the noise as a solid filled band, because a
    2 M-point spectrum has far more points than the figure has pixels and each
    pixel column ends up spanning the full noise spread. Taking the MAXIMUM
    per block keeps every peak at its true height (a peak is the maximum of
    its block) while the noise collapses to a thin line at its local maximum.
    The per-block median is returned as well, as the actual noise floor.
    """
    order = np.argsort(z)                  # mirrored axes run backwards
    z, db = np.asarray(z)[order], np.asarray(db)[order]
    if z.size <= 2 * n_out:
        return z, db, None
    edges = np.linspace(0, z.size, int(n_out) + 1).astype(int)
    zc, mx, md = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        seg = db[lo:hi]
        zc.append(0.5 * (z[lo] + z[hi - 1]))
        mx.append(seg.max())
        md.append(np.median(seg))
    return np.asarray(zc), np.asarray(mx), np.asarray(md)


def auto_xlim(peaks, fallback, margin_m=0.20, min_span_m=0.6):
    """x range that frames the detected peaks instead of the whole Nyquist span.

    A distance spectrum runs to metres while the features of interest sit in a
    few centimetres, so the default full-range view hides them. This brackets
    the peaks with a margin and enforces a minimum span so a single peak does
    not end up absurdly zoomed.
    """
    zs = [p["z_axis_m"] for p in peaks] if peaks else []
    if not zs:
        return fallback
    lo, hi = min(zs) - margin_m, max(zs) + margin_m
    if hi - lo < min_span_m:
        mid = 0.5 * (lo + hi)
        lo, hi = mid - min_span_m / 2, mid + min_span_m / 2
    return (max(lo, 0.0), hi)


def plot_two_spectra(spec_aux, peaks_aux, spec_meas, peaks_meas,
                     aux_xlim, meas_xlim, meta, out_png):
    """Just the two distance spectra, aux on top, measurement below.

    Each is computed with the OTHER channel pair as its phase reference, so the
    two panels are independent measurements that cross-check each other: the
    fringe rate of one pair predicts the peak position of the other.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(14, 8))
    for ax, res, peaks, xlim, title in (
            (axes[0], spec_aux, peaks_aux, aux_xlim,
             "reference / aux MZI  (Ch2, Ch4)  —  its own 4 m line at 4.16 m"),
            (axes[1], spec_meas, peaks_meas, meas_xlim,
             "measurement MZI  (Ch1, Ch3)  —  axis reads TEST FIBRE length")):
        z, db = res["z_m"], res["db"]
        if xlim is None:
            xlim = auto_xlim(peaks, (float(np.min(z)), float(np.max(z))))
        vis = (z >= xlim[0]) & (z <= xlim[1])
        # Reference the dB scale to the strongest peak INSIDE the window.
        top = float(np.max(db[vis])) if vis.any() else 0.0
        zc, mx, md = decimate_max(z, db - top)
        if md is not None:
            ax.fill_between(zc, -200, md, color="0.75", alpha=0.5, lw=0,
                            label="noise floor (median)")
        ax.plot(zc, mx, lw=0.8, label="peak hold")
        for p in peaks[:3]:
            if xlim[0] <= p["z_axis_m"] <= xlim[1]:
                ax.plot(p["z_axis_m"], p["db"] - top, "v", ms=6, color="tab:red")
                ax.annotate(f"{p['z_axis_m']:.4f} m", (p["z_axis_m"], p["db"] - top),
                            xytext=(0, 14), textcoords="offset points",
                            ha="center", fontsize=9, fontweight="bold")
        ax.set_xlim(*xlim)
        ax.set_ylim(-90, 10)
        ax.set_xlabel("test fibre length (m)" if res.get("axis_note")
                      else "MZI path imbalance (m)")
        ax.set_ylabel("amplitude (dB rel. peak)")
        ax.set_title(f"{title}   —   dz {res['dz_m'] * 1e6:.1f} µm, "
                     f"kaiser β={res['beta']:g}, n_g {res['n_group']:.4f}")
        ax.grid(alpha=0.3)
    fig.suptitle(f"{meta.get('source', '')} — aux-referenced OFDR", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=130)
    print(f"[out] {out_png}")


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def process(traces, rate_hz, meta, a):
    meas_sig = combine_pair(traces, a.meas, a.balance, "measurement MZI")
    aux_sig = combine_pair(traces, a.ref, a.balance, "aux MZI")

    # Crop to the part of the buffer that is actually the sweep. With a
    # calibration this window is measured; without one, use the whole buffer
    # and say so (the post-sweep tail would otherwise be treated as sweep).
    cal = None
    dnu = None
    crop = None
    if not a.no_cal:
        try:
            from lina.analysis.lina_wl_cal import find_calibration, describe
            cal = find_calibration(meta.get("requested_start_nm", a.start),
                                   meta.get("requested_stop_nm", a.stop),
                                   meta.get("requested_speed_nm_s", a.speed))
            if cal is not None:
                dnu, crop = dnu_from_calibration(cal, meas_sig.size, rate_hz)
                print(f"[cal] {describe(cal)}")
                print(f"[cal] sweep window: samples {crop[0]}..{crop[1]} "
                      f"({crop[0] / rate_hz:.3f}..{crop[1] / rate_hz:.3f} s); "
                      f"Δν = {dnu / 1e12:.4f} THz")
        except Exception as e:
            print(f"[cal] unavailable: {e}")
    if crop is not None:
        meas_sig = meas_sig[crop[0]:crop[1]]
        aux_sig = aux_sig[crop[0]:crop[1]]
    if a.highpass_hz > 0:
        meas_sig = highpass(meas_sig, rate_hz, a.highpass_hz)
        aux_sig = highpass(aux_sig, rate_hz, a.highpass_hz)
        print(f"[sig] high-passed both interferograms above "
              f"{a.highpass_hz:.0f} Hz to remove the power envelope")
    if dnu is None:
        start = meta.get("actual_start_nm", a.start)
        stop = meta.get("actual_stop_nm", a.stop)
        dnu = dnu_from_wavelengths(start, stop)
        print(f"[cal] no calibration — using the commanded range "
              f"{start:.1f}-{stop:.1f} nm, Δν = {dnu / 1e12:.4f} THz. The "
              f"distance scale inherits the EXFO's speed error (1-2 %).")

    # Aux quality — a smeared or aliased aux fringe invalidates everything.
    f_aux, tone_frac = aux_quality(aux_sig, rate_hz)
    print(f"[aux] dominant fringe {f_aux / 1e3:.3f} kHz "
          f"({100 * f_aux / (rate_hz / 2):.1f} % of Nyquist), "
          f"{100 * tone_frac:.1f} % of energy in that line")
    if f_aux > 0.45 * rate_hz:
        print("[aux] WARNING: aux fringe is near Nyquist — it may be aliased. "
              "Raise the sample rate or lower the sweep speed.")
    if tone_frac < 0.2:
        print("[aux] WARNING: aux is not single-tone; the phase reference is "
              "unreliable. Check the aux MZI and its channel assignment.")

    # Aux delay implied by that fringe, as a sanity check on the aux fibre.
    dnu_dt = dnu / (meas_sig.size / rate_hz)
    tau_aux = f_aux / dnu_dt
    print(f"[aux] implied aux delay {tau_aux * 1e9:.2f} ns = "
          f"{C_VAC * tau_aux / N_GROUP_SMF28:.3f} m of fibre imbalance "
          f"(n_g {N_GROUP_SMF28})")

    phase = aux_phase(aux_sig)
    f_phase = float((phase[-1] - phase[0]) / (2 * np.pi) /
                    (aux_sig.size / rate_hz))
    print(f"[aux] fringe rate from the unwrapped phase: {f_phase / 1e3:.3f} kHz "
          f"({(phase[-1] - phase[0]) / (2 * np.pi):.0f} fringes) — this and the "
          f"spectral estimate above should agree")
    if a.no_aux:
        print("[aux] --no-aux: FFT on the raw time samples (uniform in time, "
              "NOT in frequency) — peaks will be smeared by the nonlinearity")
        uniform = meas_sig - meas_sig.mean()
        bad, span = 0.0, float(phase[-1] - phase[0])
    else:
        uniform, bad, span = resample_on_aux(meas_sig, phase)
        print(f"[aux] resampled on {span / (2 * np.pi):.0f} aux fringes"
              + (f"; {100 * bad:.2f} % non-monotonic phase samples"
                 if bad > 0 else ""))
        if bad > 0.01:
            print("[aux] WARNING: >1 % of the aux phase is non-monotonic — "
                  "the resampling is unreliable (noisy or aliased aux).")

    res = fft_to_distance(uniform, dnu_total_hz=dnu, n_group=a.n_group,
                          beta=a.beta, pad_factor=a.pad_factor,
                          convention=a.convention)
    # Read the TEST FIBRE's length directly off the axis.
    #
    # The interferometer's own path imbalance adds to every fibre-dependent
    # peak, so the raw axis puts a 1.040 m fibre at 0.574 + 1.040 = 1.614 m.
    # Measured on this setup (three states: two fibres and none):
    #   measurement pair:  peak = 0.574 m (its own imbalance) + L_fibre
    #   aux pair:          peak = 2.081 m (its own imbalance) - L_fibre
    # so subtracting 0.574 makes the measurement axis read fibre length, and
    # MIRRORING about 2.081 does the same for the aux pair (its peak moves
    # the other way). Both offsets are the peak positions measured with no
    # fibre connected, in this convention.
    off = float(getattr(a, "offset_m", 0.0) or 0.0)
    if getattr(a, "mirror", False):
        res["z_m"] = off - res["z_m"]
        res["axis_note"] = f"mirrored about {off:.5f} m"
    elif off:
        res["z_m"] = res["z_m"] - off
        res["axis_note"] = f"offset by -{off:.5f} m"
    # With an offset subtracted, the interferometer's own peak sits at z≈0 and
    # its skirt (measured -28 dB at -0.07 m) can outrank the fibre peak, so the
    # near-zero region has to be excluded or the listing reports the artefact.
    peaks = list_peaks(res, n_peaks=a.n_peaks, min_db=a.min_db,
                       exclude_below_m=(0.15 if getattr(a, "offset_m", 0.0)
                                        else 0.02))

    print(f"\n[fft] {res['n_points']} points, pad x{a.pad_factor}, "
          f"dz {res['dz_m'] * 1e6:.2f} um, Nyquist {res['z_max_m']:.3f} m")
    print(f"\n  {'axis (m)':>10} {'dB':>7} {'delay (ns)':>11} "
          f"{'refl. (m)':>10} {'fibre L (m)':>12}")
    for p in peaks:
        print(f"  {p['z_axis_m']:10.5f} {p['db']:7.1f} {p['tau_ns']:11.3f} "
              f"{p['reflection_m']:10.5f} {p['fibre_length_m']:12.5f}")
    if not peaks:
        print("  (no peaks above the threshold — lower --min-db)")
    return res, meas_sig, aux_sig, phase, peaks


def stamp():
    return datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def cmd_measure(a):
    from lina_sweep_test import open_devices, close_devices, run_sweep, save_raw
    laser = pm = tof = None
    try:
        laser, pm, tof = open_devices(a.exfo, a.coredaq, None)
        res = run_sweep(laser, pm, start_nm=a.start, stop_nm=a.stop,
                        speed_nm_s=a.speed, channels=a.channels,
                        rate_hz=a.rate, n_samples=a.samples, pad_s=a.pad,
                        park_settle_s=a.park_settle, restore_wl=True)
    finally:
        close_devices(laser, pm, tof)
    os.makedirs(a.out, exist_ok=True)
    npz = save_raw(res, a.out, f"ofdr_{a.tag}" if a.tag else "ofdr")
    traces, rate, meta = res["traces_mw"], res["rate_hz"], res
    fft, meas, aux, phase, peaks = process(traces, rate, meta, a)
    meta = dict(meta); meta["source"] = os.path.basename(npz)
    plot_ofdr(fft, meas, aux, rate, phase, meta,
              os.path.join(a.out, f"ofdr_{a.tag + '_' if a.tag else ''}{stamp()}.png"),
              peaks, xlim=a.xlim)


def cmd_process(a):
    from lina_sweep_test import load_raw
    meta = load_raw(a.raw)
    traces = meta["traces_mw"]
    rate = meta["rate_hz"]
    fft, meas, aux, phase, peaks = process(traces, rate, meta, a)
    meta["source"] = os.path.basename(a.raw)
    out = a.png or os.path.splitext(a.raw)[0] + f"_fft_{stamp()}.png"
    plot_ofdr(fft, meas, aux, rate, phase, meta, out, peaks, xlim=a.xlim)


def cmd_peaks(a):
    """Both spectra from ONE capture, stacked: aux on top, measurement below."""
    from lina_sweep_test import load_raw
    meta = load_raw(a.raw)
    traces, rate = meta["traces_mw"], meta["rate_hz"]
    meta["source"] = os.path.basename(a.raw)

    pair_a = [int(c) for c in a.meas_channels.split(",")]
    pair_b = [int(c) for c in a.aux_channels.split(",")]

    # The test fibre is traversed TWICE (measured: a 1 m patch cord adds
    # 2.079 m of path, a 3 m one adds 6.051 m), so the measurement panel uses
    # the reflection convention (2*n_g) and subtracts the interferometer's own
    # imbalance — its axis then reads test-fibre length directly.
    print("=== measurement MZI (referenced to the other pair) ===")
    a.meas, a.ref = pair_a, pair_b
    a.offset_m, a.mirror = a.meas_offset_m, False
    a.convention = "reflection"
    spec_meas, *_ , peaks_meas = process(traces, rate, meta, a)
    print("\n=== reference / aux MZI (referenced to the other pair) ===")
    a.meas, a.ref = pair_b, pair_a
    # The aux pair's peak moves DOWN as the fibre gets longer, so its axis is
    # mirrored about its own no-fibre peak to also read fibre length.
    # The aux delay line is traversed ONCE, so it keeps the transmission
    # convention (n_g): its own imbalance then reads its physical fibre length
    # (~4.16 m for the nominal 4 m line). Its fibre-dependent peak sits at
    # |aux_own - 2*L_fibre|, which FOLDS once the fibre exceeds aux_own/2 —
    # that is why this panel is not turned into a fibre-length axis.
    a.offset_m, a.mirror = a.aux_offset_m, bool(a.aux_offset_m)
    a.convention = "transmission"
    spec_aux, *_, peaks_aux = process(traces, rate, meta, a)

    out = a.png or os.path.join(os.path.dirname(a.raw) or ".",
                                f"PEAKS_{stamp()}.png")
    plot_two_spectra(spec_aux, peaks_aux, spec_meas, peaks_meas,
                     a.aux_xlim, a.meas_xlim, meta, out)


def cmd_overlay(a):
    """Several captures overlaid, one colour each — the direct comparison.

    Both channel pairs are shown: each is computed with the OTHER pair as its
    phase reference, so the two panels are independent measurements of the
    same fibre and must agree.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from lina_sweep_test import load_raw

    entries = []
    for item in a.raw:
        label, _, path = item.partition("=")
        if not path:
            label, path = os.path.basename(item).split("_")[0], item
        entries.append((label, path))

    pair_m = [int(c) for c in a.meas_channels.split(",")]
    pair_a = [int(c) for c in a.aux_channels.split(",")]

    fig, axes = plt.subplots(2, 1, figsize=(14, 9))
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    summary = []
    for k, (label, path) in enumerate(entries):
        meta = load_raw(path)
        traces, rate = meta["traces_mw"], meta["rate_hz"]
        for ax, (meas, ref, off, mirror, xlim, conv, title) in zip(axes, (
                (pair_m, pair_a, a.meas_offset_m, False, a.meas_xlim, "reflection",
                 f"measurement MZI (Ch{pair_m[0]},Ch{pair_m[1]}) — axis reads "
                 f"TEST FIBRE length"),
                (pair_a, pair_m, a.aux_offset_m, bool(a.aux_offset_m),
                 a.aux_xlim, "transmission",
                 f"reference / aux MZI (Ch{pair_a[0]},Ch{pair_a[1]}) — its own "
                 f"4 m line at 4.16 m"))):
            a.meas, a.ref = meas, ref
            a.offset_m, a.mirror = off, mirror
            a.convention = conv
            import contextlib, io as _io
            with contextlib.redirect_stdout(_io.StringIO()):
                res, *_, peaks = process(traces, rate, meta, a)
            z, db = res["z_m"], res["db"]
            vis = (z >= xlim[0]) & (z <= xlim[1])
            top = float(np.max(db[vis])) if vis.any() else 0.0
            zc, mx, _ = decimate_max(z, db - top, a.decimate)
            ax.plot(zc, mx, lw=0.9, alpha=0.9,
                    color=colours[k % len(colours)], label=label)
            ax.set_xlim(*xlim)
            ax.set_ylim(-70, 8)
            ax.set_xlabel("test fibre length (m)" if res.get("axis_note")
                          else "MZI path imbalance (m)")
            ax.set_ylabel("amplitude (dB rel. peak in window)")
            ax.set_title(f"{title}   —   dz {res['dz_m'] * 1e6:.1f} µm, "
                         f"n_g {res['n_group']:.4f}"
                         + (f", {res['axis_note']}" if res.get("axis_note") else ""))
            ax.grid(alpha=0.3)
            if ax is axes[0]:
                fib = [p for p in peaks if xlim[0] < p["z_axis_m"] < xlim[1]]
                fib = max(fib, key=lambda p: p["db"]) if fib else None
                summary.append((label, fib["z_axis_m"] if fib else None,
                                fib["db"] if fib else None))
    for ax in axes:
        ax.legend(fontsize=9, ncol=len(entries))
    fig.suptitle("OFDR: test fibre length, aux-referenced — "
                 "same processing for every capture", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(a.out, exist_ok=True)
    out = a.png or os.path.join(a.out, f"OVERLAY_{stamp()}.png")
    fig.savefig(out, dpi=130)
    print(f"[out] {out}")
    print(f"\n  {'dataset':>10} {'strongest peak (m)':>20} {'dB':>6}")
    for label, z, db in summary:
        z_txt = "—" if z is None else f"{z:.5f}"
        db_txt = "" if db is None else f"{db:.1f}"
        print(f"  {label:>10} {z_txt:>20} {db_txt:>6}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--exfo", default="exfo-1")
    p.add_argument("--coredaq", default="coredaq-1")
    p.add_argument("--out", default="lina_ofdr_data")
    p.add_argument("--tag", default="")
    # channel roles
    p.add_argument("--meas", default="1,3",
                   help="measurement-MZI channels (default 1,3)")
    p.add_argument("--ref", default="2,4",
                   help="aux/reference-MZI channels (default 2,4)")
    p.add_argument("--balance", choices=["auto", "always", "never"], default="auto",
                   help="subtract each channel pair (balanced detection). "
                        "auto: only when the pair is anti-correlated")
    # processing
    p.add_argument("--n-group", type=float, default=N_GROUP_SMF28,
                   help=f"fibre group index (default {N_GROUP_SMF28}, SMF-28 @1550)")
    p.add_argument("--beta", type=float, default=10.0, help="Kaiser beta")
    p.add_argument("--pad-factor", type=int, default=1,
                   help="zero-pad multiple (peak interpolation, not resolution)")
    p.add_argument("--convention", choices=["reflection", "transmission"],
                   default="reflection",
                   help="distance axis. transmission (default) uses n_g, so an "
                        "MZI path imbalance of L metres of fibre reads L on the "
                        "axis — the fibre length itself. reflection uses 2*n_g "
                        "(reflectometry convention), where the same fibre reads "
                        "L/2")
    p.add_argument("--no-aux", action="store_true",
                   help="skip aux resampling (to see what it is worth)")
    p.add_argument("--no-cal", action="store_true",
                   help="ignore the lambda(t) calibration")
    p.add_argument("--highpass-hz", type=float, default=500.0,
                   help="high-pass both interferograms to strip the laser's "
                        "power envelope (0 disables)")
    p.add_argument("--decimate", type=int, default=3000,
                   help="display points per spectrum: block MAXIMUM, so peaks "
                        "keep their height while the noise band collapses to a "
                        "thin line (0 = plot every bin)")
    p.add_argument("--offset-m", type=float, default=0.0,
                   help="subtract this from the distance axis, so the axis "
                        "reads TEST FIBRE length instead of total path "
                        "imbalance (use the peak measured with no fibre)")
    p.add_argument("--mirror", action="store_true",
                   help="axis = offset - z instead of z - offset, for a pair "
                        "whose peak moves DOWN as the fibre gets longer")
    p.add_argument("--xlim", default=None,
                   help="x range of the spectrum panel, e.g. 0,3 (default: the "
                        "full unambiguous range up to Nyquist)")
    p.add_argument("--n-peaks", type=int, default=8)
    p.add_argument("--min-db", type=float, default=-45.0)
    sub = p.add_subparsers(dest="cmd", required=True)

    sp = sub.add_parser("measure", help="capture 4 channels and process")
    sp.add_argument("--start", type=float, default=1520.0)
    sp.add_argument("--stop", type=float, default=1570.0)
    sp.add_argument("--speed", type=float, default=5.0)
    sp.add_argument("--rate", type=int, default=100_000)
    sp.add_argument("--channels", default="1,2,3,4")
    sp.add_argument("--samples", type=int, default=None)
    sp.add_argument("--pad", type=float, default=0.5,
                    help="extra seconds beyond the sweep (small: the tail is "
                         "cropped by the calibration anyway)")
    sp.add_argument("--park-settle", type=float, default=5.0)
    sp.set_defaults(func=cmd_measure)

    sp = sub.add_parser("peaks", help="both spectra stacked (aux above, "
                                      "measurement below) — peaks only")
    sp.add_argument("--raw", required=True)
    sp.add_argument("--png", default=None)
    sp.add_argument("--meas-channels", default="1,3")
    sp.add_argument("--aux-channels", default="2,4")
    sp.add_argument("--meas-xlim", default="auto",
                    help="x range, or 'auto' to frame the detected peaks")
    sp.add_argument("--aux-xlim", default="auto")
    sp.add_argument("--meas-offset-m", type=float, default=0.0,
                    help="measurement pair's own imbalance (subtracted)")
    sp.add_argument("--aux-offset-m", type=float, default=0.0,
                    help="aux pair's own imbalance (axis mirrored about it)")
    sp.add_argument("--start", type=float, default=1520.0)
    sp.add_argument("--stop", type=float, default=1570.0)
    sp.add_argument("--speed", type=float, default=5.0)
    sp.set_defaults(func=cmd_peaks)

    sp = sub.add_parser("overlay", help="several captures in one plot, "
                                        "one colour each")
    sp.add_argument("--raw", nargs="+", required=True,
                    help="captures as label=path (or just paths)")
    sp.add_argument("--png", default=None)
    sp.add_argument("--meas-channels", default="1,3")
    sp.add_argument("--aux-channels", default="2,4")
    sp.add_argument("--meas-xlim", default="auto",
                    help="x range, or 'auto' to frame the detected peaks")
    sp.add_argument("--aux-xlim", default="auto")
    sp.add_argument("--meas-offset-m", type=float, default=0.0)
    sp.add_argument("--aux-offset-m", type=float, default=0.0)
    sp.add_argument("--start", type=float, default=1520.0)
    sp.add_argument("--stop", type=float, default=1570.0)
    sp.add_argument("--speed", type=float, default=5.0)
    sp.set_defaults(func=cmd_overlay)

    sp = sub.add_parser("process", help="re-process a saved capture")
    sp.add_argument("--raw", required=True)
    sp.add_argument("--png", default=None)
    sp.add_argument("--start", type=float, default=1520.0)
    sp.add_argument("--stop", type=float, default=1570.0)
    sp.add_argument("--speed", type=float, default=5.0)
    sp.set_defaults(func=cmd_process)

    a = p.parse_args()
    if a.xlim:
        a.xlim = tuple(float(v) for v in str(a.xlim).split(",")[:2])
    for _k in ("meas_xlim", "aux_xlim"):
        _v = getattr(a, _k, None)
        if _v:
            setattr(a, _k, None if str(_v).lower() == "auto"
                    else tuple(float(x) for x in str(_v).split(",")[:2]))
    a.meas = [int(c) for c in str(a.meas).split(",") if c.strip()]
    a.ref = [int(c) for c in str(a.ref).split(",") if c.strip()]
    if hasattr(a, "channels"):
        a.channels = [int(c) for c in str(a.channels).split(",") if c.strip()]
    a.func(a)


if __name__ == "__main__":
    main()
