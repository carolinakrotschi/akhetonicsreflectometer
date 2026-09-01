"""
OFDR evaluation for Lina (EXFO sweep + LINEAR coreDAQ), shared by the
measurement script (``scripts/lina_ofdr_test.py``) and the GUI window
(``LabGUI/dashboards/instruments/lina_window.py``).

The physics, all of it established by measurement on this setup:

  * a LINEAR coreDAQ capture is uniform in TIME, but an FFT only resolves
    delays if the signal is uniform in OPTICAL FREQUENCY. The EXFO is not
    (1-2 % fast, slightly non-uniform — see Analysis/lina_wl_cal.py).
  * the aux MZI has a FIXED path imbalance, so its unwrapped fringe phase IS
    the optical-frequency axis: phi(t) = 2*pi*tau_aux*nu(t). Resampling the
    measurement onto a uniform phi grid makes it uniform in nu.
  * dz = c/(factor*n_g*dnu_total) and z_max = dz*N/2, where factor is 2*n_g
    when the path is traversed twice (the test fibre here) or n_g when once
    (the aux delay line, whose nominal 4 m reads 4.16 m).
  * the interferometer's own imbalance adds to every fibre-dependent peak, so
    subtracting it makes the axis read test-fibre length directly. Measured:
    measurement pair peak = 0.574 m + L_fibre (in the 2*n_g convention),
    aux pair peak = |4.162 m - 2*L_fibre| (in the n_g convention).

Verified against three fibre states: nominal 1 m patch cords read 1.0397 m
and 1.0553 m, a nominal 3 m one reads 3.0256 m, and with no fibre connected
the fibre peak is absent. Both channel pairs, processed with each other as
reference, agree.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import hilbert, find_peaks, butter, filtfilt

C_VAC = 299_792_458.0        # m/s
N_GROUP_SMF28 = 1.4682       # SMF-28 group index at 1550 nm


def combine_pair(traces: dict, channels: list, mode: str, label: str,
                 log=print):
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
        log(f"{label}: single channel {channels[0]}")
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
    log(f"{label}: Ch{channels[0]}/Ch{channels[1]} correlation "
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


def dnu_from_wavelengths(start_nm: float, stop_nm: float) -> float:
    """Optical frequency span between two wavelengths (exact, not the c*dl/l^2
    approximation)."""
    return abs(C_VAC / (min(start_nm, stop_nm) * 1e-9)
               - C_VAC / (max(start_nm, stop_nm) * 1e-9))

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


def evaluate(traces: dict, rate_hz: float, *, meas_channels, aux_channels,
             dnu_total_hz: float, n_group: float = N_GROUP_SMF28,
             beta: float = 10.0, pad_factor: int = 4,
             convention: str = "reflection", offset_m: float = 0.0,
             mirror: bool = False, highpass_hz: float = 500.0,
             use_aux: bool = True, n_peaks: int = 6, min_db: float = -35.0,
             balance: str = "auto", log=print):
    """One channel pair's distance spectrum, referenced to the other pair.

    Args:
        traces: ``{channel: 1-D array}`` — the capture, already cropped to the
            sweep window (the post-sweep tail must not be included).
        rate_hz: sample rate of the capture.
        dnu_total_hz: optical frequency span the samples cover. This alone
            sets the distance scale; take it from the measured lambda(t)
            calibration, not from the commanded sweep range, which is 1-2 % off.
        offset_m / mirror: subtract (or mirror about) the interferometer's own
            imbalance so the axis reads test-fibre length. See the module
            docstring.
        log: where diagnostics go (``print``, or a GUI status callback).

    Returns:
        ``(spectrum, peaks, diagnostics)``.
    """
    # Both pairs use "auto": whether each is a balanced pair is decided from
    # its own correlation. ("aux" is a LABEL, not a mode — passing it here
    # silently disabled balanced subtraction for the reference pair, which
    # wrecked the phase reference.)
    meas = combine_pair(traces, meas_channels, balance, "measurement", log=log)
    aux = combine_pair(traces, aux_channels, balance, "reference/aux", log=log)
    if highpass_hz > 0:
        meas = highpass(meas, rate_hz, highpass_hz)
        aux = highpass(aux, rate_hz, highpass_hz)

    f_aux, tone = aux_quality(aux, rate_hz)
    phase = aux_phase(aux)
    fringes = float((phase[-1] - phase[0]) / (2 * np.pi))
    diag = {"aux_fringe_hz": f_aux, "aux_tone_fraction": tone,
            "aux_fringes": fringes,
            "aux_fringe_hz_from_phase": fringes / (meas.size / rate_hz)}
    log(f"aux fringe {f_aux / 1e3:.2f} kHz ({100 * tone:.0f} % single-tone), "
        f"{fringes:.0f} fringes")
    if f_aux > 0.45 * rate_hz:
        log("WARNING: aux fringe near Nyquist — may be aliased")

    if use_aux:
        uniform, bad, span = resample_on_aux(meas, phase)
        diag["nonmonotonic_fraction"] = bad
        if bad > 0.01:
            log(f"WARNING: {100 * bad:.1f} % of the aux phase is "
                f"non-monotonic — resampling unreliable")
    else:
        uniform = meas - meas.mean()

    spec = fft_to_distance(uniform, dnu_total_hz=dnu_total_hz,
                           n_group=n_group, beta=beta, pad_factor=pad_factor,
                           convention=convention)
    if mirror and offset_m:
        spec["z_m"] = offset_m - spec["z_m"]
        spec["axis_note"] = f"mirrored about {offset_m:.5f} m"
    elif offset_m:
        spec["z_m"] = spec["z_m"] - offset_m
        spec["axis_note"] = f"offset by -{offset_m:.5f} m"
    # Exclude the region around z=0: with an offset subtracted, the
    # interferometer's OWN peak sits there and its skirt (measured: -31 dB at
    # -0.07 m) can outrank the fibre peak, so "strongest peak" would report
    # the artefact instead of the measurement.
    exclude = 0.15 if (offset_m or mirror) else 0.02
    peaks = list_peaks(spec, n_peaks=n_peaks, min_db=min_db,
                       exclude_below_m=exclude)
    return spec, peaks, diag
