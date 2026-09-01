# Analysis modules

The processing shared by the scripts in `../scripts/` and the Lina window in
`LabGUI/dashboards/instruments/lina_window.py`. It lives here, in one place, so
the GUI and the scripts cannot drift apart — a divergence between them would be
invisible and would silently produce different lengths from the same data.

Reasoning: `../report.md`, §3 and §6.

---

## `lina_wl_cal.py` — the λ(t) wavelength axis

A LINEAR coreDAQ capture is one free-running buffer started by a single trigger
edge, so every sample carries a **time**, not a wavelength. This module holds
the measured mapping from time to wavelength and its application.

| Function | Purpose |
|---|---|
| `fit_time_to_wavelength(times, wavelengths, degree)` | The fit. Degree 2 measured at 3.1 pm rms. |
| `find_calibration(start, stop, speed)` | Looks up the calibration for one sweep configuration in `../calibrations/`. Returns `None` rather than guessing — a calibration from a different configuration is worse than none, because it is silently wrong. |
| `wavelength_axis(cal, n_pts, rate_hz)` | Wavelength per sample. Rate-independent: the fit is in seconds, so a 5 kHz calibration applies to a 100 kHz capture (verified to −3 pm). |
| `sweep_window(cal, n_pts, rate_hz)` | Which part of the buffer is the sweep. Measured: t = 0.004 … 9.889 s of a 12.00 s capture; everything after is the return slew. |
| `apply_calibration(cal, traces, rate_hz)` | Both of the above, cropping the traces consistently. |

Tolerances for matching a stored calibration: ±0.05 nm on start/stop, ±2 % on
speed.

---

## `lina_ofdr.py` — aux-referenced OFDR

| Function | Purpose |
|---|---|
| `combine_pair(traces, channels, mode, label)` | Turns a channel pair into one interferogram. Whether it is a balanced pair is decided from the **correlation** (measured −0.88 … −0.99), not assumed. |
| `highpass(sig, rate, cutoff)` | Strips the residual power envelope. Without it the aux "fringe frequency" is estimated as 0.000 kHz and the Hilbert phase is biased. |
| `aux_phase(sig)` | Unwrapped fringe phase — this *is* the optical-frequency axis. |
| `aux_quality(sig, rate)` | Dominant fringe frequency and how single-tone the aux really is. Both are warnings, not decorations: a non-single-tone aux invalidates the phase reference. |
| `resample_on_aux(meas, phase)` | Interpolates the measurement onto a uniform phase grid, and reports what fraction of the phase is non-monotonic. |
| `fft_to_distance(sig, dnu_total_hz, n_group, ...)` | Kaiser-windowed FFT → distance. `dz = c/(factor·n_g·Δν)`. |
| `list_peaks(spec, ...)` | Peaks, reported in both conventions so an axis label cannot mislead. Excludes the region around z = 0, where a subtracted offset leaves the interferometer's own peak (its skirt can outrank the fibre peak). |
| `decimate_max(z, db, n_out)` | Block maximum for display ("peak hold") plus the block median as the true noise floor. |
| `evaluate(...)` | The whole pipeline for one channel pair, referenced to the other. |

### Conventions, and why they differ per channel pair

- The **test fibre** is traversed **twice** (measured: 1 m → 2.079 m of path,
  3 m → 6.051 m), so its axis uses `2·n_g` and subtracts the measurement
  interferometer's own imbalance (0.57396 m) — then it reads fibre length.
- The **aux delay line** is traversed **once** (nominal 4 m reads 4.16 m with
  `n_g`), so it keeps `n_g` and its own axis. Its fibre-dependent peak sits at
  `|own − 2·L|`, which folds once the fibre exceeds half its own imbalance —
  which is why that panel is deliberately not relabelled as a fibre length.
