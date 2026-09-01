# Calibrations

One file per sweep configuration. The GUI and the scripts look them up here
automatically by matching start wavelength, stop wavelength and sweep speed
(±0.05 nm, ±2 %).

Reasoning: `../report.md`, §3.

## Contents

| File | Configuration | Fit |
|---|---|---|
| `lina_wl_cal_1520.0-1570.0nm_5.00nms.json` | 1520 → 1570 nm at 5 nm/s | degree 2, **3.1 pm rms**, 5.1 pm max |

Measured on 2026-09-01 with the TOF1550 filter at 8 positions (1530…1565 nm),
`exfo-1` + `coredaq-1` at 5 kHz. Raw data: `../data/01_wavelength_calibration/`.

## What is inside

- `fit.coeffs` — λ(t) polynomial, **nm versus seconds from the trigger edge**.
  Being in seconds makes it rate-independent (verified: measured at 5 kHz,
  applied to a 100 kHz capture to −3 pm).
- `derived.buffer_time_s` — when the sweep reaches the start wavelength
  (+0.0041 s, i.e. essentially at the trigger edge).
- `derived.sweep_end_s` — when the sweep ends (9.8895 s of a 12.00 s buffer);
  everything after that is the laser's return slew and gets cropped.
- `derived.effective_speed_nm_s` — 5.1174 nm/s against 5.0000 commanded.
- `points` — the 8 measured (time, wavelength) pairs the fit came from.

## Validity

**A calibration is valid only for the configuration it was measured at.** Any
other start/stop/speed needs its own ~3-minute run:

```bash
python lina/scripts/lina_sweep_test.py --channels 2 calibrate \
    --start <nm> --stop <nm> --speed <nm/s> --rate 5000 --trust-filter \
    --filter-positions <list of nm inside the range>
```

The lookup refuses a mismatched calibration rather than guessing. Drift is not
a concern: measured at +0.0 pm/min over 18 minutes, so a calibration does not
go stale within a session.
