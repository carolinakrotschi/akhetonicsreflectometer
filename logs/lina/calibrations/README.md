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

The lookup refuses a mismatched calibration rather than guessing. It matches on
start/stop/speed only, not on the coreDAQ — the λ(t) fit is a property of the
EXFO sweep and the trigger timing. The existing file was measured with the Mk1
LINEAR coreDAQ; it should carry over to the Mk2 (a µs-scale trigger-latency
change is fm at 5 nm/s), but **one confirming calibration run on the Mk2 is
still outstanding**. Drift is not
a concern: measured at +0.0 pm/min over 18 minutes, so a calibration does not
go stale within a session.

## `aux_mzi.json` — the aux delay, as a watchdog

**Not present yet.** This file holds one number, `tau_aux_s`, the aux MZI's
delay. It is deliberately NOT what sets any axis: `aux_wavelength_axis` takes
its scale from the commanded sweep endpoints, so a stale value here cannot
corrupt a measurement. Its two jobs are:

* **staleness.** Every sweep implies an aux delay in passing (fringes divided
  by the endpoints' frequency span). Compared against this file, a changed
  bench announces itself instead of being discovered months later.
* **aliasing.** A cleanly aliased fringe is invisible in the capture — folding
  a pure tone gives another pure tone with a smooth, monotonic phase, so every
  intrinsic check passes. The implied delay is the only thing that gives it
  away, and without this file that check cannot run.

**Why it is missing rather than filled in.** The value measured on 2026-09-01
(20.3877 ns, from `lina/analysis/` on the OFDR captures) no longer describes
the bench. Measured 2026-09-15 on two fibres over two sweep windows, the
captures imply **20.635 ns** (1520–1570 nm, anchored by the λ(t) calibration
above) and **20.66 ns** (1505–1625 nm, anchored by the commanded endpoints) —
about +1.2 % on the stored figure, i.e. roughly 5 cm of fibre. An independent
sign of the same change: the measurement pair now correlates +0.4 to +0.5
where `report.md` §6.1 recorded −0.88 to −0.99, so the arm was re-patched.

A 2.07 m reference fibre backs the new figure: it reads 2.069 m with the
correction and 2.037 m without.

Writing the file is therefore a measurement decision, not a code one. Either
take 20.635 ns from the λ(t)-anchored captures, or — better — re-run the
three-minute filter calibration, which produces a fresh λ(t) *and* a fresh
aux delay for the bench as it stands:

```bash
python lina/scripts/lina_sweep_test.py --channels 2 calibrate \
    --start 1520 --stop 1570 --speed 5 --rate 5000 --trust-filter \
    --filter-positions 1530,1535,1540,1545,1550,1555,1560,1565
```
