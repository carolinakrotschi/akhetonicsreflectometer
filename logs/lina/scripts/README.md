# Scripts

Four command-line tools written on 2026-09-01. All of them run from the
repository root and expect the Lina hardware (`exfo-1` + `coredaq-1`, plus
`tof1550-1` for the calibration) to be free — close AkheLab first, since
LabDevice allows one instance per serial id and the coreDAQ's USB port is
exclusive.

Reasoning behind every design decision: `../report.md`.

---

## Why these exist as scripts and not only as GUI features

The GUI could not be debugged interactively from a terminal session, and the
problem being investigated (a wrong wavelength axis) needed dozens of captures
with varied processing. A script gives that; a Qt window does not. The
processing itself lives in `../analysis/`, so the GUI and these scripts share
one implementation and cannot drift apart.

---

## `lina_sweep_test.py` — capture and wavelength calibration

The headless version of the Lina Response Sweep (the LINEAR coreDAQ path),
plus the λ(t) calibration workflow.

| Subcommand | Purpose |
|---|---|
| `sweep` | One raw, uncropped capture + diagnostics. With `--cal`, reports where a marker lands on the calibrated axis versus the naive one. |
| `filter-cw` | Steps the laser in CW across the TOF1550 passband → the filter's true peak wavelength, independent of the filter's own setting error. |
| `calibrate` | For each filter position: optional CW ground truth + a sweep, locating the marker in time. Fits λ(t), reports effective sweep rate and dead time, and publishes the result to `../calibrations/`. |
| `repeat` | N identical sweeps with the filter parked — run-to-run jitter and drift. |
| `refit` | Re-fits a saved calibration at another polynomial degree (no hardware). |
| `plot-cal` | Plots a saved calibration: fit, residuals, naive-axis error. |
| `apply` | Re-plots a saved capture on a calibrated axis. |

```bash
python lina/scripts/lina_sweep_test.py --channels 2 calibrate \
    --start 1520 --stop 1570 --speed 5 --rate 5000 --trust-filter \
    --filter-positions 1530,1535,1540,1545,1550,1555,1560,1565
```

Keep calibration runs at a few kHz: the fit is in seconds, so it applies at
100 kHz anyway, and the smaller USB transfer is far less likely to stall.

---

## `lina_ofdr_test.py` — OFDR measurement and evaluation

Captures all four channels and turns the two MZI channel pairs into distance
spectra.

| Subcommand | Purpose |
|---|---|
| `measure` | Capture + process + plot in one go. |
| `process` | Re-process a saved capture (different window, padding, channel roles). |
| `peaks` | Both spectra of one capture, stacked: aux above, measurement below. |
| `overlay` | Several captures in one figure, one colour each. |

Options worth knowing:

- `--meas 1,3 --ref 2,4` — channel roles (defaults match this setup).
- `--offset-m` / `--meas-offset-m` — subtract the interferometer's own path
  imbalance so the axis reads test-fibre length.
- `--convention reflection|transmission` — factor 2·n_g or n_g. Reflection is
  the default because the test fibre here is traversed twice.
- `--decimate 3000` — block-maximum ("peak hold") display: peaks keep their
  exact height while the noise collapses to a thin line.
- `--n-group 1.4682` — every length scales linearly with this.

```bash
python lina/scripts/lina_ofdr_test.py --tag myrun measure \
    --start 1520 --stop 1570 --speed 5 --rate 100000
```

---

## `lina_window_test.py` — headless GUI hardware test

Builds the **real** `LinaWindow` with the real `DeviceRegistry` on Qt's
offscreen platform and drives it like a user: selects the capture device, ticks
channels, presses "Response Sweep" (and with `--ofdr`, "OFDR Plot"). The
capture runs through the GUI's own worker, so this tests the GUI rather than a
copy of it.

```bash
python lina/scripts/lina_window_test.py --channels 1,2,3,4 --ofdr
```

This is what found the `"aux"`-in-the-mode-parameter bug and the premature-read
bug (`../report.md` §4.1, §10, §11).

---

## `lina_cal_summary.py` — one overview figure

Reads only saved files, no hardware. Six panels: filter passband, raw capture
with the sweep window marked, the λ(t) fit, residuals per degree, the axis-error
comparison, and stability over time.

```bash
python lina/scripts/lina_cal_summary.py --data lina/data/01_wavelength_calibration
```

---

## Note on the data folders

The scripts write into whatever `--out` points at (defaults: `lina_cal_data`,
`lina_ofdr_data`, `lina_gui_test` relative to the working directory). The data
recorded on 2026-09-01 was moved into `../data/` afterwards; pass `--out` to
write there directly, e.g.
`--out lina/data/02_ofdr_measurements`.
