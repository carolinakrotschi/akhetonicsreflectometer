# Scripts

Seven command-line tools. All of them run from the
repository root and expect the Lina hardware (`exfo-1` + the paired coreDAQ,
plus `tof1550-1` for the calibration) to be free — close AkheLab first, since
LabDevice allows one instance per serial id and the coreDAQ's USB port is
exclusive.

Reasoning behind every design decision: `../report.md`.

---

## The coreDAQ: Mk2 LOG since 2026-09-17

Everything up to 2026-09-17 was measured with the Mk1 **LINEAR** demo unit.
Both coreDAQs are now **Mk2 InGaAs LOG** units (`coredaq-1` = SN0001,
`coredaq-2` = SN0002), told apart only by serial. What that means here:

- **Which unit:** `coredaq_name` of `lina-1` in `Interface/config/OBR_config.json`.
  The GUI preselects it and every script's `--coredaq` defaults to it, so
  re-cabling is that one edit. The driver resolves the port by serial and
  refuses a unit whose serial does not match.
- **Capture mode** is no longer inferred from the variant. The GUI has a
  *Continuous (OFDR)* / *Stepped* combo, defaulting to continuous; the scripts
  were always continuous.
- **Needs py-coredaq 2.4.0** (`requirements.txt`). 1.2.1 rejects Mk2 firmware
  v1.6 with an error advising a firmware update — do **not** flash, it may
  cost the 1 MS/s tier.
- **The signal is still mW**, linearised from the device's 128-point LUT, but:
  no autorange jumps; noise roughly constant in dB (relative) rather than
  constant in mW; never negative — dark samples sit on a clamp floor; and a
  **silent clamp at ~3 mW** (no over-range flag). The GUI and the scripts
  (`run_sweep`, so also `lina_chip_measure` / `lina_voa_series`) count samples
  on either clamp and write it to the saved header / NPZ meta (`detector`,
  built by `lina/analysis/lina_detector.py`); the scripts also print a warning.
  A clamped fringe shows up in OFDR as ghost peaks at multiples of the real
  delay.
- **Bandwidth depends on photocurrent:** 150 kHz above ~10 nW, 50 kHz at 1 nW,
  unspecified below. Fringes at 5 nm/s are ~13 kHz (aux) and ~22 kHz (3 m fibre),
  so this only bites if fringe minima drop below ~1 nW.
- **Rate:** up to 1 MS/s on the HIGH tier (GUI ceiling follows the unit). At
  that rate the ~150 kHz analog bandwidth, not Nyquist, limits the delay range.
- **Responsivity correction** (`λ-cal correction` in the GUI) is now active on
  the continuous path; the scripts stay raw (see `lina_chip_measure.py`). The
  GUI stores the factor's ingredients in `detector.responsivity_correction`,
  so an analysis can divide it back out and compare GUI captures with raw ones
  (the reflectometer repo's `tools/scan_io.py` does that by default).
- **To confirm on first use:** repeat a reference capture (1 m / 3 m patch
  cord, 1520–1570 nm, 5 nm/s, 100 kHz) against 1.0397 m / 3.0256 m, and re-run
  one TOF1550 calibration to check the λ(t) fit carried over from the Mk1.

---

## Why these exist as scripts and not only as GUI features

The GUI could not be debugged interactively from a terminal session, and the
problem being investigated (a wrong wavelength axis) needed dozens of captures
with varied processing. A script gives that; a Qt window does not. The
processing itself lives in `../analysis/`, so the GUI and these scripts share
one implementation and cannot drift apart.

---

## `lina_sweep_test.py` — capture and wavelength calibration

The headless version of the Lina Response Sweep (the continuous coreDAQ path),
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

## `lina_chip_measure.py` — a guided measurement session

The Lina GUI's Response Sweep, repeated over channels and sweep windows without
the clicking. It asks which channel(s) to measure, then per channel: prompts
for the patch change with the emission off, runs every window unattended, and
names the files itself. Afterwards it asks again, so one session covers as
many channels as you like; an empty answer ends it.

Same capture path as `lina_sweep_test sweep`, same Auto-fit sample count as
the GUI, same wavelength axis (λ(t) calibration → aux fringe → linear).

```bash
python lina/scripts/lina_chip_measure.py --dry-run          # the plan, no hardware
python lina/scripts/lina_chip_measure.py --chip 2680 --label ligentechhi
python lina/scripts/lina_chip_measure.py --fibers 44,48 --connected 44
```

One `<chip>_<label>_<YYYYmmdd-HHMMSS>_fiber<NN>_<window>.npz` per capture in
`--out` (default `lina/raw_data`), plus a quicklook PNG. The NPZ layout is the
one `load_raw` reads, with the wavelength axis in `wl_nm` and the crop it
belongs to in `meta["crop"]`, so the OFDR tools take these files directly:

```bash
python lina/scripts/lina_ofdr_test.py peaks --raw lina/raw_data/<file>.npz
```

---

## `lina_voa_series.py` — a VOA attenuation series

Ten (or N) identical wide-window sweeps on one fibre, one per setting of the
manual variable optical attenuator. It stops before every capture, shows the
live CW power so the knob can be turned to a visible level, and records what
you type there as that capture's VOA setting — so the setting ends up in the
file rather than in your memory.

Same capture path, wavelength axis and NPZ layout as `lina_chip_measure.py`,
whose helpers it imports rather than copies. Always the wide 1505–1625 nm
window, which has no λ(t) calibration and is carried by the aux fringe.

```bash
python lina/scripts/lina_voa_series.py --dry-run        # the plan, no hardware
python lina/scripts/lina_voa_series.py --fiber 48 --runs 10
python lina/scripts/lina_voa_series.py --fiber 48 --start-index 8 --runs 3
```

Files are `<chip>_<YYYYmmdd-HHMMSS>_fiber<NN>_fullnmrange_<i>.npz`, the
timestamp being the session's, so a series is one contiguous block of names
differing only in the trailing index. Two figures close it: `_stacked.png`
(one panel per run, shared axes) and `_levels.png` (mean power against run
index). The stacked figure is linear with a min/max envelope up to ~13 dB of
spread and logarithmic with peak hold past it — `--log-y` / `--linear-y`
force either.

`--replot "<glob>"` rebuilds both figures from saved NPZs without hardware, so
a session interrupted at run 7 still gets its figure.

---

## `lina_npz_to_json.py` — captures in the house JSON format

Turns the NPZs into `{"header": …, "data": [{…}]}` JSON, the format the rest of
AKHExperiment stores measurements in and the one the GUI's "Save Data" writes.
No hardware, so it can be run long after the session.

```bash
python lina/scripts/lina_npz_to_json.py                 # all of lina/raw_data
python lina/scripts/lina_npz_to_json.py --out L:/Measurements/2680 <file>.npz
```

Columns are cropped to `meta["crop"]`, so `Wavelength [nm]` and the `Ch<n> [mW]`
columns all have the same length. It is text: ~50 MB per 1 M samples, ~123 MB
per 2.4 M. Numbers are written at 1 fm / 6 significant figures, past what the
instruments resolve but well short of float64's 17 digits.

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
