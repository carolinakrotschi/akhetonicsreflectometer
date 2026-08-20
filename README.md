# OFDR Reflectometer — what's in this folder

Quick orientation. For the physics/background, read `HANDOVER.md` first, then
`docs/handover_explained.md`. For the step-by-step lab procedure, see
`docs/measurement_procedure.md`.

## Active pipeline (main directory)

The setup as of 2026-08-18 uses the **old** trigger mode: EXFO laser sweep,
1 pm wavelength grid, 2 detector channels (Ch1/Ch2).

| Script | What it does |
|---|---|
| `process_reflectogram.py` | Main processing pipeline: loads an EXFO JSON scan, balances the two channels, corrects the frequency axis for the laser's tuning error, FFTs to a reflectogram, reports peak positions/widths, writes a CSV + PNG. This is the everyday tool. |
| `diagnose_artifacts.py` | Five discriminator tests for telling a real reflection apart from an artifact: `pm-am` (additive vs. multiplicative noise), `compare` (is a feature fixed/calibratable or a random per-scan realization?), `fold` (where would a reflector beyond Nyquist alias to?), `tuning` (measure the laser's actual sweep error), `alias` (single-file peak/aliasing-order check, or the definitive two-scan pair test with different wavelength steps). Run these *before* trusting an unexplained feature. |
| `plot_tuning_diagnosis.py` | Rebuilds the two-panel "MZI peak vs. sweep tuning error" figure (the one originally saved as `george/real_data_diagnosis.png`) for any scan, reusing `process_reflectogram.py`'s functions directly. |

The following three belong to the **new**, not-yet-built setup: free-running
acquisition with an auxiliary Mach-Zehnder interferometer ("aux MZI", the
"ruler") providing the frequency reference instead of a wavelength table,
4 channels. **Ch1/Ch2 = aux MZI (the ruler), Ch3/Ch4 = measurement
interferometer with the test fiber** (confirmed 2026-08-18 against the
real wiring — double-check this against the bench before trusting a scan).

| Script | What it does |
|---|---|
| `calc_sweep_parameters.py` | Desk calculator: given point count / sweep speed / span, computes resolution, Nyquist range, and whether the aux MZI's own signal fits inside that range. Use this before touching hardware. |
| `simulate_aux_data.py` | Generates synthetic 4-channel free-running data with *known* ground-truth reflector positions, including realistic imperfections (sweep bow, tuning ripple, gain drift, noise). Used to validate `process_reflectogram_aux.py` against a known answer before trusting it on real data. |
| `check_aux_interferometer.py` | Acceptance test for the aux MZI itself: contrast, monotonicity (laser must never tune backward), calibration (tau_aux via fringe counting), and laser quality (slow bow / fast ripple), measured cleanly through the aux signal. Run this before trusting anything the new setup measures. |
| `process_reflectogram_aux.py` | The new 4-channel processing pipeline: same balanced-subtraction/window/FFT/peak-report core as `process_reflectogram.py`, but the frequency axis comes from the aux MZI's phase instead of a wavelength table. |

## `tools/` — standalone helper scripts

Not part of the day-to-day pipeline above, but useful on their own.

| Script | What it does |
|---|---|
| `plot_band_comparison.py` | Compares the 6-32cm "dirt band" across several raw scans on one plot (built 2026-08-18 for the terminator investigation, see `raw_data/README.md`). |

(2026-08-18: `check_aliasing.py` was merged into `diagnose_artifacts.py` as
its `alias` subcommand — it overlapped with `fold`/`compare`/`pm-am` there.
`simulate_scan_data.py` and `plot_fsr_vs_wavelength.py` were retired as
dead weight — the former only validated the now-legacy CSV pipeline, the
latter was a one-off tied to a folder that no longer exists. All three are
still in git history if ever needed again.)

## `legacy/` — historical record only, no scripts left

`legacy/CHANGES.md` is the only thing left here: the documented bug list
from the very first (buggy) reflectometer attempt. The scripts themselves
(`reflectometer_v1_buggy.py`, `reflectometer_v2_fixed.py`, and the generic
CSV-based `reflectometer_csv_generic.py`) were all deleted 2026-08-18 as
fully superseded by `process_reflectogram.py` — recoverable from git
history if a plain-CSV or pre-EXFO-JSON input format is ever needed again.

## Other folders

- `raw_data/` — raw EXFO JSON scans, named `<date>_<condition>.json`; see `raw_data/README.md` for what each one is and the current findings.
- `results/<date>/` — generated CSVs and plots, one subfolder per day.
- `docs/` — explanatory and planning documents written for/with Carolina.
- `george/` — frozen reference copies of the files as originally handed over by the supervisor. Not used for running anything; the active copies above are what you actually execute.
- `logs/<date>.md` — daily work journal (see `CLAUDE.md` for the rule that maintains it).
