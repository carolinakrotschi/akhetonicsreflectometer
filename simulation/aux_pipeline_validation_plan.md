# Aux-pipeline validation — scenario matrix and findings

**For:** Carolina · **As of:** 2026-08-20 · Companion to `md files/followup_questions.md`

Harness: `simulation/validate_aux_pipeline.py`. Runs `simulation/simulate_aux_data.py`
(`simulate(a)`) → `process_reflectogram_aux.py` (`process(a)`, stays at the
repo root — it also processes real bench data) across a matrix of scenarios
and writes `simulation/results/aux_validation_<date>/summary.csv` +
`summary.md`. Each scenario gets its own subfolder with `sim.npz` and the
usual `_reflectogram.csv/.png`, so any flagged case can be re-run by hand
(note the explicit channel args — the simulator writes Ch1/Ch2=aux,
Ch3/Ch4=meas, not the real-hardware interleaved default):

```bash
python process_reflectogram_aux.py simulation/results/aux_validation_2026-08-20/<name>/sim.npz \
    --dl 4 --aux-a 1 --aux-b 2 --meas-a 3 --meas-b 4
```

## Why this exists

`process_reflectogram_aux.py` already had a truth-comparison built in, but
only for one ad-hoc scenario per manual run, judged only by a single
position threshold. This matrix systematically exercises the specific risks
already flagged in `HANDOVER.md` — resolution limit, weak-signal
floor, Nyquist aliasing, sweep-nonlinearity robustness, τ_aux calibration
error, and the currently uncorrected chromatic dispersion — instead of
relying on one optimistic demo case.

## Two bugs found and fixed along the way

1. **Silent aliasing gap** (`process_reflectogram_aux.py`): a ground-truth
   reflector beyond `z_nyq` used to be skipped with `"-- outside range"`
   *without* affecting the PASSED/FAILED verdict — a reflector the pipeline
   completely lost could still print PASSED. Fixed: the true position is now
   folded via `z_apparent = fold(z_true, period=2*z_nyq)` and checked there;
   an unaccounted-for reflector now correctly fails.
2. **Narrow-window blind spot in the harness itself** (not in
   `process_reflectogram_aux.py`): the ±20-cell search window around
   `z_true` that both the existing script and my first harness draft used is
   fine for sub-resolution errors, but a τ_aux miscalibration can move a
   peak by *hundreds* of cells — far outside that window. The first run
   reported nonsense (8-19 "cells" of error, roughly window-sized, regardless
   of the actual miscalibration). Fixed by predicting the miscalibration's
   effect analytically (`z_expect = z_true * (1 + eps)`, since a wrong τ_aux
   rescales the entire distance axis, not just one peak) and searching around
   *that*. See `grade_calibration()` in the harness.

## Scenario matrix and results (this run)

Operating point: 1,000,000 pts, 60 nm/s, `--dl 4` → dz_bin ≈ 14.2 µm,
z_nyq ≈ 6.95 m, τ_aux ≈ 19.59 ns.

| Scenario | Category | Result | Interpretation |
|---|---|---|---|
| `baseline_measured_conditions` | position | **FAILED** (errors 12-18 cells, 170-250 µm) | Ripple set to the *measured* 90-110 MHz (used 100 MHz) instead of the tool's optimistic 20 MHz default. `simulate_aux_data.py` itself reports **~21% of steps running backwards** at this ripple — the aux phase is no longer monotonic, which the tool's own docstring already says "no software can rescue." This is the single most important finding: the earlier clean demo in `docs/followup_questions.md` #9 used the tool's 20 MHz default, not the measured bench value — **measure the real ripple before trusting position accuracy.** |
| `baseline_tool_default` | position | PASSED (errors < 0.4 cells) | Regression anchor — reproduces `docs/followup_questions.md` #9 exactly. Confirms the refactor didn't change behavior, but is optimistic relative to the bench. |
| `doublet_2cell` (28 µm apart) | resolution | **resolved** (8.7/5.4 dB valley) | Two reflectors ~2 resolution cells apart are cleanly separated. |
| `doublet_1cell` (14 µm apart) | resolution | **not resolved** (correctly) | One cell apart merges into a single peak, as physically expected — not a bug, it's the resolution floor. |
| `weak_near_floor` (down to -55 dB) | amplitude | found, errors < 0.3 cells | Positional recovery of weak reflectors is accurate even near/under the -45 dB *display* floor. Note: the truth-comparison itself is not gated by `--peak-floor-db` — a weak reflector found here would still not appear in the normal "Peaks above X dB" operator listing without lowering that flag. |
| `weak_near_floor_elevated_noise` (5x noise) | amplitude | found, errors < 0.4 cells | Recovery held up even under 5x detector noise; floor wasn't reached in this run. |
| `nyquist_edge_inside` (6.80 m) | range | PASSED | Clean behavior just inside the Nyquist range. |
| `nyquist_beyond_aliased` (8.00 m, 9.00 m) | range | **aliased_as_expected** | Both reflectors correctly reappear at the predicted folded position (`2*z_nyq - z_true`), confirming the HANDOVER.md §5 aliasing story reproduces exactly in the aux-referenced pipeline too — thanks to fix #1 above, this is no longer silently missed. |
| `ripple_stress_near_limit` (110 MHz) | robustness | 23% backward steps, errors up to 13 cells | Consistent with the baseline finding — ripple above ~50 MHz (this operating point's own monotonicity limit) reliably breaks position accuracy. |
| `bow_stress_2x` (118 pm) | robustness | 0% backward steps, errors < 1.6 cells | Slow bow has a lot more margin than fast ripple — 2x the measured bow barely moves the needle. |
| `tau_aux_miscal_0.5/1/2pct` | calibration | **FAILED at every level**, scales linearly | At 0.5% τ_aux error: 16.5 cells (235 µm) at 0.046 m, up to **701 cells (1 cm)** at 2.0 m. At 2%: up to **2765 cells (~4 cm)**. This confirms the resampling itself is correct to ~0.2-0.35 cells once τ_aux is known (see "vs.predicted" column) — the entire error is the calibration offset. **Practical consequence: τ_aux must be held to well under 0.1% to keep 2-cell accuracy at 2 m** — a much tighter bench requirement than "accurate to ~1%" (`docs/followup_questions.md` #9) suggested for *positions*; that 1% figure is fine for a first picture, not for a number in a report. |
| `dispersion_beta2_off/on` (0.30 m, 1.50 m) | known limitation | position unchanged, amplitude −5.4→−9.1 dB at 1.50 m | Matches the physics: a symmetric (quadratic) spectral phase broadens the peak (lower height) without shifting its centroid. HANDOVER.md §6.5 calls dispersion correction deferred; this quantifies today's actual cost as **peak-height loss, not position error** — a reflector already near the noise floor could be pushed under threshold by dispersion at ~1.5 m+, even though its *reported position* stays fine. |

## Bottom line

The pipeline's core mechanism (aux-phase resampling, aliasing fold, weak-
signal recovery) checks out cleanly whenever its two stated preconditions
hold: (1) the fast tuning ripple stays under roughly the monotonicity limit
printed by `simulate_aux_data.py` at the chosen operating point (~50 MHz
here, well below the 90-110 MHz quoted as *measured* in `HANDOVER.md`
§2/§5), and (2) τ_aux is calibrated to well under 0.1%, not the ~1% that was
previously considered good enough for "first pictures." Both are bench
tasks, not software fixes — this is exactly the numbers-before-narratives
outcome the project's working agreements (`HANDOVER.md` §7) call for:
measure the real ripple and tighten the τ_aux calibration procedure before
trusting reported positions from the aux-referenced pipeline.
