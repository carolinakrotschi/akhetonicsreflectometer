# OFDR Reflectometer — Project Rules

## Directory structure

- `raw_data/` — raw EXFO/CoreDAQ JSON scans and their trimmed `.npz`
  derivatives. Raw `.json` scans are gitignored (some exceed GitHub's 100MB
  limit) and stay local-only; only `.npz` files and everything under
  `results/` are tracked. Keep the mapping in `raw_data/README.md` up to date.
- `results/<YYYY-MM-DD>/` — generated analysis output (CSV, PNG), sorted by date.
- `simulation/` — everything that generates or validates against synthetic
  data with known ground truth: `simulate_aux_data.py`, the validation
  harness (`validate_aux_pipeline.py`), and its output
  (`results/<date>/aux_validation/`, kept there under `results/`, not
  duplicated). Distinct from `tools/`: this is specifically the "does the
  Auswertung reconstruct a known-true signal correctly" side of the project.
  `process_reflectogram_aux.py` itself stays in the main directory (it
  processes both real bench data and simulated data, not simulation-only).
- `george/` — reference copies of the files handed over by the supervisor
  (state at handover time, frozen). The active working copies of the scripts
  live in the main directory and are what actually gets run — do NOT import
  from or invoke anything in `george/`.
- `HANDOVER.md` (root) — the living, actively-updated project reference
  (distance convention, channel pairing, known artifacts, open
  investigations). `george/HANDOVER.md` is the original frozen supervisor
  handover and is no longer updated; root `HANDOVER.md` is what's current.
- `md files/` — explanatory/planning documents written for/with Carolina
  (not the supervisor handover itself). Formerly `docs/`.
- `tools/` — standalone helper scripts not part of the core measurement
  pipeline. `helper plots/` holds the plotting-specific ones.
- `VIP/` — historical change notes (`CHANGES.md`), formerly under `legacy/`.
- `logs/<YYYY-MM-DD>.md` — daily journal (see below).

## Follow these rules in every session in this directory

- **`HANDOVER.md` (root) is binding.** Read it at the start of any task and
  honor its assumptions/conventions (distance convention, channel pairing
  Ch1/Ch3=aux vs Ch2/Ch4=measurement, mode choice diagnostic/cosmetic/none,
  known artifacts such as the 46.4mm reference peak and the fixed internal
  reflections around 587mm/518mm).

- **Daily log:** For every calendar day on which work happens in this
  directory, `logs/<YYYY-MM-DD>.md` should exist. If today's doesn't exist
  yet, create it on the first command you run. Log every notable command
  (script runs, analyses) with a short description, chronologically from the
  first to the last command of the day.

- **End of day:** When the user indicates she's done for the day (e.g. "that's
  it for today", "done", a session visibly wrapping up), briefly ask in the
  chat: "Quick recap — what got done today, what are the results, how do you
  interpret them?" Append her answer verbatim as a closing section to today's
  `logs/<YYYY-MM-DD>.md`.
