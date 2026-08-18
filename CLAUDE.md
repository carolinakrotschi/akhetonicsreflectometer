# OFDR Reflectometer — Project Rules

## Directory structure

- `raw_data/` — raw EXFO JSON scans, named `<YYYY-MM-DD>_<condition>.json`.
  Keep the mapping in `raw_data/README.md` up to date.
- `results/<YYYY-MM-DD>/` — generated analysis output (CSV, PNG), sorted by date.
- `george/` — reference copies of the files handed over by the supervisor
  (state at handover time, frozen). The active working copies of the scripts
  live in the main directory and are what actually gets run — do NOT import
  from or invoke anything in `george/`.
- `docs/` — explanatory/planning documents written for/with Carolina (not the
  supervisor handover itself).
- `tools/` — standalone helper scripts not part of the core measurement
  pipeline.
- `legacy/` — superseded earlier script versions, kept for reference.
- `logs/<YYYY-MM-DD>.md` — daily journal (see below).

## Follow these rules in every session in this directory

- **`george/HANDOVER.md` is binding.** Read it at the start of any task and
  honor its assumptions/conventions (distance convention, mode choice
  diagnostic/cosmetic/none, known artifacts such as the 46.4mm reference peak).

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
