# hhi_1 — HHI chip 1: all chip measurements (from 2026-09-11 onwards)

Everything belonging to the HHI chip-1 campaign in one place: the chip design
(GDS) and the raw-data-side files of every chip scan taken on 2026-09-11 and
later. The processed output (reflectogram/chain PNG+CSV, identified.csv) stays
where it is, under `results/2026-09-11/` and `results/2026-09-14/`.

Fiber-to-coupler rule (confirmed twice, see `logs/2026-09-11.md`): **b = 32 - fiber**.

## Chip design

| File | What |
|---|---|
| `HHI_RUN_1.gds` | GDS layout of HHI chip 1 (copy of `results/2026-09-11/HHI_RUN_1.gds`). Source of all design lengths; read it with `tools/gds_waveguide_length.py` / `tools/gds_netlist_reflectors.py`. |

## Scans

| Scan | Date | Fiber -> coupler | What |
|---|---|---|---|
| `2026-09-11-12-16hhi1fiber4` | 2026-09-11 | 4 -> b28 (TOP loopback b27/b28) | First chip measurement. n_g,chip = 3.477. |
| `2026-09-11-12-16hhi1fiber32` | 2026-09-11 | 32 -> b0 (BOTTOM loopback b0/b1) | Second. n_g,chip = 3.4842, agrees with the top pair to 0.22%. |
| `2026-09-11-12-16hhi1fiber27` | 2026-09-11 | 27 -> b5 (labelled "kanal24") | Third, single-ended. Note: the 2026-09-11 write-up assigns this scan to `b2` from a path-length match, which the rule above contradicts (rule says b5) -- unresolved, see note below. |
| `2026-09-11-12-16hhi1fiber6` | 2026-09-11 | 6 -> b26 | Fourth. Circuit channel, multipath forest -- unusable for precision work. |
| `2026-09-11-12-16hhi1fiber30` | 2026-09-11 | 30 -> b2 | Fifth. Cleanest scan of the day (noise floor -71 dB). |
| `2026-09-11-12-16hhi1fiber7` | 2026-09-11 | 7 -> b25 | Sixth. STUB negative control -- passed, nothing behind facet A. |
| `2026-09-11-11-36_3mfiber` | 2026-09-11 | — | ~3 m fiber, 40 min before the chip scans, same session. Not a chip scan; kept with the day for context. |
| `2026-09-14-08-56_nichtsconnected` | 2026-09-14 | — | Nothing connected (baseline/control for the day). |
| `2026-09-14-09-03_fiber4` | 2026-09-14 | 4 -> b28 | Repeat of the top loopback. |
| `2026-09-14-09-35_fiber32` | 2026-09-14 | 32 -> b0 | Repeat of the bottom loopback. |
| `2026-09-14-09-58_fiber23` | 2026-09-14 | 23 -> b9 | |

The coupler column applies `b = 32 - fiber` throughout. For fibers 27 and 30
that disagrees with the original per-scan path-length assignments made before
the rule was confirmed (fiber 27 -> b2 there); the n_SSC/n_waveguide split that
rested on those assignments is already retracted. Treat the two single-ended
2026-09-11 scans' coupler identity as open.

Per-scan details (numbers, conclusions, retractions) are in the top-level
`raw_data/README.md` and in `logs/2026-09-11.md` / `logs/2026-09-14.md`.

## What is actually in this folder

This folder holds the chip design and, once copied over, the raw `.json` scans
of the campaign -- **nothing else**. No PNGs: the plots live with the processed
output under `results/<date>/` and on the raw-data side in `raw_data/`.

**The raw `.json` scans are NOT here yet.** They are gitignored (~150 MB each,
past GitHub's 100 MB limit, see `.gitignore`) and therefore local-only on the
measurement PC -- they never travelled through git to this machine. To complete
this folder, copy the eleven `.json` files listed above from the measurement
PC's `raw_data/` into `raw_data/hhi_1/`. The gitignore rule was widened to
`raw_data/**/*.json`, so they stay untracked once they land here.
