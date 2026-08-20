# Raw data — file index

**Not tracked in git** (`.gitignore`'d): with the current 1M-datapoint sweeps,
raw JSON scans run ~150 MB each, past GitHub's 100 MB file limit. Only the
`_trimmed1530.npz` derivative of each scan (produced during processing) is
committed. This file stays local-only; the mapping below is the record of
what exists on disk.

All files: EXFO sweep 1505.06-1624.96 nm, 119,921 points, old 1-pm trigger mode,
2 channels (Ch1/Ch2). Format: see `george/HANDOVER.md` section 2.

Note: `2026-08-20_laserdirektangeschlossen.json` uses power keys in `[W]`
(`"Ch1 [W]"` etc.), not `[mW]` like every other file listed here -- the
instrument apparently switches units depending on power range. Convert before
comparing medians against other scans.

| File | Date | Condition at the test port |
|---|---|---|
| `2026-08-18_open_end.json` | 2026-08-18 | No terminator, fiber end open (baseline, formerly `scan_s0.json`) |
| `2026-08-18_no_fiber.json` | 2026-08-18 | No test fiber connected at all (formerly `scan_s1new.json`) |
| `2026-08-18_terminator1.json` | 2026-08-18 | Terminator 1 |
| `2026-08-18_terminator2.json` | 2026-08-18 | Terminator 2 |
| `2026-08-18_10db_coupler.json` | 2026-08-18 | 10dB coupler inserted |
| `2026-08-19_firstscanwith4channels.json` | 2026-08-19 | First test scan, new 4-channel aux-referenced free-run setup (600 pts, 1505-1625 nm, Ch1/Ch2 = aux MZI, Ch3/Ch4 = measurement). Processed with `process_reflectogram_aux.py`. |
| `2026-08-19_secondscanwith4channelsmorepoints.json` | 2026-08-19 | Same setup, 60,000 pts, 1505-1625 nm. |
| `2026-08-19_1miodatapoints.json` | 2026-08-19 | Same setup, 1,000,000 pts, 1520-1570 nm (scan_S5 of `md files/measurement_procedure.md`). First 200,000 pts (wl<1530nm) are detector-saturated; drop them (`_trimmed1530.npz`) before trusting anything -- with them included, acceptance test **fails monotonicity** (2.4-25% backward steps depending on channel pairing, limit 0%). Channel-pairing investigation (see `logs/2026-08-19.md`/`logs/2026-08-20.md`): aux = Ch2/Ch4 (dL≈4.24m), measurement = Ch1/Ch3 -- opposite of the 2026-08-18 convention above. |
| `2026-08-19_deltaLwirklich1m.json` | 2026-08-19 | Same setup/format. Measurement-interferometer LO arm rebuilt to a real, physical ~1 m round-trip match (test fiber ~1.04 m one-way, confirmed by tape measure). Same first-10nm saturation as above; use `_trimmed1530.npz`. Aux calibration matches the scan above almost exactly (dL≈4.24m) -- confirms only the measurement arm changed. Reflectogram (aux=Ch2/Ch4, meas=Ch1/Ch3) shows the real fiber-end peak at ~1061.7mm (matches the measured 1.04m length within 2%) plus an unexplained dominant peak at ~587.0mm (not fiber-related -- stable across re-coiling, see next row) and two nonlinear mixing artifacts (~1649mm sum, ~2121mm 2nd harmonic of the fiber peak). |
| `2026-08-19_deltaLwirklich1mneugerollt.json` | 2026-08-19 | Same as `deltaLwirklich1m.json`, fiber only re-coiled (same physical fiber, different bend routing) -- otherwise nothing changed. Confirms: 587.0mm peak unchanged in position/strength (supports it not being the fiber); ~1062mm fiber peak same position, 11dB weaker (consistent with added bend loss); ~1649mm sum-artifact got 13dB *stronger* (unexplained -- inconsistent with simple intermodulation, possible real multipath reflection). |
| `2026-08-20_3mfiber.json` | 2026-08-20 | ~3m fiber moved onto the measurement interferometer as its only fiber (old ~1m fiber + LO-match arm removed); aux dL backfilled with 3x1m segments instead. Main 587mm-family peak persists (4th confirmation, internal artifact); nothing at naive 3000mm; unexplained peak at 3734.5mm (later shown to be fiber-related but not matching the tape-measured 3.05m length directly -- see log). |
| `2026-08-20_nofiberattheend.json` | 2026-08-20 | No-fiber control: test-port fiber removed entirely, everything else (incl. 3x1m aux setup) left as in `3mfiber.json`. Proves the ~587mm and ~518mm peaks are fixed internal reflections, unrelated to any DUT fiber (both persist with nothing connected). |
| `2026-08-20_laserdirektangeschlossen.json` | 2026-08-20 | Laser output connected directly to a detector, bypassing the interferometer network entirely (control/reference scan, not a reflectogram). Only Ch2 carries signal (~5.800 mW, flat to ~0.03% across the whole 1520-1570nm sweep, **no saturation ramp at the start** unlike every interferometer scan); Ch1/Ch3/Ch4 sit at the electrical noise floor (~0.1-0.5 uW). Confirms the laser source itself is flat and clean -- the detector saturation seen in the first 200,000 points of other scans is caused downstream (interferometer/detector combination), not by the laser. Note: uses `[W]` units, see note above. Plot: `results/2026-08-20/2026-08-20_laserdirektangeschlossen_channels.png`. |
| `2026-08-20_nofiberattheendandnotthorlabsconnector.json` | 2026-08-20 | Second no-fiber control: same as `2026-08-20_nofiberattheend.json` (test-port fiber removed) but with the Thorlabs connector swapped for a different, "normal" (non-Thorlabs) connector -- **not removed, replaced** (corrected 2026-08-20, see log). Same 1M-point/4-channel/1520-1570nm format, same first-200,000-point saturation, trimmed the same way (`_trimmed1530.npz`). Big result: the dominant fixed-internal-reflection peak moved from ~587.8mm to **467.1mm** (0 dB) -- first evidence tying that peak's position to the physical path length of whichever connector occupies that slot. See `logs/2026-08-20.md` for the full writeup. |
| `2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleaned.json` | 2026-08-20 | Same setup as the row above (no fiber at test port, same replacement "normal" connector #1 still installed), plus the fiber was cleaned. Same format/trimming. Main peak essentially unchanged (467.0mm vs 467.1mm, negligible) -- confirms the fixed reflector is unrelated to fiber cleanliness, consistent with it not being fiber-related at all. Noise floor is ~9dB higher than the previous scan (median -59.2dB vs -68.4dB in a quiet 700-1600mm window), matching a lower aux contrast (min visibility 0.776 vs 0.846) -- explains why many more small peaks cross the fixed -50dB reporting threshold; not a new physical feature. See `logs/2026-08-20.md`. |
| `2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleanedandnewconnector.json` | 2026-08-20 | Same setup as the row above, but connector #1 was swapped for a SECOND, different "normal" connector (#2) -- so far: Thorlabs -> normal connector #1 -> normal connector #2, never actually connector-free. Same 1M-point/4-channel/1520-1570nm format, but **saturation ramp lasts ~13nm this time** (vs ~10nm on every earlier 2026-08-20 scan) -- trimmed to `_trimmed1533.npz` (dropped first 260,000 pts, wl<1533nm) instead of the usual `_trimmed1530.npz`. Big result: main peak moved back to **589.7mm**, into the same ~587-590mm family seen throughout the day with the ORIGINAL Thorlabs connector in place -- confirms the peak's position tracks the physical path length of whichever connector occupies that slot (connector #2's length happens to be close to the Thorlabs one's; connector #1's was ~120mm shorter), not the identity of one specific connector. The ~518-521mm secondary cluster also reappeared. See `logs/2026-08-20.md` for the full writeup and `HANDOVER.md` section 2 for the updated standing conclusion. |
| `2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleanedandnewconnectorneueeinschwingzeit.json` | 2026-08-20 | Same nominal setup as the row above (connector #2), but recorded AFTER a fix to the akhelabs acquisition software (`exfo_window.py`): a new "Settle margin (nm)" field makes the laser start its physical sweep that many nm before the requested start wavelength, then crops the pre-roll back out after acquisition -- so the returned data no longer contains the ~10-13nm detector-saturation transient at all. 833,333 points, 1520.00-1570.00nm, no wavelength trim needed (saved whole as `_full.npz`, verified clean in 0.5nm bins from the very first sample). Resolution improved from 22.61um (previous scan, 13nm had to be trimmed) to **16.59um** (full ~50nm span used) -- see `logs/2026-08-20.md` for the comparison plot. Note: aux dL (3.593m) and the main peak position (495.6mm) both differ noticeably from the immediately preceding connector-#2 scan (4.276m / 589.7mm) -- the physical setup was evidently touched again while testing the new firmware feature; not yet reconciled, see log. |

**Finding (6-32cm band, see `results/2026-08-18/band_6_32cm_comparison.png`):**
Terminator 1, Terminator 2, and the 10dB coupler all sit at approx. -52 to
-54 dB — essentially identical to `open_end` (-52.6 dB). Only `no_fiber` is
noticeably lower (-68.8 dB, noise floor). None of the tested terminators
effectively suppresses the reflection; the suspicion falls on Rayleigh
backscatter along the whole patchcord, or the near-end connector
(circulator <-> patchcord), rather than the far end. Next test discussed:
a short patchcord (<35 cm) to distinguish between the two hypotheses.

An earlier terminator attempt (`scan_s1.json`, deleted before 2026-08-18) had
an even louder band (-35.8 dB) — cause unresolved, file no longer available.

**Test fiber length (2026-08-18):** confirmed 50 cm one-way (1 m round trip).
`python diagnose_artifacts.py fold 0.5` shows this is BEYOND the old setup's
Nyquist range (0.417 m) and folds to **0.334 m (33.4 cm)** as an undersampled
band (1.67 points/fringe) — right at the edge of the observed 6-32cm dirt
band. This means the far-end reflection itself likely IS a real contributor
to that band (aliased, not resolvable), on top of/instead of the Rayleigh
backscatter hypothesis above. With the current 1-pm trigger mode, no
terminator or connector fix can turn this into a clean peak — the fiber is
simply longer than this mode can resolve unambiguously. The Aux-referenced
setup (~7 m Nyquist range) is required to see it cleanly at its true 50 cm
position.

**Open discrepancy:** confirmed 2026-08-18 to be the same physical patchcord
as the one behind the ~1.046 m estimate in `HANDOVER.md` §5 (2026-08-12) --
but the two lengths don't reconcile by a simple convention mixup (`fold 0.5`
and `fold 1.046` predict different apparent band positions). See the
"Test fiber length" note in `HANDOVER.md`'s hardware section for the
candidate explanations and the suggested `alias` pair test to settle it.
