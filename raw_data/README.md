# Raw data — file index

All files: EXFO sweep 1505.06-1624.96 nm, 119,921 points, old 1-pm trigger mode,
2 channels (Ch1/Ch2). Format: see `george/HANDOVER.md` section 2.

| File | Date | Condition at the test port |
|---|---|---|
| `2026-08-18_open_end.json` | 2026-08-18 | No terminator, fiber end open (baseline, formerly `scan_s0.json`) |
| `2026-08-18_no_fiber.json` | 2026-08-18 | No test fiber connected at all (formerly `scan_s1new.json`) |
| `2026-08-18_terminator1.json` | 2026-08-18 | Terminator 1 |
| `2026-08-18_terminator2.json` | 2026-08-18 | Terminator 2 |
| `2026-08-18_10db_coupler.json` | 2026-08-18 | 10dB coupler inserted |

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
