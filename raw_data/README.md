# Raw data — file index

**Scans from 2026-09-11 onwards ARE in the repository, as `.npz`.** The raw
`.json` still cannot be pushed -- each one is ~155 MB, past GitHub's hard
100 MB per-file limit -- but the same data stored as a compressed `.npz`
is 11-13 MB, a factor of ~13 smaller, and the conversion is **lossless**:
the four channel arrays are kept as float64 and verified bit-identical
after the round trip, and the wavelength axis is kept in full because it
turned out not to be linear (it deviates up to 0.147 nm from a straight
line). Converted with `tools/json_to_npz.py`; `process_reflectogram_aux.py`
reads `.npz` directly, and a reflectogram computed from the `.npz` is
bit-identical to one computed from the `.json`.

**Older raw JSON: not tracked in git** (`.gitignore`'d): with the current 1M-datapoint sweeps,
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
| `2026-08-20_3mfiber.json` | 2026-08-20 | **"3mA" fiber** (label introduced 2026-08-31 to distinguish from the physically different "3mB" fiber used on 2026-08-31 -- inferred lengths differ, ~3147mm vs ~2616mm relative to connector). ~3m fiber moved onto the measurement interferometer as its only fiber (old ~1m fiber + LO-match arm removed); aux dL backfilled with 3x1m segments instead. Main 587mm-family peak persists (4th confirmation, internal artifact); nothing at naive 3000mm; unexplained peak at 3734.5mm (later shown to be fiber-related but not matching the tape-measured 3.05m length directly -- see log). |
| `2026-08-31-ohnefiberallesnormalerstemessungdestages.json` | 2026-08-31 | No fiber, first measurement of the day. Monotonicity FAILS (0.0224%, localized contrast dip). See `logs/2026-08-31.md`. |
| `2026-08-31-ohnefiberallesnormalzweitemessungdestages.json` | 2026-08-31 | No fiber, second measurement of the day (repeat). Monotonicity PASSES. |
| `2026-08-31-ohnefiberallesnormaldrittemessungdestages.json` | 2026-08-31 | No fiber, third measurement of the day (repeat). Monotonicity FAILS badly (0.93%, real sweep anomaly, not just noise) -- dL/tau_aux from this scan unreliable, do not reuse. |
| `2026-08-31-ohnefiberandereaux1mfiber.json` | 2026-08-31 | No fiber; the aux's 1m fiber segment swapped for a different physical one (dL essentially unchanged, ~3.59m; contrast/monotonicity improved). |
| `2026-08-31-ohnefiberaux5m.json` | 2026-08-31 | No fiber; aux stretched to dL≈4.489m (labeled "aux5m" target). Best contrast of the no-fiber scans that config. |
| `2026-08-31-ohnefiberaux6m.json` | 2026-08-31 | No fiber; aux stretched further to dL≈5.39m (labeled "aux6m" target). Points/fringe dropping toward the ~4 warning floor. |
| `2026-08-31-1mfiberaux5m.json` | 2026-08-31 | ~1m fiber (no separate letter label), aux5m config (dL≈4.489m). Inferred length ≈912mm relative to connector. |
| `2026-08-31-3mfiberaux5m.json` | 2026-08-31 | **"3mB" fiber** (see note on `2026-08-20_3mfiber.json` above), aux5m config (dL≈4.489m). Inferred length ≈2616.5mm relative to connector. |
| `2026-08-31-3mfiberaux4m.json` | 2026-08-31 | **"3mB" fiber** on the measurement arm, aux4m config (dL≈3.578m; **aux arm itself built from "3mA" + "1mB" segments** -- "3mA" being the old 2026-08-20 fiber repurposed as aux filler, "1mB" a separate 1m fiber, both distinct from the measurement-arm "3mB"/"1mA" naming below -- confirmed by Carolina 2026-08-31 after an initial mix-up). Same physical measurement fiber as the aux5m scan above, inferred length ≈2615.5mm, confirms length is aux-config-independent. |
| `2026-08-31-3mplus1mfiberaux4m.json` | 2026-08-31 | **"3mB" fiber + "1mA" fiber in series** on the measurement arm (the 1m component here labeled "1mA" to distinguish it from the standalone `1mfiberaux5m` fiber, which has no separate letter); aux4m config, aux arm built from "3mA"+"1mB" (see note above -- these are DIFFERENT physical fibers from the measurement arm's "3mB"/"1mA", despite the confusingly overlapping letters). Inferred combined length ≈3523.8mm, vs. 2615.5+912.0=3527.5mm predicted from the two measurement fibers measured separately -- within 4mm (this match is what confirmed the measurement arm is "3mB"+"1mA", not "3mA"+"1mB", after Carolina's initial statement suggested the latter). |
| `2026-08-31-3mplus1mfiberaux4m_2.json` | 2026-08-31 | Repeat of the scan above (same physical setup: measurement arm "3mB"+"1mA", aux arm "3mA"+"1mB"). dL=3.5778 m (tau_aux=17.5193 ns) -- matches the first scan's 3.5782 m to 0.4 mm, confirms the aux4m build is stable between scans. |
| `2026-08-20_nofiberattheend.json` | 2026-08-20 | No-fiber control: test-port fiber removed entirely, everything else (incl. 3x1m aux setup) left as in `3mfiber.json`. Proves the ~587mm and ~518mm peaks are fixed internal reflections, unrelated to any DUT fiber (both persist with nothing connected). |
| `2026-08-20_laserdirektangeschlossen.json` | 2026-08-20 | Laser output connected directly to a detector, bypassing the interferometer network entirely (control/reference scan, not a reflectogram). Only Ch2 carries signal (~5.800 mW, flat to ~0.03% across the whole 1520-1570nm sweep, **no saturation ramp at the start** unlike every interferometer scan); Ch1/Ch3/Ch4 sit at the electrical noise floor (~0.1-0.5 uW). Confirms the laser source itself is flat and clean -- the detector saturation seen in the first 200,000 points of other scans is caused downstream (interferometer/detector combination), not by the laser. Note: uses `[W]` units, see note above. Plot: `results/2026-08-20/2026-08-20_laserdirektangeschlossen_channels.png`. |
| `2026-08-20_nofiberattheendandnotthorlabsconnector.json` | 2026-08-20 | Second no-fiber control: same as `2026-08-20_nofiberattheend.json` (test-port fiber removed) but with the Thorlabs connector swapped for a different, "normal" (non-Thorlabs) connector -- **not removed, replaced** (corrected 2026-08-20, see log). Same 1M-point/4-channel/1520-1570nm format, same first-200,000-point saturation, trimmed the same way (`_trimmed1530.npz`). Big result: the dominant fixed-internal-reflection peak moved from ~587.8mm to **467.1mm** (0 dB) -- first evidence tying that peak's position to the physical path length of whichever connector occupies that slot. See `logs/2026-08-20.md` for the full writeup. |
| `2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleaned.json` | 2026-08-20 | Same setup as the row above (no fiber at test port, same replacement "normal" connector #1 still installed), plus the fiber was cleaned. Same format/trimming. Main peak essentially unchanged (467.0mm vs 467.1mm, negligible) -- confirms the fixed reflector is unrelated to fiber cleanliness, consistent with it not being fiber-related at all. Noise floor is ~9dB higher than the previous scan (median -59.2dB vs -68.4dB in a quiet 700-1600mm window), matching a lower aux contrast (min visibility 0.776 vs 0.846) -- explains why many more small peaks cross the fixed -50dB reporting threshold; not a new physical feature. See `logs/2026-08-20.md`. |
| `2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleanedandnewconnector.json` | 2026-08-20 | Same setup as the row above, but connector #1 was swapped for a SECOND, different "normal" connector (#2) -- so far: Thorlabs -> normal connector #1 -> normal connector #2, never actually connector-free. Same 1M-point/4-channel/1520-1570nm format, but **saturation ramp lasts ~13nm this time** (vs ~10nm on every earlier 2026-08-20 scan) -- trimmed to `_trimmed1533.npz` (dropped first 260,000 pts, wl<1533nm) instead of the usual `_trimmed1530.npz`. Big result: main peak moved back to **589.7mm**, into the same ~587-590mm family seen throughout the day with the ORIGINAL Thorlabs connector in place -- confirms the peak's position tracks the physical path length of whichever connector occupies that slot (connector #2's length happens to be close to the Thorlabs one's; connector #1's was ~120mm shorter), not the identity of one specific connector. The ~518-521mm secondary cluster also reappeared. See `logs/2026-08-20.md` for the full writeup and `HANDOVER.md` section 2 for the updated standing conclusion. |
| `2026-08-20_nofiberattheendandnotthorlabsconnectorandfibercleanedandnewconnectorneueeinschwingzeit.json` | 2026-08-20 | Same nominal setup as the row above (connector #2), but recorded AFTER a fix to the akhelabs acquisition software (`exfo_window.py`): a new "Settle margin (nm)" field makes the laser start its physical sweep that many nm before the requested start wavelength, then crops the pre-roll back out after acquisition -- so the returned data no longer contains the ~10-13nm detector-saturation transient at all. 833,333 points, 1520.00-1570.00nm, no wavelength trim needed (saved whole as `_full.npz`, verified clean in 0.5nm bins from the very first sample). Resolution improved from 22.61um (previous scan, 13nm had to be trimmed) to **16.59um** (full ~50nm span used) -- see `logs/2026-08-20.md` for the comparison plot. Note: aux dL (3.593m) and the main peak position (495.6mm) both differ noticeably from the immediately preceding connector-#2 scan (4.276m / 589.7mm) -- the physical setup was evidently touched again while testing the new firmware feature; not yet reconciled, see log. |

| `2026-09-11-11-36_3mfiber.json` | 2026-09-11 | ~3 m fiber on the measurement arm, LinaWindow/coreDAQ-LIN, 1520-1570 nm. **Not yet processed** -- recorded 40 min before the chip scan below, same session. |
| `2026-09-11-12-16hhi1fiber7.json` | 2026-09-11 | **Sixth HHI chip measurement -- fiber 7 = `b25`, the STUB negative control. Passed.** Facet A at 1673.3662 mm (0 dB, -3 dB width 33.2 um vs 43.1 um window limit = isolated single reflector), pigtail 1099.492 mm, noise floor -69.2 dB, and **nothing behind facet A** beyond the satellite skirt (max -40 dB in +5..+12 mm, -44 dB in +20..+60 mm). Three conclusions: (1) the rule **b = 32 - fiber** is confirmed a second time, independently of any path-length assignment; (2) **there is no "SSC-internal reflector" at +3.1..3.7 mm** -- a stub contains a full SSC and shows nothing there, so that earlier interpretation is retracted; (3) **"strongest peak = facet A" does not hold generally** -- facet A is the fiber/chip interface reflection, set by glue/gap/index match alone, and varies by >45 dB across channels (+36 dB vs the internal reference here, -9 dB in the top loopback), so facet A must be identified via the pigtail-length family (1094.4..1099.5 mm over 6 fibers), not via peak height. tau_aux = 20.3811 ns, monotonicity passes. |
| `2026-09-11-12-16hhi1fiber30.json` | 2026-09-11 | **Fifth HHI chip measurement -- fiber 30.** Cleanest scan of the day (noise floor -71 dB, single reflector with -3 dB width 49.8 um vs 43.1 um window limit). Facet A at 1670.2383 mm (pigtail 1096.312 mm), one dominant reflector at **+9.2069 mm = 2685.8 um of on-chip waveguide -- which matches NO channel in the layout** (nearest path end is b3 at 2489.6 um, 28 resolution cells away). **This refutes the "fiber 27 was a misread 30" hypothesis**: the two scans have different pigtails (1095.32 vs 1096.31 mm) and different dominant reflectors (+14.6148 vs +9.2069 mm), so they are different couplers. **Consequence: the n_SSC/n_waveguide split reported from the 3-point fit is RETRACTED** -- it relied on assigning the +14.6148 mm peak to b2's MMI input, which no longer holds; that 2.3 um agreement was most likely coincidence. What stands: n_g,chip = 3.4803 +- 0.008 from the two loopbacks alone. tau_aux = 20.3808 ns, monotonicity passes. |
| `2026-09-11-12-16hhi1fiber6.json` | 2026-09-11 | **Fourth HHI chip measurement -- fiber 6.** Confirms the fiber-to-coupler rule **b = 32 - fiber** (position from top = fiber - 3, 29 couplers on fibers 4..32): without the +3 offset fiber 6 would be `b23`, a stub with 18 um of waveguide, which cannot produce the observed -1.31 dB peak 15 mm into the chip plus a 70 mm forest of 448 peaks. Quantitatively consistent with `b26` (run 5126.8 um, implied pigtail 1095.84 mm) but **not separable from b2** -- they differ by 0.378 mm while the multipath forest (mean peak spacing 0.156 mm) limits facet-A anchoring to +-0.4 mm. Methodological lesson: **circuit channels are unusable for precision work** (decoder has many MMIs + 2 ring resonators -> recirculation); only the two loopback pairs give clean geometry. tau_aux = 20.3815 ns, monotonicity passes. |
| `2026-09-11-12-16hhi1fiber27.json` | 2026-09-11 | **Third HHI chip measurement -- single-ended channel, fiber 27.** Named "kanal24" but the measurement identifies it as **GDS coupler `b2`**, not b24: facet A at 1669.3593 mm (fiber 1095.316 mm), dominant peak (0 dB, the strongest in the whole scan) at +14.6148 mm = the input facet of an MMI 1x2 splitter at the end of b2's 4967.278 um waveguide run -- predicted 14.6214 mm, i.e. **6.6 um off over 14.6 mm (0.045%)**. tau_aux = 20.3814 ns, monotonicity passes, best contrast of the day (0.640). Combined with the two loopback scans this gives 3 measurements for 2 unknowns: **n_g,SSC = 3.4808 +- 0.0121 and n_g,waveguide = 3.4782 +- 0.0078, agreeing to 0.07%** -- so SSC and waveguide have the same group index and the uniform-n_g assumption holds. Fiber-to-channel mapping still unconfirmed; the only simple rule matching fiber 27 -> b2 is fiber n -> b(29-n). |
| `2026-09-11-12-16hhi1fiber32.json` | 2026-09-11 | **Second HHI chip measurement -- the BOTTOM loopback pair `b0`/`b1`** (design loop 556.5494 um, chip path 2956.55 um), same setup as the scan below, only re-plugged to the other coupler pair. 988,536 pts, 1520-1570 nm. tau_aux = 20.3807 ns (stable to 0.01% vs the top-loop scan), monotonicity passes. Facet A 1668.3752 mm (-16.6 dB), facet B 1675.3923 mm (-4.1 dB) -> **chip path dz = 7.0171 mm -> n_g,chip = 3.4842**, which agrees with the top pair's 3.4765 to **0.22%**. Differential vs the top pair: dz shrinks by 0.1825 mm for an 83.577 um shorter design loop = **77.0 um of on-chip waveguide measured SSC-free** (one resolution cell = 7.0 um physical). The middle reflectors stay put in both scans -> SSC-internal, not the loop ends. Pigtails here: 1094.4005 / 1095.1139 mm (per-channel, glued). Couples better than the top pair. |
| `2026-09-11-12-16hhi1fiber4.json` | 2026-09-11 | **First HHI chip measurement -- the TOP loopback pair `b27`/`b28`.** Measurement arm = MZI -> ~1 m fiber (SQS Vlaknova Optika, FC/APC, SM, no. 0103254561) -> HHI chip 1, **both** edge couplers of the top-left pair are pigtailed with ~1 m fibers (measurement arm -> fiber 1 -> chip loop -> fiber 2 -> open end, which provides the back-reflection); on-chip link is the 640.13 um loop (see `logs/2026-09-11.md` for the GDS length extraction). 988,536 pts, 1520-1570 nm. tau_aux = 20.3809 ns (dL = 4.1622 m), monotonicity passes (0.0000%). Chip input facet at 1670.23 mm (-9.1 dB) = +1096.30 mm past the fixed internal 573.93 mm reference; chip end facet at +7.200 mm, two intermediate reflectors at +3.135/+4.064 mm (2x3.135+0.929 = 7.199, symmetric 3-section path = SSC + loop + SSC). Full path vs. design 3040.13 um gives n_g,chip = 3.477. Open end of fiber 2 found at 2773.625 mm (-36.8 dB) = +1096.199 mm past the chip, i.e. **the two ~1 m fibers match to 0.099 mm** (1096.298 vs 1096.199 mm) -- independent confirmation of the whole chain. The +-0.70/+-1.53 mm satellites are instrument artifacts (ring the 573.93 peak identically). |

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

## Ligentec-SiN-Chip (Faserarray an der unteren Kante, 95 Kanaele)

**Chip-Praefix-Konvention (ab 2026-09-16):** alle Scans bis einschliesslich
2026-09-15 stammen vom selben physischen Chip, **MAP2672** -- ihre
rawdata- und results-Dateien wurden deshalb nachtraeglich mit
`2672_ligentechhi_` praefigiert. Ab 2026-09-16 wird ein neuer Chip,
**MAP2680**, vermessen; neue Scans dieses Chips sollten entsprechend mit
`2680_ligentechhi_` beginnen, damit rawdata/results-Dateinamen allein
erkennen lassen, von welchem physischen Chip sie stammen. Ausgenommen von
der 2672-Umbenennung: `raw_data/hhi_ligentec/` (gemeinsame GDS-Designdatei,
nicht chipspezifisch) und die `..._2m07cmfiber...`-Scans (laut Tabelle
unten explizit **kein Chip** angeschlossen).

**MAP2680 traegt elektrische Kontakte (Bondpads/Leiterbahnen), MAP2672
nicht.** Auch unangeschlossen koennen sie die Reflektogramm-Struktur nahe
der Facette beeinflussen (zusaetzliche Streuzentren/Reflexe durch
Metallisierung). Erster Hinweis 2026-09-16: Faser 1 auf beiden Chips zeigt
in der auf den eigenen Peak zentrierten Ansicht rechts vom Peak aehnliche,
links vom Peak abweichende Struktur -- siehe `logs/2026-09-16.md`. Noch
nicht verifiziert (z.B. gegen GDS-Kontaktlayout oder einen Kontrollscan);
offene Hypothese, beim Interpretieren chipuebergreifender Strukturvergleiche
mitdenken.

Alle Scans: neuer PM-Aufbau, freilaufend, 4 Kanaele (Ch2/Ch4 = aux,
Ch1/Ch3 = Messung), **1520.00-1570.00 nm, 988,536 Punkte**, tau_aux
~20.633 ns (aux dL 4.21 m). Aufloesung 16.6 um faseraequivalent,
Nyquistbereich 8.0 m. Kanalnummer = 96 - Fasernummer.

Gemeinsamer interner Anker fuer alle diese Scans: **529.43 mm**; auf ihn
werden die Pegel normiert, weil die alten Anker (574/587 mm) nach dem
PM-Umbau weg sind.

| File | Datum | Faser -> Kanal | Was dahinter liegt | Facette A |
|---|---|---|---|---|
| `2672_ligentechhi_2026-09-14-16-19_ligentecnofiber.json` | 2026-09-14 | — | **Kontrolle**, nichts angeschlossen | — |
| `2672_ligentechhi_2026-09-14-16-02_Ligentecfiber1.json` | 2026-09-14 | 1 -> 95 | Loopback 199.494 um zu Kanal 94 | 1582.818 mm |
| `2672_ligentechhi_2026-09-14-16-36_Ligentecfiber44.json` | 2026-09-14 | 44 -> 52 | Schaltung, 1x2-MMI nach 8912 um | 1580.416 mm |
| `2672_ligentechhi_2026-09-14-16-41_Ligentecfiber45.json` | 2026-09-14 | 45 -> 51 | Schaltung, wie Kanal 52 an comp265 | — |
| `2672_ligentechhi_2026-09-14-16-50_Ligentecfiber45neu.json` | 2026-09-14 | 45 -> 51 | Wiederholung von 16:41 | — |
| `2672_ligentechhi_2026-09-15-08-19_ligentecfiber45wavelenght1520to1570nm.json` | 2026-09-15 | 45 -> 51 | Wiederholung, Standardspanne | — |
| `2672_ligentechhi_2026-09-15-08-19_ligentecfiber45wavelenght1505to1625nm.json` | 2026-09-15 | 45 -> 51 | **Ausnahme: 1505-1625 nm**, volle Spanne | — |
| `2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48.json` | 2026-09-15 | 48 -> 48 | Schaltung; der Kanal, der in **Transmission** Signal gibt (Laser 48 -> PM 49) | 1582.141 mm |
| `2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48fnotfullwavelenghtremeasure.json` | 2026-09-15 | 48 -> 48 | Wiederholung, Standardspanne 1520-1570 nm | Facette 1585.712 mm |
| `2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48fullwavelenght.json` | 2026-09-15 | 48 -> 48 | **Ausnahme: 1505-1625 nm**, volle Spanne | Facette 1585.669 mm |
| `2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48fullwavelenghtremeasure.json` | 2026-09-15 | 48 -> 48 | **defekt** -- JSON bricht bei Zeichen 356 457 004 mitten in der Zahlenliste ab (Aufnahme unterbrochen) | -- |
| `2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48fullwavelenght_ipa.json` | 2026-09-15 | 48 -> 48 | volle Spanne, **mit IPA** (Bedeutung noch nicht dokumentiert -- vermutlich Indexanpassung an der Facette) | Facette 1585.670 mm |
| `2026-09-15-08-19_ligentec2m07cmfiberwavelenght1520to1570nm.json` | 2026-09-15 | — | **kein Chip**, nur eine 2.07-m-Faser am Port | Faserende 2642.975 mm |
| `2026-09-15-15-14meisusmfiberarrzchannel4.json` | 2026-09-15 | Meisu-Array, Kanal 4 | **anderes Faserarray**, kein Ligentec-Chip | Arrayfacette 1632.68 mm |
| `2026-09-15-15-14meisusmfiberarrzchannel4remeasure.json` | 2026-09-15 | Meisu-Array, Kanal 4 | Wiederholung von Kanal 4, neu gesteckt | Arrayfacette 1632.52 mm |
| `2026-09-15-15-14meisusmfiberarrzchannel8.json` | 2026-09-15 | Meisu-Array, Kanal 8 | dito | Arrayfacette 1633.31 mm |
| `2026-09-15-15-14meisusmfiberarrzchannel9.json` | 2026-09-15 | Meisu-Array, Kanal 9 | dito | Arrayfacette 1633.42 mm |
| `2026-09-15-15-14meisusmfiberarrzchannel12.json` | 2026-09-15 | Meisu-Array, Kanal 12 | dito | Arrayfacette 1632.44 mm |
| `2026-09-15-15-14meisusmfiberarrzchannel12remeasure.json` | 2026-09-15 | Meisu-Array, Kanal 12 | Wiederholung von Kanal 12 | Arrayfacette 1632.56 mm |
| `2026-09-15-15-14phixsmfiberarrzchannel1.json` | 2026-09-15 | PHIX-Array, Kanal 1 | **drittes Faserarray**, kein Ligentec-Chip | Arrayfacette 1624.55 mm |
| `2026-09-15-15-14phixsmfiberarrzchannel8.json` | 2026-09-15 | PHIX-Array, Kanal 8 | dito | Arrayfacette 1626.19 mm |
| `2026-09-16-10-36_ligentec2680_fiber1.json` | 2026-09-16 | **MAP2680**, 1 -> 95 | erster Scan neuer Chip, Standardspanne | Facette 1581.697 mm |
| `2026-09-16-10-36_ligentec2680_fiber1full.json` | 2026-09-16 | **MAP2680**, 1 -> 95 | volle Spanne 1505-1625 nm | Facette 1581.696 mm |
| `2680_ligentechhi_20260916-105948_fiber44_1520-1570nm.json` | 2026-09-16 | **MAP2680**, 44 -> 52 | Standardspanne | Facette 1579.112 mm |
| `2680_ligentechhi_20260916-110010_fiber44_1505-1625nm.json` | 2026-09-16 | **MAP2680**, 44 -> 52 | volle Spanne | Facette 1579.129 mm |
| `2680_ligentechhi_20260916-110152_fiber48_1520-1570nm.json` | 2026-09-16 | **MAP2680**, 48 -> 48 | Standardspanne | Facette 1579.172 mm |
| `2680_ligentechhi_20260916-110214_fiber48_1505-1625nm.json` | 2026-09-16 | **MAP2680**, 48 -> 48 | volle Spanne | Facette 1579.127 mm |

`.npz` vorhanden bisher nur fuer `2672_ligentechhi_2026-09-15-09-34_Ligentecfiber48`; die
2026-09-14er Ligentec-Scans und die acht Array-Scans (Meisu 4/8/9/12 plus
Wdh. von 4 und 12, PHIX 1/8) sind noch nicht konvertiert und existieren
nur lokal als `.json`.

Befund ueber alle Kanaele: ausser der Facettenreflexion selbst kommt
**nichts aus dem Chip zurueck** -- siehe `logs/2026-09-14.md` und
`logs/2026-09-15.md`. Pruefwerkzeuge:
`tools/check_ligentec_channel.py` (Vorhersagefenster + zwei Kontrollen)
und `tools/compare_ligentec_channels.py` (Zentrierung auf Facette A +
Nullverteilung).
