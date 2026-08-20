# Reproducing today's real bench problems in simulation

**As of:** 2026-08-20. Uses the new `--internal-reflectors`, `--ghosts`,
`--aux-leak-db`, `--saturation-frac`, `--saturation-ceiling` flags added to
`simulation/simulate_aux_data.py` specifically to test the two biggest
mysteries from today's real scans (`logs/2026-08-20.md`): the persistent
587mm/518mm peaks and the detector-saturation-driven monotonicity failures.

## 1. The 587mm/518mm "ghost" peaks

**Mechanism modeled:** two independent, always-present fixed reflectors
inside the instrument (587mm at ~0 dB, 518mm at ~-20 dB — matching what six
real scans converged on), plus **multipath/double-bounce ghosting**: light
that partially reflects off tone *i*, continues, partially reflects off tone
*j*, and returns, travels the delay **sum** τ_i+τ_j — the standard mechanism
for "sum-position" ghosts in OFDR/OTDR. Amplitude ≈ product of the two field
amplitudes. The aux tone is included in the ghost pool too (at a tunable,
reduced "leak" coupling), since the real data showed aux-involving sum peaks
as well, implying imperfect isolation between the aux and measurement paths.

**Command:**
```bash
python simulation/simulate_aux_data.py --points 1000000 --speed 50 --lam0 1520 \
    --dl 4.24 --reflectors 1.062:-13 --internal-reflectors 0.587:0 0.518:-20 \
    --ghosts --aux-leak-db -8 --out sim_ghosts.npz
python process_reflectogram_aux.py sim_ghosts.npz --dl 4.24 \
    --aux-a 1 --aux-b 2 --meas-a 3 --meas-b 4 --zmax 4.5 --peak-floor-db -50
```

**Result — every peak found matches a fundamental or a pairwise sum, at the
predicted dB level:**

| Found (mm) | dB | Identity | Predicted |
|---|---|---|---|
| 518.0 | -19.8 | fundamental: internal reflector B | -20.0 |
| 587.0 | 0.0 | fundamental: internal reflector A | 0.0 |
| 1036.0 | -39.9 | 518+518 (2nd harmonic of B) | -40.0 |
| 1062.0 | -13.1 | fundamental: DUT fiber | -13.0 |
| 1105.0 | -20.2 | 587+518 (A+B) | -20.0 |
| 1174.0 | -0.6 | 587+587 (2nd harmonic of A) | 0.0 |
| 1580.0 | -33.0 | 518+1062 (B+fiber) | -33.0 |
| 1649.0 | -12.8 | 587+1062 (A+fiber) | -13.0 |
| 2124.0 | -26.4 | 1062+1062 (2nd harmonic of fiber) | -26.0 |
| 2638.0 | -28.4 | 518+2120 (B+aux) | -28.0 |
| 2707.0 | -8.2 | 587+2120 (A+aux) | -8.0 |
| 3182.0 | -21.4 | 1062+2120 (fiber+aux) | -21.0 |
| 4240.0 | -16.7 | 2120+2120 (2nd harmonic of aux) | -15.9 |

This is qualitatively and quantitatively the same pattern as the real data:
main+fiber sums, aux+main sums, 2nd harmonics, all present, all at
roughly the product-of-amplitudes dB level.

**Why this occurs:** a real, discrete partial reflector fixed inside the
instrument (not yet physically localized — candidates per the hardware:
a connector interface in the LO arm, the circulator, or the eVOA) sits in
series with the DUT and the return path. Because it's a *partial* reflector
(not total), a fraction of the light double-bounces between it and every
other reflector in the system (including the aux's own return path, if
isolation there isn't perfect), producing the sum-delay ghosts.

**How to remove it:**
- **Kill (or weaken) the fixed reflector itself.** Ghost amplitude scales
  with the *product* of the two reflectivities — since the 587mm reflector
  participates in nearly every ghost family, suppressing it (index-matching
  gel, a better-polished connector, an isolator at the suspect interface)
  collapses most of the ghost forest at once, not just one peak.
- **Localize it the same way the 6-32cm band was killed in the old setup**
  (`HANDOVER.md` §5's method): swap/bypass one component at a time between
  the circulator and the coupler, re-scan, and see which change makes the
  587mm peak vanish. The already-run no-fiber control test proved it's
  internal, not DUT-related — the next step is component-by-component
  isolation, not another full-fiber scan.
- **Improve aux/measurement isolation** if the aux-involving ghosts
  (2707mm, 2638mm, 3182mm here) turn out to matter in practice — these
  point at optical crosstalk between the two interferometers, not the
  587mm reflector itself.

## 2. Detector saturation -> monotonicity failure

**Mechanism modeled:** the first `--saturation-frac` of samples get each
channel hard-clipped to `--saturation-ceiling` of its own (unsaturated)
range — a plausible model of an ADC/detector ceiling being hit at the start
of the sweep, exactly as found 2026-08-19 (first ~10nm / ~20% of every real
scan).

**Command:**
```bash
python simulation/simulate_aux_data.py --points 1000000 --speed 50 --lam0 1520 \
    --dl 4.24 --reflectors 0.587:0 1.062:-13 --saturation-frac 0.2 \
    --saturation-ceiling 0.12 --out sim_saturation.npz
```

**Result, replicating Carolina's own diagnostic method (balanced Ch1/Ch2 ->
analytic signal -> unwrap -> backward-step fraction) directly on the
simulated aux channels:**

```
full scan backward-step fraction: 5.312%
  saturated (first 20%) region:   26.561%
  clean (remaining 80%) region:   0.000%
```

Backward phase steps appear **exclusively** in the saturated region and
**zero** elsewhere — the identical spatial signature found in the real data
("almost all backward steps fall in the first ~16% of the sweep... zero
afterward"). This confirms the mechanism: saturation-induced amplitude
clipping distorts the Hilbert-transform-derived instantaneous phase enough
to make it momentarily run backwards, exactly where the detector is
saturated and nowhere else.

Two side-effects also reproduced:
- `process_reflectogram_aux.py`'s own contrast diagnostic misreads this as
  **polarization fading** ("aux contrast drops somewhere") — a real,
  slightly misleading overlap worth knowing about: a low `amp_min` doesn't
  automatically mean polarization: check WHERE in the sweep it happens
  first.
- Despite the locally severe phase corruption, final reported positions
  stayed accurate in this run (`resample_on_aux()`'s
  `np.maximum.accumulate` clamp absorbs the backward excursions into the
  first few resampled points rather than smearing them across the whole
  axis) — meaning the corruption can be present and briefly invisible in
  the final PASSED/FAILED verdict while still being real and worth fixing
  before trusting a bow/ripple measurement from the same region.

**Why this occurs:** confirmed already in the real investigation as
detector/interferometer behavior at sweep start, not the laser (a direct
laser-only reference scan showed no equivalent ramp) — likely input-range
autoranging or settling at the very start of acquisition.

**How to remove it:**
- **Software (already applied, costs ~24% resolution):** trim the affected
  samples, exactly as done 2026-08-19 (`_trimmed1530.npz`). Confirmed here:
  the "clean" region alone has 0% backward steps.
  `simulate_aux_data.py --saturation-frac` lets you size how much margin a
  given trim fraction actually buys back, before touching the bench.
  - **Hardware (real fix, not yet tried):** per `HANDOVER.md`'s existing
  "fixed power range on every channel" rule — check whether the CoreDAQ's
  input range is set to accommodate whatever power spike happens at sweep
  start (or add a brief settling delay before recording), rather than
  discovering it after the fact via trimming every scan.
