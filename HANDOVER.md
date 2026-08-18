# OFDR Reflectometer — Project Handover

**Date:** 2026-08-12 · **For:** incoming student + any Claude session taking over
**Scripts:** `process_reflectogram.py` (main pipeline), `diagnose_artifacts.py`
(discriminator toolkit), `tools/simulate_scan_data.py` (synthetic validator, CSV format)

Read this top to bottom once. The single most important lesson of the project
so far is at the end of §5: three plausible explanations for one artifact were
wrong before the fourth was proven, and each was eliminated by a *measurement*,
not a better story. Work that way.

---

## 1. What is being built

An optical frequency-domain reflectometer (OFDR), functionally a home-built
Luna OBR: sweep a tunable laser, interfere the light returned from a fiber
under test with a local oscillator copy, FFT the interferogram, and read out
every reflection along the fiber as a peak — position, magnitude, micron-scale
resolution.

**Core physics in one equation.** With the laser at optical frequency ν, each
reflector at round-trip delay τᵢ contributes to detected power:

    I(ν) = Σᵢ 2·√(P_LO·Pᵢ)·cos(2π·τᵢ·ν + φᵢ)

The interferogram *as a function of ν* is a sum of sinusoids; each reflector
is one tone whose frequency is its delay. FFT over ν → reflectogram.

**The numbers that govern everything:**

| Quantity | Formula | Current value |
|---|---|---|
| Resolution | Δz = c/(2·n_g·Δν_span) | 7.0 µm (120 nm span) |
| Nyquist range | z_max = c/(4·n_g·δν_step) | 0.42 m (1 pm steps) |
| Sample spacing in ν | δν = c·δλ/λ² | 132→114 MHz across sweep (18% chirp!) |
| Fringe period of a path at z | λ²/(2·n_g·z) | 0.78 pm at z = 1.05 m |
| n_g (SMF-28, 1550 nm) | — | 1.468 |

Range·×·resolution is fixed by the point count: range/resolution = N/2 cells.
Sweep *speed* cancels out of everything except the required sample rate.

**Distance-axis convention (source of a real confusion — read twice):** all
scripts use the reflection convention z = c·τ/(2·n_g). A *transmission* MZI
with arm imbalance ΔL appears at z = ΔL/2. The dominant peak at **46.4 mm**
in all current data is the interferometer's own beat (τ = 454.5 ps): either
a 92.8 mm fiber imbalance (transmission reading) or a reflector 46.4 mm past
the LO-match point (reflection reading). **A tape measure on the bench
resolves which — this has NOT been verified yet.** Mid-sweep visibility is
0.745, meaning the two beating fields are within ~7 dB of each other; check
that this is consistent with whatever the topology turns out to be.

---

## 2. Hardware, as of handover

- **Laser:** EXFO tunable, swept 1505 → 1625 nm at 10 nm/s. Emits a trigger
  every 1 pm; after the sweep it can be queried for a **measured** wavelength
  table (one entry per trigger). The table is quantized to 1 pm — it carries
  real information (using it beats the nominal grid) but not sub-step detail.
- **Detection:** CoreDAQ multichannel power meter (dev unit). Records one
  sample per external trigger, ≤50 kHz, ~130k-point trigger-mode buffer, 1M
  free-run buffer, 1 µs minimum averaging. Ch1/Ch2 = the two complementary
  outputs of the recombining coupler.
- **Topology:** all-SM for now (PM planned later), all connectors FC/APC.
  A circulator gives the reflective path; the fiber under test is a **~1 m
  patchcord** with a transparent cap on the far end (cap currently ON).
  An eVOA exists (originally for arm balancing; destined for LO trim).
- **Not yet built:** the dedicated auxiliary interferometer (see §6),
  balanced photodetectors, polarization-diverse receiver, any absolute-λ
  reference. All deferred deliberately; none blocks current work.

**Data format** (what all scripts expect):

```json
{"header": {"instrument": "EXFO", "device_name": "exfo-1"},
 "data": [{"Device Description": {"info": "..."},
           "Wavelength [nm]": [...],
           "Ch1 [mW]": [...],
           "Ch2 [mW]": [...]}]}
```

Known quirks, all handled by the loaders: the **last sample is a retrace
point** (laser parks below the start — drop anything breaking monotonicity);
~75 triggers are missing at the sweep start and ~33 at the end (arming during
laser acceleration — harmless, but **never align two scans by index; align by
wavelength value**); the wavelength table contains ~4.7k duplicate values
(1 pm quantization).

---

## 3. The three datasets analysed

| Scan | File timestamp | Condition |
|---|---|---|
| A | 2026-08-12-17-47 | cap on, patchcord straight |
| B | 2026-08-12-18-13 | cap REMOVED, straight |
| C | 2026-08-12-18-21 | cap back on, patchcord coiled at r ≈ 2–3 cm |

All: 119,921 good points, 1505.07–1624.97 nm, visibility 0.70–0.88
(declining toward long λ — unexplained, likely polarization walk; benign so
far), channels complementary at −0.94 AC correlation.

---

## 4. Processing pipeline (ofdr_process.py) and why each step exists

1. **Load + cleanup** — retrace drop by value, not position.
2. **Balanced subtraction in software** — P = Ch1 − g(λ)·Ch2 with g fitted
   per segment (it drifts 0.84→1.14 across the span; the two photodiodes
   have different responsivity slopes). Kills the DC pedestal and common-mode
   intensity noise; the sum Ch1+Ch2 is a free power monitor.
3. **λ → ν on the nominal grid** — samples uniform in λ are NOT uniform in ν
   (18% spacing chirp). A raw FFT of the unresampled buffer smears every
   peak; this step is mandatory, not cosmetic.
4. **Slow-bow correction** — the laser's sweep deviates from its grid by
   **~59 pm p-p** (a slow bow, measured every scan; `ofdr_diagnose.py
   tuning`). Uncorrected, the main peak fragments over ±0.8 mm. The pipeline
   extracts the dominant tone's Hilbert phase, smooths it heavily (so only
   the slow bow enters), and warps the ν axis. Result: main peak 21.0 µm
   wide vs an 18.2 µm window limit — transform-limited to 15%.
   **Peak width is the standing health metric**; the script warns if it
   regresses past 2× the limit.
5. **Window + FFT** — Kaiser β = 12 default (sidelobes < −80 dB at ~2.6 bins
   width). Hann is sharper but its −31 dB sidelobes will be mistaken for
   reflectors next to strong peaks.

**The `--mode` flag encodes a hard-won distinction:**
- `diagnostic` (default): correction tracks only the slow bow → real
  features stay sharp *at true positions*. Use for identifying anything.
- `cosmetic`: correction tracks the raw phase → sharpest main peak, but it
  ABSORBS weak additive features (measured: −8 dB on real artifacts) and
  shifts their apparent positions. Never diagnose from this mode.
- `none`: shows the uncorrected smearing, useful for demonstrations.

**Known harmonics:** balanced subtraction cancels even harmonics of the
fringe and passes odd ones; the sum channel does the opposite. Peaks at
exactly 2×, 3×… the main delay (e.g. 92.8, 139.2 mm) are detector/processing
distortion, NOT reflectors. The peak lister tags them.

---

## 5. The artifact investigation — what the "noisy peaks" were

**Symptom:** dense structure inside each fringe; in the delay domain, a
smeared band at 6–32 cm, −38 to −53 dB, envelope stable between scans but
fine structure completely random (envelope corr +0.979, fine corr +0.014
between A and B).

**Hypotheses tried and KILLED, with the test that killed each:**

1. *Laser FM noise* — killed by PM/AM discriminator: additive paths modulate
   the main tone's amplitude and phase equally; phase noise gives PM≫AM.
   Measured **PM/AM = 1.04** → additive. (`ofdr_diagnose.py pm-am`)
2. *Discrete connector multipath* — killed by the smear + wavelength-drift of
   the band centroid (fixed paths hold delay to ~0.1%; this walked 26%).
3. *Polarization multipath in the circulator* — killed by scan C: coiling
   (attenuation) reduced the band; SOP scrambling alone (scan B handling)
   had not reduced the envelope. Loss-sensitive + SOP-insensitive ≠
   polarization.
4. *Fabry–Pérot in the end cap* (user hypothesis) — killed by scan B: cap
   removed, band unchanged (median −52.9 → −51.0 dB). Also no comb signature.

**The verdict, proven by simulation:** the band is the **far end of the ~1 m
patchcord, aliased**. It sits at ~1.05 m true, beyond the 0.42 m Nyquist
range; it folds to |1.05 − 0.83| ≈ 21 cm. Its true fringe period is 0.78 pm —
undersampled at 1 pm — so it cannot resample coherently, and the laser's
tuning noise shreds it into a speckled band whose fine structure is a new
random realization every sweep (hence not calibratable). A synthetic tone at
main + 1.00 m pushed through the identical pipeline with the measured tuning
noise landed at **centroid 21.4 cm, 90% of energy in 19.2–23.5 cm** — on top
of the measured band. The CoreDAQ's boxcar integration suppresses it −14.8 dB,
explaining the weak apparent level. `ofdr_diagnose.py fold 1.046` reproduces
the arithmetic.

Every observation fits: circulator-only (a transmission MZI gives the far
end no path back) ✓, PM/AM = 1 ✓, envelope-stable/fine-random ✓ (noise
*spectrum* fixed, *realization* random), cap-indifferent ✓ (bare and capped
ends reflect similarly), coil-diminished ✓ (bend loss on the round trip,
−5 dB band, ripple 30–38% → 17–20%), centroid walk ✓ (fold position depends
on δν, which chirps 18% across the sweep).

**Consequences:**
- **Frequency-axis correction cannot un-alias.** Aliasing destroys
  information at sampling; there is no software fix downstream, and no
  optical anti-alias filter exists. Only defenses: no reflectors beyond
  Nyquist (terminate, shorten, index-match), or finer sampling.
- The band is an MPI **noise floor** (~−40 dB rel. main across 6–32 cm) as
  long as the far end stays where it is. Reflectors weaker than that, in
  that zone, are currently unmeasurable.
- The cap-removal scan also disproved a previous recommendation: this
  artifact **cannot** be subtracted via a reference scan (fine structure is
  per-scan random). `ofdr_diagnose.py compare` renders this verdict
  automatically for any future mystery band.

**Secondary known artifacts:**
- mm-scale common-mode etalons at 1.5 and 3.2 mm (positively correlated
  between channels → common path; likely detector air gap / eVOA / laser
  internals). Cause near-symmetric ±1.0/±1.3 mm sidebands at −21 dB around
  every strong peak. Benign now; on the list.
- Fast tuning ripple: ~90–110 MHz rms, strongest period ~4.2–4.7 pm. Within
  Nyquist it is corrected; its main damage is smearing *aliased* content.

---

## 6. Where the design is heading (decided, with rationale)

**Immediate next step — kill the far-end reflection**, any of: index-matching
gel on the (APC) end face, a tight mandrel wrap at the last few cm, or a
<35 cm test patchcord. Prediction: the 6–32 cm band drops into the noise.
This is the first thing the student should do and verify with
`ofdr_diagnose.py compare` (before/after).

**Architecture decision made in-session:** abandon trigger-mode acquisition
for **free-running internal-clock sampling** — the Luna architecture.
At 10 nm/s with 12 µs averaging: 0.12 pm effective spacing, 600k–1M points,
**3.4 m Nyquist at full 120 nm span, 6.9 µm resolution**. The 1 m patchcord
end then appears as a sharp honest peak instead of an aliased smear.
Constraints found: don't sweep much below ~10 nm/s (fringes descend into the
acoustic band; scene must also stay still for the whole sweep), and the
range/resolution product is fixed at N/2 cells regardless of speed.

**Non-negotiable consequence:** without triggers there is NO wavelength
information — the frequency axis must come from a co-recorded **auxiliary
interferometer**. The current 92.8 mm-imbalance interferometer IS a working
aux (17 pts/fringe in trigger mode; 31+ in free-run). Calibrate τ_aux once
by running trigger mode on the aux alone: total unwrapped phase / (2π·span
from the table) → τ_aux = 454.5 ps, already measured. When processing with a
real aux channel: resample measurement channels on the aux's **raw**
unwrapped phase — the diagnostic/cosmetic tension of §4 disappears because
the correction no longer comes from the measurement channel itself.
Do NOT smooth the aux phase (it must track the ~4.7 pm ripple).

**Deferred upgrades, in order, each bolts on without rework:**
1. Real balanced photodetectors + fast digitizer — needed only for Rayleigh
   backscatter (software subtraction of Ch1/Ch2 mathematically cannot reach
   it: RIN eats the ADC range ~700× below the crossover at −73 dB sig/LO).
2. Polarization-diverse receiver (PBS + two balanced pairs, LO split
   equally) — removes fading and amplitude errors. Note for the all-fiber
   build: split the SM return FIRST on the PBS, give each branch its own PM
   LO copy; sum the two reflectograms in MAGNITUDE after the FFT, never the
   fields before.
3. PM fiber for everything the instrument owns (aux, LO path). The DUT arm
   stays SM by nature; never splice SM→PM in the return (PANDA DGD ≈ 1.5
   ps/m splits every peak into doublets 0.3 mm apart per metre).
4. Wavelength anchor: skip the gas cell; laser's absolute readout is ~100×
   better than distance accuracy needs (dn_g/dλ ≈ 5×10⁻⁶/nm → 3.5 ppm/nm).
   An etalon for τ_aux calibration if traceability is ever needed.
5. Dispersion autofocus (quadratic spectral phase, β ∝ distance) — a few-
   cell effect at 40 cm range, dominant beyond a metre. Single global β
   corrects one distance at a time; segment for quantitative return loss.

---

## 7. Working agreements for the next session

- **Numbers before narratives.** Every hypothesis about an unexpected
  feature gets a discriminating measurement (`ofdr_diagnose.py` implements
  four). The session record: 3 wrong theories (2 from Claude, 1 from the
  user) eliminated by tests before the right one survived simulation.
- **The peak-width check is the standing arbiter** of frequency-axis health:
  46.4 mm peak at ≤ ~21 µm with Kaiser-12, or stop and investigate.
- **Whiteboard numbers:** Nyquist range for the current sampling mode, and
  "no reflector beyond it, ever, unless terminated."
- Fixed power range on every CoreDAQ channel (autoranging mid-sweep = phase
  discontinuity = unwrap corruption to end of scan).
- Validate any pipeline change on synthetic data with known ground truth
  before real data (`tools/simulate_scan_data.py`; note it emits the older
  two-CSV format consumed by the legacy `legacy/reflectometer_csv_generic.py`,
  kept for reference).
- Suggested first tasks for the student, in order: (1) tape-measure the
  46.4 mm topology question (§1); (2) terminate the far end and verify the
  band dies; (3) one free-run acquisition end-to-end with aux-referenced
  processing; (4) reproduce the §5 investigation on their own scans as a
  learning exercise — every test is one command.
