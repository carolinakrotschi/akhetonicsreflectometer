# OFDR Reflectometer — Handover Explained

**For:** Carolina · **Date:** 2026-08-18
**Refers to:** `HANDOVER.md`, `process_reflectogram.py`, `diagnose_artifacts.py`, and
your supervisor's email about the new MZI

This document explains (A) the handover — every formula, every term, from
first principles — and (B) your supervisor's proposal, translated and worked
through, with the numbers you'll need at the table tomorrow.

---

# PART A — Understanding the Handover

## A.1 What the device actually does

An **OFDR** (Optical Frequency-Domain Reflectometer) answers the question:
*"Where along this fiber is there a reflection, and how strong is it?"*
So essentially a map of the fiber: a connector here, a splice there,
the fiber end-face at the end.

A classic OTDR does this with a short pulse and a stopwatch
(time domain). Resolution ~cm to m. An OFDR does it differently and much
more finely: it sends **no** pulse, but a laser that continuously sweeps
its wavelength (*sweep*), and measures **interference**. Resolution:
micrometers.

The setup in one sentence: light is split. One part goes through a
**circulator** into the test fiber and comes back reflected. The other part
travels alongside untouched as a reference (**local oscillator, LO**). Both
are recombined and interfere on a detector.

## A.2 The one fundamental equation — derived

The handover states:

```
I(ν) = Σᵢ 2·√(P_LO · Pᵢ) · cos(2π · τᵢ · ν + φᵢ)
```

Where does this come from?

**Step 1 — Two fields superimpose.** At the detector, the LO field E_LO and
a reflected field Eᵢ arrive together. A detector measures power, i.e. the
squared magnitude:

```
P = |E_LO + Eᵢ|² = |E_LO|² + |Eᵢ|² + 2·Re(E_LO* · Eᵢ)
                   └──DC──┘  └small┘  └── the interesting part ──┘
```

The third term is the **interference term**. It is proportional to
`2·√(P_LO·Pᵢ)` — this is why an OFDR is so sensitive: an extremely weak
return signal Pᵢ is *amplified* by the square root and by the strong LO.
This is called **coherent detection**, or heterodyne gain.

**Step 2 — Where the cosine comes from.** The reflected light returns later
than the LO, by the round-trip delay **τᵢ** (Greek tau, *round-trip delay*).
A delay of τ corresponds, at optical frequency ν, to a phase difference of
`2π·τ·ν`. Hence `cos(2π·τᵢ·ν + φᵢ)`.

**Step 3 — The trick.** Now the laser sweeps ν. Because τᵢ is **fixed**
while ν changes, the phase runs on continuously — the detector sees a
**beat** (*fringes*, interference fringes). And here's the crux:

> If you plot the detector signal against the **optical frequency ν**, it is
> a sum of sinusoids. The "frequency" of each sinusoid (measured in *per Hz
> of optical frequency*, i.e. in seconds) **is** exactly the round-trip delay
> τᵢ of the corresponding reflector.

So: **an FFT over ν → every reflector becomes a peak, and its position on
the FFT axis is its round-trip delay, i.e. its distance.** That's the entire
principle. The rest of the scripts is care taken around this core idea.

A picture to go with it: a nearby reflector produces a slow beat (few
fringes over the whole sweep), a distant one a fast beat. Same as radar
FMCW, or OCT in medicine — the same math.

## A.3 The table of formulas, line by line

### Resolution: `Δz = c / (2 · n_g · Δν_span)` → **7.0 µm**

*Δν_span* = the total optical frequency span swept through. From 1505 nm to
1625 nm that is **14.7 THz**.

**Why this formula?** Basic Fourier rule: if you observe a signal over a
window of width B, you can distinguish frequencies down to 1/B. Here the
window is the frequency span Δν_span, so distinguishable round-trip delays
are `Δτ = 1/Δν_span`. Converted to distance: `Δz = c·Δτ/(2·n_g)`.

**Rule of thumb:** *Resolution depends ONLY on the width of the sweep.*
Sweeping more nm = finer resolution. Nothing else helps.

### Nyquist range: `z_max = c / (4 · n_g · δν_step)` → **0.42 m**

*δν_step* = the spacing **between two measurement points** in optical
frequency (at a 1 pm step: approx. 122 MHz).

**Why?** Nyquist sampling theorem: an oscillation needs at least 2 points
per period, otherwise it gets reconstructed incorrectly. A reflector at τ
produces an oscillation with period `1/τ` in ν. Condition: `1/τ >
2·δν_step`, so `τ_max = 1/(2·δν_step)` and in distance `z_max =
c/(4·n_g·δν_step)`.

**Rule of thumb:** *Range depends ONLY on how closely spaced your
measurement points are.* Smaller step = larger range.

### The most important consequence of both

```
Range / Resolution = N/2       (N = number of measurement points)
```

Currently: 120,000 points → 60,000 "cells" → 0.42 m ÷ 7 µm ≈ 60,000. ✔

**That's the entire budget.** You can trade range against resolution, but
the product is fixed by the number of points. If you want *both* better:
more points. Full stop. And **sweep speed cancels out everywhere** — it
only determines whether your digitizer is fast enough, or whether your
memory is sufficient.

### Point spacing in ν: `δν = c·δλ/λ²` → **132 → 114 MHz (18% chirp!)**

This is just the wavelength ↔ frequency conversion: `ν = c/λ`,
differentiated `|dν| = c·dλ/λ²`.

**The catch:** The laser triggers every **1 pm in λ** — but because λ² is
in the denominator, 1 pm at 1505 nm is a frequency step of 132 MHz, while at
1625 nm it's only 114 MHz. So your points are evenly spaced in λ, but
**unevenly spaced in ν** — an 18% difference over the sweep. The handover
calls this *chirp*.

An FFT assumes evenly spaced samples. If you simply FFT the raw data, every
peak smears. That's why step 3 of the pipeline (λ→ν conversion with
resampling) is **mandatory, not cosmetic**.

### Fringe period: `λ²/(2·n_g·z)` → **0.78 pm at z = 1.05 m**

The inverse of the same calculation: how many picometers the laser must
sweep for a reflector at z to complete one full interference period.

At z = 1.05 m the period is **0.78 pm** — but you only measure every 1 pm.
**Fewer than 2 points per period = undersampled.** This is exactly where the
big problem in §5 of the handover comes from (see A.6).

### Group index `n_g = 1.468`

Light in glass fiber is slower than in vacuum: `v = c/n_g`. For round-trip
delays, what counts is the **group** index (speed of the signal/wave
packet), not the phase index (~1.4682 vs. ~1.468). For SMF-28 at 1550 nm:
1.468.

In practice: **1 m of fiber ≈ 4.9 ns one-way, i.e. 9.8 ns round trip.**

## A.4 The distance convention — source of confusion

All scripts compute `z = c·τ / (2·n_g)`. The factor of 2 means: *"τ is a
round trip"* — the reflection convention. z is the one-way distance to the
reflector.

**But:** A **transmission MZI** (light passes through only once, no
reflection) with arm-length difference ΔL produces τ = n_g·ΔL/c — only
*one* pass. The script still divides by 2 and therefore displays `z =
ΔL/2`.

The dominant peak at **46.4 mm** is therefore either
- an MZI arm-length difference of **92.8 mm** (transmission interpretation), **or**
- a genuine reflector **46.4 mm** behind the LO balance point.

The handover explicitly says: **this is not yet resolved — a tape measure
on the table will settle it.** That is your task (1). Do not skip it — the
entire interpretation of the topology hinges on it.

## A.5 The pipeline (`process_reflectogram.py`) step by step

**1. Loading + cleanup.** The last measurement point is the *retrace*: the
laser jumps back below the start wavelength after the sweep. It is
discarded by **value** (not by position), because in addition ~75 triggers
are missing at the start and ~33 at the end (laser ramp-up phase).

> **Rule from the handover: never compare two scans by index, always by
> wavelength value.** Otherwise you're comparing offset data points.

**2. Balanced subtraction in software.** `P = Ch1 − g(λ)·Ch2`.

The final coupler has **two complementary outputs**: where one is bright,
the other is dark (energy conservation). Both outputs carry the same
interference signal, but in antiphase — and the same DC pedestal and the
same intensity noise **in phase**. So:

- Difference Ch1 − Ch2 → interference **doubled**, DC & common-mode noise **removed**
- Sum Ch1 + Ch2 → interference gone, pure power monitoring for free

`g` is a correction factor because the two photodiodes have different
sensitivities, drifting over 120 nm (0.84 → 1.14). It's fitted piecewise by
median.

*(This is the software version. Real balanced photodetectors do this in
hardware and are much better — see handover §6, item 1.)*

**3. Convert λ → ν and resample.** Because of the 18% chirp. Mandatory.

**4. "Slow-bow" correction.** The laser doesn't sweep exactly as
advertised. The actual deviation from the nominal 1 pm grid is **~59 pm
peak-to-peak** in a slow bow over the sweep (see the right half of
`george/real_data_diagnosis.png` — the bow from −25 pm through +18 pm back to
−40 pm).

Without correction the main peak smears out over ±0.8 mm.

How the correction is done: take the dominant tone, extract its
instantaneous phase via the **Hilbert transform** (in the code:
`analytic_signal`), smooth it heavily, and use it to warp the ν axis.

> **Hilbert / analytic signal in one sentence:** a real cosine is converted
> into a rotating phasor whose angle can be read directly as phase — you do
> this simply by deleting the negative frequencies in the spectrum. That is
> exactly what the three lines of `analytic_signal()` do.

Result: main peak **21.0 µm** wide against a theoretical limit of 18.2 µm →
**15% above ideal, very good.**

> **Peak width is the system's health indicator.** 46.4 mm peak ≤ ~21 µm
> with Kaiser-12 → everything fine. Larger → investigate first, then keep
> measuring. The script warns above 2×.

**The `--mode` switch (important, not a cosmetic detail):**

| Mode | What it does | When to use |
|---|---|---|
| `diagnostic` (default) | Correction follows only the *slow* bow | **Always, when you want to identify something.** Real features stay sharp and at the correct location |
| `cosmetic` | Correction follows the raw phase | For display only. **Eats weak features (−8 dB measured) and shifts them.** Never diagnose from this |
| `none` | No correction | Demo of how bad it would be without correction |

You can see this directly in the left half of `george/real_data_diagnosis.png`:
blue = raw EXFO grid (peak drowning in grass at −10 dB), orange =
phase-corrected (clean peak, floor at −45 dB).

**5. Windowing + FFT.** A finite-length dataset produces **sidelobes** in
the FFT — false, small peaks next to every real one. A window (soft
weighting toward the edges) suppresses these.

- **Hann**: narrow main peak, but sidelobes only −31 dB → in an OFDR these
  get mistaken for reflectors. Dangerous.
- **Kaiser β = 12** (default here): main peak 2.6 bins wide, sidelobes
  **< −80 dB**. The right choice when hunting for weak real reflections.

**Known harmonics:** peaks at exactly 2×, 3× the main delay (92.8 mm,
139.2 mm) are **not reflectors**, but distortion products of the
detector/processing. Balanced subtraction removes even-order and lets
odd-order through; the sum does the opposite. The peak lister flags these
automatically.

## A.6 The artifact story (§5) — and why it's the most important page

**Symptom:** a diffuse band of peaks between 6 and 32 cm, −38 to −53 dB.
Between two scans, the **envelope** stayed identical (correlation +0.979),
but the **fine structure** was completely re-randomized each time (+0.014).

**Four explanations were tested, three died:**

1. *Laser FM noise* → killed by the **PM/AM test**. Logic: a real extra
   light path modulates the amplitude and phase of the main tone equally
   strongly; pure phase noise modulates almost only the phase. Measured:
   **PM/AM = 1.04** → it's real light on a real path (**additive**).
2. *Multiple reflection at connectors* → killed, because fixed paths keep
   their round-trip delay to within 0.1%; this band moved by 26%.
3. *Polarization effect in the circulator* → killed by scan C: coiling
   (= attenuation) weakened the band, pure polarization twisting did not.
   Loss-sensitive + polarization-insensitive = not polarization.
4. *Fabry-Pérot in the end cap* → killed by scan B: cap removed, band
   unchanged.

**The truth: it's the fiber end at ~1.05 m — aliased.**

**Aliasing, illustrated:** Like the wagon wheel in a film that appears to
spin backward. The camera films too slowly for the actual rotation, so the
eye invents a false, slower one. Same thing here: the fiber end at 1.05 m
produces a fringe period of **0.78 pm**, but you only measure every **1
pm** — undersampled. So it appears at a **wrong** position:

```
Nyquist range 0.42 m; fold period 2 × 0.42 = 0.83 m
1.05 m − 0.83 m = 0.22 m  →  appears at ~21 cm ✔ (measured: 21 cm)
```

`python diagnose_artifacts.py fold 1.046` computes this for you.

And because it's undersampled, the resampling cannot handle it coherently —
the laser's tuning noise tears it into a noise band that gets
**re-randomized on every sweep**.

**The three consequences you must remember:**

1. **Aliasing cannot be fixed in software.** The information is lost
   during sampling, not afterward. There is no optical anti-alias filter.
   Only two ways out: (a) don't allow any reflectors beyond Nyquist
   (terminate, shorten, index-matching), or (b) **sample more finely** —
   and that's exactly what your supervisor now wants to build.
2. As long as the fiber end sits there, a **noise floor of ~−40 dB** lies
   across 6–32 cm. Anything weaker there is unmeasurable.
3. **Cannot be calibrated away** with a reference scan — the fine
   structure is different every time.

**And the working method you should adopt:** three plausible stories were
wrong. Each was killed by a **measurement**, not by a better argument.
`diagnose_artifacts.py` is exactly these four tests as commands:

```bash
python diagnose_artifacts.py pm-am  scan.json          # additive or phase noise?
python diagnose_artifacts.py compare scanA.json scanB.json  # fixed or random?
python diagnose_artifacts.py fold   1.046              # what folds where?
python diagnose_artifacts.py tuning scan.json          # how well is the laser sweeping today?
```

## A.7 Glossary

| Term | Meaning |
|---|---|
| **MZI** | Mach-Zehnder interferometer: split, two arms of different length, recombine. Produces a clean sinusoidal fringe as a function of ν |
| **LO** (Local Oscillator) | The reference beam that does not pass through the test fiber and interferes with the return signal |
| **Circulator** | 3-port component; light goes 1→2→3 but never backward. Separates outgoing from returning light |
| **eVOA** | Electronically variable optical attenuator. Used here to adjust LO power |
| **Balanced detection** | Measuring and subtracting both complementary coupler outputs |
| **Visibility** (0.70–0.88) | Fringe contrast, `(max−min)/(max+min)`. 1 = perfectly matched field strengths. 0.745 means the two fields differ by ~7 dB |
| **SOP** | State of Polarization. Changes in normal SM fiber with every touch |
| **SM / PM fiber** | Single-Mode (polarization is free to evolve) / Polarization-Maintaining (polarization stays fixed, more expensive) |
| **MPI** | Multi-Path Interference: unwanted light via a second path |
| **Rayleigh backscatter** | The extremely weak light scattered everywhere in the fiber (~−100 dB/mm). This is how a Luna OBR measures strain/temperature — requires true balanced detectors |
| **Fringe** | An interference fringe = one full period of the cosine |
| **Chirp** | Here: that 1 pm doesn't mean the same frequency change everywhere |
| **Slow bow** | Slow bow-shaped tuning error of the laser over the sweep (~59 pm p-p) |
| **Retrace** | The laser's return sweep after the main sweep |
| **Transform-limited** | As narrow as physics allows at maximum |
| **Aux interferometer** | A second, co-running MZI that serves purely as a **frequency ruler**. **← what your supervisor wants built** |
| **Dispersion / autofocus** | Different wavelengths travel at different speeds → distant peaks smear; correctable with a quadratic phase term |
| **DGD** | Differential Group Delay: in PM fiber the two polarization axes travel at different speeds (~1.5 ps/m). Never splice SM→PM in the return path for this reason |

---

# PART B — Your supervisor's proposal, translated

## B.1 What he's saying, in one paragraph

The problem: behind the circulator there is **~2 m of test fiber** (1 m
patchcord + 1 m from the component). But your Nyquist range is **0.42 m**.
So everything beyond 0.42 m aliases into your measurement window and
creates exactly the mess described in §5 of the handover.

The solution: **sample much more finely than every 1 pm.** But the laser
can only trigger every 1 pm. So: **drop the trigger entirely** and instead
sample on the power meter's *internal clock* (every 1 µs). Then you'll have
far more points — but **no idea anymore at which wavelength each point was
taken.** And you absolutely need that, because the FFT runs over ν.

Hence the new MZI: **it's your ruler.** It measures, at every instant, how
far the laser has swept so far. From its fringe phase you reconstruct the
ν axis.

This is exactly the architecture already recorded as a decision in
handover §6 ("free-running internal-clock sampling — the Luna
architecture"). Your supervisor is saying the same thing, just with
concrete components.

## B.2 ⚠ Terminology trap: there are two "resolutions"

He writes *"we'd need a much better resolution than 1pm"*. The handover
writes *"Resolution 7.0 µm"*. These are **two different things** and you
will definitely mix them up otherwise:

| | What is meant | Formula | Depends on |
|---|---|---|---|
| **Sampling step** (this is what he means) | Spacing between two measurement points, 1 pm | δν = c·δλ/λ² | Trigger / clock rate |
| **Spatial resolution** (handover) | How close two reflectors are allowed to be, 7 µm | Δz = c/(2n_g·Δν_span) | **only** sweep width |

The sampling step determines **range**. The sweep width determines
**resolution**. He wants to improve range, not spatial resolution.

If you keep this distinction clean in conversation, you'll immediately
understand why the calculations below look the way they do.

## B.3 The setup you're supposed to build

```
        EXFO Laser
             │
        ┌────┴────┐  90/10 Splitter
        │90%      │10%
        │         │
        │      ┌──┴──┐  50/50 Splitter
        │      │     │
        │   short arm  long arm  (arm difference ΔL — see B.5!)
        │      │     │
        │      └──┬──┘  50/50 Coupler (recombine)
        │         ├────────► Ch1   ┐ AUX-MZI = your ruler
        │         └────────► Ch2   ┘ (complementary outputs)
        │
   20 dB Coupler ──────► LO arm ──┐
        │                          │
    Circulator ──► TEST FIBER      │
        │  ◄── Reflections ────────┤
        │                          │
        └──────────────► Coupler ──┘
                          ├────────► Ch3   ┐ MEASUREMENT interferometer
                          └────────► Ch4   ┘ (as before)
```

*(Corrected 2026-08-18: the aux lives on Ch1/Ch2 and the measurement on
Ch3/Ch4 — an earlier draft of this diagram had the two pairs swapped.
Check this against the physical wiring, not just against this diagram.)*

**Why 4 channels:** Ch1/Ch2 = the two complementary outputs of the new
Aux-MZI. Ch3/Ch4 = the two complementary outputs of your existing
measurement. For both pairs you perform balanced subtraction (see A.5, step
2) — for the Aux this is especially important because you need its phase
very clean.

**Why 90/10 and not 50/50:** The Aux-MZI has almost no loss (no
circulator, no test fiber), so it needs much less light. The 90% go where
the losses are.

**Important from handover §6 (don't overlook this):** When you resample
using the Aux, you use its **raw, unsmoothed** phase — not a smoothed one
as in `--mode diagnostic`. The entire diagnostic/cosmetic conflict then
disappears, because the correction no longer comes from the measurement
signal itself. **Never smooth the Aux phase** — it must carry the ~4.7 pm
ripple.

## B.4 The calculation: buffer vs. sweep vs. span

The black CoreDAQ fills its buffer with **1 point every 1 µs**. From this
it necessarily follows:

```
Acquisition time T = N_points × 1 µs
Sweep span         = sweep speed × T
δλ per point       = sweep speed × 1 µs
```

And from that, the two quantities you want. Here's the table — **you
should print this out**:

| Points/channel | Sweep | T | Span | δλ | **Range z_max** | **Resolution Δz** |
|---|---|---|---|---|---|---|
| 1,000,000 | 120 nm/s | 1.00 s | 120 nm | 0.12 pm | **3.48 m** | **7.0 µm** |
| 1,000,000 | 60 nm/s | 1.00 s | 60 nm | 0.06 pm | **6.95 m** | **13.9 µm** |
| 1,000,000 | 30 nm/s | 1.00 s | 30 nm | 0.03 pm | 13.90 m | 27.8 µm |
| 250,000 | 120 nm/s | 0.25 s | 30 nm | 0.12 pm | 3.48 m | 27.8 µm |
| 250,000 | 60 nm/s | 0.25 s | 15 nm | 0.06 pm | 6.95 m | 55.6 µm |
| *current (trigger)* | *10 nm/s* | *12 s* | *120 nm* | *1 pm* | *0.42 m* | *7.0 µm* |

**How to read the table:**

- **All rows with 1M points are good enough** for 2 m of test fiber. The
  first row is even identical in quality to now — just with **8× more
  range**.
- The jump from 0.42 m to 3.5 m fully solves your problem: the fiber end
  at 1.05 m, and everything up to 2 m, then lies **within** the range and
  appears as a sharp, honest peak instead of alias garbage.
- **At 250k points per channel you have to sacrifice resolution** — that's
  why your supervisor's question ("1M per channel or split?") is the
  **decisive** question before building. Clarify that first with Giulio.

**On the buffer overflow he mentions:** He phrased it slightly backwards,
but the logic is: the buffer fills up after `N × 1 µs`. If your sweep is
**longer** than that, you lose the end (or it overwrites the beginning).
So always **calculate beforehand**: `sweep duration = span / speed` must
be `< N × 1 µs` — plus some margin for acceleration/trigger latency.

That's why the table shows 60–120 nm/s instead of the previous 10 nm/s:
with a 1 µs clock and a 1 s buffer you **must** sweep faster, otherwise
the span won't fit. This is, incidentally, not a concern acoustically (the
warning in the handover was about sweeping too *slowly*).

## B.5 The most important open number: how long does the long arm need to be?

Your supervisor only says *"one arm should be longer than the other"* —
without a number. But the number isn't arbitrary, and it's **larger than
you think**.

**Two opposing constraints:**

**(1) Upper bound** — the Aux fringes must be sampled cleanly at a 1 µs
clock rate. The fringe frequency is `f = τ_aux × dν/dt`. At a 1 MHz sample
rate, Nyquist is 500 kHz; a comfortable margin is ~7–10 points per fringe.

**(2) Lower bound** — standard OFDR rule: the Aux interferometer should
have a delay **at least as large as the most distant reflector** you want
to measure, otherwise it no longer fully corrects the tuning noise for
distant peaks (the residual error grows roughly with τ_meas/τ_aux).

**The numbers:**

| | Value |
|---|---|
| Test fiber 2 m → round-trip delay τ_meas | 19.6 ns |
| Required Aux arm difference (rule: τ_aux ≥ τ_meas) | **ΔL ≈ 4 m of fiber** |
| Aux fringe frequency at 60 nm/s, ΔL = 4 m | 144 kHz → ~7 points/fringe ✔ |
| Aux fringe frequency at 120 nm/s, ΔL = 4 m | 288 kHz → ~3.5 points/fringe (tight) |
| Your **existing** MZI (92.8 mm) corresponds to | z ≈ 46 mm — **~40× too short** |

> **This is the point you absolutely must clarify with your supervisor:**
> "one arm longer" here means **several meters of fiber** (a small spool),
> not a few centimeters. Ask him whether that's what he means, or whether
> he's counting on sub-fringe interpolation (in which case the Aux can be
> shorter, but the noise suppression gets worse).

**Recommended starting point for discussion:** 1M points/channel, sweep 60
nm/s, span 60 nm, Aux ΔL ≈ 4 m → range 6.95 m, resolution 13.9 µm, 7
points per Aux fringe. Comfortable in every dimension, and once it's
working you can go up to 120 nm/s / 120 nm for the full 7 µm.

**Practical points about the Aux-MZI:**

- The long arm should be **thermally quiet** (spooled into a box, not
  lying loosely across the table). Its length *is* your ruler — if it
  drifts, your distance axis drifts.
- In SM fiber the two arms' polarizations can drift apart → fringe
  contrast collapses (**polarization fading**). If this happens: put a
  polarization controller in one arm, or switch to PM fiber medium-term
  (handover §6, item 3).
- **Calibrate τ_aux once**: run a sweep in the *previous* trigger mode,
  where you get the wavelength table, and compute `total accumulated phase
  / (2π × frequency span) = τ_aux`. This is exactly how the 454.5 ps of
  the existing MZI were determined (handover §6).

## B.6 The trigger changeover

| | before | new |
|---|---|---|
| EXFO | Trigger every 1 pm | **Single trigger**: one pulse at sweep start |
| CoreDAQ | red, 1 sample per trigger, ~130k buffer | **black**: 1 trigger → then free-running every 1 µs, 1M buffer |
| ν axis comes from | EXFO wavelength table | **Aux-MZI phase** |

His argument for using the black CoreDAQ instead of reconfiguring the red
one is pragmatic: the black one **can only do exactly this** out of the
box, which is what you need. Swapping devices is simpler than fighting
firmware settings. **Giulio knows this device** — ask him early, not only
once things get stuck.

**What you lose in the process and must replace:** without a trigger there
is **no wavelength table anymore**. Your ν axis only exists if the Aux
works. So: **build the Aux first and verify it in the previous trigger
mode**, only then switch to single-trigger. Otherwise you're debugging two
new things at once.

You'll still read the absolute start wavelength off the laser — per
handover §6, item 4, its absolute accuracy is ~100× better than needed.

## B.7 Your checklist

**First — costs nothing, clarifies a lot (from handover §7):**

1. **Put a tape measure on the table.** Is the 46.4 mm peak a 92.8 mm arm
   difference, or a reflector at 46.4 mm? Unresolved, blocks understanding
   of the topology.
2. **Terminate the fiber end** (index-matching gel, tight winding, or a
   < 35 cm patchcord) and verify before/after with `diagnose_artifacts.py
   compare`. Prediction: the 6–32 cm band drops into the noise. This
   confirms the aliasing diagnosis with your own hands and gives you clean
   data immediately.

**To clarify before building:**

3. **1M points per channel, or split across 4 channels?** (→ Giulio).
   Determines the entire table in B.4.
4. **What is the maximum sweep speed the EXFO can do?** You will likely
   need 60–120 nm/s.
5. **How long should the long Aux arm be?** (→ Supervisor, armed with the
   numbers from B.5). My calculation says ~4 m.
6. Does the black CoreDAQ really have 4 channels and fixed measurement
   ranges? **Autoranging mid-sweep is fatal** (phase jump → unwrap broken
   to the end of the scan, handover §7).

**During the build:**

7. Build the Aux-MZI, calibrate τ_aux **in the old trigger mode**.
8. Only then switch to single-trigger + free-run.
9. Before every pipeline change: test on synthetic data with known ground
   truth (`simulate_aux_data.py`).
10. After every scan: **check peak width.** ≤ ~21 µm with Kaiser-12,
    otherwise stop and investigate.

---

## If you take away one thing from this document

```
Range / Resolution = N/2
```

Resolution comes from the **sweep width**. Range comes from the **point
spacing**. Everything your supervisor is proposing serves exactly one
purpose: getting more points into the same sweep — and the Aux-MZI is the
price you pay for the laser no longer telling you where it currently is.
