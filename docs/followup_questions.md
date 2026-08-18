# OFDR — Answers to your follow-up questions

**For:** Carolina · **As of:** 2026-08-18 · Supplement to `docs/handover_explained.md`

Included (in the same folder):

- `calc_sweep_parameters.py` — works through every parameter combination for you, including the aux check
- `simulate_aux_data.py` — generates synthetic 4-channel data with known ground truth
- `process_reflectogram_aux.py` — the 4-channel pipeline (answer to question 10, already working)
- `vergleich_alt_neu.png` — old vs. new, on simulated data

---

## 1. Why can't you just drop the `/2` in the script?

Because the `/2` isn't a setting — it's the **definition of the axis**. And
because an axis can only carry one definition, while the detector
simultaneously carries both kinds of signal.

What the FFT gives you is always just **τ**, the delay. That part is
unambiguous — there's nothing to decide there. Only the step from τ to a
distance needs an assumption, namely *"how many times did the light travel
this stretch?"*:

| Signal type | Path | Conversion |
|---|---|---|
| **Reflection** (connector, fiber end) | there **and** back | z = c·τ/(2·n_g) |
| **Transmission** (MZI arm) | one-way only | ΔL = c·τ/n_g |

In a reflectometer, *most* signals are reflections. So `/2` is correct. If you
removed it, **all real reflectors would show up at twice their true
distance** — you'd have traded one problem for a bigger one.

The MZI signal is the one exception, and exceptions can't be handled by
changing the axis — only by knowledge. That's why the script does exactly the
right thing: it **prints both readings** for the main peak
(`process_reflectogram.py`, line 222–224):

```
delay 454.50 ps | as transmission-MZI imbalance: 92.82 mm fiber
                | as reflection: 46.41 mm fiber one-way
```

So the information is there. What's missing is the tape-measure reading that
tells you which of the two lines is correct.

**And the good news:** once the aux MZI sits on its own channels (Ch3/Ch4),
the problem disappears structurally. Then the transmission signal is no
longer on the measurement channels at all, and the axis on Ch1/Ch2 is
unambiguously a reflection axis. That's another, often overlooked, argument
for the rebuild.

---

## 2. Where does the 122 MHz per 1 pm step come from?

Pure unit conversion. Two routes, both giving the same answer:

**Route A — differentiate.** ν = c/λ, so

```
|dν/dλ| = c/λ²

c/λ²  at λ = 1565 nm  =  2.998·10⁸ / (1.565·10⁻⁶)²
                      =  2.998·10⁸ / 2.449·10⁻¹²
                      =  1.224·10²⁰  Hz per meter
```

That's the conversion factor. Now multiply by the step size:

```
δν = 1.224·10²⁰ Hz/m  ×  1 pm
   = 1.224·10²⁰  ×  1·10⁻¹² m
   = 1.224·10⁸ Hz
   = 122 MHz
```

**Route B — via relative changes (I find this more intuitive).** For ν = c/λ,
the *relative* change in λ and in ν is the same size (just with opposite
sign):

```
1 pm out of 1565 nm      =  1·10⁻¹² / 1.565·10⁻⁶  =  6.39·10⁻⁷   (i.e. 0.64 ppm)
optical frequency        =  c/1565 nm  =  191.6 THz
0.64 ppm of 191.6 THz    =  122 MHz  ✔
```

Rule of thumb: **1 pm ≈ 122 MHz at 1550 nm.** And because λ² sits in the
denominator, the factor is larger at the blue end: 132 MHz at 1505 nm, 114 MHz
at 1625 nm. That's exactly the 18% chirp.

---

## 3. What does "the reflector at τ creates an oscillation with period 1/τ" mean?

It means: **by how much does the laser have to detune its frequency for the
signal to go once from bright to dark and back to bright?**

From the base equation: the phase is `2π·τ·ν`. One full period is reached
when this phase changes by 2π:

```
2π · τ · Δν = 2π      →      Δν = 1/τ
```

The units work out: τ in seconds → 1/τ in hertz. So it's a period
**measured in hertz of optical frequency** — not in seconds. Unusual, but
that's exactly the point: we're working in the frequency domain, ν is our
x-axis.

**Concretely, for your setup:**

| Reflector | τ | 1/τ | in wavelength | Points per fringe at 1 pm |
|---|---|---|---|---|
| MZI peak (46.4 mm) | 454.5 ps | 2.20 GHz | **18.0 pm** | **18.0** ✔ |
| Fiber end at 1.05 m | 10.28 ns | 97.3 MHz | **0.79 pm** | **0.79** ✘ |
| Component at 2.0 m | 19.59 ns | 51.1 MHz | **0.42 pm** | **0.42** ✘ |

(The handover doc quotes "17 pts/fringe" for the MZI — close to my 18; the
difference comes from the choice of λ used in the conversion.)

**Read the last column:** the setup can sample the MZI signal without any
trouble (18 points per fringe, plenty). The fiber end at 1.05 m gets **less
than one** point per fringe. Nyquist requires at least 2. That's why it
aliases — and why the cutoff sits at 0.42 m, where you get exactly 2 points
per fringe.

**Intuitively:** farther away = faster fringes. Too far away = fringes faster
than your sampling = you see a false, slower oscillation. Wagon-wheel effect
in film.

---

## 4. What does the 59 pm peak-to-peak mean?

You can see the bow shape — the question is what's bad about it. Answer in
three steps.

**Step 1: what the number describes.** The laser claims that at sample point
k it's at `λ_start + k · 1 pm`. The right half of `george/real_data_diagnosis.png`
shows where it **actually** is: up to +18 pm above that (in the middle of the
sweep) and up to −40 pm below it (at the end). Total swing ≈ 59 pm.

So your x-axis is wrong. Not randomly noisy, but smoothly and systematically
warped.

**Step 2: why 59 pm out of 120,000 pm is still a lot.** In relative terms
that's 0.05% — sounds harmless. But the FFT doesn't care about relative
errors in the axis, it cares about **phase errors**. And that one is large:

```
59 pm  ×  122 MHz/pm             =  7.2 GHz frequency error
Phase error on the 454.5-ps tone:
  Δφ = 2π · τ · Δν
     = 2π · 454.5 ps · 7.2 GHz
     = 20.6 rad
     = 3.3 full periods
```

So the tone that *should* be a perfectly uniform oscillation drifts by
**3.3 full periods** out of step over the course of the sweep. To the FFT,
that's no longer a single frequency but a frequency mixture → the peak smears
out.

**Step 3: why that turns into ±0.8 mm.** A peak sits at the *slope* of the
phase (phase per frequency = delay). Because the error is a smooth bow, it
has a different slope in the middle than at the edge:

```
Bow amplitude A ≈ 10 rad, sweep span S = 14.7 THz
maximum slope deviation       dτ ≈ A·π/S  ≈ 2.2 ps
                               dz  = c·dτ/(2·n_g)  ≈ 0.22 mm
```

Same order of magnitude as the ±0.8 mm quoted in the handover doc (my bow
model is a smooth half-sine; the real curve is sharper-edged and yields
correspondingly more).

**Compare to the target width: 21 µm.** Without correction, the peak ends up
roughly 30–40 times broader than it should be. That's the whole point.

**And this is exactly why the aux MZI isn't a refinement but a necessity.**
Right now the bow is estimated from the measurement signal itself, which
works but also warps weak real features along with it (hence `diagnostic` vs.
`cosmetic`). The aux measures it **directly and independently**.

---

## 5. What is `cosmetic` for?

Short answer: **not for measuring. As an upper bound, and as a test for
whether the axis is at fault.** Three uses:

**(a) Proving the hardware can do it.** If you get a 21 µm peak in `cosmetic`
mode, the reflector *is* physically sharp and your optics deliver full
coherence. That's a statement about the setup, not about the reflector
positions.

**(b) The actual diagnostic value: "is the frequency axis the problem?"**

| `diagnostic` | `cosmetic` | Diagnosis |
|---|---|---|
| broad | **sharp** | The axis is at fault. Better correction → aux MZI |
| broad | **also broad** | The axis is *not* at fault. Look in the physics: dispersion, polarization, or the reflector really is extended |

That's a genuine discriminator, in the same family as the four tests in
`diagnose_artifacts.py`. It stops you from spending days optimizing the
software while the problem is sitting on the bench — or vice versa.

**(c) Pictures for the supervisor.** "We achieve 21 µm, 15% above the
transform limit" is a legitimate figure, and it comes out of this mode.

**Why you must not use it for diagnosis:** the correction follows the
*raw phase* of the main signal. Anything else present in the signal — which
is exactly the weak reflectors you're looking for — warps that phase too, and
gets partly **cancelled out** by the correction. Measured effects were −8 dB
loss on real features plus position shifts. The mode makes the picture
prettier by erasing the things you actually care about.

**And:** after the rebuild, `cosmetic` becomes unnecessary. Aux referencing
gives you sharp peaks *and* correct positions at the same time, because the
correction no longer comes from the measurement signal itself. The whole
mode is a symptom of the current stopgap setup.

---

## 6. Why is the signal in the PNG in the negative dB range?

Because the scale is **relative to the largest peak**, not absolute. In the
code (`process_reflectogram.py`, line 204):

```python
db = 20 * np.log10(R / R.max() + 1e-15)
```

It's divided by `R.max()`. The strongest peak is therefore, by definition,
`20·log10(1) = 0 dB`, and everything else is smaller than 1 → the logarithm
is negative. There are no positive values in this plot at all — that's
structural, by design.

**What the numbers mean** (note: `20·log10`, because `R` is an *amplitude*,
not a power):

| dB | Amplitude relative to the main peak |
|---|---|
| 0 | 1 : 1 (this is the main peak itself) |
| −21 | 1 : 11 |
| −45 | 1 : 178 |
| −53 | 1 : 447 |
| −80 | 1 : 10,000 |

So: the artifact band at −45 dB is ~180× weaker than the main peak. The
Kaiser side lobes at −80 dB are 10,000× weaker — which is why Kaiser is the
right choice when you're looking for weak real reflections.

**Why it's done this way:** what you care about is the **ratio** (is this
reflector strong or weak *compared to* the main signal), not the absolute
power in mW. For absolute values (return loss in dB) you'd need a calibration
against a known reference reflection — you don't have that yet, but it's also
not on the list.

---

## 7. 250k points per channel with a 120 nm span — does that work?

**Arithmetically yes, practically no.** And the reason is interesting. First
the calculation path, then the numbers.

### The calculation path (always the same, five lines)

The fixed constraint is the clock rate: 1 point every 1 µs. Everything
follows from that:

```
1)  T        = N × 1 µs                     acquisition duration
2)  v        = Span / T                     required sweep speed
3)  δλ       = v × 1 µs                     point spacing in wavelength
4)  δν       = δλ × 122 MHz/pm              point spacing in frequency
5a) z_max    = c / (4·n_g·δν)               RANGE
5b) Δz       = c / (2·n_g·Δν_span)          RESOLUTION
```

### Plugged in for 250,000 points and 120 nm

```
1)  T   = 250,000 × 1 µs                    = 0.25 s
2)  v   = 120 nm / 0.25 s                   = 480 nm/s
3)  δλ  = 480 nm/s × 1 µs                   = 0.48 pm
4)  δν  = 0.48 × 122.4 MHz                  = 58.8 MHz
5a) z_max = 3·10⁸ / (4 × 1.468 × 58.8·10⁶)  = 0.87 m
5b) Δν_span for 120 nm = 14.7 THz
    Δz  = 3·10⁸ / (2 × 1.468 × 14.7·10¹²)   = 7.0 µm
```

Check: 0.87 m ÷ 7.0 µm = 125,000 = N/2 ✔

### And now the three reasons why it doesn't work anyway

1. **0.87 m of range isn't enough** for the 2 m test fiber. Better than the
   current 0.42 m, but the component at 2 m still aliases. Problem not
   solved, just shrunk.
2. **480 nm/s** — can the EXFO even do that? You need to check. Many
   tunable lasers top out at 100–200 nm/s.
3. **The killer: the aux MZI no longer fits.** At 480 nm/s and ΔL = 4 m, the
   aux beats at **1.15 MHz**. Your sampling rate is 1 MHz, so Nyquist is
   500 kHz. **The ruler itself aliases.** That breaks the frequency axis,
   and without a frequency axis there's no result at all.

The script tells you this directly:

```
$ python calc_sweep_parameters.py --points 250000 --span 120 --dl 4
  RANGE       z_max = 0.869 m
  RESOLUTION  dz    = 6.95 um
  --- Aux MZI with dL = 4.00 m ---
  appears at          z = 2.000 m (must be < z_max = 0.869 m)
  beat frequency      1150.8 kHz  = 0.9 points per fringe
    ERROR: aux is beyond Nyquist -- it aliases itself.
    ERROR: fewer than 2 points per aux fringe -- phase cannot be reconstructed.
```

### What does work with 250k

| Span | v | z_max | Δz | Aux ΔL=4 m |
|---|---|---|---|---|
| 120 nm | 480 nm/s | 0.87 m | 7.0 µm | ✘ aliases |
| 60 nm | 240 nm/s | 1.74 m | 13.9 µm | ✘ aliases |
| 30 nm | 120 nm/s | 3.48 m | 27.8 µm | marginal (3.5 pts/fringe) |
| **15 nm** | **60 nm/s** | **6.95 m** | **55.6 µm** | ✔ 7.0 pts/fringe |

Compared with 1M points:

| Span | v | z_max | Δz | Aux ΔL=4 m |
|---|---|---|---|---|
| 120 nm | 120 nm/s | 3.48 m | **7.0 µm** | marginal (3.5 pts/fringe) |
| **60 nm** | **60 nm/s** | **6.95 m** | **13.9 µm** | ✔ 7.0 pts/fringe |

**Bottom line: the buffer difference is a factor of 4 in resolution** (55.6 µm
vs. 13.9 µm at the same comfort level). That's why "1M per channel, or split
across 4?" is the first question for Giulio, before any soldering happens.

Play with it yourself:

```bash
python calc_sweep_parameters.py --table
python calc_sweep_parameters.py --points 1000000 --span 60 --dl 4
```

---

## 8. How do you arrive at 4 m arm-length difference? Step by step

### Lower bound — accuracy

**Step 1.** What's the farthest thing you want to measure? → 2 m test fiber.

**Step 2.** Its round-trip delay (there **and** back, since it's a
reflection):

```
τ_meas = 2 · n_g · z / c = 2 · 1.468 · 2 m / 2.998·10⁸ = 19.59 ns
```

**Step 3.** The rule: **τ_aux ≥ τ_meas.**

Why? The aux is your ruler, and its tick marks are its fringes. The spacing
of the tick marks in frequency is `1/τ_aux` (question 3). The reflector you
want to measure oscillates at `1/τ_meas`. A ruler whose tick marks are
coarser than what you're trying to measure can no longer resolve the
distortion there — the residual error grows roughly as `τ_meas/τ_aux`.

**Step 4.** Convert. The aux is a **transmission** MZI — the light passes
through the longer arm only once — **no factor of 2**:

```
τ_aux = n_g · ΔL / c        →      ΔL = c · τ_aux / n_g
                                      = 2.998·10⁸ · 19.59 ns / 1.468
                                      = 4.00 m
```

**The shortcut worth remembering:**

```
τ_aux = τ_meas
n_g·ΔL/c = 2·n_g·z/c
        ΔL = 2 · z_max
```

> **ΔL = 2 × the farthest reflector.** Measuring 2 m → 4 m arm difference.
> Measuring 3 m → 6 m. It's that simple.

Equivalent, and more practical to check: the aux **appears in the
reflectogram at z = ΔL/2**. This position must be at least as far as your
farthest reflector. ΔL = 4 m → aux peak at 2.0 m ✔

### Upper bound — sampling

The aux is itself a signal on the detector. So it too must lie **within the
Nyquist range**:

```
ΔL/2  <  z_max          i.e.        ΔL < 2 · z_max
```

And expressed in frequency — the aux's fringe rate:

```
dν/dt  = (c/λ²) · dλ/dt = 122.4 MHz/pm · sweep speed
f_aux  = τ_aux · dν/dt
Points per fringe = 1 MHz / f_aux
```

| ΔL | Sweep | f_aux | Points/fringe | Verdict |
|---|---|---|---|---|
| 4 m | 30 nm/s | 72 kHz | 13.9 | very comfortable |
| **4 m** | **60 nm/s** | **144 kHz** | **7.0** | **good** |
| 4 m | 120 nm/s | 288 kHz | 3.5 | above Nyquist, but tight |
| 4 m | 480 nm/s | 1151 kHz | 0.9 | broken |
| 8 m | 60 nm/s | 288 kHz | 3.5 | tight |

### So, put together

```
2 · z_meas   ≤   ΔL   <   2 · z_max
   4 m       ≤   ΔL   <   13.9 m        (at 1M points, 60 nm/s)
```

**ΔL = 4 m is both the lower bound and a good choice** — enough ruler for a
2 m measurement range, 7 points per fringe, and plenty of headroom if the
test fiber gets longer later on.

**If your supervisor wants to be more conservative:** the literature often
quotes `τ_aux ≥ 2·τ_meas` (i.e. ΔL = 8 m). That applies when you **clock
hard on fringe edges** (one sample point per aux fringe) — in that case the
aux itself acts as the sampler and must satisfy Nyquist for the measurement
signal. We instead sample at 7 points per fringe and interpolate the phase,
so ΔL = 4 m is sufficient. **Clarify with him which of the two variants he
means** — it's the difference between a 4 m and an 8 m spool.

---

## 9. Do I need to write new software for the 4 channels?

**No — and I've already written and tested it for you.** `process_reflectogram_aux.py`
is in the folder and works. But read through what changes, so you understand
it and can adapt it.

### What changes in the code

| | |
|---|---|
| **Unchanged** | balanced subtraction (now called twice), windowing, FFT, peak list, width check, plot |
| **Removed** | `to_uniform_nu()` — there's no longer a wavelength table |
| **Removed** | `phase_correct()` with `diagnostic`/`cosmetic` — no more self-referencing, so the whole conflict disappears |
| **New** | `resample_on_aux()` — ~25 lines |

Net result: **less** code than now. It's a change, not a rewrite from
scratch.

### What `resample_on_aux()` does, in four lines of logic

```python
phi = np.unwrap(np.angle(analytic(aux)))      # aux phase, NOT smoothed
phi = np.maximum.accumulate(phi)              # enforce monotonicity
phi_u = np.linspace(phi[0], phi[-1], m)       # uniform phase grid
y = PchipInterpolator(phi, meas)(phi_u)       # resample the measurement onto it
```

The core idea in one sentence: the aux phase is `φ(t) = 2π·τ_aux·ν(t) + const`.
Equal steps in φ are therefore **equal steps in ν** — regardless of how
non-uniformly the laser actually swept. After that, the axis is done and the
FFT is valid.

And `δν = Δφ / (2π·τ_aux)` gives you both range and resolution, without ever
needing a wavelength anywhere.

### The proof that it works

I generated synthetic 4-channel data with **known ground truth**
(1M points, 60 nm/s, aux ΔL = 4 m, including the measured 59-pm bow and
tuning ripple) and ran it through the pipeline:

```
   z_true  0.0460 m (  0.0 dB) -> found  0.0460 m, error  +4.1 um (+0.30 cells)
   z_true  0.3000 m (-25.0 dB) -> found  0.3000 m, error  -4.8 um (-0.35 cells)
   z_true  1.0500 m (-30.0 dB) -> found  1.0500 m, error  +3.6 um (+0.26 cells)
   z_true  2.0000 m (-35.0 dB) -> found  2.0000 m, error  +4.2 um (+0.31 cells)
   largest error: 4.8 um (one cell = 13.6 um)   -> PASSED

Main peak 46.0041 mm, -3 dB width 40.9 um (window limit ~35.5 um)
```

All four reflectors within better than half a resolution cell, amplitudes
within 0.5 dB, peak width 15% above the transform limit — **the same health
metric as now, but with 6.7 m range instead of 0.42 m.**

`vergleich_alt_neu.png` shows both overlaid: on top, the same set of
reflectors in the old trigger mode (1.05 m and 2.0 m fold back to 0.22 m and
0.33 m and fill the whole range with junk), on the bottom, aux-referenced
(four clean peaks). This is a simulation, not a measurement — but it's
exactly the artifact from §5 of the handover doc, reproduced and then fixed.

### What you still need to do yourself

1. **Adapt the file format from the black CoreDAQ.** `process_reflectogram_aux.py`
   reads `.npz` and the known JSON structure with `Ch1..Ch4 [mW]`. What the
   black device actually exports, you'll only know once it's connected
   (→ Giulio). It's the `load()` function, right at the top, ~15 lines.
2. **Check which channels are aux and which are measurement.** It's stated
   as an assumption in the file header: Ch1/Ch2 measurement, Ch3/Ch4 aux.
3. **Calibrate τ_aux** and pass it via `--tau-aux-ns`. Without that, the
   script estimates it from `--dl`, which is accurate to about 1% for
   positions — good enough for first pictures, not for numbers in a report.

### One side finding you should know about

While simulating, I noticed: the fast tuning ripple quoted in the handover
doc, **90–110 MHz rms at ~4.5 pm period, sits dangerously close to a hard
limit.** If the slope of this ripple exceeds the nominal sweep slope, the
laser momentarily runs **backwards** in ν. Then the aux phase is no longer
monotonic, and **no** method can resolve that anymore — not even the aux.

```
Limit:  A_max = (c/λ²) · period / (2π)
             = 122.4 MHz/pm · 4.5 pm / (2π)
             ≈ 88 MHz
```

The handover doc itself warns that this number was "measured via the main
tone" and is therefore contaminated by additive noise sources — so it's
probably overestimated. **But it's an assumption the entire rebuild rests
on.** The good news: the aux MZI measures this ripple cleanly and directly.
So: as soon as the aux is set up, **measure the ripple first** — before even
switching to free-running mode. `simulate_aux_data.py --ripple-mhz 100`
shows you what the error pattern looks like if it turns out to be too much.

---

## 10. What exactly you should do now

### Phase 0 — today/tomorrow, without touching any screws

**0.1** Tape measure on the 46.4 mm peak. Is it a 92.8 mm arm difference, or a
reflector at 46.4 mm? Count up the fiber and connector lengths.
→ *Result goes in a lab notebook — this has been an open question since the
handover doc.*

**0.2** Record a reference scan in the current state and save it. Then:

```bash
python diagnose_artifacts.py tuning  scan_today.json
python process_reflectogram.py       scan_today.json --mode diagnostic
```

Note down: peak width (should be ≤ ~21 µm), bow in pm, fast ripple in MHz.
→ *This is your baseline. Without it you can't later tell whether anything
has improved.*

**0.3** Terminate the fiber end (index-matching gel on the APC end face, or a
tight wind over the last few centimeters, or a patch cord < 35 cm) and
record a new scan. Then:

```bash
python diagnose_artifacts.py compare scan_today.json scan_terminated.json
```

Expectation: the band at 6–32 cm drops into the noise.
→ *This confirms the aliasing diagnosis yourself, and from then on you have
clean data to work with. Highest payoff per effort of anything on this
list.*

**0.4** As a warm-up: run the simulation so you're familiar with the pipeline
before touching hardware.

```bash
python simulate_aux_data.py --points 400000 --speed 60 --dl 4 --out sim.npz
python process_reflectogram_aux.py sim.npz --dl 4 --zmax 2.5
```

### Phase 1 — clarify questions before anything is bought or soldered

**1.1 → Giulio:** Does the black CoreDAQ have **1M points per channel, or 1M
split across 4**? *(Factor of 4 in resolution. Most important number.)*

**1.2 → Giulio:** All 4 channels simultaneously? Fixed measurement ranges
settable? *(Autoranging mid-sweep ruins the scan.)* What export format?

**1.3 → EXFO manual:** maximum sweep speed? *(You'll likely need
60–120 nm/s.)* And: how do you enable single-trigger mode?

**1.4 → Supervisor:** ΔL = 4 m or 8 m? Show him section 8 with the
calculation. *(He hasn't given a number, and his existing MZI at 92.8 mm is
~40× too short.)*

**1.5 → Supervisor:** confirm the starting configuration. My suggestion:

| | |
|---|---|
| Points | 1M per channel |
| Sweep | 60 nm/s |
| Span | 60 nm |
| Aux ΔL | 4 m |
| ⇒ Range | 6.95 m |
| ⇒ Resolution | 13.9 µm |
| ⇒ Aux | 7 points per fringe |

*Only once this is working, move to 120 nm/s and 120 nm span for the full
7 µm.*

### Phase 2 — build, in this order

**2.1** 90/10 splitter right after the laser. 90% into the existing setup.
*Before/after scan: the 10% loss shouldn't break anything.*

**2.2** Build the aux MZI: 50/50 → two arms with ΔL ≈ 4 m → 50/50 back
together, both outputs onto Ch3/Ch4. Lay the long spool somewhere thermally
quiet (in a box, not loose on the bench — its length *is* your ruler).

**2.3** Check the aux **in the old trigger mode**: clean sine wave? Contrast
stable across the whole sweep? *(If it drops out somewhere → polarization
fading → put a polarization controller in one arm.)*

**2.4** Calibrate τ_aux, in the old trigger mode with the wavelength table:

```
τ_aux = total accumulated aux phase / (2π · frequency span)
```

**2.5** Measure the fast ripple with the aux (section 9, side finding). If
it's under ~50–80 MHz rms → green light.

**2.6** **Only now** switch to single-trigger + free-running mode. One new
thing at a time, otherwise you're debugging two unknowns at once.

**2.7** Buffer check on every scan: `Span / speed < N × 1 µs`.
If the fringes cut off in the last part of the buffer, the sweep was too
long.

**2.8** First real measurement through `process_reflectogram_aux.py`. Check:
peak width ≲ 1.5× the window limit, aux ≥ 4 points/fringe, and — this is the
moment of truth — **the fiber end at ~1 m shows up as a sharp peak at its
true position, not as a band around 21 cm anymore.**

### Standing rules

- After every scan: **check the peak width.** Regression → stop and
  investigate.
- Before every pipeline change: **test on synthetic data with known ground
  truth.** You now have the generator for that.
- Every hypothesis about an unexpected feature gets a **measurement**, not a
  story. The four tests in `diagnose_artifacts.py` are each one command. The
  handover doc records that three plausible theories were wrong, and each
  was killed by a measurement.
