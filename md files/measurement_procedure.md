# Procedure Plan — which program, when, with which setup

**For:** Carolina · 18.08.2026

---

## ⚠ First: a correction to my own plan

In the list you quoted, steps 2.3–2.5 were listed as *"check and calibrate
the aux in the old trigger mode"*. **That doesn't work.** I ran the numbers:

```
Old trigger mode, 1 pm step        ->  range 0.42 m
An aux with 4 m arm difference     ->  appears at 2.00 m
```

The aux would fall far outside the range and would **alias against itself**.
In trigger mode, at most an aux with **ΔL < 0.83 m** could be checked — but
for a 2 m measurement path you need 4 m.

**The solution is simpler than the original plan:** τ_aux can be calibrated
in free-run mode by **counting fringes**. For this you don't need a
wavelength table at all — just the start and end wavelength you set on the
laser:

```
τ_aux = number of fringes counted / frequency span of the sweep
```

Accuracy: ~7 ppm. Completely sufficient. There's now a dedicated program for
this, `check_aux_interferometer.py`.

The order below has been corrected accordingly.

---

## The six programs at a glance

| Program | Instrument | Trigger | Channels | What it's for |
|---|---|---|---|---|
| `process_reflectogram.py` | red | old mode | Ch1, Ch2 | the existing analysis |
| `diagnose_artifacts.py` | red | old mode | Ch1, Ch2 | the four diagnostic tests |
| `calc_sweep_parameters.py` | — | — | — | desk-based planning |
| `simulate_aux_data.py` | — | — | — | simulation for practice and testing |
| **`check_aux_interferometer.py`** | **black** | **free-run** | **Ch1, Ch2** | **ruler acceptance test** |
| **`process_reflectogram_aux.py`** | **black** | **free-run** | **Ch1–Ch4** | **the new analysis** |

---

# Step by step

## S0 · Baseline · *old setup, nothing changed*

**Setup:** as it is now. **Instrument:** red CoreDAQ, trigger mode (1 pm),
2 channels, fixed measurement ranges.

```bash
python process_reflectogram.py  scan_S0.json
python diagnose_artifacts.py tuning scan_S0.json
```

**Record:**

| | expected |
|---|---|
| Main peak position | ~46.4 mm |
| Main peak width | ≤ ~21 µm |
| Slow drift | ~59 pm p-p |
| Fast jitter | ~90–110 MHz rms |

**Why:** This is your before picture. Without these four numbers you'll have
nothing to compare against later. Keep the file.

---

## S1 · Kill the fiber end · *only the test fiber touched*

**Setup:** index-matching gel on the end face, or tightly coil the last few
centimeters, or a patch cord < 35 cm. **Instrument:** unchanged.

```bash
python process_reflectogram.py  scan_S1.json
python diagnose_artifacts.py compare scan_S0.json scan_S1.json
```

**Expectation:** The median in the 6–32 cm range drops significantly. The
dirt band disappears into the noise.

**Proceed if:** The band has gotten smaller. This confirms the aliasing
diagnosis yourself — and from here on you have clean data.

*(If it does NOT get smaller: stop. Then the explanation in the handover
isn't the whole truth after all, and you need to know that before
rebuilding.)*

---

## S2 · Fix parameters · *desk work, no instrument*

```bash
python calc_sweep_parameters.py --table
python calc_sweep_parameters.py --points 1000000 --span 60 --dl 4
```

Enter what Giulio told you about the buffer question. If 1 million **per
channel**:

| | |
|---|---|
| Points | 1,000,000 |
| Sweep speed | 60 nm/s |
| Sweep from / to | 1505 → 1565 nm |
| Aux arm difference | 4 m |
| → Range | 6.95 m |
| → Resolution | 13.9 µm |
| → Aux appears at | 2.00 m |
| → Points per aux fringe | 7.0 |

If 1 million **divided by 4**, recalculate with `--points 250000` — then
you'll need to drop to 15 nm span and end up at 55.6 µm resolution.

**Print this out and hang it at the setup.**

---

## S3 · Practice on the simulation · *desk work, no instrument*

Do this **before** touching hardware, so you know the outputs.

```bash
python simulate_aux_data.py --points 400000 --speed 60 --dl 4 \
       --lam0 1505 --out sim.npz

python check_aux_interferometer.py sim.npz --lam-start 1505 --lam-stop 1529

python process_reflectogram_aux.py sim.npz --tau-aux-ns 19.587 --zmax 2.5
```

*(The 1529 is start + 400,000 × 1 µs × 60 nm/s = 1505 + 24 nm.)*

**Expectation:** `check_aux_interferometer` reports good contrast, good
monotonicity, τ_aux ≈ 19.59 ns, ΔL ≈ 4.00 m. `process_reflectogram_aux`
finds four reflectors and ends by printing `-> PASSED`.

To see what a **broken** case looks like:

```bash
python simulate_aux_data.py --points 400000 --ripple-mhz 100 --out sim_bad.npz
python check_aux_interferometer.py sim_bad.npz --lam-start 1505 --lam-stop 1529
```

There, the monotonicity check will raise an alarm. That's what it looks like
when the laser jitters too much.

---

## S4 · Install the 90/10 splitter · *a single change*

**Setup:** 90/10 right after the laser. 90% into the existing setup, 10%
left open for now (or terminated). **Instrument:** still red CoreDAQ, still
trigger mode.

```bash
python process_reflectogram.py scan_S4.json
```

**Proceed if:**

- Main peak is at the same position as in S0 (± a few µm)
- Peak width unchanged
- Everything ~0.5 dB darker — that's the inserted loss, normal

**Why:** You're checking that the splitter hasn't broken anything, **using
tools you already know**. One new thing at a time.

---

## S5 · Build and accept the aux · *now the switchover happens*

**Setup:** 10% branch → 50/50 splitter → two arms, one 4 m longer → 50/50
coupler → both outputs to **Ch1 and Ch2**.
Put the 4 m coil in a box, not loose on the table.

**Instrument:** black CoreDAQ. Single trigger on the EXFO, free-run with
1 µs. 4 channels, **fixed measurement ranges** (autoranging mid-sweep
destroys the scan). Sweep 1505 → 1565 nm at 60 nm/s.

```bash
python check_aux_interferometer.py scan_S5.json --lam-start 1505 --lam-stop 1565
```

**The program answers four questions. Proceed only if all four pass:**

| Check | Criterion | If not |
|---|---|---|
| **1 Contrast** | Visibility > 0.4 everywhere | Polarization fading → add a polarization controller in one arm |
| **2 Monotonicity** | 0% backward-running steps | **STOP.** The laser is jittering too much. Find the cause before continuing |
| **3 Calibration** | τ_aux ≈ 19.6 ns, ΔL ≈ 4 m, ≥ 4 points/fringe | ΔL doesn't match the installed fiber → remeasure |
| **4 Laser quality** | Jitter clearly below the stated limit | see check 2 |

At the end the program prints the line you need:

```
--> USE THIS NUMBER GOING FORWARD:  process_reflectogram_aux.py ... --tau-aux-ns 19.5871
```

**Write it in the lab notebook.** It's your reference standard.

**Buffer check for this scan:**

```
Sweep duration  = 60 nm / 60 nm/s     = 1.00 s
Buffer duration = 1,000,000 × 1 µs    = 1.00 s     -> exact match, no margin
```

Better to build in some margin: sweep 55 nm instead of 60, or 65 nm/s
instead of 60. If the flicker cuts off near the end of the file, the sweep
was too long.

> **Why check 4 matters more than it looks:** The handover cites 90–110 MHz
> jitter, measured *through the measurement signal* and therefore distorted
> by stray light. At ~88 MHz with a 4.5 pm period, the laser occasionally
> runs backward, and at that point **no** method can help anymore. The aux
> measures this number cleanly for the first time. This is the assumption
> the whole rebuild rests on — which is why it gets checked before anyone
> trusts a reflectogram.

---

## S6 · First real measurement · *same scan, now the measurement channels*

You don't need a new scan — `scan_S5.json` already contains all four
channels.

```bash
python process_reflectogram_aux.py scan_S5.json --tau-aux-ns 19.5871 --zmax 3.0
```

**Three checks, in this order:**

**Check 1 — is the aux where it should be?**
In the reflectogram there must be a peak at **z = ΔL/2 = 2.00 m**. That's
the aux itself, coupling in through the 20 dB coupler. If it sits elsewhere,
τ_aux is wrong.

**Check 2 — the best test you have:**
The old, known peak at **46.4 mm** must be at **exactly the same position**
as in S0 and S4. You know the answer in advance, measured with a completely
different method. If both agree, the new chain works end-to-end.

Tolerance: one resolution cell, i.e. ~14 µm.

**Check 3 — the moment of truth:**
The fiber end at ~1 m must now appear as a **sharp peak at its true
position** — no longer as a dirt band at 21 cm.

**And the ongoing check:** peak width ≲ 1.5 × the window limit (the program
prints both numbers side by side).

---

## Afterward

Once S6 checks out, you can push the configuration further:

```bash
python calc_sweep_parameters.py --points 1000000 --speed 120 --dl 4
```

120 nm/s over a 120 nm span → range 3.48 m at **7.0 µm** resolution —
today's full sharpness with 8× the range. The catch: the aux then has only
3.5 points per fringe. That's above Nyquist, but with no margin —
`check_aux_interferometer.py` will tell you. So do this afterward, not
before.

---

## If something doesn't fit

| Symptom | First suspect | Test |
|---|---|---|
| Contrast collapses in places | Polarization in the aux arms | `check_aux_interferometer.py`, top-right panel |
| Monotonicity alarm | Laser jittering too much | `check_aux_interferometer.py`, check 2 — do **not** continue |
| τ_aux doesn't match the installed fiber | Start/end wavelength entered incorrectly, or sweep not fully within the buffer | Redo the buffer calculation |
| 46.4 mm peak at the wrong position | τ_aux wrong, or channels swapped | Check the Ch1–Ch4 assignment |
| Peaks twice as far away as expected | Reflection/transmission confused | see the "÷2" explanation |
| Flicker cuts off at end of file | Sweep longer than the buffer | Smaller span or sweep faster |
| Everything looks wrong, no idea | measure first, then explain | `diagnose_artifacts.py` has four ready-made tests |
