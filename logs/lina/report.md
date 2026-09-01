# Lina OFDR — work report, 2026-09-01

Full record of one measurement day on **Lina** (EXFO T200S swept laser +
coreDAQ-1, LINEAR variant): every consideration, every measurement, every
number, every dead end and every correction.

Everything produced that day lives in this folder. The rest of the repository
is unchanged, with two justified exceptions (§10).

---

## Contents

| § | Topic |
|---|---|
| 1 | The original problem, and how it turned out to be a different one |
| 2 | The TOF1550 as a wavelength marker (driver + first measurement) |
| 3 | What the capture really contains — the λ(t) calibration |
| 4 | Two capture bugs that made large captures unusable |
| 5 | Stability: repeatability and drift |
| 6 | OFDR: from a capture to a fibre length |
| 7 | The fibre-length measurements (1mA, 1mB, 3mB, 3mB+1mB, none) |
| 8 | Presentation: convention, offset, axes, noise |
| 9 | What was tested and **rejected** |
| 10 | Changes outside this folder |
| 11 | Open points |

---

## 1. The original problem, and how it turned out to be a different one

**Starting hypothesis (yours):** the EXFO sweeps and triggers coreDAQ-1; the
trigger is not synchronised, so the capture contains *the laser's settling
period, the sweep, and the return to the start wavelength*.

**What held up:** that the wavelength axis was wrong, and that the capture
contains more than the sweep.

**What did not hold up:** the cause. There is no settling period in the
capture. The EXFO does not raise its sync output until it has finished
settling. Proven by control experiment: laser parked at the start wavelength
for 5 s versus not parked at all — the marker moved by **< 1 pm** (5.9067 s
versus 5.9064 s). The dead time is at the **end**.

From that followed the actual finding: the GUI *constructed* the axis wrongly,
it did not trigger wrongly.

---

## 2. The TOF1550 as a wavelength marker

### 2.1 The idea

A narrowband filter in the path transmits only near its centre wavelength. Its
transmission peak therefore marks a **known wavelength at a measurable time**
in the capture. That makes the time → wavelength mapping directly measurable
instead of inferred from setpoints.

### 2.2 Identifying the device

The filter did not exist in software. Finding its port:

| Port | Device | Result |
|---|---|---|
| COM17 | "Thorlabs MTD Series Evaluation Board", VID 1313 | **not** the filter — replied in mK/ms (TEC controller) |
| COM11 / COM12 | STM32 CDC (VID 0483) | the two coreDAQs |
| COM14 | SiLabs CP210x | something else |
| **COM6** | FTDI (VID 0403), "USB Serial Port" | **the TOF1550** |

The first guess, COM17, was wrong — the Thorlabs VID was misleading. The
TOF1550 uses the generic FTDI VCP driver, is **not** identifiable from its USB
descriptor, and has no `*IDN?`. Its port is therefore pinned in configuration.

### 2.3 Protocol (User Guide TTN356097-D02 §4.2) and two deviations from it

115200 8N1, abbreviated SCPI: `WAVElength?`, `WAVElength X` (nm),
`WAVElength:STep` (pm), `:INcrease`, `:DEcrease`, `CHANnel`.

Found on the unit itself, differing from the manual:

1. Replies are **LF**-terminated, not CR as documented.
2. Replies are prose with a prefix: `OK. The center wavelength is 1550.000 nm.`
   → the driver extracts the first number from the reply.

### 2.4 Driver test

| Set | Read back | Error |
|---|---|---|
| 1545.000 nm | 1545.000 nm | 0 pm |
| 1550.500 nm | 1550.500 nm | 0 pm |
| 1560.000 nm | 1560.000 nm | 0 pm |

Step/increase/decrease also verified (100 pm step → 1560.0 → 1560.1 → 1560.0).

### 2.5 Passband measured in CW (ground truth)

Laser stepped in CW in 25 pm steps over 1.5 nm, all four channels:

| Channel | Peak | FWHM | max | baseline | contrast |
|---|---|---|---|---|---|
| Ch1 | 1549.9846 nm | 0.225 nm | 0.70 µW | 1.3 nW | 549× |
| Ch2 | 1549.9981 nm | 0.150 nm | 154 µW | 0.30 µW | 517× |
| Ch3 | 1549.9860 nm | 0.175 nm | 0.60 µW | 0.5 nW | 1173× |
| Ch4 | 1549.9643 nm | 0.175 nm | 151 µW | 0.34 µW | 447× |

FWHM 0.15–0.23 nm against a datasheet 0.21 nm at −3 dB: consistent. Filter set
to 1550.000 nm versus a true peak of ~1549.99 nm → **absolute error ~15 pm**,
inside specification (< 32 pm). So the filter is usable as a marker.

**Side finding:** both strong channels show a **notch right on the peak** —
probably coreDAQ autorange switching at full scale. Irrelevant for the
calibration (the centroid is robust), possibly not for spectra. Unresolved
(§11).

**Correction to myself:** the first channel selection ("best channel Ch3") was
wrong because it chose by *contrast* — an unconnected channel has enormous
contrast on nanowatt noise. Changed to absolute peak power.

---

## 3. What the capture really contains — the λ(t) calibration

### 3.1 Raw capture, one sweep, filter at 1550 nm

1520 → 1570 nm at 5 nm/s, 12.00 s armed window, 10 kHz:

- **one** peak at **5.9065 s**, width 59 ms
- a **second, narrow** peak at **10.1195 s**, width 2.8 ms

Interpretation: the broad peak is the sweep crossing the passband (0.18 nm at
5 nm/s = 36 ms); the narrow one is the laser's **return slew** crossing the
same wavelength 20× faster. This was the first direct sighting of the return
inside the buffer.

The naive linear model puts 1549.986 nm at 5.9972 s → **−91 ms = −454 pm**.

### 3.2 Eight filter positions → λ(t)

1530 … 1565 nm in 5 nm steps, 5 kHz, `--trust-filter`:

| Filter | Time in the buffer | Peak width |
|---|---|---|
| 1530 nm | 1.9622 s | 60.5 ms |
| 1535 nm | 2.9455 s | 61.4 ms |
| 1540 nm | 3.9314 s | 60.4 ms |
| 1545 nm | 4.9185 s | 58.4 ms |
| 1550 nm | 5.9064 s | 60.8 ms |
| 1555 nm | 6.8990 s | 60.2 ms |
| 1560 nm | 7.8934 s | 61.6 ms |
| 1565 nm | 8.8909 s | 60.4 ms |

### 3.3 Fit degree: why quadratic

| Degree | Residual |
|---|---|
| 1 (linear) | **27.2 pm rms**, 43.9 pm max |
| **2** | **3.1 pm rms**, 5.1 pm max |
| 3 | 1.9 pm rms, 3.9 pm max |
| 4 | 1.9 pm rms, 3.7 pm max |

Degree 1 leaves a **systematic arch** (−39, −6, +13, +26, **+35**, +20, −4,
−44 pm) — not noise, real curvature. Degree 2 removes it, degree 3 adds
nothing. So **degree 2**.

### 3.4 What the fit says about the laser

| Quantity | Value |
|---|---|
| effective sweep rate (t=0) | **5.1174 nm/s** versus 5.0000 commanded (+2.3 %) |
| trigger dead time | **+0.0041 s** — essentially zero |
| sweep ends at | **9.8895 s** of 12.00 s |
| sweep occupies | **82.4 %** of the armed buffer |

So the sweep starts at the trigger edge and finishes 0.13 s before its nominal
end because it runs too fast. The remaining 2.1 s are tail, including the
return slew.

### 3.5 What that cost the old evaluation

Same 100 kHz capture, marker demonstrably at 1550.000 nm:

| Axis construction | Marker appears at | Error |
|---|---|---|
| linspace over the whole buffer | 1544.612 nm | **−5.4 nm** |
| + `Buffer time` 1 s | 1542.304 nm | −7.7 nm |
| + `Buffer time` 2 s (your default) | 1539.534 nm | **−10.5 nm** |
| + `Buffer time` 3 s | 1536.149 nm | −13.9 nm |
| **measured λ(t), degree 2** | **1549.997 nm** | **−3 pm** |

Three separate error sources:

1. **Axis stretched over the whole buffer** instead of the 9.89 s sweep —
   dominant, −5.4 nm.
2. **`Buffer time` crops at the front** — corrects the wrong end, deletes real
   1520–1530 nm data and makes the error monotonically worse.
3. **Commanded versus actual sweep speed** — 1–2 %, ~0.5 nm over 50 nm,
   remains even after 1 and 2 are fixed.

### 3.6 Rate independence (the decisive test)

The fit is in **seconds**, not samples. Calibration measured at **5 kHz**,
applied to a **100 kHz** capture of 1.2 M samples: marker at **1549.9970 nm**,
i.e. **−3 pm**. So: calibrate once at a low, reliably readable rate and apply
it at the production rate.

---

## 4. Two capture bugs that made large captures unusable

### 4.1 Read too early (deterministic, this was the main cause)

A LINEAR capture free-runs `n_frames / rate` seconds **from the trigger edge**
— independent of how long `start_continuous_sweep()` blocks. With 1.2 M frames
at 100 kHz = **12.00 s** window against a **10.81 s** sweep, the buffer was
still filling when it was read. The USB transfer then stalls (the remaining
frames do not exist yet), the transport's 30 s budget expires, and
py_coreDAQ's reset-retry **destroys the arming** → the follow-on error "no
capture was armed by this session".

Observed: stalls at 4,234,228 / 4,800,000 bytes (88 %) and at
424,308 / 480,000 (88 %).

**Fix:** wait out the armed window before reading. After that, 1.05 M and
1.2 M frames read cleanly.

### 4.2 Sporadic packet loss (statistical, still present)

Even with a full buffer the read occasionally dies, and always exactly
**512 bytes** short of the end — the USB packet size. Exactly one packet is
lost.

**Rejected hypothesis:** that this is an alignment effect (byte count an exact
multiple of 512). I declared it confirmed after *one* success and *one*
failure — too early. The control experiment refuted it:

| Frames | Bytes | % 512 | Result |
|---|---|---|---|
| 1,200,000 × 2 ch | 4,800,000 | 0 | stalled at 4,799,488 |
| 1,199,999 × 2 ch | 4,799,996 | 508 | first **succeeded**, two minutes later **stalled** at 4,799,484 |

Same frame count, both outcomes. So **no** alignment effect; the workaround
that had been added was removed again.

Measured failure rate:

| Transfer | Channels × frames | Failures |
|---|---|---|
| 4.8 MB | 2 × 1.2 M | 3 of 6 attempts |
| 8.4 MB | 4 × 1.05 M | both runs affected, one exhausted 2 retries |

**Fix:** the whole park → arm → sweep → read sequence is repeated as a unit (a
failed read consumes the arming, so re-reading is impossible). Budget: 4
attempts. Plus a 5 s wait before retrying: after an aborted sweep the EXFO
stays `busy` for a few seconds and rejects the trigger reconfiguration with
`Metrino.T200S.BusyException` — with a 1 s wait the retry was burned on the
laser instead of on the failure being retried.

Verified: 3 runs of 1.2 M × 2 ch → one failed twice and succeeded on the third
attempt, one failed once, one succeeded immediately — **all three ended with
usable data**.

### 4.3 Side observation

The coreDAQ can wedge such that the port no longer opens
(`PermissionError: Access is denied`) and only a physical replug helps. I first
attributed this to two leftover Python processes — wrong: they were still
running long after the port worked again.

---

## 5. Stability: repeatability and drift

Filter fixed at 1550 nm, identical sweeps:

| Test | Runs | Span | Scatter | Drift |
|---|---|---|---|---|
| short | 5 | 1.7 min | 0.8 pm rms, 2.0 pm p-p | −0.9 pm/min (not resolvable) |
| long | 6 | 18.0 min | **0.7 pm rms**, 1.7 pm p-p | **+0.0 pm/min** |

The scatter does **not** grow over the 10× longer span → white noise at the
level of the fit residuals, no time-dependent component.

**Consequence: one calibration per sweep configuration is enough.** No
re-calibration during a session.

---

## 6. OFDR: from a capture to a fibre length

### 6.1 Setup

- **Ch1, Ch3** — measurement MZI, contains the test fibre
- **Ch2, Ch4** — reference MZI ("aux"), ΔL nominally 4 m

Each MZI is read on two channels. Whether a pair is balanced is decided from
its **correlation**, not assumed: measured −0.88 to −0.99 → balanced
difference.

### 6.2 Why the aux MZI is needed

The coreDAQ samples uniformly in **time**; the FFT needs uniformity in
**optical frequency**. The aux MZI has a fixed delay, so its fringe phase is
strictly proportional to optical frequency: `φ(t) = 2π·τ_aux·ν(t)`. Resampling
the measurement onto a uniform φ grid makes it uniform in ν.

**Independent confirmation of the calibration:** the unwrapped aux phase
deviates from a straight line by **1.08 % of the frequency span** — as a
smooth arch, i.e. the same quadratic nonlinearity the filter fit found as its
degree-2 term. Two entirely different methods, same result.

### 6.3 Scale

```
dz    = c / (factor · n_g · Δν)          resolution
z_max = dz · N/2 = c·f_s / (2·factor·n_g·dν/dt)
```

With 1520–1570 nm: Δν = 6.2813 THz, n_g = 1.4682 → dz = 16.3 µm (reflection
convention), Nyquist 8.03 m. Δν comes from the **calibrated** wavelength axis,
not the commanded range — otherwise the 1–2 % speed error would go straight
into every length.

### 6.4 Validation on synthetic data

Known 1.000 m path imbalance, 5 m aux, deliberately nonlinear sweep, noise:
**1.00000 m** recovered (error < 0.01 mm), plus a planted −28 dB secondary
reflection at 2.5 m found. Implied aux length 5.048 m against 5.000 m true.

### 6.5 High-pass before the phase extraction

Balanced subtraction does not fully remove the power envelope (unequal
detector responsivity). The remainder dominates the spectrum at DC: the aux
"fringe frequency" was estimated as **0.000 kHz**, i.e. useless. With a 500 Hz
high-pass it jumped to 11.1 kHz. Fringes are at kHz, the envelope at Hz — they
separate cleanly.

---

## 7. The fibre-length measurements

### 7.1 First assignment (before the signs were known)

Both channel pairs, each using the *other* as reference — two independent
measurements that check each other:

| Pair | Peak | Expected from that pair's fringe frequency |
|---|---|---|
| Ch1/Ch3 | 1.14795 m | 1.148 m (3.57 kHz) ✓ |
| Ch2/Ch4 | 4.16179 m | 4.16 m (12.95 kHz) ✓ |

Repeatability over 3 sweeps: 1.14795 / 1.14800 / 1.14797 m → **±0.03 mm**.

### 7.2 The fibre swap as an experiment

1mB → 1mA swapped:

| Pair | Peak | 1mB | 1mA | Shift |
|---|---|---|---|---|
| meas | main | 1.14797 | 1.14793 | −0.05 mm |
| meas | | 1.01266 | 1.01273 | +0.07 mm |
| meas | | 1.04986 | 1.04992 | +0.07 mm |
| meas | | **3.22729** | **3.25856** | **+31.26 mm** |
| aux | | **2.08219** | **2.05117** | **−31.02 mm** |
| aux | main | 4.16197 | 4.16196 | −0.01 mm |

Two peaks move by the same amount with **opposite sign**. The strong peaks do
**not** move — so they are not the fibre.

### 7.3 Measurement with no fibre (this settled the assignment)

| Peak | 1mB | 1mA | no fibre | Interpretation |
|---|---|---|---|---|
| 1.0127 m | −24.9 | −25.6 | −31.0 | fibre-independent |
| 1.0499 m | −24.8 | −26.2 | −31.2 | fibre-independent |
| **1.1480 m** | 0.0 | 0.0 | 0.0 | measurement MZI's own ΔL |
| 3.2273 / 3.2586 m | −22.4 | −18.4 | **gone** | **only with a fibre** |
| 2.0822 / 2.0512 m (aux) | −10.0 | −7.9 | **gone** | **only with a fibre** |
| **4.1620 m** (aux) | 0.0 | 0.0 | 0.0 | aux MZI's own ΔL |

Exactly two peaks vanish without a fibre. The strongest peak is **not** the
fibre — that was the most important misinterpretation this measurement ruled
out.

### 7.4 Model and its self-check

```
measurement pair peak = 1.1480 m + L        (test fibre traversed TWICE)
aux pair peak         = |4.1620 m − L|      (aux traversed ONCE)
```

| | L from the meas pair | aux peak predicted | aux peak measured | Deviation |
|---|---|---|---|---|
| 1mB | 2.0793 m | 2.08265 m | 2.08219 m | 0.46 mm |
| 1mA | 2.1106 m | 2.05133 m | 2.05117 m | 0.16 mm |

The cross-sum must be constant (1.1480 + 4.1620 = 5.3099 m): measured 5.30948
(1mB) and 5.30973 (1mA) → both within **0.5 mm**, with one free parameter.

The factor of 2 was confirmed by the 3 m fibre: 1 m → 2.079 m of path,
3 m → 6.051 m. And the aux ΔL of nominally 4 m reads **4.162 m** with factor
n_g → traversed once. (If it were exactly 4.000 m, n_g would have to be 1.528
— unrealistic for SMF-28, so it really is ~4.16 m.)

### 7.5 Final results (axis reads fibre length)

| Fibre | Peak | Prediction | dB |
|---|---|---|---|
| **1mB** | **1.0397 m** | — | −22.4 |
| **1mA** | **1.0553 m** | — | −18.4 |
| **3mB** | **3.0256 m** | — | −31.2 |
| **3mB + 1mB** | **4.06334 m** | 3.02564 + 1.03969 = 4.06533 | −31.9 |
| no fibre | no peak | — | — |

**The lengths add up to within 2.0 mm (0.05 %).** That confirms the scale
independently from 1 m to 4 m with no re-adjustment. The prediction made before
the measurement was 4.065 m — hit to 1.9 mm.

1mA is **15.6 mm longer** than 1mB.

Other real signals: at 1.51 m half the 3 m fibre (multiple reflection); at
3.03 m, in the series connection, the connector between the two fibres; a broad
bump at 2.5–2.7 m (distributed scatter, not a point reflection — its median is
raised too).

---

## 8. Presentation

### 8.1 Convention (changed twice, as the physics became clear)

| Stage | Factor | A 1 m fibre appears at |
|---|---|---|
| initially | 2·n_g | 0.5 m |
| after your request "1 m should be at 1 m" | n_g | 1.0 m — but only coincidentally right |
| **final, once the double pass was proven** | 2·n_g **+ offset subtraction** | 1.0397 m, physically correct |

The two panels need **different** conventions, because the test fibre is
traversed twice and the aux delay line once.

### 8.2 Offset subtraction, and why the aux panel gets no fibre-length axis

`--offset-m 0.57396` subtracts the measurement MZI's own ΔL → the axis reads
fibre length. The equivalent for the aux panel would be a **mirror** about
2.08098 m — but that **folds** as soon as the fibre exceeds half the aux ΔL:
for the 3 m fibre, \|4.162 − 6.051\| is negative, the peak folds, and the
mirrored axis wrongly showed 1.136 m instead of 3.026 m. The aux panel
therefore keeps its own axis.

### 8.3 Noise display

The problem was not the scaling (dB is already logarithmic) but that ~2 M FFT
bins are drawn and the noise fills in as a solid block. Solution: **block-maximum
decimation** ("peak hold", 3000 points) — peaks keep their exact height (a peak
*is* the maximum of its block) and the noise collapses to a thin line. Plus the
block median as a grey area = the actual noise floor.

### 8.4 Axis windows

`auto` frames the detected peaks (±0.2 m, minimum 0.6 m span) instead of the
whole Nyquist range.

---

## 9. What was tested and rejected

This section is deliberately complete — rejected ideas are results.

| Idea | Test | Outcome |
|---|---|---|
| The trigger captures the settling phase | 5 s park versus no park | **refuted**, < 1 pm difference |
| `Buffer time` (front crop) fixes the axis | marker of known wavelength | **refuted**, makes it monotonically worse (−5.4 → −13.9 nm) |
| 512-byte alignment causes the read stall | non-aligned byte count | **refuted**, fails identically; workaround removed |
| Leftover Python processes blocked COM12 | processes still running while the port worked | **refuted**, the device had wedged |
| Pick "best channel" by contrast | an unconnected channel wins | **wrong**, changed to peak power |
| The first 10 nm are bad and must be excluded | fringe amplitude per 0.5 nm; crop scan 0/2/5/10/15/20 nm | **no effect** (below) — the option was removed again |
| Degree 3 or 4 for λ(t) | residuals | **no gain** over degree 2 |

### 9.1 The crop scan in detail

| Crop | Samples | dz | Peak | Shift | SNR | FWHM |
|---|---|---|---|---|---|---|
| 0 nm | 988,536 | 16.25 µm | 4.06334 m | — | 41.8 dB | 0.05 mm |
| 2 nm | 949,435 | 16.95 µm | 4.06415 m | +0.81 mm | 41.9 dB | 0.06 mm |
| 5 nm | 890,717 | 18.12 µm | 4.06435 m | +1.01 mm | 41.7 dB | 0.05 mm |
| 10 nm | 792,672 | 20.45 µm | 4.06520 m | +1.86 mm | 40.6 dB | 0.05 mm |
| 15 nm | 694,401 | 23.45 µm | 4.06530 m | +1.96 mm | 41.3 dB | 0.09 mm |
| 20 nm | 595,900 | 27.45 µm | 4.06591 m | +2.57 mm | 42.2 dB | 0.08 mm |

SNR stays at 41–42 dB, and FWHM gets worse at heavy crops (dz grows as Δν
shrinks). The fringe amplitude is constant across the whole sweep (0.78–1.05 ×
median in the aux, 0.96–1.03 in the measurement channel, nothing below 60 %).

**One side finding remains:** the peak moves systematically by up to +2.6 mm.
That is not a length change but the peak's own structure (2–3 sub-lobes over
~50 mm) with the maximum jumping between them. So the **systematic**
uncertainty is ~±2 mm while the **repeatability** is ±0.03 mm. Irrelevant for
comparisons, relevant for absolute values.

---

## 10. Changes outside this folder

Everything new lives in `lina/`. Two existing files had to change because you
wanted the functionality in the GUI:

### 10.1 `LabGUI/dashboards/instruments/lina_window.py`

| Change | Reason |
|---|---|
| "Use λ(t) calibration" (default on) + status line | shows **before** a sweep whether a calibration exists for Start/Stop/Speed |
| λ(t) axis + tail cropping in `_CoreDAQSweepWorker` | §3 |
| **`Buffer time` removed entirely** (field, front crop, arming arithmetic) | §3.5 — its premise was refuted |
| `Settle margin` kept, default 10 → **0 nm** | a non-zero value would silently disable the calibration; its own justification (saturation transient) is **not** refuted, only untested |
| wait for the armed window | §4.1 |
| capture retry (4 attempts, 5 s pause) | §4.2 |
| Samples/ch default 1,000,000 → 1,050,000, auto-fit on parameter change | matches sweep + 5 % head room, no more ⚠ warning |
| all four coreDAQ channels by default | OFDR needs both MZI pairs from the **same** sweep |
| "restore λ" default on | matches the scripted measurement runs |
| button **"OFDR Plot"** + n_g / MZI-offset fields | the two spectra from the last capture |

### 10.2 `Interface/PowerMeters.py`

| Change | Reason |
|---|---|
| `arm_sweep_capture` returns the armed frame count | the GUI needs the value actually armed |
| `read_sweep_capture` clamps to the armed count | avoids a misleading validation error on FW v4.2, where `captured_frames()` returns −1 |
| docstring documents the 512-byte packet loss | §4.2 |

**Reverted:** `.gitignore` and `CAROLINA_OFDR_CHANGES.md` are back to their
pre-today state (`519e18f`).

**The TOF1550 driver lives in `lina/interface/`**, not `Interface/`, so that
`Interface/` stays unchanged. If the filter should become a regular GUI
instrument type, it belongs in `Interface/` and in the device registry — its
configuration would then move to `Interface/config/` as well.

---

## 11. Open points

1. **The GUI button does not yet give the right length.** The evaluation is
   demonstrably identical to the script (the same `evaluate()` function,
   verified on script data), but captures taken through the GUI have an aux
   phase with ~5 % too few fringes (121,636 instead of 128,011), which
   stretches the distance axis by the same factor: 2.846 m instead of 3.026 m.
   Cause not yet found. **The script is the reliable path.**
2. **Autorange notch** on strong peaks (§2.5) — probably TIA range switching.
   Testable with a fixed range instead of autorange.
3. **`Settle margin`** — the saturation in the first ~10 nm was never
   cross-checked, because with the filter in the path the start of the sweep is
   dark. A sweep without the filter, looking for clipping in the first few nm,
   would settle it.
4. **Absolute scale** depends linearly on n_g = 1.4682. The differences are
   sound; absolute values would need a fibre of known length or an OSA
   comparison. Calibrating with `--trust-filter` also means: relative accuracy
   ~3–5 pm, absolute ±32 pm (the filter's specification).
5. **The calibration is valid only for 1520–1570 nm at 5 nm/s.** Any other
   configuration needs a ~3-minute run; the lookup deliberately refuses a
   mismatched calibration.
6. **Large reads remain statistically unreliable** (§4.2). At four channels,
   8.4 MB often fails; the retry absorbs it but costs one sweep each time.
