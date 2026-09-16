# The OFDR reflectometer — how I built it, what I measured, what I still don't know

Carolina Krotsch · 16 September 2026

This is a general write-up for the company. First how the setup works, then how a raw scan turns into a length, then everything I have measured so far (fibres and chips) with the results, and at the end the list of things I don't understand and would like help with.

---

## 1. What the instrument is

It's an **optical frequency-domain reflectometer (OFDR)** — basically a home-built version of a Luna OBR. You sweep a tunable laser in wavelength, mix the light coming back from the device under test with a reference copy of the same light, record the interference, and Fourier-transform it. Every reflecting interface along the fibre then shows up as one peak, at its own distance.

Why that works: with the laser at optical frequency ν, a reflector sitting at round-trip delay τᵢ gives

```
I(ν) = Σᵢ 2·√(P_LO·Pᵢ)·cos(2π·τᵢ·ν + φᵢ)
```

So if you look at the signal as a function of ν, it's just a sum of sine waves — one per reflector — and the frequency of each sine wave *is* its delay. That's why an FFT over ν gives you the reflectogram directly.

Two numbers describe what the instrument can do:

| quantity | formula | what I have |
| --- | --- | --- |
| resolution | Δz = c / (2·n_g·Δν) | **16.3 µm** at 1520–1570 nm, **6.9 µm** at 1505–1625 nm |
| max. distance | z_max = c / (4·n_g·δν) | **8.0 m** at 1520–1570 nm, **3.4 m** at 1505–1625 nm |

with n_g = 1.468 for SMF-28 at 1550 nm. Note that resolution and range trade off against each other: range/resolution = N/2, where N is the number of points. So a wider sweep gives better resolution but less range. How *fast* you sweep doesn't matter for either.

### 1.1 The hardware

- **Laser:** EXFO tunable, 1505–1625 nm, swept at 10 nm/s.
- **Detector:** CoreDAQ multichannel power meter (dev unit), 4 channels, running free on its own internal clock. A normal scan is **988,536 points** over 1520–1570 nm, which is about 150 MB of JSON.
- **Measurement interferometer:** built around a circulator, so I reach the device in reflection. The device under test hangs on a ~1 m patchcord.
- **Auxiliary interferometer ("the ruler"):** a Mach–Zehnder with ΔL ≈ 4.2 m, made from two ~2.07 m fibres. This is the part that makes everything work — see section 2.
- **Which channel is which:** Ch2/Ch4 = aux (the strong pair), Ch1/Ch3 = measurement (the weak pair). The pairs are interleaved, not next to each other. I got this wrong for two weeks; it was settled by checking the wiring by hand on 31 August.
- Everything the instrument owns was rebuilt in **PM fibre** in September. The device arm stays SM because that's what the chips have. All connectors FC/APC.

### 1.2 The one decision that changed everything

At the start I sampled on the laser's own trigger output, one trigger per 1 pm. That limits the range to **0.42 m** — shorter than my own patchcord. So the far end of the test fibre was outside the range and folded back into the middle of the trace as an ugly smeared band at about −40 dB. I tested and killed four different explanations for that band before a simulation showed it was simply aliasing. Once something aliases there is no software fix, because the information is already gone at the sampling stage.

The fix was to **stop using triggers and sample free-running**, like a Luna does. That gives about 20× more range at better resolution. The catch is that without triggers there is *no wavelength information at all* in the data — and that's exactly what the auxiliary interferometer provides.

---

## 2. How a raw JSON scan becomes a length

This is the part worth understanding, because every number further down comes out of it. The script is `process_reflectogram_aux.py`.

**What I start with.** A JSON file with four arrays in mW, one sample per clock tick, and no wavelength axis at all:

```json
{"data": [{"Ch1 [mW]": [...], "Ch2 [mW]": [...],
           "Ch3 [mW]": [...], "Ch4 [mW]": [...]}]}
```

**Step 1 — balanced subtraction.** The two channels of each pair are the two complementary outputs of the same coupler, so P = Ch_a − g·Ch_b cancels the DC offset and the common-mode intensity noise. g is fitted piecewise because the two photodiodes don't have the same responsivity slope — it drifts from about 0.84 to 1.14 across the sweep. Done separately for the aux pair and the measurement pair.

**Step 2 — get the frequency axis back from the aux.** The aux is a fixed delay τ_aux, so its phase is

```
φ(t) = 2π · τ_aux · ν(t) + const
```

I take the Hilbert transform of the balanced aux signal and unwrap the phase. Then: **equal steps in φ are equal steps in ν.** It doesn't matter how uneven or jittery the laser sweep actually was — the aux measures it. That's the whole trick.

**Step 3 — resample.** I interpolate the measurement signal (PCHIP) onto a uniform φ grid. Now it's evenly spaced in optical frequency and the FFT is allowed. Important: I do **not** smooth the aux phase. It has to follow the laser's fast jitter, otherwise it can't correct it.

**Step 4 — the scale factor.** τ_aux is what converts phase into frequency. I get it by counting fringes over a sweep whose span I know:

```
τ_aux = (total unwrapped aux phase) / (2π · frequency span of the sweep)
```

Right now it's about **20.61 ns** (ΔL ≈ 4.21 m). It has to be redone every time the aux arm gets re-plugged. From it:

```
δν      = Δφ / (2π · τ_aux)         spacing in optical frequency
span_ν  = Δφ_total / (2π · τ_aux)   total span
```

**Step 5 — window and FFT.** Subtract a 5th-order polynomial baseline, apply a Kaiser β = 12 window, then a real FFT. Kaiser because its sidelobes are below −80 dB. Hann gives a sharper peak, but its −31 dB sidelobes look like reflectors when they sit next to a strong peak, and I'd rather not fool myself.

**Step 6 — the distance axis.** FFT bin k sits at

```
z_k = k · c / (2 · n_g · δν · m)
```

with m the number of resampled points. I use the **reflection convention**, z = c·τ/(2·n_g), so z is the one-way distance to the reflector, not the round trip.

**What comes out** is a peak list in mm with amplitudes in dB relative to the strongest peak, plus a CSV and PNG of the whole reflectogram. A worked example for a standard scan:

| quantity | value |
| --- | --- |
| N = 988,536 points, 1520.000 → 1570.000 nm | mean step 0.0506 pm |
| δν | 6.08 … 6.64 MHz (9 % chirp over the sweep) |
| Δν | 6.281 THz |
| resolution Δz = c/(2 n_g Δν) | **16.26 µm** per bin |
| …times the Kaiser factor 2.6 | **≈ 42 µm** actual −3 dB width |
| Nyquist range | **8.03 m** |
| dynamic range I actually get | 63 … 72 dB |

Two checks run automatically on every scan, and I quote them whenever someone doubts a number: the **−3 dB width of the main peak** compared to the window limit (if the frequency axis went wrong, this is where you see it first), and the **minimum aux contrast** (if the aux fades due to polarisation, its phase — and therefore every position — is unreliable).

**For lengths on a chip** there's one more step. The instrument measures delay and reports it in *fibre-equivalent* mm. A waveguide of length s on the chip appears at

```
dz = s · n_g,chip / n_g,fibre
```

So to get s I need the chip's group index. Not having that number is the main problem in all the chip work (sections 4.2 and 5).

---

## 3. Results — fibres

### 3.1 Do the lengths add up?

> 📈 **PLOT — insert `fiber_length_series_august.png` here**

Each trace is one scan, all shifted so their own connector peak sits at x = 0. Grey = nothing connected (control), red = a ~1 m fibre, blue = a ~3 m fibre, purple = the 3 m and the 1 m fibre connected in series.

What you can see:

- Each fibre gives **one clean end reflection**, about 30 dB above everything around it, and the control trace has nothing there.
- The four different "1 m" fibres land at 1002–1019 mm, so they really are different fibres and the instrument can tell them apart.
- The series combination lands at **3939.3 mm**, against 2923.9 + 1015.7 = 3939.6 mm predicted from the same two fibres measured separately. So the instrument **adds up correctly**.

One caveat: this series is from the end of August, and back then the fringe counting for τ_aux had a scale error of roughly 10–15 % that I never fully explained. This plot is the version where the aux length was set to its nominal value instead, which is why the numbers look so round. The *additivity* result doesn't depend on the scale being right, which is why I still like the plot. The current setup is verified much better — next section.

### 3.2 Three 2.07 m PM fibres — the best accuracy test I have

Three nominally identical 2.07 m PM fibres (a, b, c). Two of them always sit in the aux and the third is the device under test, so measuring all three means swapping them around. That turns into a free consistency check, because the aux length enters differently each time.

| fibre | length from the connector plane | vs. 2070 mm from the tape measure |
| --- | --- | --- |
| a | **2085.93 mm** | +15.93 mm |
| b | **2093.89 mm** | +23.89 mm |
| c | **2095.35 mm** | +25.35 mm |

The extra couple of cm is the connector bodies, which you don't include when you lay the bare fibre next to a tape measure.

Because I have three scans for fewer unknowns, two things fall out:

1. **The aux has 23.7 ± 0.1 mm of extra fixed path** on top of the two fibres in it — coupler pigtails and connector bodies. Three independent estimates, one per scan, agree to 0.19 mm. I had assumed the aux was just the two fibres, and that was wrong by 24 mm.
2. **The pair test closes to one resolution cell.** If the aux holds the two fibres that are *not* being measured, then `aux_a − aux_b` has to equal `L(b) − L(a)`. The left side comes from fringe counting, the right side from peak positions, and nothing in the analysis links them. Measured: **+7.943 mm vs +7.960 mm — 17 µm apart.**

That 17 µm is my honest accuracy statement for the instrument on a 2 m path: about **8 ppm**, and limited by the resolution cell rather than by the calibration.

### 3.3 Something I can't explain: every fibre end is a doublet

All three fibre ends show up as **two peaks**, not one, separated by 5.26 / 3.18 / 3.76 mm. The separation is different for each fibre, so it's not a fixed instrument artifact. And the stronger of the two isn't always the end: for fibres a and b the first peak is stronger, for c the second one. So "just take the strongest peak" gives the wrong length for fibre c.

I'm confident the *first* peak is the fibre end, because that assignment makes the three estimates of the aux extra path agree to 0.19 mm while the other one leaves 3.58 mm of scatter. But what the second peak actually *is*, I don't know. See section 5.

---

## 4. Results — chips

### 4.1 Everything on one axis

> 📈 **PLOT — insert `overview_all_scans_from_connector_plane.png` here**

One row per scan, all on the same distance axis with zero at the connector plane, each row normalised to its own strongest peak. The grey underneath each row is the control scan with nothing connected, so anything sticking out above the grey is real.

From top to bottom:

- **Rows 1–3, blue — the Ligentec SiN chip through its fibre array.** One big peak at ~1023–1026 mm, which is the chip input facet. Behind it, nothing above the noise for the next 1.5 m.
- **Rows 4–6, orange — the three bare 2.07 m PM fibres.** A single sharp end reflection at ~2086–2095 mm and a flat floor everywhere else. This is what a clean measurement looks like.
- **Rows 7–14, green — bare fibre arrays with no chip on them** (Meisu and PHIX). Array facet at ~1067–1076 mm, again one peak and nothing behind it.

The reason for putting them all on one axis is the comparison between the blue rows and the green rows: **the chip scans look the same as the scans with no chip at all.** Whatever the chip adds, I can't see it.

The small peak at ~1550 mm in every row is the aux ghost. It's an artifact of the aux itself, sits at ΔL/2, and is there in every scan including the ones with nothing connected. It moves when the aux gets re-plugged while real reflectors don't — that's how I identify it.

### 4.2 The HHI InP chip — the one chip where I got real geometry out

This is the chip measurement that worked. The structure is a **loopback**: two neighbouring edge couplers connected to each other by a short waveguide loop on the chip, and both are permanently pigtailed with ~1 m fibres.

> 📈 **PLOT — insert `hhi_loopback_chain_single_scan.png` here**

Row 1 is what's physically there: reflectometer → fibre 1 → through the chip → fibre 2 → open fibre end, and that open end is what reflects back. Row 2 is the same thing as measured — four peaks in a row, with the coloured stretches in between being the two fibres. Row 3 zooms in on the chip, which is **7.2 mm out of a 2.2 m path**.

| station | z | section in front of it |
| --- | --- | --- |
| internal reference reflection | 573.93 mm | — |
| chip facet A | 1670.23 mm | fibre 1 = **1096.30 mm** |
| chip facet B | 1677.43 mm | **chip = 7.200 mm** |
| open end of fibre 2 | 2773.62 mm | fibre 2 = **1096.20 mm** |

The bit I like most: fibre 1 and fibre 2 are found **completely independently of each other** and come out 1096.298 and 1096.199 mm, so the same to 0.099 mm. Nothing in the analysis forces that, so it tells me all four stations are correctly identified.

> 📈 **PLOT — insert `hhi_loops_measured_explained.png` here**

Same measurement, but for **both** loopback pairs on the chip (design lengths 640.1 µm and 556.5 µm). Row 3 has the design model drawn on top, row 4 shows how far each measured peak lands from where the design says it should be.

Row 3 also carries an important warning: **the loop edges are not in the data.** The waveguide goes smoothly into the spot-size converter, so its ends don't reflect. The only two things the instrument knows about this chip are facet A and facet B — the hatched bars are model, not measurement. And the dips in the curve are window sidelobes plus noise, not structure.

**What this gives me, and what it doesn't.** The reflectometer measures delay, so to turn 7.200 mm into a length I need the chip's group index. It's in none of the documents I have — not the Jeppix validation report (that one has waveguide loss 1.8 dB/cm, butt-joint loss 1.2 dB, SSC coupling loss 1.5 dB typ., but no index), not the device mapping, not the pin mapping — and I don't have the HHI PDK.

So I did the obvious thing and fitted n_g so that the two known design lengths come out right: **n_g = 3.4803 ± 0.008**. I want to flag clearly that **this is circular.** If I fit n_g so that both loops come out right, then of course both loops come out right. I nearly reported it as a result.

What is genuinely *not* circular:

- **The ratio test.** Two loops, one n_g, so one degree of freedom is left over. The measured ratio of the two chip paths is 1.026008 and the design ratio is 1.028268 — they agree to **0.22 %**, against a resolution limit of 0.47 % on that ratio. If the geometry or my peak identification were wrong, this wouldn't close.
- **Both loops reproduce their design path** to within one and two resolution cells (one cell = 7.0 µm on the chip). That's row 4 of the plot.

So what I can defend is: *the measurement agrees with the design geometry at the resolution limit, and the instrument resolves about 7 µm of on-chip waveguide.* What I can't defend is an absolute length, because that needs n_g from outside.

I also tried to get n_g from the single-ended channels instead. For that I built an optical netlist out of the GDS — every waveguide polygon contributes its centre-line length, every component contributes the distance between its own ports, and Dijkstra from a coupler gives the distance to every interface. It found what looked like a clean answer: three predicted interfaces matching within 0.7 of a resolution cell and giving the same n_g to 0.08 %.

**It doesn't hold, and I want to be explicit about why.** The netlist has 690 interfaces across 29 channels, so at a tolerance of two resolution cells there is a candidate near almost any position you predict. Those matches were *selected* by the n_g I had assumed, so they can't then be evidence for it. The test that settles it: sweep n_g from 3.0 to 4.0 and see how well the best assignment fits. If the data contained the information there would be a clear maximum near 3.48. Instead the curve is flat across the whole range, with its maximum at n_g = 3.02.

Two side results from that chip. The two coupler pairs couple very differently — the bottom pair returns 13.5 dB more than the top one, which puts the top pair at about 4.9 dB per transit against a spec of 1.5 dB typical / 4 dB max. I think that's the fibre-array glueing rather than the SSCs. And I had no map of which fibre goes to which coupler, so I worked it out from the measurements: fibre n goes to coupler b(32−n). I confirmed it twice, once on a channel that should run into the circuit (it does, with a whole forest of peaks) and once on a channel that should be a dead-end stub (it is — one sharp reflector and nothing behind it).

### 4.3 The Ligentec SiN chip — the facet, and nothing else

The current chip is two chips in one die: a **Ligentec SiN** part with a 95-channel fibre array along its bottom edge (127 µm pitch, lensed edge couplers), and an **HHI InP** part above it, joined by a 132-port chip-to-chip interface at 92 µm pitch. Chip channel number = 96 − patchcord fibre number. I've measured two of them: **MAP2672**, and since today **MAP2680**, which has metallisation (bondpads and traces) where MAP2672 has none.

> 📈 **PLOT — insert `ligentec_facet_zoom_four_channels.png` here**

**Top panel:** where the facets of the four measured channels sit, counted from the connector plane. They spread over 2.40 mm between channel 52 and channel 95 — and that spread is the **length difference of the array fibres**, not anything about the chip.

**Bottom panel — this is the strongest single piece of evidence I have.** Same four channels, but each one re-centred on its own facet, so the x axis means "distance behind your own facet". These four channels are physically very different: channel 95 is a pure 199.5 µm loopback with no components at all, channel 52 runs into a 1×2 MMI at 8.9 mm, channel 48 has a structure at 6.5 mm. If those components reflected anything, the curves would have to differ exactly there.

**They lie on top of each other**, bump for bump, including the features at +7.5 mm and +11.5 mm. So what I see behind the facet is the instrument, not the chip.

> 📈 **PLOT — insert `ligentec2680_channel52_component_search.png` here**

The same conclusion, but done as a statistical test, on the new chip MAP2680, channel 52, with both sweep spans. The bottom two panels have the design positions of the real components marked as grey bands (wide because they span n_g,chip = 1.80–2.00), taken from the GDS: edge coupler at 0.5 mm, an X1 crossing at 3.5–3.9 mm, a pitch split and a 1×2 MMI at 10.8–12.1 mm behind facet A.

The test: take the maximum inside each predicted window, and compare it against equally wide control windows slid across 1–40 mm behind the facet. If a component is really there, its window should stand out from that null distribution.

| channel | component | window | 1520–1570 nm | 1505–1625 nm |
| --- | --- | --- | --- | --- |
| 52 | X1 crossing | 3.53–3.93 mm | +31.6 dB, **7 %** | +34.5 dB, **7 %** |
| 52 | pitch split | 10.80–12.01 mm | +28.1 dB, 34 % | +31.1 dB, 24 % |
| 52 | MMI 1×2 | 10.93–12.14 mm | +28.1 dB, 28 % | +31.1 dB, 21 % |
| 48 | pitch split / MZM | 7.90–8.78 mm | +29.0 dB, 24 % | +30.8 dB, 22 % |
| 48 | MMI 1×2 | 7.93–8.81 mm | +29.0 dB, 27 % | +29.3 dB, 24 % |
| 48 | straight heater | 8.16–9.07 mm | +29.0 dB, 28 % | +29.3 dB, 28 % |
| 48 | MMI 2×2 | 8.69–9.65 mm | +25.8 dB, 32 % | +25.2 dB, 39 % |

(The percentage is the fraction of control windows that reach or beat the predicted one. Below 5 % would be interesting.)

**Not a single detection.** The peak near 3.85 mm looks like something in the plot, but at 7 % it isn't separable from the multipath forest — with 92 control windows you expect about 5 hits purely by chance. The full sweep span buys 3–4 dB more dynamic range, but the percentages don't improve, because the forest rises along with the signal.

**Why that's actually the expected answer.** The only things on a photonic chip that reflect strongly are air facets. Taper tips are adiabatic on purpose, MMI and crossing interfaces only reflect through mode mismatch (−40…−50 dB), and a heater is optically invisible. On top of that, everything except the input facet has to get back out through the fibre-chip coupling a second time.

That also explains something that confused me: channels that give a clear signal in plain transmission (laser in one side, power meter on the other) give nothing in reflection. Three effects stack up. First, the return trip doubles every loss in dB — a path with 30 dB insertion loss is still fine in transmission but comes back with 60 dB in reflection, and I only have ~71 dB of dynamic range. Second, the SiN-to-InP edge is deliberately **angled** to suppress reflections, and it adds a second die-to-die coupling that gets crossed twice. Third, the InP part contains SOAs, and an **unbiased** SOA is strongly absorbing at 1550 nm — easily tens of dB over a 1150 µm section, one way. I think that third one dominates.

> 📈 **PLOT — insert `facet_comparison_both_ligentec_chips.png` here**

Three fibres × two chips, all six in one figure: overview on top, the facet region at absolute position in the middle, and each trace re-centred on its own facet at the bottom. The dotted lines in the overview are each scan's RMS noise floor, so you can read the dynamic range straight off.

Two results:

1. **The new chip sits ~1.7 mm closer.** Measured relative to the internal anchor (which sits *before* the chip, so most systematics cancel), fibre 1 and fibre 44 agree: −1.66 and −1.82 mm. That's just a different chip mounting. Fibre 48 disagrees — but fibre 48 is also the channel where the two MAP2672 sessions contradict *each other* by 2.4 mm, so the odd one out is the old data, not the new chip.
2. **The −3 dB peak width belongs to the fibre, not to the chip:**

| fibre | MAP2672 | MAP2680 |
| --- | --- | --- |
| f1 | 33.0 µm | 31.0 µm |
| f44 | 35.4 µm | 43.8 µm |
| f48 | **81.9 µm** | **82.8 µm** |

Fibre 48 is twice as wide as the others on **both** chips, and the two numbers agree to 1 µm across two different physical chips. That can't be a chip effect.

### 4.4 Two smaller results

**Index matching with IPA works, and by the amount you'd expect.** I compared the raw channels with and without isopropanol on the facet, using the aux as a control:

| quantity | without IPA | with IPA | Δ |
| --- | --- | --- | --- |
| Ch1 AC (fringe amplitude) | 0.006938 mW | 0.002954 mW | **−7.42 dB** |
| Ch1 DC | 0.075185 mW | 0.015552 mW | **−13.69 dB** |
| Ch2 AC (aux, control) | 0.149320 mW | 0.148551 mW | **−0.04 dB** |

The aux is identical to 0.04 dB, so the laser, the LO and the aux path were unchanged and everything that moved in Ch1/Ch3 came from the measurement path. And the two numbers agree with each other: AC ∝ √(P_LO·P_back) while DC ∝ P_back, so the AC change should be exactly half the DC change in dB. −13.69/2 = −6.85 against a measured −7.42 — that's 0.6 dB, two independent quantities telling the same story. So **IPA dropped the back-reflection by about 14 dB**, which sits between the Fresnel predictions for glass/air → glass/IPA (~15 dB) and SiN/air → SiN/IPA (~6.6 dB).

One problem: with IPA the facet peak got **wider** (16.3 → 26.9 µm). My guess is that the isopropanol evaporates during the 2–4 minute sweep, so the index match changes *while* I'm measuring. Non-volatile gel would be the fix.

**Bare fibre arrays without any chip.** Meisu and PHIX arrays measured on their own:

| array | channels | facet position from the connector plane |
| --- | --- | --- |
| Meisu | 4, 8, 9, 12 | 1075.40 … 1076.37 mm (spread 0.97 mm) |
| PHIX | 1, 8 | 1067.51, 1069.14 mm |

The ~8 mm between the two manufacturers is just array length. But the **coupling varies by 13 dB between channels of the same array**, and re-plugging a bad channel only recovers 1–2 dB of that — so it's the array, not the connection.

---

## 5. What I don't know, and where I'd like ideas

Roughly in the order of how much they're blocking me.

**1. The chip group index — my biggest blocker.** Every on-chip length I report is only as good as an n_g I fitted against design lengths. Could someone ask **HHI (or Jeppix)** for the design manual / PDK for the E1700 waveguide (the GDS references HHI_InP_PDK 6-15-0), and **LIGENTEC** for the AN350 component port definitions and their n_g? With a real n_g I could measure unknown lengths directly instead of bootstrapping. The LIGENTEC port definitions would also let me build a proper map of the 95-channel array, which I couldn't do from the GDS black boxes alone.

**2. A request for the next tapeout.** Two loopback structures with **very different** lengths, e.g. 0.5 mm and 5 mm, instead of the 640 µm and 557 µm on the HHI chip. The couplers cancel exactly in the difference, so with a long lever the waveguide length comes out almost error-free. On the Ligentec chip it's worse: **both** loopbacks are exactly 199.5 µm, so the ratio test — the only non-circular evidence I had on the HHI chip — isn't even possible there. And the expected echo lands at 1.27–1.42 mm, right inside the instrument's own sidelobe skirt (±0.70 and ±1.53 mm around every strong peak). A useful target needs **more than ~2 mm of waveguide behind the facet.**

**3. The fibre-end doublet.** Every fibre end shows two peaks, 3–5 mm apart, with the separation different per fibre. I have a clean test for it — put index-matching gel on the end face, the end reflection dies, and whatever is left is the other feature — but I don't have a mechanism in mind at all. Ideas very welcome.

**4. The internal anchor moves between sessions.** There's a fixed reflector sitting before the chip that I use as my zero. Within one session it's stable to 70 µm over six scans; **between** sessions it jumps by up to 1.2 mm. That's either a τ_aux error of ~0.1 % or an internal connector pair that got re-plugged. Until it's sorted, absolute positions carry ±1 mm systematically between sessions, so all cross-session comparisons have to be done relative to the anchor.

**5. A broad raised band at 2.0–2.6 m in every high-dynamic-range scan.** About 30 dB above the control, ~600 mm wide, with no individual peaks in it. My first guess — second harmonic of the very strong facet — doesn't work: the maximum is 60–70 mm away from where the harmonic would be. It shows up in *every* scan whose noise floor is below about −58 dB, which makes me think it's always there and just becomes visible when the dynamic range is good enough. Still unexplained.

**6. The coupling spread across array channels.** 13 dB between channels of the same Meisu array, and it survives a re-plug; 13.5 dB between the two HHI coupler pairs. I believe it's the fibre-array glueing rather than the couplers themselves, but I haven't proven it, and it costs me dynamic range on the bad channels.

**7. A calibration bias I never got to the bottom of.** In August the fringe counting for τ_aux had a roughly constant 10–15 % scale error. It went away with the aux rebuild, and the current setup is verified to 17 µm on a 2 m path (section 3.2), but I never found the cause — so I can't promise it won't come back. The cheapest next check is to verify the aux fibre's real physical length with a tape measure.

**8. Upgrades that are planned but not built.** Each of these bolts on without redoing anything: real balanced photodetectors plus a fast digitiser (only needed if we ever want Rayleigh backscatter — software subtraction mathematically can't reach it), a polarisation-diverse receiver (removes fading and amplitude errors), and dispersion autofocus (a few-cell effect at 40 cm, dominant beyond a metre). I'd like input on whether any of these is worth doing now, or whether the chip work should stay the priority.

---

## One thing I learned about method

The pattern that has held every single time: **results where the geometry forces the answer survive; results that rest on "a peak in a dense spectrum matches a prediction" don't.** Three separate results this month looked convincing and then fell over on exactly that — an n_g determination from the netlist, an SSC/waveguide index split, and a "third loopback pair" that turned out not to exist. Each one was killed by a test, not by a better argument. So the rule on the bench is that every guess about an unexpected feature gets a discriminating *measurement* before it gets written down.
