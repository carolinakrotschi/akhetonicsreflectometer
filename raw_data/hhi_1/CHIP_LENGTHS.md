# On-chip waveguide length per scan, from the raw data

Produced by `analyze_hhi1_chiplength.py` in this folder, which reads only the
`*.json` scans here. Design numbers come from `HHI_RUN_1.gds` in this folder.
Fiber-to-coupler mapping used throughout: **b = 32 − fiber**.

## How the number is obtained

The reflectogram axis is the usual reflection convention
`z = c·τ / (2·n_g)` with the **fibre** index n_g = 1.468. A path of physical
length `L` inside the chip therefore appears at

    Δz = L · n_chip / n_fiber        →        L = Δz · n_fiber / n_chip

So every on-chip length is a **distance between two peaks**, not an absolute
position: facet A (the fibre → InP interface where the light enters the chip)
and whatever reflects further in.

`τ_aux` is calibrated per scan from that scan's own EXFO wavelength table
(total unwrapped aux phase / 2π·frequency span). It came out at
**20.3797 … 20.3816 ns** over all eleven scans — a spread of 4.5 ppm, i.e.
the z axis is consistent between scans to ±8 µm at z = 1.7 m.

Instrument figures, identical for all scans: 988,536 points, 1520–1570 nm,
resolution **16.59 µm** in z (≈ 7.0 µm of chip), Nyquist range 8.04 m,
7.7 points per aux fringe.

### Why the facets are located by cross-correlation

Facet A and facet B are not single peaks but clusters of sub-peaks ~80–100 µm
apart, and which sub-peak is the highest flips between scans. Picking "the"
peak gave 7.2003 vs 7.2782 mm for the same channel on two days — a 78 µm
disagreement that is an artefact of that choice.

Facet A and facet B are the same kind of Fresnel reflection seen through the
same instrument, so the whole cluster has the same shape at both ends.
Cross-correlating the two clusters instead of picking peaks gives
**7.1999 vs 7.1997 mm** for the same pair — reproducible to 0.2 µm across
three days and a re-plug. That is the method of record.

## Result 1 — the loopback channels (solid)

For fibers 4 and 32 the light leaves the chip again through the *neighbouring*
coupler, so the far reflection is a real chip facet and the topology forces
the identification. Nothing has to be assumed.

| scan | fiber | coupler | Δz [mm] | L_chip [µm] | GDS design [µm] | diff |
|---|---|---|---|---|---|---|
| `2026-09-11-12-16hhi1fiber4` | 4 | b28 | 7.1999 | 3039.0 | 3038.8 | +0.2 µm |
| `2026-09-14-09-03_fiber4` | 4 | b28 | 7.1997 | 3038.9 | 3038.8 | +0.1 µm |
| `2026-09-11-12-16hhi1fiber32` | 32 | b0 | 7.0012 | 2955.1 | 2955.3 | −0.2 µm |
| `2026-09-14-09-35_fiber32` | 32 | b0 | 7.0012 | 2955.1 | 2955.3 | −0.2 µm |

Design path = SSC + loop + SSC = 2 × 1199.4 µm + 640.0 µm (b28) and
+ 556.5 µm (b0).

**Chip group index**, fitted from these four measurements:

    n_chip = 3.47794      spread between the two channels 0.00039 (0.011 %)

This is the one fitted quantity, so the four rows above are not an independent
test on their own. The **non-circular** statement is the cross-prediction:

| calibrate on | predict | predicted L [µm] | GDS [µm] | error |
|---|---|---|---|---|
| b0 (n = 3.47774) | b28 | 3039.1 | 3038.8 | **+0.3 µm** |
| b28 (n = 3.47812) | b0 | 2955.0 | 2955.3 | **−0.3 µm** |

i.e. two loops differing by 83.5 µm in design agree with each other and with
the GDS to **1 part in 10⁴** on a ~3 mm path.

Confirmed independently by the double-bounce: a second echo appears at
2 × Δz (14.44 mm for fiber 4, 14.04 mm for fiber 32) at −43 to −52 dB, which
is what a Fabry–Pérot between facet A and facet B must produce.

### What is *not* resolved: the SSC / waveguide index split

The two loops differ only in plain waveguide (83.5 µm), the SSCs being
identical. Solving the two equations separately gives

    n_waveguide = 3.4916,   n_SSC = 3.4745

but the lever arm is only 83.5 µm, so **1 µm of error in Δz moves n_waveguide
by 0.025**. With a realistic ±2 µm systematic that is 3.49 ± 0.05 — consistent
with the lumped 3.4779 and not actually a measurement. Use the lumped index;
the split needs a channel pair with a much larger design difference.

## Result 2 — the one-ended channels (provisional)

These have no second facet, and the chip cluster has no clean onset, so
**facet A cannot be measured from the scan itself**. Two independent
approaches failed here: the leading-edge detector picks up a long low-level
skirt, and a design-comb fit is under-determined (the rich channels have
80–170 interfaces, so almost any offset "fits").

What can be done is to ask, for each design interface of the assigned channel,
what facet-A position it would imply, and check that against the two
directly measured values (fiber 4: 1670.15 / 1670.23 mm; fiber 32: 1668.31 /
1668.55 mm).

| scan | fiber | coupler | strongest chip peak | best design interface | implied facet A | verdict |
|---|---|---|---|---|---|---|
| `...hhi1fiber7` | 7 | b25 | 1673.298 mm | 1199.4 µm (SSC end) | 1670.456 mm | **consistent** — b25 is a stub, the design has exactly one interface and the data has exactly one strong peak |
| `...hhi1fiber27` | 27 | b5 | 1683.893 mm | 5671.1 µm (far end of b5) | 1670.457 mm | **consistent** — agrees with fiber 7 to 1 µm |
| `...hhi1fiber6` | 6 | b26 | 1697.255 mm | 11377.8 µm | 1670.299 mm | **consistent** |
| `...hhi1fiber30` | 30 | b2 | 1679.373 mm | — | — | **fails**: b2's design jumps from 1199.4 µm straight to 6166.7 µm, and the measured peak falls in that gap (3830 µm for facet A = 1670.3) |
| `2026-09-14-09-58_fiber23` | 23 | b9 | 1685.537 mm | — | — | **fails**: b9's design only reaches 2637 µm; the peak needs ~6430 µm |

Fibers 7, 27 and 6 independently imply facet A = 1670.46 / 1670.46 / 1670.30 mm,
agreeing with each other to 0.16 mm and with fiber 4's *measured* 1670.15–1670.23 mm.
That is a genuine side result: **array fibers 4, 6, 7 and 27 have the same
pigtail length to within ~0.3 mm**. Fiber 32 is 1.8 mm shorter.

Caveats on the two failures:

- For **fiber 23 / b9** the design model is known to be incomplete: the netlist
  deliberately refuses to carry the path through the `WGMETxE1700twin`
  crossings (they carry two *independent* guides side by side; connecting them
  produced the phantom b9↔b10 loopback that was retracted). So b9's design comb
  stopping at 2637 µm is a limitation of the extraction, not necessarily of the
  chip. Do not read this as disproving b = 32 − fiber.
- For **fiber 30 / b2** the gap in the design comb is real, so this one is an
  open discrepancy.
- Either strongest peak could also be multipath rather than a design interface.
  A single peak cannot identify a channel: with 80–170 interfaces, almost every
  rich channel has *some* interface that would place facet A in the plausible
  window.

**To settle a one-ended channel you need a second known reference in the same
scan.** The loopbacks have one by construction; that is exactly why they are
the only solid rows here.

## Files

| file | what |
|---|---|
| `analyze_hhi1_chiplength.py` | the whole analysis, reads only this folder |
| `chip_lengths.csv` | the table above in machine-readable form |
| `loopbacks.png` | the quantitative result: facet A → facet B for fibers 4 and 32 |
| `reflectograms_chip.png` | chip region (1660–1702 mm) of all eleven scans, facets marked |
| `reflectograms_full.png` | full 0–3 m reflectogram of all eleven scans |

The `2026-09-14-08-56_nichtsconnected` control shows **nothing at all** between
1660 and 1702 mm, which is what makes everything in that window chip signal.
