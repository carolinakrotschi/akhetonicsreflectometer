Hi — update for this week. This one is about the HHI chip, not fibres.

## What I set out to do

I wanted to measure a waveguide length on the chip with the reflectometer. The structure I picked is the little loop at the top left of HHI_RUN_1 — the one that connects two neighbouring edge couplers to each other. I first extracted its length from the GDS by summing up the arcs and straights along the centre line:

| # | type | R (µm) | angle | L (µm) |
|---|---|---|---|---|
| 0, 3 | arc | 93.544 | 91.306° | 149.072 each |
| 1, 2, 6, 7 | arc (S-bends) | 58.967 | 57.466° | 59.143 each |
| 8 | straight | – | – | 45.513 |
| 4, 5, 9, 10 | port S-bends | – | – | 11.611 / 18.338 each |
| | **total** | | | **640.1266** |

The waveguide is on layer 12/0 and is exactly 2.300 µm wide. I got that width out of the layout as a by-product: for a ribbon of constant width, area and perimeter give you both the centre-line length and the width, and it comes out 2.3000 µm on every single polygon, which is how I knew I had the right layer.

The second loop pair, at the bottom left, is 556.5494 µm. Those two are the only loopback pairs on the chip — the other 27 couplers all run into the decoder circuit.

## How the measurement is set up

Measurement MZI → about 1 m of fibre → into the chip through coupler A → around the loop → out through coupler B → about 1 m of fibre → open fibre end, which is what reflects back. Both couplers are permanently pigtailed.

So in the reflectogram I see four things in a row: the internal reflection I already know about (around 574 mm), then facet A, then facet B, then the open end of the second fibre. For the top loop:

| station | z (mm) | section before it |
|---|---|---|
| internal reflection | 573.928 | – |
| facet A | 1670.227 | fibre 1 = 1096.298 mm |
| facet B | 1677.426 | **chip = 7.200 mm** |
| open end of fibre 2 | 2773.625 | fibre 2 = 1096.199 mm |

The bit I liked here: I found fibre 1 and fibre 2 completely independently of each other, and they came out 1096.298 and 1096.199 mm — the same to 0.099 mm. Nothing in the analysis forces that, so it tells me the four stations are correctly identified.

## The problem I ran into

The reflectometer measures time, not length. To turn the 7.200 mm into a length I need the group index of the chip, and it is in none of the documents I have — not the Jeppix validation report (that one has waveguide loss 1.8 dB/cm, butt-joint loss 1.2 dB, SSC coupling loss typ. 1.5 dB, but no index), not the device mapping, not the pin mapping. And I don't have the HHI PDK.

So at first I did the obvious thing: I took the two loops, used their known design lengths, and solved for n_g. That gives 3.4803 ± 0.008. Then I used that n_g to "measure" the loops and got 636.8 and 559.8 µm against 640.1 and 556.5 design.

**That is circular** and I want to flag it clearly, because I nearly reported it as a result. If I fit n_g so that both loops come out right, then of course they come out right. The only genuinely non-circular thing in there is the *ratio*: two measurements, one n_g, so one degree of freedom is left over. The ratio of the two measured chip paths is 1.026008 and the ratio of the two design paths is 1.028268 — they agree to 0.22 %, and the resolution limit on that ratio is 0.47 %. That part is a real test, and it passes. But it is not a length measurement.

## How I got out of it

I needed n_g from a structure that is *not* one of the loops. The single-ended channels should do that — they run from a coupler into the circuit, and every component interface along the way is a potential reflector at a known design distance.

The first attempts failed, and it took me a while to understand why:

- The waveguide trace I had only followed the waveguide up to the *first* component. But the light carries on *through* the components, so I was predicting one reflector when there are dozens.
- Channels that run deep into the decoder give a dense forest of peaks — 448 peaks above −45 dB over 70 mm for one of them. That is multipath: the decoder has 23 MMI splitters (20 of them 1×2) and six 2R units with two rings each, so twelve ring resonators in which light recirculates. In a forest that dense there is always some peak within 0.17 mm of any prediction, so matching a single peak to a single prediction proves nothing. I fell into exactly that trap once.

So I built a proper optical netlist out of the GDS instead: every waveguide polygon contributes its centre-line length, every HHI component contributes the distance between its own optical ports, and coincident ports are connected. Then Dijkstra from a coupler gives the path length to *every* interface in the chip. I validated it on the top loop, where I know the answer: it returns 640.0 µm from one coupler port to the other.

With that, channel b2 gives three interfaces that I can find in the measurement:

| measured Δz | what it is | design path from facet A | → n_g | miss |
|---|---|---|---|---|
| 14.6148 mm | MMI 1×2 input | 6166.7 µm | 3.47909 | −6.4 µm |
| 14.9300 mm | phase shifter PMTO500 | 6296.9 µm | 3.48064 | +0.1 µm |
| 15.2618 mm | metal crossing | 6442.0 µm | 3.47785 | −12.1 µm |

All three within 0.7 of a resolution cell, and the three n_g values agree with each other to 0.08 %. From the strongest peak, with the longest lever:

**n_g = 3.4791 ± 0.0059 (0.17 %), and this does not use the loops at all.**

It agrees with the circular loop value 3.4803 to **0.035 %**. There is a nice self-check in there too: if I had picked the wrong peak for facet A, off by one satellite spacing (0.73 mm), n_g would come out 3.65, i.e. 5 % away. The fact that it lands on 0.035 % confirms the facet identification as well.

## The actual result

Applying the b2 group index to the loops — now not circular:

| loop | Δz measured | loop measured | design | error |
|---|---|---|---|---|
| top (b27/b28) | 7.1996 mm | **639.1 µm** | 640.1266 µm | **−1.1 µm** |
| bottom (b0/b1) | 7.0171 mm | **562.1 µm** | 556.5494 µm | **+5.5 µm** |

One resolution cell is 7.0 µm on the chip, so both are inside one cell. I would say the reflectometer measures on-chip waveguide length to a few µm.

Two side results:

- The fibre-to-coupler wiring. I had no map for this, so I worked it out from the measurements: fibre n goes to coupler b(32−n), i.e. the channel number is the position counted from the top. I confirmed it twice — once with a channel that should run into the circuit (it does, forest and all) and once with a channel that should be a dead-end stub (it is: a single sharp reflector 33.2 µm wide, which is *narrower* than the window limit of 43.1 µm, and nothing at all behind it). 13 of the 29 couplers are stubs.
- The two coupler pairs couple very differently. The return signal through the bottom pair is 13.5 dB stronger than through the top pair. The path lengths differ by 83 µm, which is 0.03 dB of propagation loss, so that is all coupling. Against the validation report spec (SSC coupling loss typ. 1.5 dB, max 4 dB), the top pair is around 4.9 dB per transit, i.e. out of spec. I think that is the fibre-array glueing on those channels rather than the SSCs themselves. I will use the bottom pair from now on.

## What is still open

- One channel's scan has a single strong reflector at +9.2069 mm behind facet A, which is 2686 µm of on-chip waveguide, and that matches **no interface of any channel** in the full netlist. I have no explanation.
- Both loop scans show two peaks at +3.135 and +4.064 mm behind facet A that I cannot place. I thought for a while they were inside the SSC, but the stub channel has a complete SSC in its path and shows nothing there, so that idea is dead.
- There is a bookkeeping conflict I could not resolve: the channel that physically behaves like b2 should be fibre 30 under the wiring rule, but I had recorded fibre 27 for that scan. The physics is unambiguous (three matched interfaces), so I trust the geometry over my label, but I don't know where the discrepancy comes from.

## What would help most

Could you ask HHI (or Jeppix) for the design manual / PDK documentation for the E1700 waveguide? The GDS references HHI_InP_PDK 6-15-0. If I had their n_g value, I would no longer need any of this bootstrapping — I could just measure unknown lengths directly. As it is, my n_g is only as good as the assumption that the fabricated lengths match the drawn ones.

And a design request for the next tapeout: two loopback structures with **very** different lengths, e.g. 0.5 mm and 5 mm, instead of 640 and 557 µm. The couplers cancel exactly in the difference, so with a long lever the waveguide length comes out almost error-free. With the 83 µm difference this chip offers, the lever is too short — that was the main limitation all week.

## Plots attached

- **measured_loops_explained.png** — the whole measurement in four rows: the physical chain, the same thing as measured, a zoom on the chip, and a detail on how close facet B lands to its predicted position. The hatched bars in row 3 are model, not data — the loop edges are genuinely not visible in the reflectogram, because the waveguide passes smoothly into the SSC and its ends do not reflect.
- **b2_netlist_match.png** — the three component interfaces of channel b2 and the three n_g values they give, against the loop value.
- **measured_loops_comparison.png** — both loops aligned on facet A; the chip path shrinks by exactly the expected amount for the shorter loop.
- **gds_loop_top.png** — the loop as drawn, with the length of every segment.
- **why_circuit_channels_fail.png** — why I can only use loopback channels: loopback gives two peaks and nothing else, a circuit channel gives a forest.

— Carolina
