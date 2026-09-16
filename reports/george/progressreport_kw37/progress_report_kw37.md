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

## How I tried to get out of it, and why it did not work

I needed n_g from a structure that is *not* one of the loops. The single-ended channels looked like the way to do it: they run from a coupler into the circuit, and every component interface along the way sits at a design distance I can compute.

Two things had to be fixed first:

- My waveguide trace only followed the waveguide up to the *first* component. The light carries on *through* the components, so I was predicting one reflector where there are dozens. I therefore built a proper optical netlist out of the GDS: every waveguide polygon contributes its centre-line length, every HHI component contributes the distance between its own optical ports, coincident ports are joined, and Dijkstra from a coupler then gives the path length to every interface in the chip. That part works, and I validated it against known geometry — starting from one coupler of the top loop pair it returns 640.0 µm to the other coupler's port, which is the loop length to 0.13 µm.
- Channels running deep into the decoder give a dense forest of peaks — 448 peaks above −45 dB over 70 mm for one of them. That is multipath: the decoder has 23 MMI splitters (20 of them 1×2) and six 2R units with two rings each, so twelve ring resonators in which light recirculates.

With the netlist I then found what looked like a clean result: one channel where three predicted interfaces (MMI input, phase shifter, metal crossing) all showed up in the measurement within 0.7 of a resolution cell and gave the same n_g to 0.08 %. I was about to report that as an independent determination of n_g.

**It does not hold, and I want to be explicit about why.** The netlist has 690 interfaces spread over 29 channels. At a matching tolerance of two resolution cells there is a candidate interface near almost any predicted position — so those matches were *selected* by the n_g I had assumed, and cannot then be evidence for it. The test that settles it: sweep n_g from 3.0 to 4.0 and ask how well the best channel assignment fits at each value. If n_g were determinable there would be a clear maximum near 3.48. Instead the score is flat — it stays between 7.3 and 19.3 across the whole range, with its maximum at n_g = 3.02. A flat curve means the data does not contain the information.

The same applies to identifying *which* channel was plugged: the two best candidates come out within 1 % of each other in score, so a reflectogram alone cannot tell me which coupler I am on.

## Where that leaves the result

Honestly: **the circularity is still there.** n_g = 3.4803 ± 0.008 is a calibration against the design lengths of the two loops, not an independent measurement, and I could not verify it from the chip itself.

What does hold, and does not depend on any peak-hunting:

- **The ratio test.** Two loop measurements, one n_g, so one degree of freedom is left over. The measured ratio of the two chip paths is 1.026008 and the design ratio is 1.028268 — they agree to 0.22 %, against a resolution limit of 0.47 % on that ratio. If the geometry or the peak identification were wrong, this would not close.
- **Each loop reproduces its design path.** Given n_g, the top pair lands 2.3 µm from its predicted facet-to-facet distance and the bottom pair 13.1 µm, i.e. within one and two resolution cells (one cell is 7.0 µm on the chip). The loopback topology is what makes this trustworthy: there are exactly two facets and nothing to interpret.

So what I can defend is: *the measurement is consistent with the design geometry at the resolution limit, and the instrument resolves about 7 µm of on-chip waveguide.* What I cannot yet defend is an absolute length measurement, because that needs n_g from outside.

The methodological lesson, which held three times this week: every result that rested on "a peak in a dense spectrum matches a prediction" fell over. Every result where the geometry forces the assignment survived.

I also wrote a script that takes a raw scan, calibrates itself, and checks the reflectogram against a *named* channel from the netlist — it prints which component each peak corresponds to and what n_g each one implies, and it has the n_g sweep built in as a check so the failure mode above cannot hide. Without being told the channel it only prints a ranking, with a warning when the top candidates are too close to separate.

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
- **netlist_fit_that_failed.png** — the fit that looked like an independent confirmation: three component interfaces of one candidate channel and the three consistent n_g values they give. I am attaching it as the cautionary example, not as evidence — the n_g sweep shows equally good fits exist at any n_g between 3.0 and 4.0, so this picture is convincing and wrong.
- **loopback_verified.png** — what a trustworthy check looks like instead: the top loop pair against its predicted facet-to-facet distance, where the topology leaves nothing to interpret.
- **measured_loops_comparison.png** — both loops aligned on facet A; the chip path shrinks by exactly the expected amount for the shorter loop.
- **gds_loop_top.png** — the loop as drawn, with the length of every segment.
- **why_circuit_channels_fail.png** — why I can only use loopback channels: loopback gives two peaks and nothing else, a circuit channel gives a forest.

— Carolina
