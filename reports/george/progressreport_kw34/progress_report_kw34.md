

Hi — since I got your claude handover docs,
I've got the  setup with the aux interferometer built and
measured with it and also tuned the akhelab software so that I could measure in this new way. 

---

## 1. How I compute the distance axis now

The laser doesn't sweep perfectly linearly (you already knew this). Instead
of correcting it with a wavelength table, I use a second reference
interferometer (the **aux MZI**, fixed arm-length difference ΔL) as a clock
for the actual laser frequency:

- The aux's phase is proportional to optical frequency. Equal phase steps
  mean equal frequency steps, no matter how non-linear the sweep actually is.
- I resample the measurement signal onto uniform aux phase steps, then FFT.
- Distance: `z = c · τ / (2 · n_g)`, where τ comes from the FFT frequency
  and the calibrated aux delay `τ_aux` 

 how accurate a position reading
is depends almost entirely on that one number, τ_aux — not on the rest of
the code. I ran this against simulated data with known ground truth, and a
0.5% error in τ_aux already costs about 1cm at 2m distance, and it scales
linearly . The annoying part: the plot looks
completely fine while this happens...

---

## 2. measurement results
With the new setup I tested 5 different fibers (the ones that were available in the lab), 4 one-meter fibers, one 3-meter fiber and one configurationwith no fiber connected at all. It stood out, that there was a big peak at  (~587mm) shows up in every single scan, even with no
fiber connected. And also the 3 1m fibers had a clustered peak at around 1,60m and the 3m fiber had a peak at 3.6m.

Regarding the big peak in every scan:  it can't be coming from the test fiber — it's a fixed
internal partial reflection. I tracked it down with connector-swap tests:

| Connector in the slot | Peak position |
| Original (Thorlabs) | ~587–590mm |
| random green Connector from the box on the table #1 | ~467mm |
| random green Connector from the box on the table #2 | ~590mm |

 so this peak is just the
back-reflection of whatever connector is there but I couldnt resolve it entirely - even by cleaning the connectors, fibers,...

but very annoyingly:  light that partially reflects
off the connector, keeps going, reflects off the fiber end, and comes back
travels the **sum** of both path lengths:

```
so my theory, why the 1m fibers all arive at 1.6m and the 3m fiber at 3.6: peak position you see  =  connector position (~587mm)  +  true fiber length (this calculated value is also what the legend of the plot says)
```

Checked on two fibers: the 1m fiber (1.04m by tape measure) shows up at
1649mm (587 + 1040 ≈ 1627, off by ~1%). The 3m fiber (3.05m by tape) shows
up at 3735mm (587 + 3050 ≈ 3637, off by ~3%). What gives it away is that
it's the same **fixed offset** both times.

So I land at this current hypothesis:  subtract the
connector's position from whatever peak you see, and that gets you the real
fiber length.

All the other peaks in the signal showed no structure with all the fibers plotted together, so maxbe they are real back reflections on the way in the fiber or peaks from the aux,...

---

## 3. Problems

First problem I faced and how I fixed it: 

      Every scan had the first ~10-13nm saturated (detector clipping), which threw
      off the monotonicity/bow numbers pretty badly. Turned out it's not a
      software range issue — it still saturates even at the lowest gain setting —
      it's a real optical settling transient right at the start of the sweep.

      **Fix, in the AkheLab software:** added a *"Settle margin (nm)"* field — the
      physical sweep now starts a bit earlier than requested, so the transient is
      already over by the time the laser reaches your actual measurement range,
      then the extra lead-in samples get cropped back out automatically. leading to  resolution is back to its full value (**16.6µm instead
      of 22.6µm with trimming, about 27% better**).

Second Problem: Akhelabs Software was not taylored to the new setup yet

      To run the free-running setup at 1M samples/channel, I extended
      `exfo_window.py`:

      - **Samples/ch + Rate (Hz):** the coreDAQ capture is one single free-running
        trigger — there's no point-by-point sync with the laser. So the capture
        duration (`Samples/ch ÷ Rate`) and the sweep duration (`span ÷ speed`) have
        to line up, otherwise the buffer fills before the sweep is done, or sits
        there waiting. Added a live readout, `≈X.XXs capture (sweep ≈X.Xs)`, that
        warns you when the two are off by more than 15%.
      - **Settle margin (nm):** see above — stretches the sweep a bit at the start,
        then trims it back down automatically. If you set samples/rate manually,
        you might need to bump Samples/ch up slightly so enough points still land
        in your actual target range despite the longer physical sweep.
      - **Anti-saturation range:** forces the lowest gain on every channel before
        the sweep starts — my first attempt at fixing the saturation, but it
        wasn't enough on its own (see point 3).

---

## Still open

- I've found the connector reflection, but haven't gotten rid of it yet —
  next step is probably index-matching or a lower-back-reflection connector.
- The τ_aux calibration needs to get tighter before I'd trust a fiber length
  as an actual number in a report (see point 1).
- Also maybe buy a proper fiber-optic retroreflector from Thorlabs so we
  have a known, strong, well-defined reflection to test/calibrate
  against instead of relying on ambiguous connector-family reflections:
  https://www.thorlabs.com/fiber-optic-retroreflectors?tabName=Overview

Happy to walk through any of this in more detail, or show you the raw data/plots
if you want to look yourself.

— Carolina
