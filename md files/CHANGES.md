# Findings and Corrections

## Measurement: aliased

| Quantity | Value |
|---|---|
| Step size δλ | 0.9998 pm |
| Bandwidth | 49.908 nm = 6.268 THz |
| Resolution δz = c/(2 n_g B) | 16.3 µm |
| **Nyquist range** λ²/(4 n_g δλ) | **40.66 cm reflector distance** (81.3 cm OPD) |
| Reflectogram peak | 34.8 cm — not a peak, a broad hump |
| Energy in the top 20% before Nyquist | **68%** |

Ch1 and Ch2 correlate at **−0.93** (out-of-phase coupler ports), the sum
Ch1+Ch2 is nearly constant (σ/µ = 0.10), the difference carries the full
fringe (σ/µ = 0.52). So the signal is **real interference with full
modulation** — not noise. It's simply beyond Nyquist.

Best single-target hypothesis from three independent methods (first-order
folding, chirp slope across sub-bands, dechirped NUDFT scan):
**OPD ≈ 92 cm, corresponding to a 46 cm reflector distance.** That's ~13%
beyond the limit. The number only becomes trustworthy with a sweep at
δλ ≤ 0.5 pm.

## Bugs in reflectometer_v1_buggy.py (formerly fft_reflectometer.py)

*(2026-08-18: this script and its corrected sibling `reflectometer_v2_fixed.py`
were deleted — the bug list below is the reason they existed in the first
place, and it's preserved here along with git history; keeping the actual
buggy code around had no further value.)*

1. **Last sample is garbage.** `Wavelength[-1] = 1519.95` (jumps back to the
   sweep start). This turns `np.mean(np.diff(wavelength))` into **2.40e-6 nm**
   instead of 1e-3 nm — the frequency axis ends up scaled wrong by a factor
   of **416**. → discard non-monotonic points, use `np.median` instead of
   `np.mean`.

2. **Wrong conversion in the λ branch (default).** `length_nm = peak_frequency
   / refractive_index` only holds in the wavenumber branch. In the λ domain,
   f = n_g·ΔL/λ², so **ΔL = f·λ²/n_g**. The λ² ≈ 2.39e6 nm² factor is missing.

3. **Missing factor of 2.** For a reflection measurement, the reflector
   distance is OPD/2.

4. **No DC/envelope removal.** `argmax` over `freq > 0` picks up bin 1–2 (the
   laser power envelope) instead of the actual signal. Result from the old
   version: 5.7e-9 m and 1.3e-2 m — both artifacts.

5. **No window.** Without a Hann window, sidelobes smear weak reflections.

6. **Ch2 unused.** `(Ch1−Ch2)/(Ch1+Ch2)` removes the power envelope exactly
   and doubles the fringe contrast. Free improvement.

7. **No Nyquist check.** The script never warned that the measurement was
   outside the valid range.

8. `plt.xlim(0, 3)` suggests a 3 m measurement range. The real range is
   40.7 cm.

9. The reported wavelengths are **rounded** to 1 pm (step sequence 0/1/2 pm).
   For the phase, an ideally uniform grid (`linspace(λ_0, λ_N, N)`) is the
   better assumption — the rounded values add up to ±0.5 pm of jitter, which
   at OPD ≈ 1 m already means ~1.9 rad of phase error.

## Next step

δλ = 0.5 pm → range 81 cm, 100,000 points.
δλ = 0.2 pm → range 2.03 m, 250,000 points (safety factor 2 for 1 m).

Resolution stays at 16.3 µm in all cases — it only depends on bandwidth.
