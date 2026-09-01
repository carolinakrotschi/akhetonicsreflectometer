# Data: OFDR fibre-length measurements (2026-09-01, 17:03–17:57)

Recorded with the **TOF1550 filter removed** — these are the real
measurements. Each capture is one sweep, all four channels, 100 kHz,
1,050,000 samples per channel.

Setup: `exfo-1` + `coredaq-1`, 1520 → 1570 nm at 5 nm/s.

- **Ch1, Ch3** — measurement MZI, contains the test fibre
- **Ch2, Ch4** — reference ("aux") MZI, own path imbalance 4.16 m

Produced with `lina/scripts/lina_ofdr_test.py`. Full reasoning: `../../report.md`,
§6–§9.

---

## The consideration

An FFT of the interferogram gives a delay, and therefore a length — but only if
the signal is uniform in **optical frequency**. The coreDAQ samples uniformly
in **time**, and the EXFO sweeps 1–2 % fast and slightly non-uniformly. The aux
MZI's fixed delay makes its fringe phase the optical-frequency axis, so the
measurement is resampled onto a uniform aux-phase grid before the FFT.

Both channel pairs are evaluated, each using the *other* pair as its phase
reference. That makes them two independent measurements of the same fibre,
which is the strongest cross-check available without an external standard.

---

## Naming

`<state>_...` — which fibre was connected when the capture was taken:

| Prefix | State |
|---|---|
| `1mA_` | test fibre "1mA" (nominally 1 m) |
| `1mB_` | test fibre "1mB" (nominally 1 m) |
| `3mB_` | test fibre "3mB" (nominally 3 m) |
| `3mB1mB_` | 3mB and 1mB connected in series |
| `nofiber_` | no test fibre connected |

`sweep_ofdr_*.npz` = raw capture (4 channels, cropped to the sweep window by
the λ(t) calibration). `ofdr_*.png` = the diagnostic figure produced with the
capture. `*_PEAKS_*.png` / `OVERLAY_*.png` = evaluation figures, regenerable
from the NPZ at any time.

---

## Files

| File | What it shows |
|---|---|
| `OVERLAY_all.png` | **The main result.** All five states in one figure, one colour each. Upper panel: axis reads test-fibre length. Lower panel: the aux pair in its own axis. |
| `1mB_sweep_ofdr_run1/2/3_*.npz` | Three repeats of the same state → peak reproducibility ±0.03 mm (1.14795 / 1.14800 / 1.14797 m in the raw axis). |
| `1mA_sweep_ofdr_*.npz` | The fibre swap. Two peaks move by ±31 mm, the strong peaks do not — this is what showed that the strongest peak is *not* the fibre. |
| `nofiber_sweep_ofdr_*.npz` | **The decisive measurement.** Exactly the two peaks that depend on the fibre disappear; everything else stays. Without this, the peak assignment would be guesswork. |
| `3mB_sweep_ofdr_*.npz` | Confirms the factor of two: 3 m fibre → 6.051 m of path, so the fibre is traversed twice. |
| `3mB1mB_sweep_ofdr_*.npz` | Additivity test: predicted 4.06533 m from the two single measurements, measured **4.06334 m** → 2.0 mm. Also shows the connector between the fibres at 3.03 m. |
| `sweep_ofdr_3mBonly_check_*.npz` | Control capture taken after the GUI button produced a wrong length — the script found 3.02522 m on the same hardware minutes later, which localised the problem to the GUI capture path (report §11.1). |
| `CROPSCAN_3mB1mB.png` | The crop scan that **refuted** the idea that the first 10 nm must be discarded: SNR unchanged at 41–42 dB, FWHM worse at heavy crops. |
| `*_PEAKS_stacked.png` | Aux above, measurement below, in the raw axis (earlier presentation stage). |
| `*_PEAKS_fibrelength.png` | Same, but the measurement axis reads fibre length (offset subtracted). |
| `1mB_PEAKS_testfiber_0-3m.png`, `1mB_PEAKS_aux_1-5m.png` | Intermediate stage: separate zoomed views before the stacked figure existed. |
| `1mB_sweep_..._fft_*.png` | Re-evaluations of the same capture with different processing (with/without high-pass, swapped channel roles) — the steps that identified which pair is which. |

---

## Results

Axis reads test-fibre length (reflection convention, MZI offset 0.57396 m
subtracted):

| Fibre | Peak | Prediction | dB |
|---|---|---|---|
| 1mB | **1.0397 m** | — | −22.4 |
| 1mA | **1.0553 m** | — | −18.4 |
| 3mB | **3.0256 m** | — | −31.2 |
| 3mB + 1mB | **4.06334 m** | 4.06533 m (sum of the two) | −31.9 |
| no fibre | no peak | — | — |

- 1mA is **15.6 mm longer** than 1mB.
- Lengths are additive to **2.0 mm (0.05 %)**.
- Repeatability ±0.03 mm; systematic uncertainty ~±2 mm from the peak's
  sub-lobe structure (report §9.1).
- All lengths scale linearly with `--n-group` (1.4682, SMF-28 at 1550 nm).

---

## Reproducing

```bash
# New measurement (filter must be OUT of the path)
python lina/scripts/lina_ofdr_test.py --tag myrun measure \
    --start 1520 --stop 1570 --speed 5 --rate 100000

# Both spectra of one capture, stacked, axis reading fibre length
python lina/scripts/lina_ofdr_test.py --pad-factor 4 --n-peaks 3 peaks \
    --raw lina/data/02_ofdr_measurements/3mB_sweep_ofdr_run1_*.npz \
    --meas-offset-m 0.57396

# The overlay of all states
python lina/scripts/lina_ofdr_test.py --pad-factor 4 overlay \
    --raw 1mA=<npz> 1mB=<npz> 3mB=<npz> 3mB+1mB=<npz> nofiber=<npz> \
    --meas-offset-m 0.57396 --meas-xlim 0,4.5 --aux-xlim 1.5,4.5
```

The `--meas-offset-m 0.57396` value is the measurement MZI's own path
imbalance, taken from the `nofiber_` capture. Re-measure it whenever the
interferometer itself is changed.
