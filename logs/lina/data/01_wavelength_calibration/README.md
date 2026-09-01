# Data: λ(t) wavelength calibration (2026-09-01, 15:27–16:21)

Everything in this folder was recorded with the **TOF1550 filter in the optical
path**. The filter was purely an instrument: it transmits only near its centre
wavelength and therefore marks a known wavelength at a measurable time in the
capture. For the actual OFDR measurements (`../02_ofdr_measurements`) it was
removed again.

Devices: `exfo-1` (T200S) + `coredaq-1` (LINEAR, FW v4.2) + `tof1550-1` (COM6).
Produced with `lina/scripts/lina_sweep_test.py`.

Full reasoning and all results: `../../report.md`, §2–§5.

---

## The consideration

The GUI mapped the capture onto wavelengths with
`linspace(start, stop, n_samples)`. That assumes the buffer contains exactly
the sweep, and that the sweep runs at the commanded speed. Both are false.
Rather than infer the mapping, it was measured.

---

## Files

| File | Produced by | What it shows |
|---|---|---|
| `filter_cw_t1_*.csv/.png` | `filter-cw --filter-nm 1550 --span 1.5 --step 0.025 --channels 1,2,3,4` | The passband measured in CW: peak 1549.98–1550.00 nm, FWHM 0.15–0.23 nm (datasheet 0.21 nm), contrast 447–1173×. Ground truth **independent** of the filter's own scale. Also shows the autorange notch on the peak. |
| `sweep_t2_*.npz/.csv/.png` | `sweep --rate 10000 --filter-nm 1550` | First raw capture: **one** broad peak at 5.9065 s (sweep crossing, 59 ms) and **one** narrow peak at 10.1195 s (2.8 ms) = the laser's **return slew**. First direct evidence that the tail is inside the buffer. |
| `sweep_cal_15*_2026-09-01_15-39*` … `15-41*.npz` | `calibrate --filter-positions 1530,…,1565 --rate 5000` | The 8 calibration points. Peak times 1.9622 s (1530 nm) … 8.8909 s (1565 nm). |
| `sweep_cal_1530_…15-36-02.npz` | first, aborted calibration run | That run failed at its second point (USB read stall, report §4.2). Only the 1530 nm point exists. **Not** used in the calibration. |
| `lina_wl_cal.json` | result of the calibration run | The fit. Identical to `../../calibrations/lina_wl_cal_1520.0-1570.0nm_5.00nms.json`, which is the copy the GUI and the scripts look up automatically. |
| `calibration_cal1_*.png` | `plot-cal` | Fit, residuals per degree, error of the uncorrected axis, residual-versus-degree, stability. This is the figure that justifies degree 2. |
| `SUMMARY_lina_wl_cal.png` | `lina_cal_summary.py` | Six-panel overview of all calibration results. Regenerable. |
| `sweep_t100k_*.npz/.png` | `sweep --rate 100000 --cal …` | **The rate-independence test:** calibration measured at 5 kHz, applied here to 100 kHz / 1.2 M samples → marker at 1549.9970 nm, **−3 pm**. The same file also shows the naive axis at −5388 pm. |
| `sweep_nopark_*.npz/.png` | `sweep --park-settle 0` | **Refutes the settling hypothesis:** without parking/settling the marker sits at 5.9067 s versus 5.9064 s with a 5 s park → < 1 pm difference. |
| `repeat_rep1_*.csv` | `repeat --runs 5` | Repeatability over 1.7 min: 0.8 pm rms, 2.0 pm p-p. |
| `repeat_drift_*.csv` | `repeat --runs 6 --interval-s 180` | Drift over 18 min: 0.7 pm rms, **+0.0 pm/min**. No measurable drift → one calibration is enough. |

The `.csv` files are text copies of the `.npz` captures (one row per sample).
They are no longer written automatically (`--csv` enables them), because a
1.2 M capture is ~25 MB as CSV versus ~30 kB as compressed NPZ.

---

## Key results

| Quantity | Value |
|---|---|
| fit | degree 2, **3.1 pm rms** (degree 1: 27.2 pm with a systematic arch) |
| effective sweep rate | 5.1174 nm/s versus 5.0000 commanded |
| trigger dead time | +0.0041 s (essentially zero) |
| sweep in the buffer | t = 0.004 … 9.889 s of 12.00 s (82.4 %) |
| error of the old linear axis | **−5.4 nm** at 1550 nm, **−10.5 nm** with `Buffer time = 2 s` |
| with the calibration | **−3 pm** |
| repeatability / drift | 0.7 pm rms / +0.0 pm/min over 18 min |

---

## Reproducing

```bash
# Passband (the filter must be in the path)
python lina/scripts/lina_sweep_test.py --channels 1,2,3,4 \
    filter-cw --filter-nm 1550 --span 1.5 --step 0.025

# Calibration (~3 min); publishes the result to lina/calibrations/
python lina/scripts/lina_sweep_test.py --channels 2 calibrate \
    --start 1520 --stop 1570 --speed 5 --rate 5000 --trust-filter \
    --filter-positions 1530,1535,1540,1545,1550,1555,1560,1565

# Regenerate the overview figure (no hardware)
python lina/scripts/lina_cal_summary.py --data lina/data/01_wavelength_calibration
```
