# Data: headless GUI hardware tests (2026-09-01, 16:43–17:58)

Output of `lina/scripts/lina_window_test.py`. Each file is the result of the
**real** Lina window from AkheLab being driven on real hardware without
clicking: the actual `LinaWindow` with the actual `DeviceRegistry` is built on
Qt's offscreen platform, the capture device is selected in its combo box, the
channels are ticked, and "Response Sweep" is pressed.

Purpose: verify a GUI change against hardware, and see exactly what the window
would show. The capture therefore runs through the GUI's own
`_CoreDAQSweepWorker` — the same acquire path, arming arithmetic and λ(t) axis
code a user gets.

Full reasoning: `../../report.md`, §4 and §10–§11.

---

## The consideration

The scripts were written first because the GUI could not be debugged
interactively from here. Once the GUI carried the same features, it needed to
be tested on hardware too — but a Qt window cannot be clicked from a terminal
session. Driving the real window offscreen tests the GUI's own code rather than
a parallel re-implementation, which is the only way a GUI test proves anything.

---

## Files

| File | What it shows |
|---|---|
| `gui_response_sweep_*.png/.npz` | The sweep the window produced: power versus the calibrated wavelength axis, plus the per-channel arrays. |
| `gui_ofdr_*.png` | The figure the new **"OFDR Plot"** button draws (aux above, measurement below). |

---

## What these runs established

1. **The λ(t) axis works in the GUI.** 988,536 points spanning exactly
   1520.0000 → 1570.0000 nm, monotonic, with the status line reporting
   `✓ λ(t) deg 2, 3 pm rms`.

2. **The premature-read bug and its fix.** With 1.2 M frames at 100 kHz
   (12.00 s window against a 10.81 s sweep) the read stalled at
   4,799,488 / 4,800,000 bytes and destroyed the arming. After adding the wait
   for the armed window, the same setting read cleanly.

3. **The retry works.** Three runs at 1.2 M × 2 channels: one failed twice and
   succeeded on the third attempt, one failed once, one succeeded immediately —
   all three ended with usable data. At four channels (8.4 MB) one run
   exhausted the original 2-retry budget, which is why it is now 4.

4. **A real bug in my own module wiring, found by this test.** The reference
   channel pair was passed the string `"aux"` in `combine_pair`'s *mode*
   parameter (which only accepts `auto`/`always`/`never`), so it was silently
   **not** balanced-subtracted — only Ch2 was used. The aux phase was
   consequently unusable (3 % single-tone, 18 % non-monotonic). Fixed.

5. **An unresolved problem (report §11.1).** Even after that fix, captures
   taken through the GUI yield an aux phase with ~5 % too few fringes
   (121,636 versus 128,011 in a script capture minutes later on the same
   hardware), which stretches the distance axis by the same factor: the
   3 m fibre reads 2.846 m instead of 3.026 m. The evaluation itself was
   verified identical to the script's by running the same `evaluate()` on
   script data. Cause not yet found — **use the script for real measurements.**

---

## Reproducing

```bash
# Defaults as shipped, all four channels, plus the OFDR button
python lina/scripts/lina_window_test.py --channels 1,2,3,4 --ofdr

# Force a specific arming to reproduce the read behaviour
python lina/scripts/lina_window_test.py --channels 1,2 --samples 1200000
```

Requires the hardware to be free — close AkheLab first (LabDevice allows one
instance per serial id, and the coreDAQ's USB port is exclusive).
