# Changes on `carolina-ofdr` vs. upstream `main`

This document describes everything that has been changed on the
`carolina-ofdr` branch relative to the original AkheLab codebase
(`Akhetonics/AKHExperiment`, `main`). It exists because this account does
not have push access to the upstream repo, so these changes currently only
live on this branch / on a personal fork.

**Relocated:** the plain EXFO laser window (`exfo_window.py`) has been
reverted to its original, pre-this-branch state. Everything described below
(Response Sweep's coreDAQ-LIN support, the Y-axis toggle, manual LIN
arming, settle margin) now lives instead in a new, separate "OBR"
instrument type — `Interface/OBR.py` (`OBR` superclass, `Lina` concrete
class: an EXFO + coreDAQ-LIN pairing) and
`LabGUI/dashboards/instruments/lina_window.py` (`LinaWindow`, a copy of the
old EXFO window preselected to that pairing). The plain "Lasers" category
now shows the unmodified upstream EXFO window again; the new "OBR" category
shows Lina. The technical write-up below still describes the original work
accurately — read `exfo_window.py`/`Interface/Lasers.py` as
`lina_window.py`/`Interface/OBR.py` throughout. `Interface/Lasers.py` itself
was **not** reverted (§4's `trig_out()` fix is a generic, backward-compatible
driver improvement that Lina also depends on).

All changes so far touch the EXFO tunable-laser "Response Sweep" workflow
and its coreDAQ power-meter integration, plus one supporting fix in the
EXFO SCPI driver. Three files were originally affected:

- `Interface/Lasers.py` (kept, unchanged by the relocation)
- `LabGUI/dashboards/instruments/exfo_window.py` (reverted; logic moved to
  `lina_window.py`)
- `.vscode/settings.json`

Commits (oldest first): `ad7b8f0`, `f4adc55`, `4eaf141`, `965f64f` (the last
one removed the anti-saturation-range checkbox mentioned in §5, keeping
only the settle-margin field).

---

## 1. Made coreDAQ-1 (LINEAR variant) work in Response Sweep at all

**Problem:** Response Sweep's coreDAQ capture path only supported the
custom-stepper coreDAQ variant (e.g. coredaq-2, LOG frontend), which arms
the device in STEPPED trigger mode — one ADC sample per EXFO trigger edge,
re-arming after each one. That mode requires coreDAQ firmware v4.3+ and
only exists on the stepper-capable units. The standard LINEAR coreDAQ (e.g.
coredaq-1) doesn't support per-edge STEPPED arming and doesn't need
firmware v4.3, but the code always used STEPPED mode regardless of which
variant was selected.

**Fix (`exfo_window.py`, `_CoreDAQSweepWorker`):** the worker now branches
on the coreDAQ's PM_config `variant` field (`LINEAR` vs `LOG`) via a new
`stepped` flag:

- `stepped=True` (LOG/custom-stepper units, unchanged behavior): one EXFO
  trigger edge per wavelength step, coreDAQ re-arms after each sample, the
  real per-edge wavelength table is read back from the laser afterward
  (`get_pulse_wavelengths`).
- `stepped=False` (LINEAR units, e.g. coredaq-1, new): a single trigger
  edge starts the capture; the coreDAQ then free-runs continuously at its
  own sample rate for the whole sweep. There is no per-edge wavelength
  table in this mode, so the wavelength axis is reconstructed linearly over
  the sweep's actual (readback) start/stop wavelength.

### Two firmware-interaction bugs found and fixed while wiring this up

1. **`finalize_capture()` wiping a completed capture.** `finalize_capture()`
   sends an ACQ STOP to the device. When the capture had already reached
   the firmware's DONE state (4) on its own, sending STOP anyway reset the
   stored frame count to 0 instead of being a harmless no-op — silently
   destroying a buffer that had fully filled. Fixed by querying the
   acquisition state first and only calling `finalize_capture()` when the
   capture is *not* already DONE.
2. **`captured_frames()` silently reading as 0 on older firmware.** This
   query (`FRAMES?`) is itself gated behind firmware v4.3+. On older
   firmware it fails and the PowerMeters.py wrapper swallows the exception,
   returning `-1` — which looked identical to "0 frames captured" even
   though the SDRAM buffer was genuinely full. Fixed for the LINEAR/
   continuous path: since a DONE state in that mode means "all `n_frames`
   armed were collected" by construction, the code now uses the originally
   armed `n_frames` directly instead of re-querying the device.

---

## 2. Y-axis linear/log toggle for the Response Sweep plot

**`exfo_window.py`:** added two buttons under the Response Sweep plot to
switch the Y-axis between linear and logarithmic scale. The selection
persists across re-plots (re-running a sweep, or replaying the same sweep
after a Y-axis change, keeps the chosen scale instead of resetting to
linear).

---

## 3. Manual sample-count / sample-rate arming for LINEAR coreDAQ

**Problem:** in continuous/LINEAR mode there's no per-edge trigger table,
so the frame count has to be derived some other way. Originally it was
always derived from the EXFO's reported sweep duration × sample rate —
with no way to request an exact, fixed number of samples per channel (e.g.
for a fixed ~1,000,000-point-per-channel OFDR capture, independent of how
long the sweep happens to take).

**Fix (`exfo_window.py`):** when a LINEAR coreDAQ (e.g. coredaq-1) is
selected as the capture device, the UI now exposes direct **Samples/channel**
and **Sample-rate (Hz)** controls, plus a live capture-duration estimate
that warns if it looks mismatched against the actual configured sweep time.
Device buffer sizing note surfaced in the tooltip: 32 MB SDRAM ≈ 16.7M
points on 1 channel, ≈ 4.2M points/channel across 4 channels — so the
default of ~1,000,000 points/channel is comfortably within range.

(A short-lived standalone script, `carolina/ofdr.py`, was written first to
prototype this single-trigger-edge + continuous-capture approach outside
the GUI, then removed once the same functionality was integrated directly
into Response Sweep's manual-arming path above.)

---

## 4. Real SCPI trigger in/out control for SWEEP mode (`Interface/Lasers.py`)

**Problem:** `EXFO.trig_out()` only supported `mode="PULSE"` in a way that
sent SCPI commands (`SOUR:SWE:SYNC:*`, `SOUR:WAV:SWE:SYNC:SAMP`,
`SOUR:OPT:BCOM:STAT`) that do not appear anywhere in the T200S SCPI Command
Reference and almost certainly originate from a different EXFO laser
model. `mode="SWEEP"` existed but only toggled a trigger-type SCPI value
with no way to actually configure trigger source (external vs. immediate)
or the trigger-out BNC signal independently.

**Fix:** `trig_out()` gained explicit parameters for the SWEEP path, each
mapped to a command that *is* documented in the T200S SCPI Command
Reference:

- `trigger_in: bool` → `TRIGger[:SEQuence]:SOURce EXTernal|IMMediate`
  (p.133) — `True` arms the sweep to start from the external Trig In BNC,
  `False` starts it manually/by software.
- `trigger_out: bool` → `OUTPut:SYNChronization[:STATe] ON|OFF` (p.95) —
  the Trig Out BNC goes high for the duration of the sweep.
- `out_type` / `sampling_pm` (Window vs. Pulse trigger-out shape, and
  per-step sampling interval): **not supported remotely** — no matching
  command exists in the SCPI reference for this instrument. These
  parameters are accepted for API compatibility but are logged as a
  warning and not sent to the device; they must be set on the laser's own
  touchscreen.
- After configuring, `trig_out()` now reads back the device's error queue
  (`SYST:ERR?`, via a new `get_error()` helper) and raises if the laser
  rejected the configuration, instead of failing silently.

The original `mode="PULSE"` code path (used by the STEPPED coreDAQ sweep
workers) was left functionally unchanged — its SCPI commands are flagged
as unverified/likely wrong in a docstring, but not touched, so as not to
break the existing STEPPED-mode callers without a documented replacement.

**`exfo_window.py` UI side:** the "Trigger Output" panel gained real
**Trigger in** / **Trigger out** checkboxes and **Window** / **Pulse**
radio buttons (mirroring the laser's own front-panel "Sweep triggers"
options), wired to the corrected `trig_out()` call above; the Window/Pulse
choice is visibly marked as "set on the laser touchscreen — not remotely
controllable" wherever it can't actually be sent over SCPI.

---

## 5. coreDAQ-1 saturation in the first ~10 nm of a sweep

**Symptom:** on coreDAQ-1 (LINEAR), the first ~10 nm of every Response
Sweep showed the channel pinned in over-range / saturation.

**Investigation:**

- Added an **"Anti-saturation range"** checkbox (LINEAR coreDAQ only):
  when checked, disables autorange and forces every armed channel to
  manual TIA range index 0 (5 mW full-scale, the least sensitive available
  range) before arming.
- Tested with the checkbox on — the channel still saturated. This ruled
  out a range/gain misconfiguration: range index 0 is already the most
  power-tolerant range the LINEAR frontend has, so if that still trips, no
  range/gain setting can fix it — the over-range trip point is confirmed
  (via the `py-coreDAQ` documentation) to be a **fixed hardware threshold**
  on the ADC voltage (`|V| > 4.2 V`), not adjustable in Watts or Volts.
- With no EDFA or variable optical attenuator in the optical path (direct
  laser → coreDAQ), the real cause was identified as the **laser's own
  output-power settling transient** right as a continuous sweep starts
  moving — a fixed *wavelength* distance of ringing/overshoot, not a
  power-level or gain problem.

**Fix — "Settle margin (nm)" (LINEAR coreDAQ only, default 10 nm):**
the *physical* EXFO sweep is started this many nm before the user's
requested start wavelength (in whichever direction is "before" for that
sweep's direction), so the laser's settling transient plays out before it
reaches the wavelength range the user actually asked for. The coreDAQ
still captures through that pre-roll — sized automatically into the
existing frame-count calculation — and after capture, the pre-roll samples
(anything before the originally requested start wavelength) are cropped
back out of the result. The returned sweep therefore still spans exactly
the user's requested start/stop, but without the settling artifact at the
front.

Both the anti-saturation-range checkbox and the settle-margin field are
independent, opt-in controls (both default to their current/no-op state
except settle margin, which defaults to 10 nm) — they don't change
behavior for the STEPPED (LOG-variant) coreDAQ path at all.

---

## 6. Misc

- `.vscode/settings.json`: default Python environment manager switched
  from `ms-python.python:venv` to `ms-python.python:system` (editor-local
  preference, not functional code).
