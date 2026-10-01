"""
Guided Lina measurement session: several fibres x several sweep windows.

WHY THIS EXISTS
---------------
The Lina GUI runs ONE Response Sweep per button press and saves it under a
name typed into a dialog. Measuring a chip is that same capture repeated over
fibres and sweep windows, where the only manual step is re-patching the fibre.
Done by hand that is a dozen clicks per capture and a filename typed each
time; this runs the whole sequence instead, pauses only for the fibre change,
and names every file the same way.

What it reproduces from the GUI (LabGUI/dashboards/instruments/lina_window.py,
``_CoreDAQSweepWorker``, continuous path):

  * capture      -- lina_sweep_test.run_sweep: park + 5 s settle, arm the
                   coreDAQ on the EXFO sweep trigger, wait the armed window
                   out, read, with the re-arm retries the USB transfer needs.
  * samples/ch   -- the GUI's Auto-fit: ceil(sweep_s * 1.05 * rate), i.e.
                   1 050 000 for 1520-1570 nm and 2 520 000 for 1505-1625 nm
                   at 100 kHz and 5 nm/s.
  * axis         -- a measured lambda(t) calibration when one exists for the
                   configuration, otherwise the aux fringe, otherwise an
                   honestly-labelled linear axis (see build_axis).
  * restore      -- the laser returns to its pre-sweep CW wavelength.

The GUI's "lambda-cal correction" box (responsivity R(lambda_set)/R(lambda))
was a no-op on the Mk1 LINEAR coreDAQ. On the Mk2 LOG units it is ACTIVE in the
GUI, but this script deliberately stays raw: what is saved is the power the
meter reported at its set wavelength, so a few-% responsivity tilt across the
band is NOT removed here. Apply it in analysis if absolute spectra matter
(``CoreDAQ.responsivity_a_per_w``).

Emission is switched OFF around every fibre change, so no fibre is ever
unplugged with light in it.

OUTPUT
------
One NPZ per capture in ``--out``::

    <chip>_<label>_<YYYYmmdd-HHMMSS>_fiber<NN>_<window>.npz

holding the RAW uncropped traces (``ch1``..``ch4``, mW), the wavelength axis
of the sweep (``wl_nm``), and a JSON ``meta`` with the sweep parameters, the
crop indices the axis corresponds to (``meta["crop"] = [lo, hi]``, so
``ch1[lo:hi]`` lines up with ``wl_nm``) and how the axis was built. The layout
is the one ``lina_sweep_test.load_raw`` reads, so the OFDR tools work on these
files unchanged::

    python lina/scripts/lina_ofdr_test.py peaks --raw lina/raw_data/<file>.npz

A quicklook PNG is written next to each capture unless ``--no-plot``.

HOW A SESSION RUNS
------------------
It asks which channel(s) to measure (one number, or several), then for each of
them in turn: prompts for the patch change with the emission off, runs every
sweep window, and writes the files. When the batch is done it asks again, so a
session continues for as long as there are channels to measure; an empty
answer ends it. ``--fibers 44,48`` skips the first question, ``--connected 44``
the first patch prompt.

EXAMPLES
--------
  python lina/scripts/lina_chip_measure.py --dry-run
  python lina/scripts/lina_chip_measure.py --chip 2680 --label ligentechhi
  python lina/scripts/lina_chip_measure.py --fibers 44,48 --connected 44

NOTE: close the LabGUI Lina/coreDAQ windows first -- LabDevice allows only one
open instance per serial id, and the coreDAQ USB-CDC port is exclusive.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from lina_sweep_test import open_devices, close_devices, run_sweep, lina_coredaq_name
from lina.analysis.lina_wl_cal import (
    find_calibration, wavelength_axis, sweep_window, describe as describe_cal,
    aux_wavelength_axis, describe_aux,
)


# The two sweep windows of the Lina GUI's presets ("Short" and "Full"), at the
# 5 nm/s the calibration was measured at. Only 1520-1570 has a lambda(t)
# calibration; 1505-1625 is outside the TOF1550 filter's range and is carried
# by the aux fringe instead (lina/calibrations/README.md).
WINDOWS = {
    "1520-1570nm": {"start": 1520.0, "stop": 1570.0, "speed": 5.0},
    "1505-1625nm": {"start": 1505.0, "stop": 1625.0, "speed": 5.0},
}


def autofit_samples(win: dict, rate_hz: int) -> int:
    """Samples/channel the GUI's Auto-fit button would put in the box."""
    sweep_s = abs(win["stop"] - win["start"]) / win["speed"]
    return int(math.ceil(sweep_s * 1.05 * rate_hz))


def build_axis(res, log=print):
    """Wavelength axis for a raw capture, as the GUI's Response Sweep builds it.

    Returns ``(wl_nm, lo, hi, provenance)`` -- ``wl_nm`` covers the sweep only
    and lines up with ``trace[lo:hi]``. The provenance dict records which of
    the three sources produced it, because all three look alike once saved.
    """
    traces, rate, n_pts = res["traces_mw"], res["rate_hz"], res["n_pts"]

    cal = find_calibration(res["requested_start_nm"], res["requested_stop_nm"],
                           res["requested_speed_nm_s"])
    if cal is not None:
        lo, hi = sweep_window(cal, n_pts, rate)
        wl = wavelength_axis(cal, n_pts, rate)[lo:hi]
        log(f"[axis] lambda(t) calibration: {describe_cal(cal)}")
        return wl, lo, hi, {
            "wavelength_axis_source": "lambda_t",
            "wavelength_axis_note": describe_cal(cal),
            "wavelength_axis_calibration": os.path.basename(cal.get("_path", "")),
        }

    try:
        wl, _cropped, diag = aux_wavelength_axis(
            traces, res["actual_start_nm"], res["actual_stop_nm"], rate,
            log=lambda m: log("[axis] " + m))
        lo, hi = (int(diag["window"][0]), int(diag["window"][1]))
        log(f"[axis] {describe_aux(diag)}")
        return wl, lo, hi, {
            "wavelength_axis_source": "aux",
            "wavelength_axis_note": describe_aux(diag),
            "wavelength_axis_aux": {
                k: diag[k] for k in ("fringes", "tau_aux_implied_s",
                                     "tau_aux_reference_s", "tau_aux_deviation",
                                     "start_nm", "stop_nm", "tail_note")
                if k in diag},
        }
    except Exception as e:
        log(f"[axis] aux axis unavailable ({e}) -- falling back to the linear axis")
        wl = np.linspace(res["actual_start_nm"], res["actual_stop_nm"], n_pts)
        return wl, 0, n_pts, {
            "wavelength_axis_source": "linear",
            "wavelength_axis_note": f"uncalibrated linear axis ({e})",
        }


def save_capture(res, wl_nm, lo, hi, provenance, path, extra=None):
    meta = {k: v for k, v in res.items() if k != "traces_mw"}
    meta.update(provenance)
    meta.update(extra or {})
    meta["crop"] = [int(lo), int(hi)]
    np.savez_compressed(path, meta=json.dumps(meta), wl_nm=np.asarray(wl_nm),
                        **{f"ch{ch}": tr for ch, tr in res["traces_mw"].items()})
    print(f"[out] {path}  ({os.path.getsize(path) / 1e6:.1f} MB)")
    return path


def plot_capture(res, wl_nm, lo, hi, title, png, max_pts=3000):
    """Quicklook: per-block min/max envelope, NOT every nth sample.

    The aux channels fringe at ~13 kHz, so taking every nth sample aliases
    them into a band of hash that hides whether the capture is any good --
    which is the only question this figure exists to answer.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    step = max(int(np.ceil(wl_nm.size / max_pts)), 1)
    n = (wl_nm.size // step) * step
    wl_blocks = wl_nm[:n].reshape(-1, step)[:, 0]
    chans = sorted(res["traces_mw"])
    fig, axes = plt.subplots(len(chans), 1, figsize=(9, 2 + 1.6 * len(chans)),
                             sharex=True, squeeze=False)
    for ax, ch in zip(axes[:, 0], chans):
        blocks = np.asarray(res["traces_mw"][ch])[lo:hi][:n].reshape(-1, step)
        ax.fill_between(wl_blocks, blocks.min(axis=1), blocks.max(axis=1),
                        lw=0.4, color=f"C{ch - 1}")
        ax.set_ylabel(f"Ch{ch} (mW)", fontsize=8)
        ax.grid(alpha=0.3)
    axes[-1, 0].set_xlabel("Wavelength (nm)")
    axes[0, 0].set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)
    print(f"[out] {png}")


def summarize(res, wl_nm, lo, hi):
    print(f"[cap] axis {wl_nm[0]:.3f} -> {wl_nm[-1]:.3f} nm over "
          f"{wl_nm.size} samples (of {res['n_pts']} captured)")
    for ch, tr in sorted(res["traces_mw"].items()):
        seg = np.asarray(tr)[lo:hi]
        print(f"      Ch{ch}: {seg.min():.4g} ... {seg.max():.4g} mW "
              f"(mean {seg.mean():.4g}, p-p {seg.max() - seg.min():.4g})")


def emission_off(laser):
    """No light in a fibre that is about to be unplugged."""
    if laser is None:
        return
    try:
        laser.enable("OFF")
    except Exception as e:
        print(f"[safe] could not switch emission off: {e}")


def wait_for_fiber(laser, fiber):
    """Emission off, ask for the patch change, wait for the user."""
    emission_off(laser)
    print("\n" + "=" * 66)
    print(f"  CONNECT CHANNEL {fiber}".ljust(64))
    print("  Laser emission is OFF. Press Enter when it is patched "
          "(Ctrl-C to stop).")
    print("=" * 66)
    try:
        input()
    except EOFError:
        raise KeyboardInterrupt


def ask_fibers():
    """Which channel(s) to measure next. An empty answer ends the session."""
    print("\n" + "=" * 66)
    print("  WHICH CHANNEL(S)?  e.g. 44   or   44,48")
    print("  Enter alone ends the session.")
    print("=" * 66)
    while True:
        try:
            answer = input("  channels: ")
        except EOFError:
            return []
        fibers = [f for f in answer.replace(";", ",").replace(" ", ",").split(",")
                  if f]
        if all(re.fullmatch(r"[\w.-]+", f) for f in fibers):
            return fibers
        print("  [?] numbers, comma separated -- they go into the file name")


def main():
    # The shared describe()/describe_aux() helpers return the GUI's labels,
    # which contain lambda and en-dash; a cp1252 console raises on those and
    # would abort a running session mid-capture.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except Exception:
            pass

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chip", default="2680", help="chip id, first filename field")
    p.add_argument("--label", default="ligentechhi", help="second filename field")
    p.add_argument("--fibers", default="",
                   help="channel numbers to start with; the script asks when "
                        "this is empty, and asks again after every batch")
    p.add_argument("--connected", default=None,
                   help="channel already patched -- skips its connect prompt")
    p.add_argument("--windows", default=",".join(WINDOWS),
                   help="sweep windows per fibre, in order (default: both)")
    p.add_argument("--out", default=os.path.join("lina", "raw_data"))
    p.add_argument("--channels", default="1,2,3,4",
                   help="coreDAQ channels; the aux axis needs 2 and 4")
    p.add_argument("--rate", type=int, default=100_000, help="coreDAQ sample rate (Hz)")
    p.add_argument("--exfo", default="exfo-1")
    p.add_argument("--coredaq", default=lina_coredaq_name(),
                   help="default: the coreDAQ paired with lina-1 in OBR_config.json")
    p.add_argument("--park-settle", type=float, default=5.0,
                   help="seconds at the start wavelength before arming")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan and the filenames, touch no hardware")
    a = p.parse_args()

    fibers = [f.strip() for f in a.fibers.split(",") if f.strip()]
    windows = [w.strip() for w in a.windows.split(",") if w.strip()]
    unknown = [w for w in windows if w not in WINDOWS]
    if unknown:
        p.error(f"unknown window(s) {unknown}; known: {list(WINDOWS)}")
    channels = [int(c) for c in a.channels.split(",")]

    print(f"\nSession: chip {a.chip} / {a.label}")
    print(f"  channels : {', '.join(fibers) if fibers else '(asked for below)'}"
          + (f"  ({a.connected} already patched)" if a.connected else ""))
    print(f"  out      : {os.path.abspath(a.out)}")
    print(f"  coreDAQ  : {channels} @ {a.rate} Hz on {a.coredaq}")
    total_s = 0.0
    for w in windows:
        win = WINDOWS[w]
        n = autofit_samples(win, a.rate)
        sweep_s = abs(win["stop"] - win["start"]) / win["speed"]
        total_s += a.park_settle + n / a.rate
        print(f"  window {w}: {win['start']}->{win['stop']} nm @ "
              f"{win['speed']} nm/s, sweep {sweep_s:.1f} s, arm {n} samples/ch "
              f"= {n / a.rate:.2f} s, {n * 2 * len(channels) / 1e6:.1f} MB to read")
    print(f"  ~ {total_s / 60:.1f} min of capture per channel, plus the buffer "
          f"reads and the patch change\n")

    if a.dry_run:
        when = datetime.now()
        for fiber in fibers or ["<NN>"]:
            for w in windows:
                print("  " + os.path.join(
                    a.out, f"{a.chip}_{a.label}_{when:%Y%m%d-%H%M%S}"
                           f"_fiber{fiber}_{w}.npz"))
        print("\n(dry run -- no hardware was touched)")
        return

    os.makedirs(a.out, exist_ok=True)
    laser = pm = None
    written = []
    patched = str(a.connected) if a.connected else None
    try:
        laser, pm, _ = open_devices(a.exfo, a.coredaq, None)
        while True:
            if not fibers:
                # Between batches, so the session can go on without restarting
                # the script (and without re-opening the devices).
                emission_off(laser)
                fibers = ask_fibers()
                if not fibers:
                    break

            for fiber in fibers:
                if fiber != patched:
                    wait_for_fiber(laser, fiber)
                else:
                    print(f"\n[chan] {fiber} is already patched -- starting")
                patched = fiber

                for w in windows:
                    win = WINDOWS[w]
                    stem = (f"{a.chip}_{a.label}_{datetime.now():%Y%m%d-%H%M%S}"
                            f"_fiber{fiber}_{w}")
                    print(f"\n--- channel {fiber}, {w} -> {stem}")
                    res = run_sweep(laser, pm,
                                    start_nm=win["start"], stop_nm=win["stop"],
                                    speed_nm_s=win["speed"], channels=channels,
                                    rate_hz=a.rate,
                                    n_samples=autofit_samples(win, a.rate),
                                    park_settle_s=a.park_settle, restore_wl=True)
                    wl_nm, lo, hi, prov = build_axis(res)
                    summarize(res, wl_nm, lo, hi)
                    path = save_capture(res, wl_nm, lo, hi, prov,
                                        os.path.join(a.out, stem + ".npz"),
                                        extra={"chip": a.chip, "label": a.label,
                                               "fiber": str(fiber), "window": w})
                    written.append(path)
                    if not a.no_plot:
                        plot_capture(res, wl_nm, lo, hi,
                                     f"{stem}  [{prov['wavelength_axis_source']}]",
                                     os.path.join(a.out, stem + ".png"))
                print(f"\n[chan] {fiber} done ({len(windows)} capture(s))")
            fibers = []
    except KeyboardInterrupt:
        print("\n[stop] interrupted -- leaving the laser off")
    finally:
        if laser is not None:
            try:
                laser.enable("OFF")
            except Exception:
                pass
        close_devices(laser, pm, None)

    print(f"\n{len(written)} capture(s) written to {os.path.abspath(a.out)}")
    for path in written:
        print("  " + os.path.basename(path))


if __name__ == "__main__":
    main()
