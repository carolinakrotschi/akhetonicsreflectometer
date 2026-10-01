"""
VOA attenuation series: N identical Lina sweeps, one per attenuator setting.

WHY THIS EXISTS
---------------
A variable optical attenuator in the path is set by hand -- it has a knob, no
bus, and nothing reads its position back. So a series over attenuation is the
same capture repeated with a human turning the knob in between, which is
exactly the part a script cannot do and exactly the part that ruins a series
when it is done from memory: a capture taken before the knob was turned, or a
file named after the setting it did not have.

This runs the sequence instead. Before every capture it stops, says which
number in the series is coming, shows the CW power the coreDAQ sees right now
(so the knob can be turned to a visible level), and waits. Whatever is typed
at that prompt is recorded with the capture as the VOA setting and becomes the
label of that trace in the final figure. Then it sweeps, saves, and moves on.
At the end it plots all captures stacked on a shared axis.

The capture itself is the one ``lina_chip_measure.py`` runs, imported from it
rather than copied: same park + settle, same trigger-armed coreDAQ read with
the USB retries, same Auto-fit sample count as the GUI, same wavelength axis
(lambda(t) calibration -> aux fringe -> honest linear), same NPZ layout. So
the OFDR and analysis tools read these files unchanged.

Every capture is the WIDE window, 1505-1625 nm at 5 nm/s. That window is
outside the TOF1550 filter's range, so it has no lambda(t) calibration and its
wavelength axis comes from the aux fringe -- the same axis the GUI builds for
it, and ``meta["wavelength_axis_source"]`` says so per file.

The laser is NOT switched off between runs. Nothing is unplugged in a VOA
series -- the knob is turned in-line -- and leaving emission on is what makes
the live power readout at the prompt possible. Emission goes off at the end
and on any interrupt.

OUTPUT
------
One NPZ + quicklook PNG per capture in ``--out``::

    <chip>_<YYYYmmdd-HHMMSS>_fiber<NN>_fullnmrange_<i>.npz      i = 1..N

The timestamp is the SESSION's, shared by every file of the series, so the ten
files form one visibly contiguous block that differs only in the trailing
index. ``meta`` carries the sweep parameters, the crop the wavelength axis
belongs to, and the series fields ``series_id``, ``run_index``, ``run_total``
and ``voa_setting`` (what was typed at the prompt).

Two figures close the session, named after the series:

    <series>_stacked.png   one panel per run, shared wavelength AND power axis
    <series>_levels.png    mean power per run, per channel, against run index

The stacked figure shares its y axis on purpose: with each panel autoscaled
the attenuation -- the whole point of the series -- disappears. ``--free-y``
turns that off. Its power scale follows the series: up to ~13 dB of spread it
is linear with a min/max envelope; past that it goes logarithmic and switches
to peak hold, because on a shared linear axis the attenuated runs are flat
lines on zero. ``--log-y`` / ``--linear-y`` force it, and the figure title and
the console both say which was used.

EXAMPLES
--------
  python lina/scripts/lina_voa_series.py --dry-run
  python lina/scripts/lina_voa_series.py --fiber 48 --runs 10
  python lina/scripts/lina_voa_series.py --fiber 48 --start-index 8 --runs 3
  python lina/scripts/lina_voa_series.py --replot "lina/raw_data/2680_20260917-*_fiber48_fullnmrange_*.npz"

``--replot`` touches no hardware: it rebuilds both figures from saved NPZs, so
a session interrupted at run 7 still gets its figure out of the seven captures
it did take.

NOTE: close the LabGUI Lina/coreDAQ windows first -- LabDevice allows only one
open instance per serial id, and the coreDAQ USB-CDC port is exclusive.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time
from datetime import datetime

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from lina_sweep_test import (open_devices, close_devices, run_sweep, read_power_mw,
                             lina_coredaq_name)
# The capture-side helpers are lina_chip_measure's; importing them keeps the
# two scripts on one wavelength axis and one file layout by construction.
from lina_chip_measure import (
    WINDOWS, autofit_samples, build_axis, save_capture, plot_capture,
    summarize, emission_off,
)


# The wide preset, and the name the files carry for it.
WINDOW_KEY = "1505-1625nm"
RANGE_NAME = "fullnmrange"


def monitor_power(laser, pm, channels, wl_nm, settle_s=1.0):
    """CW power per channel at ``wl_nm``, so the knob can be turned to a level.

    Returns None if it cannot be read -- a missing readout must not stop a
    series, it only costs the operator a number on screen.
    """
    try:
        laser.enable("ON")
        laser.set_wavelength(wl_nm)
        time.sleep(settle_s)
        return read_power_mw(pm, channels)
    except Exception as e:
        print(f"[mon] no live power readout: {e}")
        return None


def prompt_voa(index, last, powers, wl_nm, first):
    """Stop for the knob. What is typed becomes this capture's VOA setting."""
    print("\n" + "=" * 66)
    if first:
        print(f"  RUN {index}/{last} -- set the VOA to the FIRST position")
    else:
        print(f"  RUN {index}/{last} -- TURN THE VOA KNOB to the next position")
    if powers:
        levels = "  ".join(f"Ch{ch} {p:.4g} mW" for ch, p in sorted(powers.items()))
        print(f"  live CW power @ {wl_nm:.1f} nm:  {levels}")
    print("  Laser emission is ON (nothing is unplugged in a VOA series).")
    print("  Type the setting (e.g. 3.0dB, or the knob position) and press")
    print("  Enter; Enter alone records no setting. Ctrl-C stops the session.")
    print("=" * 66)
    try:
        return input("  VOA setting: ").strip()
    except EOFError:
        raise KeyboardInterrupt


def envelope(y, wl, max_pts=3000):
    """Block min/max of ``y`` against ``wl`` -- see plot_capture's comment.

    Every nth sample aliases the ~13 kHz aux fringe into a band of hash; the
    per-block envelope keeps the extremes a capture is judged by.
    """
    step = max(int(np.ceil(wl.size / max_pts)), 1)
    n = (wl.size // step) * step
    blocks = np.asarray(y)[:n].reshape(-1, step)
    return wl[:n].reshape(-1, step)[:, 0], blocks.min(axis=1), blocks.max(axis=1)


def load_series(paths):
    """Read the saved captures of a series, in run order, cropped to the axis."""
    runs = []
    for path in paths:
        with np.load(path, allow_pickle=False) as z:
            meta = json.loads(str(z["meta"]))
            lo, hi = meta.get("crop", [0, None])
            runs.append({"path": path, "meta": meta, "wl_nm": z["wl_nm"],
                         "traces": {int(k[2:]): z[k][lo:hi] for k in z.files
                                    if k.startswith("ch")}})
    runs.sort(key=lambda r: (r["meta"].get("run_index", 0), r["path"]))
    return runs


def trace_label(run):
    i = run["meta"].get("run_index", "?")
    setting = run["meta"].get("voa_setting") or ""
    return f"#{i}" + (f"   VOA {setting}" if setting else "")


#: Span of run levels above which a shared LINEAR power axis stops being
#: readable -- past ~13 dB the attenuated traces sit on the zero line and the
#: shape that is being compared is gone. Measured against the figure, not a
#: rule of thumb: at 3.5 dB per knob step, run 5 onwards is already a flat line.
AUTO_LOG_RATIO = 20.0


def plot_stacked(runs, channel, png, free_y=False, log_y=None):
    """All captures underneath each other, on one shared wavelength axis.

    ``log_y=None`` picks the scale from the data: a series that spans more
    than ~13 dB gets a log axis, because on a shared linear one everything
    below the first couple of settings is a flat line at zero. The choice is
    printed and written into the title, since the two figures look alike.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = [r for r in runs if channel in r["traces"]]
    if not runs:
        print(f"[plot] no capture has Ch{channel} -- stacked figure skipped")
        return None

    means = [abs(float(np.mean(r["traces"][channel]))) for r in runs]
    span = max(means) / max(min(means), 1e-12)
    if log_y is None:
        log_y = (not free_y) and span > AUTO_LOG_RATIO
        if log_y:
            print(f"[plot] runs span {10 * np.log10(span):.1f} dB -- log power "
                  f"axis, peak hold (--linear-y overrides)")

    envs = [envelope(r["traces"][channel], r["wl_nm"]) for r in runs]
    if log_y:
        # On a LINEAR coreDAQ an attenuated trace straddles zero -- its zero
        # has noise on both sides -- and a log axis cannot draw the lower
        # envelope at all. A LOG coreDAQ (both Mk2 units) never goes negative,
        # but pins dark samples to its clamp floor, so the lower envelope is
        # then a flat line at the floor rather than the trace's own minimum. Clipping it to a floor does not fail visibly, it fills the panel
        # with a solid block that hides the trace. So on log, show the block
        # MAXIMUM as a line: peak hold, the same display lina_ofdr_test uses,
        # where peaks keep their exact height and the noise collapses to a line.
        tops = np.concatenate([hi for _, _, hi in envs])
        tops = tops[tops > 0]
        ylim = (float(tops.min()) * 0.7, float(tops.max()) * 1.4)

    fig, axes = plt.subplots(len(runs), 1, figsize=(9, 1.2 + 1.25 * len(runs)),
                             sharex=True, sharey=not free_y, squeeze=False)
    for ax, run, (wl, lo_env, hi_env) in zip(axes[:, 0], runs, envs):
        if log_y:
            ax.plot(wl, hi_env, lw=0.5, color="C0")
            ax.set_yscale("log")
            ax.set_ylim(*ylim)
        else:
            ax.fill_between(wl, lo_env, hi_env, lw=0.4, color="C0")
        ax.set_ylabel("mW", fontsize=8)
        ax.grid(alpha=0.3)          # decades only; the minor log grid is hash
        ax.tick_params(labelsize=7)
        ax.text(0.995, 0.92, trace_label(run), transform=ax.transAxes,
                ha="right", va="top", fontsize=8,
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
    axes[-1, 0].set_xlabel("Wavelength (nm)")
    axes[0, 0].set_title(
        f"{runs[0]['meta'].get('series_id', '')}   Ch{channel}   "
        f"{len(runs)} VOA settings, "
        f"{'independent' if free_y else 'shared'} "
        + ("log power axis, peak hold" if log_y
           else "linear power axis, min/max envelope"), fontsize=9)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)
    print(f"[out] {png}")
    return png


def plot_levels(runs, png):
    """Mean power per run: what the knob actually did, as a number per step."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    chans = sorted({ch for r in runs for ch in r["traces"]})
    idx = [r["meta"].get("run_index", i + 1) for i, r in enumerate(runs)]
    fig, ax = plt.subplots(figsize=(8, 4))
    for ch in chans:
        means = [float(np.mean(r["traces"][ch])) if ch in r["traces"] else np.nan
                 for r in runs]
        ax.plot(idx, means, "o-", label=f"Ch{ch}", color=f"C{ch - 1}")
    ax.set_yscale("log")
    ax.set_xticks(idx)
    settings = [r["meta"].get("voa_setting") or "" for r in runs]
    if any(settings):
        ax.set_xticklabels([f"{i}\n{s}" for i, s in zip(idx, settings)], fontsize=7)
    ax.set_xlabel("Run index")
    ax.set_ylabel("Mean power over the sweep (mW)")
    ax.set_title(f"{runs[0]['meta'].get('series_id', '')} -- level per VOA setting",
                 fontsize=9)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)
    print(f"[out] {png}")
    return png


def report_table(runs, channel):
    """Per-run level, and its ratio to the first run -- the attenuation steps."""
    rows = [r for r in runs if channel in r["traces"]]
    if not rows:
        return
    print(f"\n  Ch{channel}")
    print("  run  VOA setting       mean (mW)       min        max   vs #1")
    print("  " + "-" * 62)
    ref = None
    for run in rows:
        seg = run["traces"][channel]
        mean = float(np.mean(seg))
        if ref is None:
            ref = mean
        rel = (10 * np.log10(mean / ref) if ref > 0 and mean > 0 else float("nan"))
        print(f"  {str(run['meta'].get('run_index', '?')):>3}  "
              f"{(run['meta'].get('voa_setting') or '-'):<15} "
              f"{mean:>11.4g} {seg.min():>10.4g} {seg.max():>10.4g}  "
              f"{rel:+7.2f} dB")


def main():
    # describe()/describe_aux() return the GUI's labels, which contain lambda
    # and an en-dash; a cp1252 console raises on those and would abort a
    # running session mid-capture.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(errors="replace")
        except Exception:
            pass

    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--chip", default="2680", help="chip id, first filename field")
    p.add_argument("--fiber", default="48", help="channel number in the filename")
    p.add_argument("--runs", type=int, default=10, help="captures in the series")
    p.add_argument("--start-index", type=int, default=1,
                   help="first run number -- to continue an interrupted series")
    p.add_argument("--out", default=os.path.join("lina", "raw_data"))
    p.add_argument("--channels", default="1,2,3,4",
                   help="coreDAQ channels; the aux axis needs 2 and 4")
    p.add_argument("--plot-channel", type=int, default=1,
                   help="channel shown in the stacked figure")
    p.add_argument("--free-y", action="store_true",
                   help="autoscale each stacked panel (hides the attenuation)")
    scale = p.add_mutually_exclusive_group()
    scale.add_argument("--log-y", dest="log_y", action="store_true", default=None,
                       help="force a log power axis on the stacked figure")
    scale.add_argument("--linear-y", dest="log_y", action="store_false",
                       help="force a linear power axis (the default is chosen "
                            "from how far the series spans)")
    p.add_argument("--rate", type=int, default=100_000, help="coreDAQ sample rate (Hz)")
    p.add_argument("--exfo", default="exfo-1")
    p.add_argument("--coredaq", default=lina_coredaq_name(),
                   help="default: the coreDAQ paired with lina-1 in OBR_config.json")
    p.add_argument("--park-settle", type=float, default=5.0,
                   help="seconds at the start wavelength before arming")
    p.add_argument("--monitor-nm", type=float, default=1550.0,
                   help="CW wavelength for the live power shown at the prompt")
    p.add_argument("--no-monitor", action="store_true",
                   help="do not read CW power before each prompt")
    p.add_argument("--no-plot", action="store_true", help="no per-capture quicklook")
    p.add_argument("--replot", default=None,
                   help="glob of saved NPZs -- rebuild the figures, no hardware")
    p.add_argument("--dry-run", action="store_true",
                   help="print the plan and the filenames, touch no hardware")
    a = p.parse_args()

    channels = [int(c) for c in a.channels.split(",")]
    win = WINDOWS[WINDOW_KEY]

    if a.replot:
        paths = sorted(glob.glob(a.replot))
        if not paths:
            p.error(f"--replot matched no files: {a.replot}")
        runs = load_series(paths)
        print(f"[replot] {len(runs)} capture(s)")
        stem = (runs[0]["meta"].get("series_id")
                or os.path.splitext(os.path.basename(paths[0]))[0])
        out_dir = os.path.dirname(os.path.abspath(paths[0]))
        report_table(runs, a.plot_channel)
        plot_stacked(runs, a.plot_channel,
                     os.path.join(out_dir, stem + "_stacked.png"), a.free_y,
                     a.log_y)
        plot_levels(runs, os.path.join(out_dir, stem + "_levels.png"))
        return

    n = autofit_samples(win, a.rate)
    sweep_s = abs(win["stop"] - win["start"]) / win["speed"]
    indices = list(range(a.start_index, a.start_index + a.runs))
    session = f"{datetime.now():%Y%m%d-%H%M%S}"
    series_id = f"{a.chip}_{session}_fiber{a.fiber}_{RANGE_NAME}"

    print(f"\nVOA series: chip {a.chip}, fibre {a.fiber}, {RANGE_NAME}")
    print(f"  window   : {win['start']}->{win['stop']} nm @ {win['speed']} nm/s, "
          f"sweep {sweep_s:.1f} s")
    print(f"  coreDAQ  : {channels} @ {a.rate} Hz on {a.coredaq}, arm {n} "
          f"samples/ch = {n / a.rate:.2f} s, "
          f"{n * 2 * len(channels) / 1e6:.1f} MB per read")
    print(f"  runs     : {a.runs} (index {indices[0]}..{indices[-1]})")
    print(f"  out      : {os.path.abspath(a.out)}")
    print(f"  ~ {a.runs * (a.park_settle + n / a.rate) / 60:.1f} min of capture "
          f"in total, plus the buffer reads and your knob turns\n")

    if a.dry_run:
        for i in indices:
            print("  " + os.path.join(a.out, f"{series_id}_{i}.npz"))
        print("  " + os.path.join(a.out, f"{series_id}_stacked.png"))
        print("  " + os.path.join(a.out, f"{series_id}_levels.png"))
        print("\n(dry run -- no hardware was touched)")
        return

    os.makedirs(a.out, exist_ok=True)
    laser = pm = None
    written = []
    try:
        laser, pm, _ = open_devices(a.exfo, a.coredaq, None)
        for pos, i in enumerate(indices):
            powers = None if a.no_monitor else monitor_power(
                laser, pm, channels, a.monitor_nm)
            setting = prompt_voa(i, indices[-1], powers, a.monitor_nm, pos == 0)

            stem = f"{series_id}_{i}"
            print(f"\n--- run {i}/{indices[-1]} -> {stem}"
                  + (f"   [VOA {setting}]" if setting else ""))
            res = run_sweep(laser, pm,
                            start_nm=win["start"], stop_nm=win["stop"],
                            speed_nm_s=win["speed"], channels=channels,
                            rate_hz=a.rate, n_samples=n,
                            park_settle_s=a.park_settle, restore_wl=True)
            wl_nm, lo, hi, prov = build_axis(res)
            summarize(res, wl_nm, lo, hi)
            path = save_capture(res, wl_nm, lo, hi, prov,
                                os.path.join(a.out, stem + ".npz"),
                                extra={"chip": a.chip, "label": RANGE_NAME,
                                       "fiber": str(a.fiber), "window": WINDOW_KEY,
                                       "series_id": series_id, "run_index": i,
                                       "run_total": indices[-1],
                                       "voa_setting": setting,
                                       "voa_monitor_nm": a.monitor_nm,
                                       "voa_monitor_mw": powers})
            written.append(path)
            if not a.no_plot:
                plot_capture(res, wl_nm, lo, hi,
                             f"{stem}   VOA {setting or '-'}   "
                             f"[{prov['wavelength_axis_source']}]",
                             os.path.join(a.out, stem + ".png"))
    except KeyboardInterrupt:
        print("\n[stop] interrupted -- plotting what was captured")
    finally:
        if laser is not None:
            emission_off(laser)
        close_devices(laser, pm, None)

    print(f"\n{len(written)} capture(s) written to {os.path.abspath(a.out)}")
    if not written:
        return
    runs = load_series(written)
    report_table(runs, a.plot_channel)
    plot_stacked(runs, a.plot_channel,
                 os.path.join(a.out, series_id + "_stacked.png"), a.free_y,
                 a.log_y)
    plot_levels(runs, os.path.join(a.out, series_id + "_levels.png"))


if __name__ == "__main__":
    main()
