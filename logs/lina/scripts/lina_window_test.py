"""
Headless hardware test of the REAL Lina window from AkheLab.

This does not reimplement anything. It builds the actual
``LabGUI.dashboards.instruments.lina_window.LinaWindow`` with the actual
``DeviceRegistry``, on Qt's offscreen platform, then drives it exactly as a
user would: pick the capture device in the combo box, tick the channels, and
press "Response Sweep". The capture therefore runs through the GUI's own
``_CoreDAQSweepWorker`` — the same acquire path, the same arming arithmetic,
the same lambda(t) axis code — and the result dict the GUI would plot is
printed and plotted here instead.

Use it to verify a GUI change against hardware without clicking, and to see
exactly what the window would show.

  python lina/scripts/lina_window_test.py                       # defaults as shipped
  python lina/scripts/lina_window_test.py --start 1530 --stop 1560 --speed 5
  python lina/scripts/lina_window_test.py --channels 1,2,3,4 --no-cal

Requires the real Lina + coreDAQ hardware, and that no other process holds
them (close AkheLab first — one instance per serial id, and the coreDAQ's USB
port is exclusive).
"""

from __future__ import annotations

import argparse
import os
import sys
import time

# Qt must be told to run windowless BEFORE QApplication is created.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLBACKEND", "Agg")

sys.path.insert(0, os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from PyQt5.QtCore import QTimer, QEventLoop
from PyQt5.QtWidgets import QApplication


def wait_until(predicate, timeout_s, app, poll_ms=100, label=""):
    """Spin the Qt event loop until predicate() is true (worker threads need it
    to deliver their signals) or the timeout expires."""
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        app.processEvents(QEventLoop.AllEvents, poll_ms)
        if predicate():
            return True
        time.sleep(poll_ms / 1000.0)
    print(f"[test] TIMEOUT after {timeout_s:.0f} s waiting for {label}")
    return False


def select_capture_device(win, name):
    """Pick an entry in the Response Sweep capture-device combo by device name."""
    cmb = win._cmb_rs_osc
    for i in range(cmb.count()):
        data = cmb.itemData(i)
        if data and data.get("name") == name:
            cmb.setCurrentIndex(i)
            print(f"[test] capture device: '{cmb.itemText(i)}' "
                  f"(type={data.get('_capture_type')}, "
                  f"variant={data.get('variant')})")
            return data
    listing = [f"{cmb.itemText(i)}" for i in range(cmb.count())]
    raise SystemExit(f"'{name}' not in the capture combo. Available: {listing}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--device", default="lina-1")
    p.add_argument("--capture", default="coredaq-1",
                   help="capture device to select in the Response Sweep combo")
    p.add_argument("--channels", default="1",
                   help="coreDAQ channels to tick (1-based)")
    p.add_argument("--start", type=float, default=None, help="override Start (nm)")
    p.add_argument("--stop", type=float, default=None, help="override Stop (nm)")
    p.add_argument("--speed", type=float, default=None, help="override Speed (nm/s)")
    p.add_argument("--rate", type=int, default=None, help="override Rate (Hz)")
    p.add_argument("--samples", type=int, default=None,
                   help="override Samples/ch (default: whatever Auto-fit set)")
    p.add_argument("--no-cal", action="store_true",
                   help="untick 'Use λ(t) calibration' to see the old behaviour")
    p.add_argument("--timeout", type=float, default=300.0,
                   help="seconds to wait for the sweep to finish")
    p.add_argument("--ofdr", action="store_true",
                   help="also press 'OFDR Plot' and save the figure it draws")
    p.add_argument("--out", default="lina_gui_test")
    a = p.parse_args()
    channels = [int(c) for c in a.channels.split(",") if c.strip()]

    app = QApplication(sys.argv)
    from LabGUI.core.device_registry import DeviceRegistry
    from LabGUI.core.window_manager import WindowManager
    from LabGUI.dashboards.instruments.lina_window import LinaWindow

    registry = DeviceRegistry()
    wm = WindowManager()

    print(f"[test] building LinaWindow('{a.device}') offscreen…")
    win = LinaWindow(a.device, registry=registry, wm=wm)
    # Offscreen still needs show(): the window's own code gates the channel
    # selection on QCheckBox.isVisible(), which stays False until shown.
    win.show()

    # The window acquires its device in a worker thread; _device is set only on
    # success, so this also tests the real acquire path.
    if not wait_until(lambda: win._device is not None, 30, app,
                      label="device acquire"):
        raise SystemExit("[test] FAILED: device never acquired — is AkheLab or "
                         "another script holding the EXFO/coreDAQ?")
    print(f"[test] acquired: {win._device} (λ = {win._device.get_wavelength()} nm)")

    cap = select_capture_device(win, a.capture)
    app.processEvents()      # let _on_capture_device_changed rebuild the panel

    # Sweep parameters: leave the shipped defaults unless overridden, so this
    # tests what a user actually gets on opening the window.
    if a.start is not None:
        win._spin_sw_start.setValue(a.start)
    if a.stop is not None:
        win._spin_sw_stop.setValue(a.stop)
    if a.speed is not None:
        win._spin_sw_speed.setValue(a.speed)
    if a.rate is not None:
        win._spin_rs_coredaq_rate.setValue(a.rate)
    if a.samples is not None:
        win._spin_rs_coredaq_samples.setValue(a.samples)
    if a.no_cal:
        win._chk_rs_wl_axis_cal.setChecked(False)

    for i, chk in enumerate(win._chk_rs_pm_ch):
        chk.setChecked((i + 1) in channels)
    app.processEvents()

    print("\n[test] window state as the user would see it:")
    print(f"  sweep            : {win._spin_sw_start.value()} -> "
          f"{win._spin_sw_stop.value()} nm @ {win._spin_sw_speed.value()} nm/s, "
          f"{win._spin_sw_cycles.value()} cycle(s)")
    print(f"  rate / samples   : {win._spin_rs_coredaq_rate.value()} Hz / "
          f"{win._spin_rs_coredaq_samples.value()} per channel")
    print(f"  settle margin    : {win._spin_rs_coredaq_settle.value()} nm")
    print(f"  channels ticked  : "
          f"{[i + 1 for i, c in enumerate(win._chk_rs_pm_ch) if c.isChecked()]}")
    print(f"  λ(t) calibration : {win._chk_rs_wl_axis_cal.isChecked()}  "
          f"[{win._lbl_rs_wl_cal.text()}]")
    print(f"  arm estimate     : {win._lbl_rs_coredaq_duration.text()}")
    print(f"  restore λ after  : {win._chk_rs_restore_wl.isChecked()}")

    # Press the button. Everything after this is the GUI's own code path.
    # Hook the window's OWN handlers rather than one worker's signals: the
    # window may build a further worker to repeat a capture whose USB buffer
    # read failed, and signals connected to the first worker would miss it.
    results = []
    _orig_progress = win._on_rs_progress
    _orig_finished = win._on_rs_finished

    def progress_hook(msg):
        print(f"  [gui] {msg}")
        _orig_progress(msg)

    def finished_hook(r):
        results.append(dict(r))
        _orig_finished(r)

    win._on_rs_progress = progress_hook
    win._on_rs_finished = finished_hook

    print("\n[test] pressing '▶ Response Sweep'…")
    win._start_response_sweep()
    if win._rs_worker is None:
        raise SystemExit(f"[test] FAILED: no worker started — status says "
                         f"'{win._lbl_rs_status.text()}'")

    def settled():
        """A retry may follow an error, so only stop on a FINAL outcome."""
        if win._rs_result is not None:
            return True
        return bool(results) and "error" in results[-1] \
            and win._rs_retry_data is None

    ok = wait_until(settled, a.timeout, app, label="sweep to finish")
    if not ok:
        raise SystemExit("[test] FAILED: sweep did not finish in time")
    result = results[-1] if results else {}
    if len(results) > 1:
        print(f"\n[test] the window repeated the capture {len(results) - 1}x "
              f"after a failed buffer read")

    print("\n[test] result the GUI received:")
    if "error" in result:
        print(f"  ERROR: {result['error']}")
        print(f"  status label: {win._lbl_rs_status.text()}")
        raise SystemExit(1)

    wl = np.asarray(result.get("wavelengths_at_edges", []), dtype=float)
    ds = result.get("downsampled_converted", {})
    print(f"  points           : {wl.size}")
    print(f"  wavelength range : {wl.min():.4f} -> {wl.max():.4f} nm "
          f"(span {wl.max() - wl.min():.4f} nm)")
    print(f"  monotonic        : {bool(np.all(np.diff(wl) > 0) or np.all(np.diff(wl) < 0))}")
    print(f"  y label          : {result.get('y_label')}")
    print(f"  axis note        : {result.get('axis_note', '(none)')}")
    print(f"  λ(t) calibrated  : {result.get('wl_axis_calibrated')}")
    for ch, arr in sorted(ds.items()):
        arr = np.asarray(arr, dtype=float)
        print(f"  {ch:>6}: {arr.size} pts, min {arr.min():.4g}, "
              f"max {arr.max():.4g} {result.get('y_unit', '')}")

    # The window keeps the result for its own plot/save buttons — confirm that
    # the GUI-side handler ran too, not just the worker.
    print(f"  window._rs_result set: {win._rs_result is not None}")
    print(f"  status label         : {win._lbl_rs_status.text()}")

    # Plot exactly what the window plots.
    os.makedirs(a.out, exist_ok=True)
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 5))
    for ch, arr in sorted(ds.items()):
        ax.plot(wl, np.asarray(arr, dtype=float), lw=0.6, label=str(ch))
    ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel(result.get("y_label", ""))
    ax.set_title(f"LinaWindow Response Sweep (headless) — "
                 f"{win._spin_sw_start.value():.0f}→{win._spin_sw_stop.value():.0f} nm "
                 f"@ {win._spin_sw_speed.value():.0f} nm/s, "
                 f"{win._spin_rs_coredaq_rate.value()} Hz"
                 f"{result.get('axis_note', '')}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    stamp = time.strftime("%Y-%m-%d_%H-%M-%S")
    png = os.path.join(a.out, f"gui_response_sweep_{stamp}.png")
    fig.savefig(png, dpi=130)
    print(f"\n[out] {png}")

    npz = os.path.join(a.out, f"gui_response_sweep_{stamp}.npz")
    np.savez_compressed(npz, wavelengths_nm=wl,
                        **{str(k): np.asarray(v, dtype=float) for k, v in ds.items()})
    print(f"[out] {npz}")

    if a.ofdr:
        print("\n[test] pressing 'OFDR Plot'…")
        if not win._btn_rs_ofdr.isEnabled():
            raise SystemExit("[test] FAILED: the OFDR button is disabled — "
                             "capture all four channels (--channels 1,2,3,4)")
        win._show_ofdr_plot()
        app.processEvents()
        ofdr_win = getattr(win, "_ofdr_window", None)
        if ofdr_win is None:
            raise SystemExit(f"[test] FAILED: no OFDR window opened — status "
                             f"says '{win._lbl_rs_status.text()}'")
        print(f"[test] OFDR status: {win._lbl_rs_status.text()}")
        # Save the figure the GUI itself drew, so the result is inspectable.
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
        canvases = ofdr_win.findChildren(FigureCanvasQTAgg)
        if canvases:
            png = os.path.join(a.out, f"gui_ofdr_{stamp}.png")
            canvases[0].figure.savefig(png, dpi=130)
            print(f"[out] {png}")
        else:
            raise SystemExit("[test] FAILED: the OFDR window has no plot canvas")
        ofdr_win.close()

    win.close()
    print("\n[test] PASSED — the GUI path produced a usable sweep.")


if __name__ == "__main__":
    main()
