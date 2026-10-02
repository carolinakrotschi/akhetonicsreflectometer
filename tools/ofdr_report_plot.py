#!/usr/bin/env python3
"""
Prototype of the new OFDR report for the Lina / Odyssey "OFDR Plot" button.

Layout (top to bottom):
  0. text box: DUT port with uncertainty, fibre end and length, then every
     peak found, position in mm and level in dB
  1. full distance axis, 0 .. Nyquist
  2. same spectrum with 0 at the DUT port, so the axis reads fibre length;
     fibre end marked
  3. zoom onto the fibre end

The fibre length is NOT known beforehand. It is found as:
  - DUT port P   = strongest peak in --port-window (0.3-0.8 m; since the
                   20 dB coupler came out on 2026-10-02 the port is the
                   strongest reflection, no offset to calibrate)
  - fibre end E  = strongest peak beyond P + 10 mm
  - L            = E - P, measured inside the same scan (the common scale
                   drift between runs cancels, logs/2026-10-01.md)

The spectrum comes from exactly the Odyssey pipeline (lina_ofdr_test.process:
lambda(t) calibration crop, aux-phase resampling, Kaiser beta 10), with no
MZI offset, so the axis is the absolute one.

Usage:
    <AKHExperiment>/.venv/Scripts/python.exe tools/ofdr_report_plot.py <npz> \
        --out results/2026-10-02/ofdr_layout_example/run9
"""

import argparse
import json
import os
import sys
from types import SimpleNamespace

import numpy as np

AKH = os.environ.get("AKHEXPERIMENT",
                     r"C:\Users\Carolina.krotsch\Documents\Akhelab\AKHExperiment")
sys.path.insert(0, AKH)
sys.path.insert(0, os.path.join(AKH, "lina", "scripts"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lina_ofdr_test import process                                # noqa: E402
from lina.analysis.lina_ofdr import N_GROUP_SMF28                 # noqa: E402
from lina.analysis.lina_ofdr_report import (                      # noqa: E402
    analyse, draw, summary_lines, PORT_WINDOW_M)


def load_npz(path):
    d = np.load(path, allow_pickle=True)
    meta = d["meta"].item() if "meta" in d.files else {}
    if isinstance(meta, (str, bytes)):
        meta = json.loads(meta)
    traces = {int(k[2:]): np.asarray(d[k], float)
              for k in d.files if k.startswith("ch")}
    rate = float(meta.get("rate_hz") or meta["detector"]["sample_rate_hz"])
    return traces, rate, meta


def spectrum(traces, rate, meta, n_group):
    a = SimpleNamespace(meas=[1, 3], ref=[2, 4], balance="auto", no_cal=False,
                        start=1520.0, stop=1570.0, speed=5.0, highpass_hz=500.0,
                        no_aux=False, n_group=n_group, beta=10.0, pad_factor=4,
                        convention="reflection", offset_m=0.0, mirror=False,
                        n_peaks=8, min_db=-45.0)
    res, *_ = process(traces, rate, meta, a)
    return res


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("npz")
    p.add_argument("--out", required=True, help="output prefix (no extension)")
    p.add_argument("--n-group", type=float, default=N_GROUP_SMF28)
    p.add_argument("--port-window", nargs=2, type=float, default=list(PORT_WINDOW_M),
                   metavar=("LO", "HI"), help="window for the DUT port peak, m")
    a = p.parse_args()

    traces, rate, meta = load_npz(a.npz)
    spec = spectrum(traces, rate, meta, a.n_group)
    rep = analyse(spec, port_window_m=tuple(a.port_window))
    fig = plt.figure(figsize=(12, 16.5))
    draw(fig, rep, title=f"OFDR report  —  {os.path.basename(a.npz)}",
         subtitle=(f"Δz {rep['dz']*1e6:.1f} µm, range {rep['zmax']:.3f} m, "
                   f"{rate/1e3:.0f} kHz, dB rel. strongest peak; ± = fit ⊕ "
                   f"scale 1e-4·z (n_g uncertainty not included)"))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    fig.savefig(a.out + ".png", dpi=130, bbox_inches="tight")
    print("\n".join(summary_lines(rep)))
    print("wrote", a.out + ".png")


if __name__ == "__main__":
    main()
