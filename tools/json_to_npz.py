"""Convert a raw EXFO/coreDAQ JSON scan to a compressed .npz.

Why: the raw JSON scans run ~155 MB each, past GitHub's hard 100 MB
per-file limit, so they cannot be pushed. Stored as a compressed .npz the
same data is ~7 MB -- a factor of 22 -- and small enough to keep in the
repository alongside the analysis.

The conversion is LOSSLESS: the four channel arrays are kept as float64,
exactly the values the JSON carries. Only the JSON text overhead goes
away (every sample is written as a decimal string in the JSON, which is
what makes those files so large).

The wavelength axis is stored as start/stop/count when it is a linear ramp
(the uncalibrated fallback axis), and in full otherwise -- which is the
normal case now: the lambda(t)-calibrated and the aux-referenced axes are
not linear. The script checks rather than assumes.

The power unit is kept as written (the Lina GUI auto-scales its columns to
[uW] / [mW] / [W]) and recorded under "unit", so the round trip stays
bit-identical; scan_io.load_scan converts to mW on reading. The header goes
along whole, so a Mk2 capture keeps its "detector" block (frontend, serial,
rate, clamp fractions, R(lambda) correction) -- the analysis needs it to stay
comparable with the Mk1 scans.

The result loads with scan_io.load_scan, which every analysis tool uses
(process_reflectogram_aux.load and plot_voa_series.reflectogram included).

Usage
    python tools/json_to_npz.py raw_data/2026-09-14-09-03_fiber4.json
    python tools/json_to_npz.py raw_data/2026-09-1*.json --verify
"""

import argparse
import glob
import json
import os
import re
import sys

import numpy as np

CHANNELS = ["Ch1 [mW]", "Ch2 [mW]", "Ch3 [mW]", "Ch4 [mW]"]


def convert(src, out=None, verify=False):
    with open(src) as fh:
        doc = json.load(fh)
    e = doc["data"][0]
    # "ChN [unit]" in any power unit -- older scans use [W], the GUI writes
    # whatever it auto-scaled to. All four must share one unit.
    found = {}
    for k in e:
        m = re.match(r"^Ch(\d) \[([^\]]+)\]$", k)
        if m:
            found[int(m.group(1))] = (k, m.group(2))
    missing = [c for c in CHANNELS if int(c[2]) not in found]
    if missing:
        sys.exit("%s: missing channels %s" % (src, missing))
    units = {found[n][1] for n in (1, 2, 3, 4)}
    if len(units) != 1:
        sys.exit("%s: channels in different units %s" % (src, sorted(units)))
    keys, unit = [found[n][0] for n in (1, 2, 3, 4)], units.pop()

    ch = [np.asarray(e[k], float) for k in keys]
    n = len(ch[0])
    meta = {}

    wl = e.get("Wavelength [nm]")
    if wl is not None:
        wl = np.asarray(wl, float)
        lin = np.linspace(wl[0], wl[-1], len(wl))
        dev = float(np.max(np.abs(wl - lin)))
        meta["lam_start_nm"] = np.array(wl[0])
        meta["lam_stop_nm"] = np.array(wl[-1])
        meta["lam_max_dev_nm"] = np.array(dev)
        if dev > 1e-6:
            meta["wavelength_nm"] = wl        # nicht linear -> voll behalten
            print("   wavelength axis is not linear (max dev %.2e nm) -> kept in full" % dev)

    hdr = doc.get("header", {})
    for k, v in hdr.items():
        if isinstance(v, str):
            meta["hdr_" + k] = np.array(v)
    # Der ganze Header zusaetzlich als JSON-String unter "meta" -- so wie ihn
    # der Lina-Rekorder in seinen .npz ablegt. Die hdr_*-Schluessel oben
    # verlieren jede geschachtelte Angabe, und genau dort steht
    # wavelength_axis_aux.tau_aux_implied_s, ohne das die aux-referenzierten
    # Werkzeuge (process_reflectogram_aux, plot_voa_series) ein --tau-aux-ns
    # von Hand brauchen.
    if hdr:
        meta["meta"] = np.array(json.dumps(hdr))
    meta["unit"] = np.array(unit)
    meta["n_points"] = np.array(n)
    meta["source_json"] = np.array(os.path.basename(src))

    out = out or os.path.splitext(src)[0] + ".npz"
    np.savez_compressed(out, ch1=ch[0], ch2=ch[1], ch3=ch[2], ch4=ch[3], **meta)

    s_in, s_out = os.path.getsize(src), os.path.getsize(out)
    print("%s -> %s   %.1f MB -> %.1f MB  (%.1f%%, factor %.1f)"
          % (os.path.basename(src), os.path.basename(out),
             s_in / 1e6, s_out / 1e6, 100 * s_out / s_in, s_in / s_out))

    if verify:
        d = np.load(out)
        for i, c in enumerate(ch):
            back = d["ch%d" % (i + 1)]
            if not np.array_equal(back, c):
                sys.exit("   VERIFY FAILED on ch%d" % (i + 1))
        print("   verified: all four channels bit-identical after round trip")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scans", nargs="+", help="one or more .json scans (globs ok)")
    ap.add_argument("--verify", action="store_true",
                    help="read the .npz back and check the arrays are identical")
    a = ap.parse_args()

    files = []
    for pat in a.scans:
        files.extend(sorted(glob.glob(pat)) or [pat])
    tot_in = tot_out = 0
    for f in files:
        out = convert(f, verify=a.verify)
        tot_in += os.path.getsize(f)
        tot_out += os.path.getsize(out)
    if len(files) > 1:
        print("\ntotal: %.1f MB -> %.1f MB  (factor %.1f)"
              % (tot_in / 1e6, tot_out / 1e6, tot_in / tot_out))


if __name__ == "__main__":
    main()
