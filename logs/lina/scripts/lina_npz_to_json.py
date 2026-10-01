"""
Convert Lina captures (NPZ) to the house JSON data format.

WHY THIS EXISTS
---------------
A session writes NPZ (lina_chip_measure.py): compact, and what the analysis
scripts read. The rest of AKHExperiment stores measurements as
``{"header": {...}, "data": [{label: value, ...}]}`` JSON, which is also what
the Lina GUI's "Save Data" button writes -- so a capture that has to be read
by anything other than these scripts needs that format. This converts one
into the other without touching the bench again.

WHAT COMES OUT
--------------
One JSON per NPZ, same stem::

    {"header": {...everything the NPZ recorded...},
     "data": [{"Device Description": {...}, "Wavelength [nm]": [...],
               "Ch1 [mW]": [...], ...}]}

The columns are cropped to ``meta["crop"]``, i.e. the part of the buffer the
wavelength axis actually covers (the trigger dead time and the laser's return
slew are dropped) -- the same extent the GUI saves, so all columns are the
same length.

SIZE. This is text: a 2.4 M-sample capture is ~115 MB of JSON against 32 MB of
NPZ, and takes a while to write. Numbers are formatted at a precision well
past what the instrument resolves (6 significant figures for power, 1 fm for
wavelength), not at float64's 17 digits, which would double the file for no
information.

EXAMPLES
--------
  python lina/scripts/lina_npz_to_json.py                  # all of lina/raw_data
  python lina/scripts/lina_npz_to_json.py lina/raw_data/2680_*_fiber44_*.npz
  python lina/scripts/lina_npz_to_json.py --out L:/Measurements/2680 <file>.npz
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

# Keys of the NPZ meta that describe the run rather than the capture; they go
# into the Device Description, where this repo keeps the bench context.
DESCRIPTION_KEYS = ("chip", "label", "fiber", "window")

WL_FMT = "%.6f"      # 1 fm; the calibration itself is good to ~3 pm
POWER_FMT = "%.6g"   # the coreDAQ is 16-bit, so this is ~2 digits to spare


def _write_array(fh, arr, fmt, chunk=200_000):
    """Stream one array as a JSON list, without building it in memory first."""
    arr = np.asarray(arr, dtype=float)
    n_bad = int(np.count_nonzero(~np.isfinite(arr)))
    fh.write("[")
    for i in range(0, arr.size, chunk):
        block = arr[i:i + chunk]
        if i:
            fh.write(",")
        # A non-finite value would be written as "nan"/"inf", which no JSON
        # parser accepts -- null keeps the file readable and the gap visible.
        if n_bad:
            fh.write(",".join(fmt % v if np.isfinite(v) else "null" for v in block))
        else:
            fh.write(",".join(fmt % v for v in block))
    fh.write("]")
    return n_bad


def write_trace_json(path, columns, header, description):
    """Write the AKHExperiment data format, streaming the big columns."""
    with open(path, "w") as fh:
        fh.write('{\n"header": ')
        json.dump(header, fh, indent=4, sort_keys=True)
        fh.write(',\n"data": [\n{\n"Device Description": ')
        json.dump(description, fh, indent=4, sort_keys=True)
        bad = {}
        for label, (arr, fmt) in columns.items():
            fh.write(',\n' + json.dumps(label) + ': ')
            n_bad = _write_array(fh, arr, fmt)
            if n_bad:
                bad[label] = n_bad
        fh.write("\n}\n]\n}\n")
    return bad


def convert(npz_path, out_dir=None, overwrite=False):
    d = np.load(npz_path, allow_pickle=False)
    meta = json.loads(str(d["meta"]))
    out_dir = out_dir or os.path.dirname(os.path.abspath(npz_path))
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir,
                       os.path.splitext(os.path.basename(npz_path))[0] + ".json")
    if os.path.exists(out) and not overwrite:
        print(f"[skip] {out} exists (--overwrite to replace)")
        return None

    lo, hi = meta.get("crop", [0, int(meta.get("n_pts", 0))])
    columns = {}
    if "wl_nm" in d.files:
        columns["Wavelength [nm]"] = (d["wl_nm"], WL_FMT)
    else:
        # A capture saved before the axis was stored with it: no wavelengths
        # rather than a guessed linear axis.
        print(f"[warn] {os.path.basename(npz_path)} has no wl_nm array")
    for key in sorted(k for k in d.files if k.startswith("ch")):
        columns[f"Ch{key[2:]} [mW]"] = (d[key][lo:hi], POWER_FMT)

    lengths = {k: int(np.asarray(v).size) for k, (v, _) in columns.items()}
    if len(set(lengths.values())) > 1:
        raise RuntimeError(f"column lengths disagree in {npz_path}: {lengths}")

    header = dict(meta)
    header["instrument"] = "Lina"
    header["source_file"] = os.path.basename(npz_path)
    description = {k: meta[k] for k in DESCRIPTION_KEYS if k in meta}
    description["capture"] = (
        f"{meta.get('actual_start_nm', '?')}-{meta.get('actual_stop_nm', '?')} nm "
        f"@ {meta.get('actual_speed_nm_s', '?')} nm/s, "
        f"{meta.get('rate_hz', '?')} Hz, axis from "
        f"{meta.get('wavelength_axis_source', 'unknown')}")

    bad = write_trace_json(out, columns, header, description)
    n = next(iter(lengths.values()), 0)
    print(f"[out] {out}  ({os.path.getsize(out) / 1e6:.1f} MB, "
          f"{len(columns)} columns x {n} rows)")
    for label, count in bad.items():
        print(f"      {label}: {count} non-finite value(s) written as null")
    return out


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paths", nargs="*",
                   default=[os.path.join("lina", "raw_data")],
                   help="NPZ files, globs, or directories (default lina/raw_data)")
    p.add_argument("--out", default=None,
                   help="output directory (default: next to each NPZ)")
    p.add_argument("--overwrite", action="store_true")
    a = p.parse_args()

    files = []
    for path in a.paths:
        if os.path.isdir(path):
            files.extend(sorted(glob.glob(os.path.join(path, "*.npz"))))
        else:
            files.extend(sorted(glob.glob(path)) or [path])
    if not files:
        p.error(f"no NPZ files found in {a.paths}")

    print(f"{len(files)} file(s) to convert\n")
    for path in files:
        print(f"--- {os.path.basename(path)}")
        convert(path, a.out, a.overwrite)


if __name__ == "__main__":
    main()
