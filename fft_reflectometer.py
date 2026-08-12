#!/usr/bin/env python3
"""Simple script to compute FFT from a JSON measurement file in the test data folder."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError as exc:
    raise ImportError("matplotlib is required for plotting. Install it with 'pip install matplotlib'.") from exc


DATA_DIR = Path(__file__).resolve().parent / "test data"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
REFRACTIVE_INDEX_DEFAULT = 1.468


def load_json_file(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def find_signal_arrays(data: dict[str, Any]) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    if not isinstance(data, dict):
        raise ValueError("Expected JSON object as a dictionary.")

    if "data" not in data or not data["data"]:
        raise ValueError("JSON does not contain a valid 'data' entry.")

    first_block = data["data"][0]
    if not isinstance(first_block, dict):
        raise ValueError("The first entry in 'data' is not a dictionary.")

    arrays: dict[str, np.ndarray] = {}
    wavelength: np.ndarray | None = None

    for key, value in first_block.items():
        if isinstance(value, list) and value and all(isinstance(x, (int, float)) for x in value):
            arrays[key] = np.array(value, dtype=float)
            if key.lower().startswith("wavelength"):
                wavelength = arrays[key]

    if wavelength is None:
        raise ValueError("No wavelength data found in JSON (e.g. 'Wavelength [nm]').")

    return wavelength, arrays


def make_uniform_wavenumber(wavelength: np.ndarray, intensity: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if wavelength.ndim != 1 or intensity.ndim != 1 or wavelength.size != intensity.size:
        raise ValueError("Wavelength and intensity arrays must be 1D and the same length.")

    k = 1.0 / wavelength
    if k[0] > k[-1]:
        k = k[::-1]
        intensity = intensity[::-1]

    k_uniform = np.linspace(k[0], k[-1], wavelength.size)
    intensity_uniform = np.interp(k_uniform, k, intensity)
    return k_uniform, intensity_uniform


def compute_fft(signal: np.ndarray, sample_spacing: float) -> tuple[np.ndarray, np.ndarray]:
    n = signal.size
    spectrum = np.fft.rfft(signal)
    freq = np.fft.rfftfreq(n, d=abs(sample_spacing))
    amplitude = np.abs(spectrum)
    return freq, amplitude


def pick_channel(arrays: dict[str, np.ndarray], name: str | None) -> tuple[str, np.ndarray]:
    if name:
        if name not in arrays:
            raise ValueError(f"Channel '{name}' not found. Available channels: {', '.join(sorted(arrays))}")
        return name, arrays[name]

    preferred = [key for key in arrays if key.lower().startswith("ch1") or key.lower().startswith("ch2")]
    if preferred:
        return preferred[0], arrays[preferred[0]]

    alternatives = [key for key in arrays if not key.lower().startswith("wavelength")]
    if not alternatives:
        raise ValueError("No suitable signal channels found in JSON.")
    return alternatives[0], arrays[alternatives[0]]


def save_fft_result(output_dir: Path, file_stem: str, freq: np.ndarray, amplitude: np.ndarray) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{file_stem}.npz"
    np.savez_compressed(output_path, freq=freq, amplitude=amplitude)
    return output_path


def plot_fft_result(output_dir: Path, file_stem: str, freq: np.ndarray, amplitude: np.ndarray, axis_label: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / f"{file_stem}_fft.png"

    plt.figure(figsize=(10, 6))
    plt.plot(freq, amplitude, color="blue", linewidth=1)
    plt.title(f"FFT of {file_stem}")
    plt.xlabel(axis_label)
    plt.ylabel("Amplitude")
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200)
    plt.close()
    return plot_path


def calculate_fiber_length(freq: np.ndarray, amplitude: np.ndarray, refractive_index: float) -> tuple[float, float]:
    if freq.size < 2:
        raise ValueError("Frequency array too small to calculate fiber length.")

    positive_mask = freq > 0
    if not np.any(positive_mask):
        raise ValueError("No positive frequency components found for length calculation.")

    search_freq = freq[positive_mask]
    search_amplitude = amplitude[positive_mask]
    peak_index = np.argmax(search_amplitude)
    peak_frequency = search_freq[peak_index]

    length_nm = peak_frequency / refractive_index
    length_m = length_nm * 1e-9
    return peak_frequency, length_m


def plot_distance_result(output_dir: Path, file_stem: str, freq: np.ndarray, amplitude: np.ndarray, refractive_index: float, peak_frequency: float) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_path = output_dir / f"{file_stem}_distance.png"

    positive_mask = freq > 0
    if not np.any(positive_mask):
        raise ValueError("No positive frequency components available for distance plotting.")

    distance_m = (freq[positive_mask] / refractive_index) * 1e-9
    amplitude_plot = amplitude[positive_mask]
    peak_distance_m = (peak_frequency / refractive_index) * 1e-9

    plt.figure(figsize=(10, 6))
    plt.plot(distance_m, amplitude_plot, color="green", linewidth=1)
    plt.axvline(peak_distance_m, color="red", linestyle="--", label="Detected reflector")
    plt.title(f"Distance-domain reflectogram of {file_stem}")
    plt.xlabel("Distance [m]")
    plt.ylabel("Amplitude")
    plt.xlim(0, 3)
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig(plot_path, dpi=200)
    plt.close()
    return plot_path


def process_file(path: Path, channel: str | None, use_wavenumber: bool, output_dir: Path, refractive_index: float) -> None:
    json_data = load_json_file(path)
    wavelength, arrays = find_signal_arrays(json_data)
    channel_name, signal = pick_channel(arrays, channel)

    if wavelength.size != signal.size:
        raise ValueError("Wavelength and signal arrays have different lengths.")

    if use_wavenumber:
        x, y = make_uniform_wavenumber(wavelength, signal)
        sample_spacing = x[1] - x[0]
        freq, amplitude = compute_fft(y, sample_spacing)
        axis_label = "Wavenumber frequency [1/nm]"
    else:
        delta_lambda = np.diff(wavelength)
        sample_spacing = float(abs(np.mean(delta_lambda)))
        freq, amplitude = compute_fft(signal, sample_spacing)
        axis_label = "cycles per nm"

    print(f"Processing: {path.name}")
    print(f"  Channel: {channel_name}")
    print(f"  Samples: {signal.size}")
    print(f"  Sample spacing: {sample_spacing:.6e}")
    print(f"  Axis label: {axis_label}")

    top_indices = np.argsort(amplitude)[-10:][::-1]
    print("  Top 10 FFT peaks:")
    for idx in top_indices:
        print(f"    freq={freq[idx]:.6e}, amplitude={amplitude[idx]:.6e}")

    peak_frequency, fiber_length_m = calculate_fiber_length(freq, amplitude, refractive_index)
    print(f"  Detected fiber peak frequency: {peak_frequency:.6e} cycles per nm")
    print(f"  Estimated fiber length: {fiber_length_m:.6e} m")
    print("  Note: 0 m appears in the distance plot because the FFT includes the DC/zero-frequency component at the origin.")

    saved_path = save_fft_result(output_dir, path.stem, freq, amplitude)
    plot_path = plot_fft_result(output_dir, path.stem, freq, amplitude, axis_label)
    distance_plot_path = plot_distance_result(output_dir, path.stem, freq, amplitude, refractive_index, peak_frequency)
    print(f"  Saved FFT result to: {saved_path}")
    print(f"  Saved FFT plot to: {plot_path}")
    print(f"  Saved distance plot to: {distance_plot_path}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute FFT for JSON measurements in the test data folder.")
    parser.add_argument("filename", nargs="?", default=None, help="JSON filename inside the test data folder. If omitted, all JSON files in test data are processed.")
    parser.add_argument("--channel", help="Data channel to use, e.g. 'Ch1 [mW]' or 'Ch2 [mW]'.")
    parser.add_argument("--use-wavenumber", action="store_true", help="Interpolate to a uniform wavenumber grid before FFT.")
    parser.add_argument("--refractive-index", type=float, default=REFRACTIVE_INDEX_DEFAULT, help=f"Refractive index of the fiber. Default: {REFRACTIVE_INDEX_DEFAULT}")
    parser.add_argument("--output-dir", default="results", help="Output folder for FFT results. Default: results")
    args = parser.parse_args()

    if not DATA_DIR.exists():
        raise FileNotFoundError(f"Data folder not found: {DATA_DIR}")

    if args.filename:
        paths = [DATA_DIR / args.filename]
    else:
        paths = sorted(DATA_DIR.glob("*.json"))

    if not paths:
        raise FileNotFoundError(f"No JSON files found in {DATA_DIR}.")

    results_folder = Path(__file__).resolve().parent / args.output_dir

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        process_file(path, args.channel, args.use_wavenumber, results_folder, args.refractive_index)


if __name__ == "__main__":
    main()
