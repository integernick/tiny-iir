#!/usr/bin/env python3

import sys
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import sosfreqz

def parse_sos_coeffs(lines):
    """
    Parse lines of the form:
      1 -0.000000000000000 -1.000000000000000 1 -0.000000000000000 0.887369786202070
    Returns a NumPy array of shape (num_sections, 6):
       [b0, b1, b2, a0, a1, a2]
    """
    sos = []
    for line in lines:
        line = line.strip()
        if not line:
            continue  # skip empty lines
        parts = [p.strip() for p in line.split()]
        # Expect 6 parts: [b0, b1, b2, a0, a1, a2]
        if len(parts) != 6:
            print(f"Skipping malformed line: {line}")
            continue
        nums = list(map(float, parts))
        sos.append(nums)
    return np.array(sos, dtype=np.float64)

def plot_frequency_response(sos, gain=1.0, fs=2.0, title="Filter Frequency Response"):
    """
    Plots the magnitude response of the given SOS filter at sampling rate fs.
    'gain' can be applied as a scalar multiplier if the filter has a separate gain factor.
    """
    w, h = sosfreqz(sos, worN=2048, fs=fs)
    h_dB = 20 * np.log10(np.abs(h) * gain + 1e-20)

    plt.figure()
    plt.plot(w, h_dB)
    plt.xlim(0.0, fs/2)
    plt.title(title, wrap=True)
    if fs > 2.0:
        plt.xlabel("Frequency [Hz]")
    else:
        plt.xlabel("Normalized Frequency")
    plt.ylabel("Magnitude [dB]")
    plt.grid(True)
    plt.show()

def run_default_demo():
    """
    Runs the "default" lines if no arguments are provided.
    This is just an example with some arbitrary filter lines & gain.
    """

    print("No arguments provided. Using a default filter demo...")

    # Example: a 4-section Butterworth LPF
    gain = 0.001233279705048
    lines = [
        "0.999999999534339 0.000000000000000 -1.000000000000000 1 -0.000000000000000 0.781860132236034",
        "0.500000000000000 0.999999999534339 0.500000000000000 1 -0.270521447993815 0.805343665648252",
                          "0.500000000000000 -1.000000000000000 0.500000000000000 1 0.270521447993815 0.805343665648252",
        "0.500000000000000 0.999999999534339 0.500000000000000 1 -0.490443787537515 0.867302551865578",
                          "0.500000000000000 -1.000000000000000 0.500000000000000 1 0.490443787537515 0.867302551865578",
        "0.500000000000000 0.999999999534339 0.500000000000000 1 -0.626166986301541 0.952287625055760",
                          "0.500000000000000 -1.000000000000000 0.500000000000000 1 0.626166986301541 0.952287625055760"
    ]

    sos = parse_sos_coeffs(lines)
    plot_frequency_response(sos, gain=gain, title="Default Demo Filter")

def main():
    # If no arguments, run the default example
    if len(sys.argv) == 1:
        run_default_demo()
        return

    # Otherwise, user must specify at least one argument = JSON file path
    json_file = sys.argv[1]
    if not os.path.isfile(json_file):
        print(f"Error: File '{json_file}' not found.")
        sys.exit(1)

    # Attempt to load JSON
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON file '{json_file}': {e}")
        sys.exit(1)

    # If ONLY the JSON file is given, list top-level keys
    if len(sys.argv) == 2:
        print(f"Loaded '{json_file}'. Top-level filter presets:")
        for key in data.keys():
            print(" ", key)
        return

    # Else, we have more arguments => each is a preset key to plot
    for preset_key in sys.argv[2:]:
        if preset_key not in data:
            print(f"Warning: Key '{preset_key}' not found in {json_file}.")
            print("Available keys:")
            for k in data.keys():
                print(" ", k)
            continue

        preset = data[preset_key]

        # Extract gain (if present) and lines (the SOS strings)
        gain = preset.get("gain", 1.0)
        lines = preset.get("sos", [])

        # Build a string of all other fields for title
        # (skip the actual 'lines'/'sos' array itself)
        meta_fields = {
            k: v
            for k, v in preset.items()
            if k not in ("sos", "gain")
        }
        # Build a short string like "order=4, freq_cutoff=0.75, ripple_dB=0.1"
        meta_str = ", ".join(f"{k}={v}" for k, v in meta_fields.items())

        # We'll use that info in the console printout AND in the plot title
        print(f"\n=== Preset: {preset_key} ===")
        if meta_str:
            print(f"  Fields: {meta_str}")
        print(f"  gain = {gain}")
        print(f"  lines (count={len(lines)}): {lines}")

        if not lines:
            print(f"Preset '{preset_key}' has no 'lines' or 'sos' field!")
            continue

        # Parse the lines into a numeric SOS array
        sos = parse_sos_coeffs(lines)

        # Compose a title for the plot
        # e.g. "Preset: butter_lpf (order=4, freq_cutoff=0.75)"
        title = f"Preset: {preset_key}"
        if meta_str:
            title += f" ({meta_str})"

        # Plot
        plot_frequency_response(sos, gain=gain, title=title)

if __name__ == "__main__":
    main()