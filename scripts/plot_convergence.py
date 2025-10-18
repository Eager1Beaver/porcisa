#!/usr/bin/env python3
"""
plot_convergence.py
===================

Visualisation of mesh-resolution convergence results from TestConvergence.hpp.

Usage:
    python plot_convergence.py --dir <output_folder>

Looks for:
    AT_error_summary.csv
    AT_per_resolution.csv

and produces:
    convergence_error.png
    activation_times.png
inside the same folder.
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def plot_error(df, out_path):
    plt.figure()
    plt.loglog(df["res_mm"], df["L2"], "o-", label="L2 error")
    plt.loglog(df["res_mm"], df["Linf"], "s--", label="L∞ error")
    plt.xlabel("Spatial step h (mm)")
    plt.ylabel("Error (ms)")
    plt.title("Mesh Convergence: Activation Time Error")
    plt.grid(True, which="both", ls=":")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def plot_activation(df, out_path):
    plt.figure()
    probes = [c for c in df.columns if c.startswith("AT_probe")]
    for p in probes:
        plt.plot(df["res_mm"], df[p], "o-", label=p)
    plt.xlabel("Spatial step h (mm)")
    plt.ylabel("Activation Time (ms)")
    plt.title("Activation Times at Probe Points")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="Output directory containing CSVs")
    args = ap.parse_args()

    err_csv = os.path.join(args.dir, "AT_error_summary.csv")
    at_csv = os.path.join(args.dir, "AT_per_resolution.csv")

    if not os.path.exists(err_csv) or not os.path.exists(at_csv):
        raise FileNotFoundError("Expected both AT_error_summary.csv and AT_per_resolution.csv in the directory.")

    df_err = pd.read_csv(err_csv)
    df_at = pd.read_csv(at_csv)

    print("Loaded:")
    print(f" - {err_csv}: {df_err.shape} entries")
    print(f" - {at_csv}: {df_at.shape} entries")

    plot_error(df_err, os.path.join(args.dir, "convergence_error.png"))
    plot_activation(df_at, os.path.join(args.dir, "activation_times.png"))

    print("Plots saved to:")
    print(f" - convergence_error.png")
    print(f" - activation_times.png")


if __name__ == "__main__":
    main()
