#!/usr/bin/env python3
"""
Plot single-cell AP trace and print biomarkers.

Usage:
  python plot_singlecell.py --dir ../Porcisa/SingleCell
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from textwrap import indent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="Directory with singlecell_*.csv")
    ap.add_argument("--show", action="store_true", help="Show interactive plot")
    args = ap.parse_args()

    trace_csv = os.path.join(args.dir, "singlecell_trace.csv")
    bio_csv   = os.path.join(args.dir, "singlecell_biomarkers.csv")

    if not os.path.exists(trace_csv) or not os.path.exists(bio_csv):
        raise FileNotFoundError("Expected singlecell_trace.csv and singlecell_biomarkers.csv in --dir")

    # Load
    df_trace = pd.read_csv(trace_csv)       # columns: time_ms, Vm_mV
    df_bio   = pd.read_csv(bio_csv)         # single row with biomarkers

    # Pretty-print biomarkers
    bio = df_bio.iloc[0].to_dict()
    ordered_cols = [
        "period_ms", "stim_start_ms", "stim_dur_ms", "stim_amp_uAcm2",
        "RMP_mV", "Vpeak_mV", "APA_mV", "dVdt_max_mV_per_ms",
        "APD50_ms", "APD90_ms", "Triang_ms", "t_up_ms"
    ]
    rows = []
    for k in ordered_cols:
        if k in bio:
            rows.append(f"{k:>22s} : {bio[k]:.6g}")
    print("\nSingle-cell biomarkers:\n" + indent("\n".join(rows), "  ") + "\n")

    # Plot
    fig, ax = plt.subplots()
    ax.plot(df_trace["time_ms"], df_trace["Vm_mV"], lw=1.5)
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Vm (mV)")
    title = f"Single-cell AP  |  APA={bio.get('APA_mV', float('nan')):.1f} mV,  APD90={bio.get('APD90_ms', float('nan')):.1f} ms"
    ax.set_title(title)
    ax.grid(True, ls=":")

    # Annotations: t_up, APD50, APD90
    t_up   = bio.get("t_up_ms", None)
    apd50  = bio.get("APD50_ms", None)
    apd90  = bio.get("APD90_ms", None)

    def vline(x, label):
        if pd.notna(x):
            ax.axvline(x, ls="--", alpha=0.8)
            ax.text(x, ax.get_ylim()[1], f" {label}", va="bottom", ha="left", rotation=90)

    if pd.notna(t_up):
        vline(t_up, "t_up")
    if pd.notna(t_up) and pd.notna(apd50):
        vline(t_up + apd50, "APD50")
    if pd.notna(t_up) and pd.notna(apd90):
        vline(t_up + apd90, "APD90")

    fig.tight_layout()
    out_png = os.path.join(args.dir, "singlecell_trace.png")
    fig.savefig(out_png, dpi=300)
    print(f"Saved plot: {out_png}")

    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
