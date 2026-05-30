"""plot_qc.py — QC figures for genome completeness/contamination and assembly statistics.

Outputs (written to --outdir):
  qc_completeness_contamination.{png,svg}
  qc_assembly_stats.{png,svg}
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO, TIER_COLORS
import utils

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

apply_style()


def _load_meta(path):
    if not path or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None
    return df


def plot_completeness_contamination(df, outdir):
    prefix = os.path.join(outdir, "qc_completeness_contamination")
    if df is None:
        placeholder(prefix, "no data: genome_metadata.tsv missing or empty")
        return

    required = {"completeness", "contamination"}
    if not required.issubset(df.columns):
        placeholder(prefix, "no data: completeness/contamination columns absent")
        return

    df = df.dropna(subset=["completeness", "contamination"])
    if df.empty:
        placeholder(prefix, "no data: all rows have NaN completeness/contamination")
        return

    fig, ax = plt.subplots(figsize=(3.5, 3.0))

    # Default tier colour if column absent
    if "quality_tier" in df.columns:
        tiers = df["quality_tier"].fillna("C")
    else:
        tiers = pd.Series(["C"] * len(df), index=df.index)

    tier_order = ["A", "B", "C", "Rejected", "fail"]
    seen = tiers.unique().tolist()
    plotted_tiers = [t for t in tier_order if t in seen] + [t for t in seen if t not in tier_order]

    for tier in plotted_tiers:
        mask = tiers == tier
        color = TIER_COLORS.get(tier, OKABE_ITO["grey"])
        ax.scatter(
            df.loc[mask, "completeness"],
            df.loc[mask, "contamination"],
            c=color, s=8, alpha=0.75, linewidths=0, label=tier, zorder=3,
        )

    # MIMAG-ish threshold lines
    for comp_thresh, ls in [(50, ":"), (90, "--")]:
        ax.axvline(comp_thresh, color="black", linewidth=0.5, linestyle=ls, zorder=2)
    for cont_thresh, ls in [(5, "--"), (10, ":")]:
        ax.axhline(cont_thresh, color="black", linewidth=0.5, linestyle=ls, zorder=2)

    ax.set_xlabel("Completeness (%)")
    ax.set_ylabel("Contamination (%)")
    ax.set_title("Genome quality (MIMAG thresholds)")

    handles = [mpatches.Patch(color=TIER_COLORS.get(t, OKABE_ITO["grey"]), label=t) for t in plotted_tiers]
    ax.legend(handles=handles, title="Tier", frameon=False)

    save_fig(fig, prefix)


def plot_assembly_stats(df, outdir):
    prefix = os.path.join(outdir, "qc_assembly_stats")
    if df is None:
        placeholder(prefix, "no data: genome_metadata.tsv missing or empty")
        return

    # Check at least one stat column is present
    stat_cols = {"genome_size", "n50", "n_contigs", "gc"}
    if not stat_cols.intersection(df.columns):
        placeholder(prefix, "no data: no assembly stat columns found")
        return

    # Colour by genome_type if available
    if "genome_type" in df.columns:
        types = df["genome_type"].fillna("unknown").astype(str)
    else:
        types = pd.Series(["all"] * len(df), index=df.index)

    type_list = sorted(types.unique())
    palette = list(OKABE_ITO.values())
    type_color = {t: palette[i % len(palette)] for i, t in enumerate(type_list)}

    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.5))
    fig.subplots_adjust(hspace=0.45, wspace=0.4)

    panels = [
        ("genome_size", "Genome size (Mb)", False, lambda x: x / 1e6 if x.max() > 1e4 else x),
        ("n50",         "N50 (bp)",          True,  lambda x: x),
        ("n_contigs",   "No. contigs",       True,  lambda x: x),
        ("gc",          "GC content",        False, lambda x: x),
    ]

    for ax, (col, xlabel, log_x, transform) in zip(axes.flat, panels):
        if col not in df.columns:
            ax.text(0.5, 0.5, f"{col}\nnot available", ha="center", va="center",
                    transform=ax.transAxes, fontsize=6)
            ax.axis("off")
            continue

        for gtype in type_list:
            vals = df.loc[types == gtype, col].dropna()
            if vals.empty:
                continue
            vals = transform(vals)
            if log_x:
                vals = np.log10(vals.clip(lower=1))
            ax.hist(vals, bins=20, alpha=0.6, color=type_color[gtype], label=gtype, edgecolor="none")

        ax.set_xlabel(("log10 " if log_x else "") + xlabel)
        ax.set_ylabel("Count")

    if len(type_list) > 1:
        handles = [mpatches.Patch(color=type_color[t], label=t) for t in type_list]
        fig.legend(handles=handles, title="Genome type", frameon=False,
                   loc="upper right", bbox_to_anchor=(1.0, 1.0))

    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description="QC figures")
    parser.add_argument("--meta",   required=True, help="qc/genome_metadata.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = _load_meta(args.meta)
    plot_completeness_contamination(df, args.outdir)
    plot_assembly_stats(df, args.outdir)


if __name__ == "__main__":
    main()
