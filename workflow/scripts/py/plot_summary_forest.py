"""plot_summary_forest.py — synthesis summary forest plot.

Outputs (written to --outdir):
  summary_forest.{png,svg}
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


def _load(path, name):
    if not path or not os.path.exists(path):
        return None, f"{name} missing"
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None, f"{name} empty"
    return df, None


def plot_summary_forest(registry, outdir):
    prefix = os.path.join(outdir, "summary_forest")

    if registry is None:
        placeholder(prefix, "no data: results_registry.tsv missing or empty")
        return

    req = {"estimate", "hypothesis"}
    if not req.issubset(registry.columns):
        placeholder(prefix, "no data: estimate/hypothesis columns missing")
        return

    df = registry[registry["hypothesis"] == "confirmatory"].copy()
    if df.empty:
        placeholder(prefix, "no data: no rows with hypothesis==confirmatory")
        return

    df = df.dropna(subset=["estimate"])
    if df.empty:
        placeholder(prefix, "no data: all confirmatory estimates are NaN")
        return

    # Two-column labels: family/term  +  p_adj
    row_labels = []
    p_adj_vals = []
    for _, row in df.iterrows():
        fam  = str(row.get("family", "")) if "family" in df.columns and pd.notna(row.get("family")) else ""
        term = str(row.get("term",   "")) if "term"   in df.columns and pd.notna(row.get("term"))   else ""
        label = " / ".join(filter(None, [fam, term])) or "?"
        row_labels.append(label)

        padj = row.get("p_adj", row.get("p_value", np.nan))
        p_adj_vals.append(padj)

    df["_label"] = row_labels
    p_adj_arr = pd.array(p_adj_vals, dtype=float)
    sig = np.array([p < 0.05 if not np.isnan(p) else False for p in p_adj_arr])
    colors = [OKABE_ITO["blue"] if s else OKABE_ITO["grey"] for s in sig]

    ci_lo = df.get("ci_low",  df["estimate"] - df.get("std_error", pd.Series(0, index=df.index)))
    ci_hi = df.get("ci_high", df["estimate"] + df.get("std_error", pd.Series(0, index=df.index)))

    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(5.0, max(2.5, len(df) * 0.30)))

    ax.scatter(df["estimate"].values, y, color=colors, s=18, zorder=4)
    ax.hlines(y, ci_lo.values, ci_hi.values, color=colors, linewidth=0.9)
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--", zorder=2)

    ax.set_yticks(y)
    ax.set_yticklabels(df["_label"].values, fontsize=5)
    ax.set_xlabel("Estimate (log-odds / statistic)")
    ax.set_title("Confirmatory results (synthesis forest)")

    # p_adj annotations on the right
    x_right = ax.get_xlim()[1]
    for yi, padj in zip(y, p_adj_arr):
        if not np.isnan(padj):
            pstr = f"{padj:.2g}"
        else:
            pstr = "—"
        ax.text(1.02, yi, pstr, transform=ax.get_yaxis_transform(),
                va="center", fontsize=5, ha="left")
    ax.text(1.02, len(df), "p_adj", transform=ax.get_yaxis_transform(),
            va="center", fontsize=5, ha="left", fontweight="bold")

    handles = [mpatches.Patch(color=OKABE_ITO["blue"], label="p_adj<0.05"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="p_adj≥0.05")]
    ax.legend(handles=handles, frameon=False, fontsize=6)

    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description="Summary forest figure")
    parser.add_argument("--registry", required=True, help="synthesis/results_registry.tsv")
    parser.add_argument("--outdir",   required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    registry, _ = _load(args.registry, "results_registry")
    plot_summary_forest(registry, args.outdir)


if __name__ == "__main__":
    main()
