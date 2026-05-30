"""plot_distribution.py — prevalence, association forest, and co-occurrence heatmap.

Outputs (written to --outdir):
  prevalence_by_taxon.{png,svg}
  association_forest.{png,svg}
  cooccurrence_heatmap.{png,svg}
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
import matplotlib.colors as mcolors

apply_style()


def _load(path, name):
    if not path or not os.path.exists(path):
        return None, f"{name} missing"
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None, f"{name} empty"
    return df, None


# ---------------------------------------------------------------------------
# Prevalence horizontal bar
# ---------------------------------------------------------------------------

def plot_prevalence(prevalence, outdir):
    prefix = os.path.join(outdir, "prevalence_by_taxon")

    if prevalence is None:
        placeholder(prefix, "no data: prevalence_by_taxon.tsv missing or empty")
        return

    req = {"rank", "taxon", "prevalence", "n_genomes"}
    if not req.issubset(prevalence.columns):
        placeholder(prefix, "no data: required columns missing from prevalence file")
        return

    df = prevalence[prevalence["rank"] == "gtdb_family"].copy()
    if df.empty:
        placeholder(prefix, "no data: no rows with rank==gtdb_family")
        return

    df = df.sort_values("n_genomes", ascending=False).head(20)
    df = df.sort_values("prevalence", ascending=True)

    # Wilson CIs
    ci_low  = df.get("ci_low",  pd.Series(np.zeros(len(df)), index=df.index))
    ci_high = df.get("ci_high", pd.Series(np.zeros(len(df)), index=df.index))
    xerr_low  = (df["prevalence"] - ci_low).clip(lower=0)
    xerr_high = (ci_high - df["prevalence"]).clip(lower=0)

    fig, ax = plt.subplots(figsize=(3.5, max(3.0, len(df) * 0.22)))
    y = np.arange(len(df))
    ax.barh(y, df["prevalence"].values, height=0.6,
            color=OKABE_ITO["blue"], alpha=0.8, edgecolor="none")
    ax.errorbar(df["prevalence"].values, y,
                xerr=[xerr_low.values, xerr_high.values],
                fmt="none", color="black", linewidth=0.8, capsize=2)
    ax.set_yticks(y)
    ax.set_yticklabels(df["taxon"].values, fontsize=5)
    ax.set_xlabel("Prevalence")
    ax.set_title("Prometheus prevalence by family (top 20)")
    ax.set_xlim(0, 1)

    for yi, (_, row) in zip(y, df.iterrows()):
        ax.text(1.01, yi, f"n={int(row['n_genomes'])}", va="center", fontsize=4,
                transform=ax.get_yaxis_transform())

    save_fig(fig, prefix)


# ---------------------------------------------------------------------------
# Association forest plot
# ---------------------------------------------------------------------------

def plot_association_forest(tests, outdir):
    prefix = os.path.join(outdir, "association_forest")

    if tests is None:
        placeholder(prefix, "no data: association_tests.tsv missing or empty")
        return

    req = {"estimate", "term"}
    if not req.issubset(tests.columns):
        placeholder(prefix, "no data: estimate/term columns missing")
        return

    df = tests.dropna(subset=["estimate"]).copy()
    if df.empty:
        placeholder(prefix, "no data: no rows with estimate values")
        return

    # Build row labels
    label_parts = []
    for _, row in df.iterrows():
        parts = []
        if "model" in df.columns and pd.notna(row.get("model")):
            parts.append(str(row["model"]))
        if "term" in df.columns and pd.notna(row.get("term")):
            parts.append(str(row["term"]))
        label_parts.append(" | ".join(parts) if parts else str(_))
    df["_label"] = label_parts

    ci_low  = df.get("ci_low",  df["estimate"] - df.get("std_error", pd.Series(0, index=df.index)))
    ci_high = df.get("ci_high", df["estimate"] + df.get("std_error", pd.Series(0, index=df.index)))

    fig, ax = plt.subplots(figsize=(4.5, max(2.5, len(df) * 0.28)))
    y = np.arange(len(df))

    p_vals = df.get("p_value", pd.Series(np.ones(len(df)), index=df.index)).fillna(1.0)
    sig = p_vals < 0.05
    colors = [OKABE_ITO["blue"] if s else OKABE_ITO["grey"] for s in sig]

    ax.scatter(df["estimate"].values, y, color=colors, s=14, zorder=4)
    ax.hlines(y,
              ci_low.values, ci_high.values,
              color=colors, linewidth=0.8)
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--", zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(df["_label"].values, fontsize=5)
    ax.set_xlabel("Estimate (log-odds / statistic)")
    ax.set_title("Association tests")

    handles = [mpatches.Patch(color=OKABE_ITO["blue"], label="p<0.05"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="p≥0.05")]
    ax.legend(handles=handles, frameon=False, fontsize=6)
    save_fig(fig, prefix)


# ---------------------------------------------------------------------------
# Co-occurrence heatmap strip
# ---------------------------------------------------------------------------

def plot_cooccurrence_heatmap(cooc, outdir):
    prefix = os.path.join(outdir, "cooccurrence_heatmap")

    if cooc is None:
        placeholder(prefix, "no data: cooccurrence.tsv missing or empty")
        return

    req = {"system_b", "p_adj"}
    if not req.issubset(cooc.columns):
        placeholder(prefix, "no data: system_b/p_adj columns missing")
        return

    df = cooc.dropna(subset=["p_adj"]).copy()
    if df.empty:
        placeholder(prefix, "no data: no rows with p_adj values")
        return

    df = df.sort_values("p_adj")

    # Signed -log10(p_adj): sign from direction column
    log_p = -np.log10(df["p_adj"].clip(lower=1e-300))
    if "direction" in df.columns:
        sign = df["direction"].apply(
            lambda d: -1 if str(d).strip().lower() in ("negative", "exclusion", "-1", "-") else 1)
    else:
        sign = pd.Series(np.ones(len(df)), index=df.index)

    signed = (log_p * sign).values.reshape(1, -1)
    labels = df["system_b"].values
    sig_mask = df["p_adj"].values < 0.05

    vmax = np.abs(signed).max() if np.abs(signed).max() > 0 else 1.0
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig_w = max(3.5, len(labels) * 0.35)
    fig, ax = plt.subplots(figsize=(fig_w, 1.8))
    im = ax.imshow(signed, aspect="auto", cmap="RdBu_r", norm=norm)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=5)
    ax.set_yticks([])
    ax.set_title("Prometheus co-occurrence (signed −log10 p_adj)")
    cbar = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.3, fraction=0.04)
    cbar.set_label("signed −log10(p_adj)", fontsize=6)

    # Annotate significant
    for j, s in enumerate(sig_mask):
        if s:
            ax.text(j, 0, "*", ha="center", va="center", fontsize=7, color="white")

    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description="Distribution figures")
    parser.add_argument("--prevalence",   required=True, help="distribution/prevalence_by_taxon.tsv")
    parser.add_argument("--tests",        required=True, help="distribution/association_tests.tsv")
    parser.add_argument("--cooccurrence", required=True, help="distribution/cooccurrence.tsv")
    parser.add_argument("--outdir",       required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    prevalence, _ = _load(args.prevalence,   "prevalence")
    tests,      _ = _load(args.tests,        "association_tests")
    cooc,       _ = _load(args.cooccurrence, "cooccurrence")

    plot_prevalence(prevalence, args.outdir)
    plot_association_forest(tests, args.outdir)
    plot_cooccurrence_heatmap(cooc, args.outdir)


if __name__ == "__main__":
    main()
