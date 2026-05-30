"""plot_ecology.py — ecology enrichment and viral/R-loop association figures.

Outputs (written to --outdir):
  ecology_enrichment.{png,svg}
  viral_rloop_association.{png,svg}
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


# ---------------------------------------------------------------------------
# Ecology forest/bar
# ---------------------------------------------------------------------------

def plot_ecology_enrichment(associations, outdir):
    prefix = os.path.join(outdir, "ecology_enrichment")

    if associations is None:
        placeholder(prefix, "no data: ecology_associations.tsv missing or empty")
        return

    if "estimate" not in associations.columns:
        placeholder(prefix, "no data: estimate column missing")
        return

    df = associations.dropna(subset=["estimate"]).copy()
    if df.empty:
        placeholder(prefix, "no data: all estimate values are NaN")
        return

    # Row labels: predictor:level
    labels = []
    for _, row in df.iterrows():
        parts = []
        if "predictor" in df.columns and pd.notna(row.get("predictor")):
            parts.append(str(row["predictor"]))
        if "level" in df.columns and pd.notna(row.get("level")):
            parts.append(str(row["level"]))
        labels.append(":".join(parts) if parts else "?")
    df["_label"] = labels

    p_adj = df.get("p_adj", df.get("p_value", pd.Series(np.ones(len(df)), index=df.index))).fillna(1.0)
    sig = (p_adj < 0.05).values
    colors = [OKABE_ITO["vermillion"] if s else OKABE_ITO["grey"] for s in sig]

    ci_lo = df.get("ci_low",  df["estimate"] - df.get("std_error", pd.Series(0, index=df.index)))
    ci_hi = df.get("ci_high", df["estimate"] + df.get("std_error", pd.Series(0, index=df.index)))

    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(4.0, max(2.5, len(df) * 0.28)))
    ax.scatter(df["estimate"].values, y, color=colors, s=14, zorder=4)
    ax.hlines(y, ci_lo.values, ci_hi.values, color=colors, linewidth=0.8)
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(df["_label"].values, fontsize=5)
    ax.set_xlabel("Log-odds estimate")
    ax.set_title("Ecology associations")
    handles = [mpatches.Patch(color=OKABE_ITO["vermillion"], label="p_adj<0.05"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="p_adj≥0.05")]
    ax.legend(handles=handles, frameon=False, fontsize=6)
    save_fig(fig, prefix)


# ---------------------------------------------------------------------------
# Viral / R-loop association
# ---------------------------------------------------------------------------

def plot_viral_rloop(viral, master, outdir):
    prefix = os.path.join(outdir, "viral_rloop_association")

    if viral is None or master is None:
        placeholder(prefix, "no data: viral_summary.tsv or master_table.tsv missing")
        return

    req_v = {"genome_id", "n_proviruses"}
    req_m = {"genome_id", "prometheus_present"}
    if not req_v.issubset(viral.columns):
        placeholder(prefix, f"no data: viral_summary missing columns {req_v - set(viral.columns)}")
        return
    if not req_m.issubset(master.columns):
        placeholder(prefix, f"no data: master_table missing columns {req_m - set(master.columns)}")
        return

    merge_cols = ["genome_id", "prometheus_present"]
    merged = viral.merge(master[merge_cols], on="genome_id", how="inner")
    if merged.empty:
        placeholder(prefix, "no data: no genome_id overlap between viral_summary and master_table")
        return

    merged["prometheus_present"] = (
        merged["prometheus_present"].astype(str).str.strip().str.lower()
        .isin(["1", "true", "yes"])
    )

    panels = []
    if "n_proviruses" in merged.columns:
        panels.append(("n_proviruses", "No. proviruses"))
    if "frac_rloop_proviruses" in merged.columns:
        panels.append(("frac_rloop_proviruses", "Frac. R-loop proviruses"))

    if not panels:
        placeholder(prefix, "no data: n_proviruses/frac_rloop_proviruses columns missing")
        return

    groups = {True: "Prometheus+", False: "Prometheus−"}
    colors = {True: OKABE_ITO["blue"], False: OKABE_ITO["grey"]}

    fig, axes = plt.subplots(1, len(panels), figsize=(3.5 * len(panels), 3.2))
    if len(panels) == 1:
        axes = [axes]

    for ax, (col, ylabel) in zip(axes, panels):
        vals_t = merged.loc[merged["prometheus_present"],  col].dropna().values
        vals_f = merged.loc[~merged["prometheus_present"], col].dropna().values

        data   = [vals_f, vals_t]
        xlabs  = [groups[False], groups[True]]
        xcols  = [colors[False], colors[True]]

        for xi, (d, c) in enumerate(zip(data, xcols)):
            if len(d) == 0:
                continue
            bp = ax.boxplot(d, positions=[xi], widths=0.5, patch_artist=True,
                            boxprops=dict(facecolor=c, alpha=0.55),
                            medianprops=dict(color="black", linewidth=0.8),
                            flierprops=dict(marker="o", markersize=2, alpha=0.3,
                                            markerfacecolor=c, markeredgewidth=0),
                            whiskerprops=dict(linewidth=0.7),
                            capprops=dict(linewidth=0.7))
            jitter = np.random.default_rng(42).uniform(-0.18, 0.18, size=len(d))
            ax.scatter(xi + jitter, d, s=4, alpha=0.4, color=c, zorder=3)

        # Mann-Whitney p
        if len(vals_t) > 0 and len(vals_f) > 0:
            try:
                from scipy.stats import mannwhitneyu
                _, pval = mannwhitneyu(vals_t, vals_f, alternative="two-sided")
                pstr = f"MW p={pval:.3g}"
                ax.set_title(pstr, fontsize=6)
            except Exception:
                pass

        ax.set_xticks([0, 1])
        ax.set_xticklabels(xlabs, fontsize=6)
        ax.set_ylabel(ylabel)

    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description="Ecology figures")
    parser.add_argument("--associations", required=True, help="ecology/ecology_associations.tsv")
    parser.add_argument("--viral",        required=True, help="ecology/viral/viral_summary.tsv")
    parser.add_argument("--master",       required=True, help="scaffolds/master_table.tsv")
    parser.add_argument("--outdir",       required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    associations, _ = _load(args.associations, "ecology_associations")
    viral,        _ = _load(args.viral,        "viral_summary")
    master,       _ = _load(args.master,       "master_table")

    plot_ecology_enrichment(associations, args.outdir)
    plot_viral_rloop(viral, master, args.outdir)


if __name__ == "__main__":
    main()
