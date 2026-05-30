"""plot_selection.py — molecular selection figures.

Outputs (written to --outdir):
  selection_by_domain.{png,svg}
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

apply_style()


def _load(path, name):
    if not path or not os.path.exists(path):
        return None, f"{name} missing"
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return None, f"{name} empty"
    return df, None


def plot_selection_by_domain(summary, sites, outdir):
    prefix = os.path.join(outdir, "selection_by_domain")

    if sites is None:
        placeholder(prefix, "no data: site_table.tsv missing or empty")
        return

    if "domain" not in sites.columns:
        placeholder(prefix, "no data: domain column missing from site_table.tsv")
        return

    # Determine method from summary
    method_label = "FEL/MEME"
    if summary is not None and "method" in summary.columns:
        methods = summary["method"].dropna().unique().tolist()
        if methods:
            method_label = "/".join(methods)

    # Positive sites per domain
    if "classification" in sites.columns:
        pos = sites[sites["classification"] == "positive"]
    elif "p_value" in sites.columns:
        pos = sites[sites["p_value"] < 0.05]
    else:
        pos = sites

    pos_counts = pos.groupby("domain").size().reset_index(name="n_positive")
    if pos_counts.empty:
        placeholder(prefix, "no data: no positively-selected sites detected")
        return

    # Check for continuous dN/dS column
    dnds_col = None
    for candidate in ["beta_alpha_or_dnds", "dNdS", "dnds", "omega"]:
        if candidate in sites.columns:
            dnds_col = candidate
            break

    has_dnds = dnds_col is not None and not sites[dnds_col].dropna().empty
    n_panels = 2 if has_dnds else 1

    fig, axes = plt.subplots(1, n_panels, figsize=(3.5 * n_panels, 3.0))
    if n_panels == 1:
        axes = [axes]

    # Panel 1: bar of positive sites per domain
    ax = axes[0]
    domains = pos_counts["domain"].values
    counts  = pos_counts["n_positive"].values
    palette = list(OKABE_ITO.values())
    colors  = [palette[i % len(palette)] for i in range(len(domains))]
    bars = ax.bar(range(len(domains)), counts, color=colors, edgecolor="none", width=0.65)
    ax.set_xticks(range(len(domains)))
    ax.set_xticklabels(domains, rotation=30, ha="right", fontsize=6)
    ax.set_ylabel("Positively-selected sites")
    ax.set_title(f"Selection per domain ({method_label})")

    # Panel 2: violin/strip of dN/dS by domain (optional)
    if has_dnds:
        ax2 = axes[1]
        domain_order = pos_counts["domain"].tolist()
        dnds_vals = sites[["domain", dnds_col]].dropna()
        dnds_vals = dnds_vals[dnds_vals["domain"].isin(domain_order)]

        data_by_domain = [dnds_vals.loc[dnds_vals["domain"] == d, dnds_col].values
                          for d in domain_order]
        data_by_domain = [v for v in data_by_domain if len(v) > 0]
        valid_domains  = [d for d, v in zip(domain_order, [
            dnds_vals.loc[dnds_vals["domain"] == d, dnds_col].values for d in domain_order
        ]) if len(v) > 0]

        if data_by_domain:
            parts = ax2.violinplot(data_by_domain, positions=range(len(valid_domains)),
                                   showmedians=True, widths=0.7)
            for i, (body, d) in enumerate(zip(parts["bodies"], valid_domains)):
                body.set_facecolor(palette[i % len(palette)])
                body.set_alpha(0.7)
            # Strip
            for i, (vals, d) in enumerate(zip(data_by_domain, valid_domains)):
                jitter = np.random.default_rng(42).uniform(-0.1, 0.1, size=len(vals))
                ax2.scatter(i + jitter, vals, s=4, alpha=0.4,
                            color=palette[i % len(palette)], zorder=3)
        ax2.set_xticks(range(len(valid_domains)))
        ax2.set_xticklabels(valid_domains, rotation=30, ha="right", fontsize=6)
        ax2.set_ylabel(dnds_col)
        ax2.set_title("dN/dS by domain")
        ax2.axhline(1.0, color="black", linewidth=0.6, linestyle="--")

    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description="Selection figures")
    parser.add_argument("--summary", required=True, help="evolution/selection/selection_summary.tsv")
    parser.add_argument("--sites",   required=True, help="evolution/selection/site_table.tsv")
    parser.add_argument("--outdir",  required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    summary, _ = _load(args.summary, "selection_summary")
    sites,   _ = _load(args.sites,   "site_table")

    plot_selection_by_domain(summary, sites, args.outdir)


if __name__ == "__main__":
    main()
