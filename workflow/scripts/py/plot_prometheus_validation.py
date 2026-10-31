"""plot_prometheus_validation.py — Prometheus detection→validation summary + architecture.

Detection is DefenseFinder-only, so there is no multi-method overlap to draw.
Instead the first panel shows the validation funnel: how DefenseFinder Prometheus
candidates distribute across confidence tiers (A/B/C/Rejected), with the number of
genomes carrying a retained (non-Rejected) locus annotated.

Outputs (written to --outdir):
  prometheus_detection_summary.{png,svg}
  prometheus_architecture.{png,svg}
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO, TIER_COLORS
import utils  # noqa: F401  (kept for genome_of if needed by future panels)

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

apply_style()

TIER_ORDER = ["A", "B", "C", "Rejected"]


def _load(path, name):
    """DataFrame, or None for a missing / header-only table. A malformed table raises
    (a crashed upstream step must not turn into a 'no data' placeholder)."""
    if not path or not os.path.exists(path):
        return None, f"{name} file missing"
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df.empty:
        return None, f"{name} empty"
    return df, None


def _truthy(series):
    return series.astype(str).str.strip().str.lower().isin(["1", "true", "yes"])


def plot_detection_summary(catalog, sensitivity, outdir):
    """Validation funnel: DefenseFinder candidates by confidence tier."""
    prefix = os.path.join(outdir, "prometheus_detection_summary")

    if catalog is None or "confidence_tier" not in catalog.columns:
        placeholder(prefix, "no data: prometheus_catalog.tsv missing or has no confidence_tier")
        return

    tiers = catalog["confidence_tier"].astype(str)
    unknown = sorted(set(tiers) - set(TIER_ORDER))
    if unknown:
        raise SystemExit(f"[plot_prometheus_validation] FATAL: unexpected confidence_tier "
                         f"value(s) {unknown}")
    counts = [int((tiers == t).sum()) for t in TIER_ORDER]

    # genomes carrying a retained (present) Prometheus locus
    n_present = ""
    if sensitivity is not None and "prometheus_present" in sensitivity.columns:
        n_present = int(_truthy(sensitivity["prometheus_present"]).sum())

    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    x = np.arange(len(TIER_ORDER))
    bars = ax.bar(x, counts, width=0.6, edgecolor="none",
                  color=[TIER_COLORS.get(t, OKABE_ITO["grey"]) for t in TIER_ORDER])
    ax.set_xticks(x)
    ax.set_xticklabels(TIER_ORDER)
    ax.set_ylabel("DefenseFinder candidate loci")
    ax.set_xlabel("Confidence tier")
    title = "Prometheus detection → validation"
    if n_present != "":
        title += f"  (genomes present: {n_present})"
    ax.set_title(title, fontsize=7)
    for bar, c in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(),
                str(c), ha="center", va="bottom", fontsize=6)
    save_fig(fig, prefix)


def plot_architecture(catalog, outdir):
    prefix = os.path.join(outdir, "prometheus_architecture")

    if catalog is None:
        placeholder(prefix, "no data: prometheus_catalog.tsv missing or empty")
        return
    if "architecture" not in catalog.columns or "confidence_tier" not in catalog.columns:
        placeholder(prefix, "no data: required columns (architecture, confidence_tier) absent")
        return

    df = catalog.copy()
    df["architecture"] = df["architecture"].astype(str).replace("", "none")
    df["confidence_tier"] = df["confidence_tier"].astype(str)
    top_arch = df["architecture"].value_counts().head(12).index.tolist()
    df_top = df[df["architecture"].isin(top_arch)].copy()

    tiers_present = [t for t in TIER_ORDER if t in df_top["confidence_tier"].unique()]
    counts = (df_top.groupby(["architecture", "confidence_tier"]).size()
              .unstack(fill_value=0).reindex(columns=tiers_present, fill_value=0)
              .loc[top_arch])

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    bottom = np.zeros(len(counts))
    for tier in tiers_present:
        vals = counts[tier].values
        ax.bar(range(len(counts)), vals, bottom=bottom,
               color=TIER_COLORS.get(tier, OKABE_ITO["grey"]),
               label=tier, edgecolor="none", width=0.7)
        bottom += vals

    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=40, ha="right", fontsize=5)
    ax.set_ylabel("No. loci")
    ax.set_xlabel("Domain architecture (N->C; DUF3320 optional)")
    ax.set_title("Prometheus architecture × confidence tier")
    handles = [mpatches.Patch(color=TIER_COLORS.get(t, OKABE_ITO["grey"]), label=t)
               for t in tiers_present]
    ax.legend(handles=handles, title="Tier", frameon=False, fontsize=6)
    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description="Prometheus validation figures")
    parser.add_argument("--catalog", required=True, help="prometheus_catalog.tsv")
    parser.add_argument("--sensitivity", required=True, help="detection_sensitivity.tsv")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    catalog, _ = _load(args.catalog, "catalog")
    sensitivity, _ = _load(args.sensitivity, "sensitivity")

    plot_detection_summary(catalog, sensitivity, args.outdir)
    plot_architecture(catalog, args.outdir)


if __name__ == "__main__":
    main()
