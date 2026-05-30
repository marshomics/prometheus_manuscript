"""Diagnostic figure for the preprocessing stage: completeness vs contamination of
all input genomes (pass/fail shaded), and the input -> passing -> species funnel.

Output (to --outdir): preprocess_quality.{png,svg}
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

apply_style()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, "preprocess_quality")

    try:
        df = pd.read_csv(args.annotated, sep="\t")
    except Exception:
        df = pd.DataFrame()
    if df.empty or "completeness" not in df.columns:
        placeholder(prefix, "no data: annotated_metadata.tsv missing/empty")
        return

    comp = pd.to_numeric(df["completeness"], errors="coerce")
    cont = pd.to_numeric(df.get("contamination"), errors="coerce")
    passing = (comp >= args.min_completeness) & (cont <= args.max_contamination)
    try:
        n_species = max(0, sum(1 for _ in open(args.species)) - 1)
    except Exception:
        n_species = 0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0))
    m = comp.notna() & cont.notna()
    ax1.scatter(comp[m & ~passing], cont[m & ~passing], s=3, alpha=0.3,
                color=OKABE_ITO["grey"], label="fail", edgecolors="none", rasterized=True)
    ax1.scatter(comp[m & passing], cont[m & passing], s=3, alpha=0.4,
                color=OKABE_ITO["blue"], label="pass (comp/cont)", edgecolors="none", rasterized=True)
    ax1.axvline(args.min_completeness, ls="--", lw=0.6, color="k")
    ax1.axhline(args.max_contamination, ls="--", lw=0.6, color="k")
    ax1.set_xlabel("Completeness (%)"); ax1.set_ylabel("Contamination (%)")
    ax1.set_ylim(0, min(50, np.nanmax(cont) if cont.notna().any() else 10))
    ax1.legend(frameon=False, fontsize=6, markerscale=2)
    ax1.set_title("Input genome quality")

    labels = ["input", "pass (comp/cont)", "species reps"]
    vals = [len(df), int(passing.sum()), n_species]
    ax2.bar(range(3), vals, color=[OKABE_ITO["grey"], OKABE_ITO["blue"], OKABE_ITO["green"]],
            width=0.6, edgecolor="none")
    ax2.set_xticks(range(3)); ax2.set_xticklabels(labels, rotation=20, ha="right")
    ax2.set_ylabel("Genomes"); ax2.set_title("Filtering funnel")
    for i, v in enumerate(vals):
        ax2.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=6)
    save_fig(fig, prefix)


if __name__ == "__main__":
    main()
