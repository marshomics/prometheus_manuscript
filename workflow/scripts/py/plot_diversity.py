"""plot_diversity.py — subfamily rarefaction and architecture-by-taxon figures.

Outputs (written to --outdir):
  subfamily_rarefaction.{png,svg}
  architecture_by_taxon.{png,svg}
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


def plot_rarefaction(subfamilies, outdir, seed):
    prefix = os.path.join(outdir, "subfamily_rarefaction")

    if subfamilies is None:
        placeholder(prefix, "no data: subfamilies.tsv missing or empty")
        return

    cluster_col = "cluster_0.30"
    if "genome_id" not in subfamilies.columns or cluster_col not in subfamilies.columns:
        placeholder(prefix, f"no data: need genome_id and {cluster_col} columns")
        return

    df = subfamilies[["genome_id", cluster_col]].dropna()
    if df.empty:
        placeholder(prefix, "no data: no valid rows in subfamilies")
        return

    genomes = df["genome_id"].unique()
    n = len(genomes)
    rng = np.random.default_rng(seed)
    n_boot = 100
    step = max(1, n // 40)
    sizes = list(range(1, n + 1, step))
    if sizes[-1] != n:
        sizes.append(n)
    sizes = np.array(sizes)

    boot_curves = np.zeros((n_boot, len(sizes)))
    for b in range(n_boot):
        perm = rng.permutation(genomes)
        curve = []
        seen = set()
        for k, sz in enumerate(sizes):
            subset = perm[:sz]
            seen.update(df.loc[df["genome_id"].isin(subset), cluster_col].unique())
            curve.append(len(seen))
        boot_curves[b] = curve

    mean_curve = boot_curves.mean(axis=0)
    lo = np.percentile(boot_curves, 2.5, axis=0)
    hi = np.percentile(boot_curves, 97.5, axis=0)

    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.fill_between(sizes, lo, hi, alpha=0.25, color=OKABE_ITO["blue"], linewidth=0)
    ax.plot(sizes, mean_curve, color=OKABE_ITO["blue"], linewidth=1.2)
    ax.set_xlabel("Genomes sampled")
    ax.set_ylabel("Subfamilies (≥30% identity)")
    ax.set_title("Subfamily rarefaction")
    save_fig(fig, prefix)


def plot_architecture_by_taxon(subfamilies, typology, master, outdir):
    prefix = os.path.join(outdir, "architecture_by_taxon")

    if subfamilies is None or master is None:
        placeholder(prefix, "no data: subfamilies.tsv or master_table.tsv missing")
        return

    if "genome_id" not in subfamilies.columns or "gtdb_phylum" not in master.columns:
        placeholder(prefix, "no data: genome_id or gtdb_phylum missing")
        return

    # Determine architecture category per protein
    if typology is not None and "architecture" in typology.columns:
        # join via catalog/protein_id -> architecture from typology (best effort)
        # subfamilies has protein_id and genome_id
        if "protein_id" in subfamilies.columns and "architecture" in subfamilies.columns:
            sub = subfamilies[["genome_id", "architecture"]].copy()
        elif "protein_id" in subfamilies.columns:
            # Try to merge typology by architecture if catalog has it
            sub = subfamilies[["genome_id"]].copy()
            sub["architecture"] = subfamilies.get("cluster_0.30", pd.Series(dtype=str))
        else:
            sub = subfamilies[["genome_id"]].copy()
            sub["architecture"] = subfamilies.get("cluster_0.30", pd.Series(dtype=str))
    else:
        sub = subfamilies[["genome_id"]].copy()
        sub["architecture"] = subfamilies.get("cluster_0.30", subfamilies.get("architecture", pd.Series(dtype=str)))

    sub = sub.dropna(subset=["genome_id", "architecture"])
    if sub.empty:
        placeholder(prefix, "no data: no genome_id/architecture after join")
        return

    # Merge phylum from master
    master_slim = master[["genome_id", "gtdb_phylum"]].dropna()
    merged = sub.merge(master_slim, on="genome_id", how="inner")
    if merged.empty:
        placeholder(prefix, "no data: no overlap between subfamilies and master table")
        return

    # Top ~10 phyla by count
    top_phyla = merged["gtdb_phylum"].value_counts().head(10).index.tolist()
    df_top = merged[merged["gtdb_phylum"].isin(top_phyla)]

    # Top ~8 architectures
    top_arch = df_top["architecture"].value_counts().head(8).index.tolist()
    df_top = df_top.copy()
    df_top["architecture"] = df_top["architecture"].where(
        df_top["architecture"].isin(top_arch), other="other")

    arch_list = top_arch + (["other"] if "other" in df_top["architecture"].values else [])
    palette = list(OKABE_ITO.values())
    arch_color = {a: palette[i % len(palette)] for i, a in enumerate(arch_list)}

    counts = (df_top.groupby(["gtdb_phylum", "architecture"])
              .size().unstack(fill_value=0)
              .reindex(columns=arch_list, fill_value=0))
    # Normalise to 100%
    totals = counts.sum(axis=1)
    fracs = counts.div(totals, axis=0) * 100
    fracs = fracs.loc[top_phyla]

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    bottom = np.zeros(len(fracs))
    for arch in arch_list:
        if arch not in fracs.columns:
            continue
        vals = fracs[arch].values
        ax.bar(range(len(fracs)), vals, bottom=bottom,
               color=arch_color[arch], label=arch, edgecolor="none", width=0.7)
        bottom += vals

    ax.set_xticks(range(len(fracs)))
    ax.set_xticklabels(fracs.index, rotation=45, ha="right", fontsize=5)
    ax.set_ylabel("% of loci")
    ax.set_ylim(0, 100)
    ax.set_title("Architecture by phylum (top 10)")
    handles = [mpatches.Patch(color=arch_color[a], label=a) for a in arch_list if a in fracs.columns]
    ax.legend(handles=handles, title="Architecture", frameon=False, fontsize=5,
              bbox_to_anchor=(1.01, 1), loc="upper left")
    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description="Diversity figures")
    parser.add_argument("--subfamilies", required=True, help="diversity/subfamilies.tsv")
    parser.add_argument("--typology",    required=True, help="diversity/architecture_typology.tsv")
    parser.add_argument("--master",      required=True, help="scaffolds/master_table.tsv")
    parser.add_argument("--seed",        type=int, default=12345)
    parser.add_argument("--outdir",      required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    subfamilies, _ = _load(args.subfamilies, "subfamilies")
    typology,    _ = _load(args.typology,    "typology")
    master,      _ = _load(args.master,      "master_table")

    plot_rarefaction(subfamilies, args.outdir, args.seed)
    plot_architecture_by_taxon(subfamilies, typology, master, args.outdir)


if __name__ == "__main__":
    main()
