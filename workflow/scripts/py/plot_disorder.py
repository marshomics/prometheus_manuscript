"""Disorder figures for the ProA/Prometheus disorder module.

Outputs (written to --outdir):
  disorder_along_protein.{png,svg}   mean consensus disorder per alignment column with a
                                     domain track (disorder should pile up in linkers/
                                     tails and drop inside folded domains -- a built-in
                                     sanity check on the prediction).
  disorder_by_group.{png,svg}        per-genome disorder split by methanogen status (and,
                                     if a bacterial table is given, by domain of life).
  disorder_linker_vs_domain.{png,svg} disorder in domain vs linker columns.

All distributions show the spread; the inferential claims are the phylogeny-controlled
phylolm effects in disorder_associations.tsv, not these splits.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

apply_style()
DOMAIN_COLORS = {"DUF4011": "#E69F00", "DUF3320": "#56B4E9", "DNA2": "#009E73",
                 "nuclease": "#CC79A7", "REase": "#882255", "NA": "#DDDDDD"}


def _load(path):
    if not path or not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    return None if df.empty else df


def _bool(s):
    return s.astype(str).str.lower().isin(("true", "1", "yes"))


def _resp(df):
    return "frac_disordered_nonlcr" if "frac_disordered_nonlcr" in df.columns else "frac_disordered"


# --------------------------------------------------------------------------- #
def plot_along_protein(bycol, outdir):
    prefix = os.path.join(outdir, "disorder_along_protein")
    if bycol is None or not {"alignment_col", "mean_disorder"}.issubset(bycol.columns):
        placeholder(prefix, "no data: disorder-by-column table missing/empty")
        return
    df = bycol.copy()
    df["x"] = pd.to_numeric(df["alignment_col"], errors="coerce")
    df["y"] = pd.to_numeric(df["mean_disorder"], errors="coerce")
    df = df.dropna(subset=["x", "y"]).sort_values("x")
    if df.empty:
        placeholder(prefix, "no data: no usable columns")
        return
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(7.0, 3.0), sharex=True,
                                  gridspec_kw={"height_ratios": [6, 1], "hspace": 0.06})
    ax.fill_between(df["x"].values, 0, df["y"].values, color=OKABE_ITO["blue"], alpha=0.35)
    ax.plot(df["x"].values, df["y"].values, color=OKABE_ITO["blue"], linewidth=0.8)
    ax.set_ylim(0, 1); ax.set_ylabel("mean consensus disorder")
    ax.set_title("Disorder along the Prometheus alignment")
    if "domain" in df.columns:
        for _, r in df.iterrows():
            axd.add_patch(plt.Rectangle((r["x"] - 0.5, 0), 1, 1,
                          color=DOMAIN_COLORS.get(str(r["domain"]), DOMAIN_COLORS["NA"]), linewidth=0))
        axd.set_xlim(df["x"].min() - 1, df["x"].max() + 1)
    axd.set_yticks([]); axd.set_ylim(0, 1); axd.set_xlabel("alignment column")
    from matplotlib.patches import Patch
    handles = [Patch(color=v, label=k) for k, v in DOMAIN_COLORS.items() if k != "NA"]
    handles.append(Patch(color=DOMAIN_COLORS["NA"], label="linker/NA"))
    ax.legend(handles=handles, frameon=False, fontsize=5, ncol=3, loc="upper right")
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_by_group(per_locus, master, bact_locus, outdir):
    prefix = os.path.join(outdir, "disorder_by_group")
    if per_locus is None or "genome_id" not in per_locus.columns:
        placeholder(prefix, "no data: disorder_per_locus.tsv missing/empty")
        return
    resp = _resp(per_locus)
    pl = per_locus.copy()
    pl[resp] = pd.to_numeric(pl[resp], errors="coerce")
    pg = pl.groupby("genome_id")[resp].mean().reset_index()
    panels = []
    # methanogen split (archaea)
    if master is not None and {"genome_id", "methanogen"}.issubset(master.columns):
        m = master[["genome_id", "methanogen"]].copy()
        m["genome_id"] = m["genome_id"].astype(str)
        pg["genome_id"] = pg["genome_id"].astype(str)
        j = pg.merge(m, on="genome_id", how="left")
        j["meth"] = _bool(j["methanogen"]).astype(int)
        if j["meth"].nunique() == 2:
            panels.append(("methanogen", [j.loc[j["meth"] == g, resp].dropna().values for g in (0, 1)],
                           ["non-methanogen", "methanogen"]))
    # domain-of-life split (cross-domain) if bacterial loci provided
    if bact_locus is not None and "genome_id" in bact_locus.columns:
        bresp = _resp(bact_locus)
        bl = bact_locus.copy(); bl[bresp] = pd.to_numeric(bl[bresp], errors="coerce")
        bpg = bl.groupby("genome_id")[bresp].mean()
        panels.append(("domain of life", [pg[resp].dropna().values, bpg.dropna().values],
                       ["Archaea", "Bacteria"]))
    if not panels:
        placeholder(prefix, "no data: no grouping (methanogen / bacterial) available")
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(2.6 * len(panels) + 0.5, 2.7))
    if len(panels) == 1:
        axes = [axes]
    for ax, (title, groups, labels) in zip(axes, panels):
        groups = [g for g in groups]
        if all(len(g) for g in groups):
            parts = ax.violinplot(groups, showmedians=True, widths=0.8)
            for b, c in zip(parts["bodies"], (OKABE_ITO["grey"], OKABE_ITO["blue"])):
                b.set_facecolor(c); b.set_alpha(0.6)
        ax.set_xticks(range(1, len(labels) + 1)); ax.set_xticklabels(labels, fontsize=6, rotation=15)
        ax.set_ylabel(f"per-genome disorder ({resp})", fontsize=6); ax.set_title(title, fontsize=7)
    fig.suptitle("ProA disorder by group (distributions; tests are phylo-controlled)", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_linker_vs_domain(bycol, outdir):
    prefix = os.path.join(outdir, "disorder_linker_vs_domain")
    if bycol is None or not {"mean_disorder", "domain"}.issubset(bycol.columns):
        placeholder(prefix, "no data: disorder-by-column table missing domain/disorder")
        return
    df = bycol.copy()
    df["y"] = pd.to_numeric(df["mean_disorder"], errors="coerce")
    df = df.dropna(subset=["y"])
    df["grp"] = np.where(df["domain"].astype(str).isin(["NA", "nan", ""]), "linker/tail", "domain")
    groups = [df.loc[df["grp"] == g, "y"].values for g in ("domain", "linker/tail")]
    if not all(len(g) for g in groups):
        placeholder(prefix, "no data: need both domain and linker columns")
        return
    fig, ax = plt.subplots(figsize=(3.2, 2.8))
    parts = ax.violinplot(groups, showmedians=True, widths=0.8)
    for b, c in zip(parts["bodies"], (OKABE_ITO["green"], OKABE_ITO["orange"])):
        b.set_facecolor(c); b.set_alpha(0.6)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["folded domain", "linker/tail"], fontsize=6)
    ax.set_ylabel("mean consensus disorder per column"); ax.set_ylim(0, 1)
    ax.set_title("Disorder: domain vs linker columns", fontsize=7)
    save_fig(fig, prefix)


def main():
    ap = argparse.ArgumentParser(description="Disorder figures")
    ap.add_argument("--by-column", required=True)
    ap.add_argument("--per-locus", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--bact-per-locus", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    bycol = _load(args.by_column); per_locus = _load(args.per_locus)
    master = _load(args.master); bact = _load(args.bact_per_locus)
    plot_along_protein(bycol, args.outdir)
    plot_by_group(per_locus, master, bact, args.outdir)
    plot_linker_vs_domain(bycol, args.outdir)


if __name__ == "__main__":
    main()
