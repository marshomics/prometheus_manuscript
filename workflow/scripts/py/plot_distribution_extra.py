"""Extra distribution figures (added alongside plot_distribution.py).

Outputs (written to --outdir):
  prevalence_caterpillar.{png,svg}   prevalence +/- Wilson CI per taxon (phylum+family),
                                     point size ~ n_genomes, taxa whose CI excludes the
                                     overall prevalence flagged enriched/depleted.
  cooccurrence_oddsratio.{png,svg}   Prometheus vs each system: log2 odds ratio (direction
                                     + magnitude) with BH significance, a complement to the
                                     signed-p heatmap.

Both show effect size AND uncertainty AND sample size, so a clade or system is never
read as enriched on the strength of a handful of genomes.
"""
import os
import re
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

apply_style()


def _load(path, name):
    if not path or not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    return None if df.empty else df


# --------------------------------------------------------------------------- #
def plot_prevalence_caterpillar(prevalence, master, outdir):
    prefix = os.path.join(outdir, "prevalence_caterpillar")
    if prevalence is None or not {"rank", "taxon", "prevalence", "n_genomes"}.issubset(prevalence.columns):
        placeholder(prefix, "no data: prevalence_by_taxon.tsv missing/empty")
        return

    # overall prevalence (denominator = pass_primary genomes, present = tier-A Prometheus)
    overall = np.nan
    if master is not None and {"pass_primary", "prometheus_present"}.issubset(master.columns):
        pp = master["pass_primary"].astype(str).str.lower().isin(("true", "1", "yes"))
        sub = master[pp]
        if len(sub):
            overall = sub["prometheus_present"].astype(str).str.lower().isin(("true", "1", "yes")).mean()

    rows = []
    for rank, n_top in (("gtdb_phylum", 15), ("gtdb_family", 20)):
        d = prevalence[prevalence["rank"] == rank].copy()
        if d.empty:
            continue
        d = d.sort_values("n_genomes", ascending=False).head(n_top)
        rows.append(d.assign(_rank=rank))
    if not rows:
        placeholder(prefix, "no data: no phylum/family prevalence rows")
        return
    df = pd.concat(rows, ignore_index=True)
    df = df.sort_values(["_rank", "prevalence"], ascending=[True, True]).reset_index(drop=True)

    ci_low = df.get("ci_low", df["prevalence"])
    ci_high = df.get("ci_high", df["prevalence"])
    # flag taxa whose Wilson CI excludes the overall prevalence
    if np.isfinite(overall):
        enriched = ci_low > overall
        depleted = ci_high < overall
    else:
        enriched = pd.Series(False, index=df.index); depleted = pd.Series(False, index=df.index)
    color = np.where(enriched, OKABE_ITO["blue"],
                     np.where(depleted, OKABE_ITO["vermillion"], OKABE_ITO["grey"]))
    sizes = 6 + 22 * (df["n_genomes"] / max(df["n_genomes"].max(), 1)) ** 0.5

    fig, ax = plt.subplots(figsize=(4.2, max(3.2, len(df) * 0.20)))
    y = np.arange(len(df))
    ax.hlines(y, ci_low.values, ci_high.values, color=color, linewidth=0.9, zorder=2)
    ax.scatter(df["prevalence"].values, y, s=sizes.values, c=color, zorder=4,
               edgecolor="white", linewidth=0.3)
    if np.isfinite(overall):
        ax.axvline(overall, color="black", linestyle="--", linewidth=0.7, zorder=1)
        ax.text(overall, len(df) - 0.5, f" overall={overall:.2f}", fontsize=5, va="top")
    # rank separators / labels
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.taxon}" for r in df.itertuples()], fontsize=4.5)
    ax.set_xlabel("Prometheus prevalence (Wilson 95% CI)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_title("Prevalence by clade (phylum + family); point size ~ n genomes")
    handles = [mpatches.Patch(color=OKABE_ITO["blue"], label="enriched (CI > overall)"),
               mpatches.Patch(color=OKABE_ITO["vermillion"], label="depleted (CI < overall)"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="not distinguishable")]
    ax.legend(handles=handles, frameon=False, fontsize=5, loc="lower right")
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def _parse_or(direction):
    m = re.search(r"OR\s*=\s*([0-9.eE+-]+)", str(direction))
    if not m:
        return np.nan
    try:
        return float(m.group(1))
    except ValueError:
        return np.nan


def plot_cooccurrence_oddsratio(cooc, outdir):
    prefix = os.path.join(outdir, "cooccurrence_oddsratio")
    if cooc is None or "system_b" not in cooc.columns:
        placeholder(prefix, "no data: cooccurrence.tsv missing/empty")
        return
    df = cooc.copy()
    df["OR"] = df["direction"].apply(_parse_or) if "direction" in df.columns else np.nan
    df = df[np.isfinite(df["OR"]) & (df["OR"] > 0)].copy()
    if df.empty:
        placeholder(prefix, "no data: no parseable odds ratios in cooccurrence.tsv")
        return
    df["log2OR"] = np.log2(df["OR"].clip(lower=1e-6))
    df = df.sort_values("log2OR")
    padj = df.get("p_adj", pd.Series(np.nan, index=df.index))
    sig = padj < 0.05
    color = np.where(df["log2OR"] >= 0, OKABE_ITO["blue"], OKABE_ITO["vermillion"])
    nb = df.get("n_both", pd.Series(1, index=df.index)).fillna(1)
    sizes = 8 + 26 * (nb / max(nb.max(), 1)) ** 0.5

    fig, ax = plt.subplots(figsize=(4.6, max(2.6, len(df) * 0.26)))
    y = np.arange(len(df))
    ax.hlines(y, 0, df["log2OR"].values, color=color, linewidth=0.9, zorder=2)
    ax.scatter(df["log2OR"].values, y, s=sizes.values, c=color, zorder=4,
               edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", linewidth=0.6)
    for yi, s in zip(y, sig.values):
        if s:
            ax.text(ax.get_xlim()[1], yi, "*", va="center", ha="left", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(df["system_b"].values, fontsize=5)
    ax.set_xlabel("log2 odds ratio  (>0 co-occur, <0 exclude)")
    ax.set_title("Prometheus vs other defense systems  (* BH p_adj<0.05; size ~ n co-present)")
    save_fig(fig, prefix)


def main():
    ap = argparse.ArgumentParser(description="Extra distribution figures")
    ap.add_argument("--prevalence", required=True)
    ap.add_argument("--cooccurrence", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plot_prevalence_caterpillar(_load(args.prevalence, "prev"), _load(args.master, "master"), args.outdir)
    plot_cooccurrence_oddsratio(_load(args.cooccurrence, "cooc"), args.outdir)


if __name__ == "__main__":
    main()
