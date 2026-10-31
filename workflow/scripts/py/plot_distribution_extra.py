"""Extra distribution figures (added alongside plot_distribution.py).

Outputs (written to --outdir, one file per --formats extension):
  prevalence_caterpillar.*   prevalence +/- Wilson 95% CI per taxon (15 best-sampled phyla +
                             20 best-sampled families), point size ~ n_genomes. Taxa whose CI
                             lies above / below the OVERALL prevalence (same tier set and
                             denominator) are coloured; this is descriptive -- not corrected for
                             multiple comparisons or phylogeny.
  cooccurrence_oddsratio.*   Prometheus vs each system: log2 odds ratio from the 2x2 counts
                             (Haldane-Anscombe +0.5 when a cell is 0, so perfect co-occurrence
                             or exclusion stays on the axis), * = fitPagel BH p_adj < 0.05.

The overall prevalence is the prevalence table's rank == "overall" row; for older tables
without it, it is recomputed from --master with the same tier rule (SCHEMA Contract v2
item 6: present if prometheus_tier is in the table's tier_set, absent if "", excluded
otherwise; denominator pass_primary genomes).
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

TRUE = {"true", "t", "1", "yes", "y", "1.0"}


def _load(path, required=(), numeric=()):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[plot_distribution_extra] {path}: missing required column(s) {missing}")
    if df.empty:
        return None
    for c in numeric:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c].replace({"": np.nan, "NA": np.nan}), errors="coerce")
    return df


def overall_prevalence(prev, master):
    """Overall prevalence on the prevalence table's own tier set and denominator."""
    ov = prev[prev["rank"] == "overall"]
    if len(ov) and np.isfinite(ov["prevalence"].iloc[0]):
        return float(ov["prevalence"].iloc[0])
    if master is None or not {"pass_primary", "prometheus_tier"}.issubset(master.columns):
        return np.nan
    tiers = set(str(prev["tier_set"].iloc[0]).split(",")) if "tier_set" in prev.columns else {"A"}
    sub = master[master["pass_primary"].astype(str).str.strip().str.lower().isin(TRUE)]
    t = sub["prometheus_tier"].astype(str).str.strip()
    keep = t.isin(tiers) | (t == "")
    return float(t[keep].isin(tiers).mean()) if keep.any() else np.nan


def plot_prevalence_caterpillar(prev, master, outdir, formats):
    prefix = os.path.join(outdir, "prevalence_caterpillar")
    if prev is None:
        placeholder(prefix, "no data: prevalence_by_taxon.tsv is empty", formats)
        return
    overall = overall_prevalence(prev, master)
    rows = []
    for rank, n_top in (("gtdb_phylum", 15), ("gtdb_family", 20)):
        d = prev[(prev["rank"] == rank) & np.isfinite(prev["prevalence"]) & (prev["n_genomes"] > 0)]
        if len(d):
            rows.append(d.sort_values(["n_genomes", "taxon"], ascending=[False, True]).head(n_top).assign(_rank=rank))
    if not rows:
        placeholder(prefix, "no data: no phylum/family rows with analysed genomes", formats)
        return
    df = pd.concat(rows, ignore_index=True)
    df["_rank_order"] = df["_rank"].map({"gtdb_family": 0, "gtdb_phylum": 1})
    df = df.sort_values(["_rank_order", "prevalence", "taxon"]).reset_index(drop=True)
    lo, hi = df["ci_low"], df["ci_high"]
    if np.isfinite(overall):
        enriched, depleted = (lo > overall).values, (hi < overall).values
    else:
        enriched = depleted = np.zeros(len(df), dtype=bool)
    color = np.where(enriched, OKABE_ITO["blue"], np.where(depleted, OKABE_ITO["vermillion"], OKABE_ITO["grey"]))
    sizes = 6 + 22 * (df["n_genomes"] / max(df["n_genomes"].max(), 1)) ** 0.5
    fig, ax = plt.subplots(figsize=(4.4, max(3.2, len(df) * 0.20)))
    y = np.arange(len(df))
    ax.hlines(y, lo.values, hi.values, color=color, linewidth=0.9, zorder=2)
    ax.scatter(df["prevalence"].values, y, s=sizes.values, c=color, zorder=4, edgecolor="white", linewidth=0.3)
    handles = [mpatches.Patch(color=OKABE_ITO["blue"], label="CI above overall (descriptive)"),
               mpatches.Patch(color=OKABE_ITO["vermillion"], label="CI below overall (descriptive)"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="CI includes overall")]
    if np.isfinite(overall):
        handles.append(ax.axvline(overall, color="black", linestyle="--", linewidth=0.7, zorder=1,
                                  label=f"overall prevalence = {overall:.2f}"))
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t} [{'P' if r == 'gtdb_phylum' else 'F'}]" for t, r in zip(df["taxon"], df["_rank"])],
                       fontsize=4.5)
    tier = str(df["tier_set"].iloc[0]) if "tier_set" in df.columns else "?"
    ax.set_xlabel(f"Prometheus prevalence, tier {tier} (Wilson 95% CI)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_title("Prevalence by clade ([P] phylum, [F] family); point size ~ n genomes", fontsize=6)
    ax.legend(handles=handles, frameon=False, fontsize=5, loc="best")
    save_fig(fig, prefix, formats)


def _parse_or(direction):
    m = re.search(r"OR\s*=\s*([0-9.eE+-]+|Inf)", str(direction))
    if not m:
        return np.nan
    return np.inf if m.group(1) == "Inf" else float(m.group(1))


def log2_or_from_counts(n, n_a, n_b, n_both):
    """Sample log2 odds ratio from 2x2 counts, +0.5 to every cell when any cell is 0."""
    a = n_both
    b = n_a - n_both
    c = n_b - n_both
    d = n - n_a - n_b + n_both
    cells = np.array([a, b, c, d], dtype=float)
    if np.any(cells < 0) or not np.all(np.isfinite(cells)):
        return np.nan
    if np.any(cells == 0):
        cells = cells + 0.5
    a, b, c, d = cells
    return float(np.log2((a * d) / (b * c)))


def plot_cooccurrence_oddsratio(cooc, outdir, formats):
    prefix = os.path.join(outdir, "cooccurrence_oddsratio")
    if cooc is None:
        placeholder(prefix, "no data: cooccurrence.tsv is empty (no testable system)", formats)
        return
    df = cooc.copy()
    if {"n", "n_a", "n_b", "n_both"}.issubset(df.columns):
        df["log2OR"] = [log2_or_from_counts(*v) for v in df[["n", "n_a", "n_b", "n_both"]].to_numpy()]
    else:  # tables written before the count columns existed
        orv = df["direction"].map(_parse_or)
        df["log2OR"] = np.where(np.isfinite(orv) & (orv > 0), np.log2(orv.clip(lower=1e-6)), np.nan)
    df = df[np.isfinite(df["log2OR"])].sort_values(["log2OR", "system_b"])
    if df.empty:
        placeholder(prefix, "no data: no odds ratio computable from cooccurrence.tsv", formats)
        return
    sig = (df["p_adj"] < 0.05).values if "p_adj" in df.columns else np.zeros(len(df), dtype=bool)
    color = np.where(df["log2OR"] >= 0, OKABE_ITO["blue"], OKABE_ITO["vermillion"])
    nb = df["n_both"].fillna(1) if "n_both" in df.columns else pd.Series(1, index=df.index)
    sizes = 8 + 26 * (nb / max(nb.max(), 1)) ** 0.5
    fig, ax = plt.subplots(figsize=(4.6, max(2.6, len(df) * 0.26)))
    y = np.arange(len(df))
    ax.hlines(y, 0, df["log2OR"].values, color=color, linewidth=0.9, zorder=2)
    ax.scatter(df["log2OR"].values, y, s=sizes.values, c=color, zorder=4, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", linewidth=0.6)
    xmax = ax.get_xlim()[1]
    for yi, s in zip(y, sig):
        if s:
            ax.text(xmax, yi, "*", va="center", ha="left", fontsize=7)
    ax.set_yticks(y)
    ax.set_yticklabels(df["system_b"].values, fontsize=5)
    ax.set_xlabel("log2 odds ratio (>0 co-occur, <0 exclude; +0.5 when a cell is 0)")
    ax.set_title("Prometheus vs other defense systems (* fitPagel BH p_adj<0.05; size ~ n co-present)",
                 fontsize=6)
    save_fig(fig, prefix, formats)


def main():
    ap = argparse.ArgumentParser(description="Extra distribution figures")
    ap.add_argument("--prevalence", required=True)
    ap.add_argument("--cooccurrence", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--formats", default="png,svg", help="comma-separated figure formats (default png,svg)")
    args = ap.parse_args()
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    os.makedirs(args.outdir, exist_ok=True)
    prev = _load(args.prevalence, ["rank", "taxon", "n_genomes", "prevalence", "ci_low", "ci_high"],
                 ["n_genomes", "n_present", "prevalence", "ci_low", "ci_high"])
    cooc = _load(args.cooccurrence, ["system_b", "p_adj", "direction"],
                 ["n", "n_a", "n_b", "n_both", "p_value", "p_adj"])
    master = _load(args.master)
    plot_prevalence_caterpillar(prev, master, args.outdir, formats)
    plot_cooccurrence_oddsratio(cooc, args.outdir, formats)


if __name__ == "__main__":
    main()
