"""plot_distribution.py -- prevalence, association forest and co-occurrence heatmap.

Outputs (written to --outdir, one file per --formats extension):
  prevalence_by_taxon.*    Prometheus prevalence (Wilson 95% CI) for the 20 best-sampled
                           GTDB families, primary tier set; n and median completeness shown.
  association_forest.*     methanogen log-odds with 95% CI: naive GLM / Fisher baselines
                           vs phylogenetic logistic regression (unadjusted and the
                           confounder-adjusted HEADLINE). fitPagel and D are not log-odds,
                           so they are printed as text, not drawn on the axis.
  cooccurrence_heatmap.*   signed -log10(BH p_adj) of the fitPagel test per system; sign
                           from the Fisher odds-ratio direction (exclusion < 0 < co-occurrence).

A header-only input table (genuinely empty upstream input) gives a labelled placeholder.
A table that lacks its required columns is an interface error (exit 1).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

apply_style()

LOGODDS_TESTS = ["naive_glm", "naive_fisher", "phyloglm_methanogen", "phyloglm_confounders"]
TEST_LABEL = {
    "naive_glm": "naive logistic GLM (independence)",
    "naive_fisher": "naive Fisher exact (independence)",
    "phyloglm_methanogen": "phylogenetic logistic regression",
    "phyloglm_confounders": "phylogenetic + confounders (HEADLINE)",
}


def load(path, required, numeric=()):
    """Read a TSV; None if it has no data rows; exit 1 if required columns are missing."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[plot_distribution] {path}: missing required column(s) {missing}")
    if df.empty:
        return None
    for c in numeric:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c].replace({"": np.nan, "NA": np.nan}), errors="coerce")
    return df


def direction_sign(direction) -> int:
    """+1 co-occurrence, -1 exclusion, 0 neutral/ambiguous (label format 'exclusion (OR=0.2)')."""
    d = str(direction).strip().lower()
    if d.startswith("co-occurrence") or d.startswith("positive"):
        return 1
    if d.startswith("exclusion") or d.startswith("negative"):
        return -1
    return 0


def plot_prevalence(prev, outdir, formats):
    prefix = os.path.join(outdir, "prevalence_by_taxon")
    if prev is None:
        placeholder(prefix, "no data: prevalence_by_taxon.tsv is empty", formats)
        return
    df = prev[(prev["rank"] == "gtdb_family") & np.isfinite(prev["prevalence"]) & (prev["n_genomes"] > 0)]
    if df.empty:
        placeholder(prefix, "no data: no gtdb_family rows with analysed genomes", formats)
        return
    df = df.sort_values(["n_genomes", "taxon"], ascending=[False, True]).head(20)
    df = df.sort_values(["prevalence", "taxon"], ascending=[True, True])
    lo = (df["prevalence"] - df["ci_low"]).clip(lower=0).fillna(0)
    hi = (df["ci_high"] - df["prevalence"]).clip(lower=0).fillna(0)
    tier = str(df["tier_set"].iloc[0]) if "tier_set" in df.columns else "?"
    fig, ax = plt.subplots(figsize=(3.6, max(3.0, len(df) * 0.22)))
    y = np.arange(len(df))
    ax.barh(y, df["prevalence"].values, height=0.6, color=OKABE_ITO["blue"], alpha=0.85, edgecolor="none")
    ax.errorbar(df["prevalence"].values, y, xerr=[lo.values, hi.values], fmt="none",
                color="black", linewidth=0.8, capsize=2)
    ax.set_yticks(y)
    ax.set_yticklabels(df["taxon"].values, fontsize=5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Prevalence (Wilson 95% CI)")
    ax.set_title(f"Prometheus prevalence by family (tier {tier}; 20 best-sampled)", fontsize=6)
    comp = df["median_completeness"] if "median_completeness" in df.columns else pd.Series(np.nan, index=df.index)
    for yi, (n_g, c) in zip(y, zip(df["n_genomes"], comp)):
        lab = f"n={int(n_g)}" + (f", {c:.0f}%" if np.isfinite(c) else "")
        ax.text(1.01, yi, lab, va="center", fontsize=4, transform=ax.get_yaxis_transform())
    save_fig(fig, prefix, formats)


def plot_association_forest(tests, outdir, formats):
    prefix = os.path.join(outdir, "association_forest")
    if tests is None:
        placeholder(prefix, "no data: association_tests.tsv is empty", formats)
        return
    df = tests[tests["test"].isin(LOGODDS_TESTS) & (tests["term"] == "methanogen")
               & np.isfinite(tests["estimate"])].copy()
    others = tests[~tests["test"].isin(LOGODDS_TESTS)]
    if df.empty:
        placeholder(prefix, "no data: no finite methanogen log-odds estimate", formats)
        return
    df["_order"] = df["test"].map({t: i for i, t in enumerate(LOGODDS_TESTS)})
    df = df.sort_values("_order", ascending=False)
    z = 1.959963984540054
    lo = df["ci_low"].where(np.isfinite(df["ci_low"]), df["estimate"] - z * df["std_error"])
    hi = df["ci_high"].where(np.isfinite(df["ci_high"]), df["estimate"] + z * df["std_error"])
    colors = [OKABE_ITO["grey"] if t.startswith("naive") else OKABE_ITO["blue"] for t in df["test"]]
    fig, ax = plt.subplots(figsize=(4.6, 1.2 + 0.35 * len(df)))
    y = np.arange(len(df))
    ax.hlines(y, lo.values, hi.values, color=colors, linewidth=1.0)
    ax.scatter(df["estimate"].values, y, color=colors, s=[30 if t == "phyloglm_confounders" else 14
                                                          for t in df["test"]], zorder=4)
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--", zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{TEST_LABEL[t]}  (n={int(n)}, p={p:.2g})" if np.isfinite(p)
                        else f"{TEST_LABEL[t]}  (n={int(n)})"
                        for t, n, p in zip(df["test"], df["n"], df["p_value"])], fontsize=5)
    ax.set_xlabel("methanogen effect on Prometheus presence (log-odds, 95% CI)")
    ax.set_title("Methanogen association: naive baselines vs phylogenetic models", fontsize=6)
    txt = ["grey: assumes independence (baseline); blue: phylogenetically controlled"]
    for _, r in others.iterrows():
        if r["test"] == "fitPagel" and np.isfinite(r["p_value"]):
            txt.append(f"Pagel correlated evolution: LR={r['statistic']:.2f}, p={r['p_value']:.2g}")
        if r["test"] == "D_statistic" and np.isfinite(r["estimate"]):
            txt.append(f"Fritz-Purvis D={r['estimate']:.2f} (p vs random={r['p_value']:.2g})")
    fig.subplots_adjust(bottom=min(0.5, 0.16 + 0.07 * len(txt)))
    fig.text(0.01, 0.01, "\n".join(txt), fontsize=5, va="bottom", ha="left")
    save_fig(fig, prefix, formats)


def plot_cooccurrence_heatmap(cooc, outdir, formats):
    prefix = os.path.join(outdir, "cooccurrence_heatmap")
    if cooc is None:
        placeholder(prefix, "no data: cooccurrence.tsv is empty (no testable system)", formats)
        return
    df = cooc[np.isfinite(cooc["p_adj"])].copy()
    if df.empty:
        placeholder(prefix, "no data: no rows with p_adj", formats)
        return
    df = df.sort_values(["p_adj", "system_b"])
    sign = df["direction"].map(direction_sign).values if "direction" in df.columns else np.ones(len(df))
    signed = (-np.log10(df["p_adj"].clip(lower=1e-300)).values * sign).reshape(1, -1)
    vmax = float(np.nanmax(np.abs(signed))) or 1.0
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    fig, ax = plt.subplots(figsize=(max(3.5, len(df) * 0.35), 1.9))
    im = ax.imshow(signed, aspect="auto", cmap="RdBu_r", norm=norm)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["system_b"].values, rotation=55, ha="right", fontsize=5)
    ax.set_yticks([])
    ax.set_title("Prometheus vs other systems: fitPagel, signed by Fisher direction (* BH p_adj<0.05)",
                 fontsize=6)
    cbar = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.35, fraction=0.05)
    cbar.set_label("signed -log10(p_adj)  (<0 exclusion, >0 co-occurrence)", fontsize=5)
    for j, p in enumerate(df["p_adj"].values):
        if p < 0.05:
            ax.text(j, 0, "*", ha="center", va="center", fontsize=7, color="black")
    save_fig(fig, prefix, formats)


def main():
    ap = argparse.ArgumentParser(description="Distribution figures")
    ap.add_argument("--prevalence", required=True, help="distribution/prevalence_by_taxon.tsv")
    ap.add_argument("--tests", required=True, help="distribution/association_tests.tsv")
    ap.add_argument("--cooccurrence", required=True, help="distribution/cooccurrence.tsv")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--formats", default="png,svg", help="comma-separated figure formats (default png,svg)")
    args = ap.parse_args()
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    os.makedirs(args.outdir, exist_ok=True)
    prev = load(args.prevalence, ["rank", "taxon", "n_genomes", "prevalence", "ci_low", "ci_high"],
                ["n_genomes", "n_present", "prevalence", "ci_low", "ci_high", "median_completeness"])
    tests = load(args.tests, ["test", "term", "estimate", "std_error", "p_value", "ci_low", "ci_high", "n"],
                 ["estimate", "std_error", "statistic", "p_value", "ci_low", "ci_high", "n"])
    cooc = load(args.cooccurrence, ["system_b", "p_adj", "direction"], ["p_value", "p_adj", "statistic"])
    plot_prevalence(prev, args.outdir, formats)
    plot_association_forest(tests, args.outdir, formats)
    plot_cooccurrence_heatmap(cooc, args.outdir, formats)


if __name__ == "__main__":
    main()
