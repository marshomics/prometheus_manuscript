"""plot_ecology.py -- Stage 6 figures.

Outputs (--outdir, one file per --formats entry):
  ecology_enrichment.*       phylogenetic-logistic (phyloglm) log-odds of every fitted
                             group_/niche_ predictor with its 95 % Wald interval (ci_low /
                             ci_high, else estimate +/- 1.96 SE), coloured by BH p_adj < 0.05.
                             At most 25 rows (smallest p-values first), sorted by estimate.
  viral_rloop_association.*  proviral load and R-loop-prone provirus fraction by Prometheus
                             presence (raw distributions only); the panel titles carry the
                             phyloglm estimate / p_adj from ecology_associations.tsv -- the
                             phylogeny-controlled inference -- not a naive two-sample test.
Presence follows SCHEMA Contract v2 item 6 for the tier set --primary-tiers (default A):
present = best tier in the set, absent = no validated locus, genomes whose validated loci
are all outside the set are EXCLUDED (and counted in the log). Only pass_primary genomes
are shown, as in the statistics. A header-only input gives a labelled placeholder; a
missing or unparsable input is an error.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO  # noqa: E402
import utils  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patches as mpatches  # noqa: E402

apply_style()

TRUE_S = {"1", "true", "yes", "t"}
MAX_ROWS = 25


def load_table(path, name):
    """DataFrame of strings; None for an empty / header-only table (genuinely empty input)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"{name} not found: {path}")
    if os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return None if df.empty else df


def numeric(s):
    return pd.to_numeric(s.replace({"": np.nan, "NA": np.nan}), errors="coerce")


def presence(master, tiers):
    """Contract v2 item 6: 1 / 0 / NaN (excluded) per genome."""
    if "prometheus_tier" in master.columns:
        t = master["prometheus_tier"].fillna("").str.strip().str.upper()
        return pd.Series(np.where(t.isin(tiers), 1.0, np.where(t == "", 0.0, np.nan)),
                         index=master.index)
    if "prometheus_present" in master.columns:
        utils.log("WARNING: master has no prometheus_tier; presence = prometheus_present")
        return master["prometheus_present"].str.strip().str.lower().isin(TRUE_S).astype(float)
    raise ValueError("master has neither prometheus_tier nor prometheus_present")


def analysis_set(master, tiers, label):
    m = master.copy()
    if "pass_primary" in m.columns:
        m = m[m["pass_primary"].str.strip().str.lower().isin(TRUE_S)]
    m["present"] = presence(m, tiers)
    n_ex = int(m["present"].isna().sum())
    m = m[m["present"].notna()].copy()
    m["present"] = m["present"].astype(int)
    utils.log(f"{label}: tiers={'+'.join(tiers)}: {int((m['present'] == 1).sum())} present, "
              f"{int((m['present'] == 0).sum())} absent, {n_ex} excluded (validated locus "
              f"outside the tier set)")
    return m, n_ex


def wald_interval(df):
    est = numeric(df["estimate"])
    lo = numeric(df["ci_low"]) if "ci_low" in df.columns else pd.Series(np.nan, index=df.index)
    hi = numeric(df["ci_high"]) if "ci_high" in df.columns else pd.Series(np.nan, index=df.index)
    se = numeric(df["std_error"]) if "std_error" in df.columns else pd.Series(np.nan, index=df.index)
    lo = lo.where(lo.notna(), est - 1.96 * se)
    hi = hi.where(hi.notna(), est + 1.96 * se)
    return est, lo, hi


def plot_ecology_enrichment(assoc, outdir, formats):
    prefix = os.path.join(outdir, "ecology_enrichment")
    if assoc is None or not {"predictor", "estimate"}.issubset(assoc.columns):
        placeholder(prefix, "no data: ecology_associations.tsv is empty", formats)
        return 0
    df = assoc[assoc["predictor"].str.match(r"^(group|niche)_")].copy()
    df["est"], df["lo"], df["hi"] = wald_interval(df)
    df = df[np.isfinite(df["est"])]
    if df.empty:
        placeholder(prefix, "no fitted group_/niche_ predictor (see ecology_associations.tsv notes)",
                    formats)
        return 0
    df["p"] = numeric(df["p_value"]) if "p_value" in df.columns else np.nan
    df["q"] = numeric(df["p_adj"]) if "p_adj" in df.columns else np.nan
    if len(df) > MAX_ROWS:
        df = df.sort_values(["p", "est"], na_position="last").head(MAX_ROWS)
    df = df.sort_values("est")
    sig = (df["q"] < 0.05).fillna(False).values
    colors = [OKABE_ITO["vermillion"] if s else OKABE_ITO["grey"] for s in sig]
    nlab = df["n"] if "n" in df.columns else pd.Series("", index=df.index)
    labels = [f"{p}  (n={n})" if n not in ("", "NA") else p for p, n in zip(df["predictor"], nlab)]
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(4.4, max(2.5, len(df) * 0.28 + 0.8)))
    ax.hlines(y, df["lo"].values, df["hi"].values, color=colors, linewidth=0.8)
    ax.scatter(df["est"].values, y, color=colors, s=14, zorder=4)
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5)
    ax.set_xlabel("phyloglm log-odds (95 % Wald CI)")
    ax.set_title("Niche / group associations (phylogeny-controlled)")
    handles = [mpatches.Patch(color=OKABE_ITO["vermillion"], label="BH p_adj < 0.05"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="n.s.")]
    ax.legend(handles=handles, frameon=False, fontsize=6, loc="best")
    save_fig(fig, prefix, formats)
    return len(df)


def _phyloglm_title(assoc, predictor):
    if assoc is None or "predictor" not in assoc.columns:
        return "phyloglm: no result"
    row = assoc[assoc["predictor"] == predictor]
    if row.empty:
        return "phyloglm: not tested"
    est = numeric(row["estimate"]).iloc[0]
    if not np.isfinite(est):
        return "phyloglm: not fitted"
    q = numeric(row["p_adj"]).iloc[0] if "p_adj" in row.columns else np.nan
    return f"phyloglm β={est:.2f}/SD" + (f", BH p={q:.2g}" if np.isfinite(q) else "")


def plot_viral_rloop(viral, master, assoc, outdir, formats, tiers):
    prefix = os.path.join(outdir, "viral_rloop_association")
    if viral is None or master is None:
        placeholder(prefix, "no data: viral_summary.tsv or master_table.tsv is empty", formats)
        return {}
    need_v = {"genome_id", "n_proviruses", "frac_rloop_proviruses"}
    if not need_v.issubset(viral.columns):
        raise ValueError(f"viral summary lacks {sorted(need_v - set(viral.columns))}")
    m, n_ex = analysis_set(master, tiers, "viral_rloop_association")
    d = m[["genome_id", "present"]].merge(viral, on="genome_id", how="inner")
    if d.empty or d["present"].nunique() < 1:
        placeholder(prefix, "no data: no genome shared by the master and viral tables", formats)
        return {}
    panels = [("n_proviruses", "proviruses per genome"),
              ("frac_rloop_proviruses", "fraction of proviruses R-loop-prone\n(genomes with >= 1 provirus)")]
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.0))
    counts = {"excluded": n_ex}
    for k, (ax, (col, ylabel)) in enumerate(zip(axes, panels)):
        vals = numeric(d[col])
        groups = [vals[d["present"] == g].dropna().values for g in (0, 1)]
        counts[col] = (len(groups[0]), len(groups[1]))
        rng = np.random.default_rng(42 + k)
        for xi, (g, c) in enumerate(zip(groups, (OKABE_ITO["grey"], OKABE_ITO["blue"]))):
            if len(g) == 0:
                continue
            ax.boxplot(g, positions=[xi], widths=0.5, patch_artist=True, showfliers=False,
                       boxprops=dict(facecolor=c, alpha=0.5),
                       medianprops=dict(color="black", linewidth=0.8),
                       whiskerprops=dict(linewidth=0.7), capprops=dict(linewidth=0.7))
            ax.scatter(xi + rng.uniform(-0.18, 0.18, size=len(g)), g, s=4, alpha=0.45,
                       color=c, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Prometheus−\n(n={len(groups[0])})",
                            f"Prometheus+\n(n={len(groups[1])})"], fontsize=6)
        ax.set_xlim(-0.6, 1.6)
        ax.set_ylabel(ylabel, fontsize=6)
        ax.set_title(_phyloglm_title(assoc, col), fontsize=6)
    fig.suptitle(f"Viral pressure and R-loop-prone proviruses (tiers {'+'.join(tiers)}; "
                 f"{n_ex} excluded)", fontsize=7)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, prefix, formats)
    utils.log(f"viral_rloop_association: n absent/present per panel {counts}")
    return counts


def main():
    ap = argparse.ArgumentParser(description="Stage 6 ecology figures")
    ap.add_argument("--associations", required=True, help="ecology/ecology_associations.tsv")
    ap.add_argument("--viral", required=True, help="ecology/viral/viral_summary.tsv")
    ap.add_argument("--master", required=True, help="scaffolds/master_table.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--primary-tiers", default="A", help="tier set T for presence (comma list)")
    ap.add_argument("--formats", default="png,svg", help="comma list of figure formats")
    args = ap.parse_args()
    tiers = [t.strip().upper() for t in args.primary_tiers.split(",") if t.strip()]
    if not tiers or not set(tiers) <= {"A", "B", "C"}:
        raise SystemExit(f"--primary-tiers must be a subset of A,B,C (got {args.primary_tiers!r})")
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    os.makedirs(args.outdir, exist_ok=True)
    assoc = load_table(args.associations, "ecology_associations")
    viral = load_table(args.viral, "viral_summary")
    master = load_table(args.master, "master_table")
    n = plot_ecology_enrichment(assoc, args.outdir, formats)
    utils.log(f"ecology_enrichment: {n} predictor(s) plotted")
    plot_viral_rloop(viral, master, assoc, args.outdir, formats, tiers)


if __name__ == "__main__":
    main()
