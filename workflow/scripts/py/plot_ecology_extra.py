"""Extra ecology figures (rule plot_ecology_extra in plots_extra.smk; CLI frozen, the
--primary-tiers / --formats flags are optional additions).

Outputs (--outdir):
  viral_pressure_prometheus.*  proviral load and R-loop-prone fraction by Prometheus
                               presence (raw distributions) + Prometheus prevalence across
                               R-loop-pressure quantile bins with Wilson 95 % CIs; the
                               phyloglm estimate (the inference) is printed in the title.
  niche_effect_forest.*        phyloglm log-odds of every fitted group_/niche_ predictor
                               with 95 % Wald CIs, coloured by BH p_adj (enriched/depleted);
                               the continuous viral predictors (per-SD effects) are not mixed in.
  prevalence_niche_phylum.*    Prometheus prevalence in each niche x phylum cell (point size
                               ~ genomes in the cell): where the system concentrates,
                               separating habitat from clade.
Presence follows SCHEMA Contract v2 item 6 (tier set --primary-tiers, default A; genomes
with validated loci only outside the set are excluded and counted); only pass_primary
genomes are used. Effects are the phylogeny-controlled ones; raw splits are shown only as
distributions, never as the inferential claim. Header-only inputs give placeholders;
missing or unparsable inputs are errors.
"""
import argparse
import math
import os
import re
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


def _load(path, name=""):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{name or 'input'} not found: {path}")
    if os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return None if df.empty else df


def _num(s):
    return pd.to_numeric(s.replace({"": np.nan, "NA": np.nan}), errors="coerce")


def _bool(s):
    return s.astype(str).str.strip().str.lower().isin(TRUE_S)


def presence(master, tiers):
    """Contract v2 item 6: 1 / 0 / NaN (excluded)."""
    if "prometheus_tier" in master.columns:
        t = master["prometheus_tier"].fillna("").str.strip().str.upper()
        return pd.Series(np.where(t.isin(tiers), 1.0, np.where(t == "", 0.0, np.nan)),
                         index=master.index)
    utils.log("WARNING: master has no prometheus_tier; presence = prometheus_present")
    return _bool(master["prometheus_present"]).astype(float)


def analysis_set(master, tiers, label):
    m = master.copy()
    if "pass_primary" in m.columns:
        m = m[_bool(m["pass_primary"])]
    m["present"] = presence(m, tiers)
    n_ex = int(m["present"].isna().sum())
    m = m[m["present"].notna()].copy()
    m["present"] = m["present"].astype(int)
    utils.log(f"{label}: tiers={'+'.join(tiers)}: {int(m['present'].sum())} present, "
              f"{int((m['present'] == 0).sum())} absent, {n_ex} excluded")
    return m, n_ex


def wilson(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def _assoc_note(assoc, predictor):
    if assoc is None or not {"predictor", "estimate"}.issubset(assoc.columns):
        return ""
    row = assoc[assoc["predictor"] == predictor]
    if row.empty:
        return ""
    est = _num(row["estimate"]).iloc[0]
    if not np.isfinite(est):
        return ""
    q = _num(row["p_adj"]).iloc[0] if "p_adj" in row.columns else np.nan
    return f"phyloglm β={est:.2f}/SD" + (f", BH p={q:.3g}" if np.isfinite(q) else "")


# --------------------------------------------------------------------------- #
def plot_viral_pressure(master, viral, assoc, outdir, tiers, formats):
    prefix = os.path.join(outdir, "viral_pressure_prometheus")
    if master is None or viral is None:
        placeholder(prefix, "no data: master table or viral summary is empty", formats)
        return {}
    if "prometheus_present" not in master.columns and "prometheus_tier" not in master.columns:
        raise ValueError("master table lacks prometheus_tier/prometheus_present")
    m, n_ex = analysis_set(master, tiers, "viral_pressure_prometheus")
    d = m[["genome_id", "present"]].merge(viral, on="genome_id", how="inner")
    for c in ("n_proviruses", "frac_rloop_proviruses"):
        d[c] = _num(d[c]) if c in d.columns else np.nan
    if d["present"].nunique() < 2 or d["n_proviruses"].notna().sum() < 4:
        placeholder(prefix, "no data: need both presence states and >= 4 genomes with viral values",
                    formats)
        return {}
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))
    counts = {"excluded": n_ex}
    for ax, col, lab in ((axes[0], "n_proviruses", "proviruses / genome"),
                         (axes[1], "frac_rloop_proviruses", "fraction R-loop-prone\n(genomes with proviruses)")):
        groups = [d.loc[d["present"] == g, col].dropna().values for g in (0, 1)]
        counts[col] = (len(groups[0]), len(groups[1]))
        pos = [i + 1 for i, g in enumerate(groups) if len(g) > 0]
        data = [g for g in groups if len(g) > 0]
        if data:
            parts = ax.violinplot(data, positions=pos, showmeans=False, showmedians=True, widths=0.8)
            for b, p in zip(parts["bodies"], pos):
                b.set_facecolor(OKABE_ITO["grey"] if p == 1 else OKABE_ITO["blue"])
                b.set_alpha(0.6)
        else:
            ax.text(0.5, 0.5, "no values", transform=ax.transAxes, ha="center", fontsize=6)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"absent\n(n={len(groups[0])})", f"present\n(n={len(groups[1])})"],
                           fontsize=6)
        ax.set_xlim(0.4, 2.6)
        ax.set_ylabel(lab, fontsize=6)
        ax.set_title(_assoc_note(assoc, col) or "phyloglm: not fitted", fontsize=6)
    ax = axes[2]
    dd = d.dropna(subset=["frac_rloop_proviruses"])
    bins = None
    if len(dd) >= 8 and dd["frac_rloop_proviruses"].nunique() >= 2:
        # value-based quantile bins: tied fractions (e.g. many 0.0) stay in one bin
        q = min(4, dd["frac_rloop_proviruses"].nunique())
        bins = pd.qcut(dd["frac_rloop_proviruses"], q=q, labels=False, duplicates="drop")
        if bins.nunique() < 2:
            bins = None
    if bins is not None:
        dd = dd.assign(_bin=bins)
        xs, ys, los, his, labs = [], [], [], [], []
        for b, g in dd.groupby("_bin"):
            p, lo, hi = wilson(int(g["present"].sum()), len(g))
            xs.append(int(b)); ys.append(p); los.append(p - lo); his.append(hi - p)
            labs.append(f"{g['frac_rloop_proviruses'].min():.2f}-{g['frac_rloop_proviruses'].max():.2f}\n(n={len(g)})")
        ax.errorbar(xs, ys, yerr=[los, his], fmt="o-", color=OKABE_ITO["vermillion"],
                    markersize=4, capsize=2, linewidth=1)
        ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=5)
        ax.set_ylim(0, 1)
        counts["bins"] = len(xs)
    else:
        ax.text(0.5, 0.5, "too few genomes with proviruses", transform=ax.transAxes,
                ha="center", fontsize=6)
    ax.set_xlabel("R-loop-prone fraction (quantile bins; ties kept together)", fontsize=6)
    ax.set_ylabel("Prometheus prevalence (Wilson 95 % CI)", fontsize=6)
    fig.suptitle(f"Viral pressure & R-loop targeting (tiers {'+'.join(tiers)}; {n_ex} excluded; "
                 "raw distributions, inference = phyloglm)", fontsize=7)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, prefix, formats)
    utils.log(f"viral_pressure_prometheus: {counts}")
    return counts


# --------------------------------------------------------------------------- #
def plot_niche_forest(assoc, outdir, formats):
    prefix = os.path.join(outdir, "niche_effect_forest")
    if assoc is None or not {"predictor", "estimate"}.issubset(assoc.columns):
        placeholder(prefix, "no data: ecology_associations.tsv is empty", formats)
        return 0
    df = assoc[assoc["predictor"].str.match(r"^(group|niche)_")].copy()
    df["estimate"] = _num(df["estimate"])
    df = df[np.isfinite(df["estimate"])].copy()
    if df.empty:
        placeholder(prefix, "no fitted niche/group effect", formats)
        return 0
    se = _num(df["std_error"]) if "std_error" in df.columns else pd.Series(np.nan, index=df.index)
    lo = _num(df["ci_low"]) if "ci_low" in df.columns else pd.Series(np.nan, index=df.index)
    hi = _num(df["ci_high"]) if "ci_high" in df.columns else pd.Series(np.nan, index=df.index)
    df["lo"] = lo.where(lo.notna(), df["estimate"] - 1.96 * se)
    df["hi"] = hi.where(hi.notna(), df["estimate"] + 1.96 * se)
    padj = _num(df["p_adj"]) if "p_adj" in df.columns else pd.Series(np.nan, index=df.index)
    df["sig"] = (padj < 0.05).fillna(False)
    df["label"] = df["predictor"].str.replace(r"^(niche|group)_", "", regex=True)
    df = df.sort_values("estimate")
    y = np.arange(len(df))
    color = np.where(df["sig"], np.where(df["estimate"] >= 0, OKABE_ITO["blue"], OKABE_ITO["vermillion"]),
                     OKABE_ITO["grey"])
    fig, ax = plt.subplots(figsize=(4.6, max(2.6, len(df) * 0.26 + 0.8)))
    ax.hlines(y, df["lo"].values, df["hi"].values, color=color, linewidth=0.9)
    ax.scatter(df["estimate"].values, y, c=color, s=16, zorder=4, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.6)
    ax.set_yticks(y); ax.set_yticklabels(df["label"].values, fontsize=5)
    if "n" in df.columns:
        xmax = np.nanmax(df["hi"].values) if np.isfinite(df["hi"]).any() else df["estimate"].max()
        for yi, nv in zip(y, df["n"].values):
            ax.text(xmax, yi, f" n={nv}", va="center", fontsize=4, color="#555")
    ax.set_xlabel("phyloglm log-odds, 95 % Wald CI (niche/group → Prometheus)")
    ax.set_title("Ecological associations (phylogeny-controlled)")
    handles = [mpatches.Patch(color=OKABE_ITO["blue"], label="enriched (BH<0.05)"),
               mpatches.Patch(color=OKABE_ITO["vermillion"], label="depleted (BH<0.05)"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="n.s.")]
    ax.legend(handles=handles, frameon=False, fontsize=5, loc="lower right")
    save_fig(fig, prefix, formats)
    return len(df)


# --------------------------------------------------------------------------- #
def plot_prevalence_niche_phylum(master, outdir, tiers, formats):
    prefix = os.path.join(outdir, "prevalence_niche_phylum")
    if master is None:
        placeholder(prefix, "no data: master table is empty", formats)
        return 0
    niche_cols = [c for c in master.columns if c.startswith("niche_") or c.startswith("group_")]
    if not niche_cols or "gtdb_phylum" not in master.columns:
        placeholder(prefix, "no data: no niche/group columns or gtdb_phylum", formats)
        return 0
    m, n_ex = analysis_set(master, tiers, "prevalence_niche_phylum")
    for c in niche_cols:
        m[c] = _num(m[c]).fillna(0).astype(int)
    top_niches = sorted(niche_cols, key=lambda c: (-m[c].sum(), c))[:12]
    top_phyla = m["gtdb_phylum"].value_counts().head(12).index.tolist()
    xs, ys, sizes, colors = [], [], [], []
    for i, niche in enumerate(top_niches):
        for j, ph in enumerate(top_phyla):
            g = m[(m[niche] == 1) & (m["gtdb_phylum"] == ph)]
            if len(g) == 0:
                continue
            xs.append(j); ys.append(i); colors.append(g["present"].mean())
            sizes.append(8 + 30 * (len(g) / max(1, len(m))) ** 0.4)
    if not xs:
        placeholder(prefix, "no data: no niche x phylum cell with genomes", formats)
        return 0
    fig, ax = plt.subplots(figsize=(max(4.5, 0.45 * len(top_phyla) + 2),
                                    max(3.0, 0.34 * len(top_niches) + 1.5)))
    sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="viridis", vmin=0, vmax=1,
                    edgecolor="#888", linewidth=0.3)
    ax.set_xticks(range(len(top_phyla)))
    ax.set_xticklabels(top_phyla, rotation=55, ha="right", fontsize=5)
    ax.set_yticks(range(len(top_niches)))
    ax.set_yticklabels([re.sub("^niche_|^group_", "", n) for n in top_niches], fontsize=5)
    ax.set_title(f"Prometheus prevalence by niche x phylum (tiers {'+'.join(tiers)}; "
                 f"{n_ex} excluded; size ~ n genomes)", fontsize=6)
    fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, label="prevalence")
    save_fig(fig, prefix, formats)
    return len(xs)


def main():
    ap = argparse.ArgumentParser(description="Extra ecology figures")
    ap.add_argument("--master", required=True)
    ap.add_argument("--viral", required=True)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--primary-tiers", default="A", help="tier set T for presence (comma list)")
    ap.add_argument("--formats", default="png,svg")
    args = ap.parse_args()
    tiers = [t.strip().upper() for t in args.primary_tiers.split(",") if t.strip()]
    if not tiers or not set(tiers) <= {"A", "B", "C"}:
        raise SystemExit(f"--primary-tiers must be a subset of A,B,C (got {args.primary_tiers!r})")
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    os.makedirs(args.outdir, exist_ok=True)
    master = _load(args.master, "master table")
    viral = _load(args.viral, "viral summary")
    assoc = _load(args.associations, "ecology associations")
    plot_viral_pressure(master, viral, assoc, args.outdir, tiers, formats)
    n = plot_niche_forest(assoc, args.outdir, formats)
    utils.log(f"niche_effect_forest: {n} predictor(s)")
    plot_prevalence_niche_phylum(master, args.outdir, tiers, formats)


if __name__ == "__main__":
    main()
