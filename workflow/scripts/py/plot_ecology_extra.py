"""Extra ecology figures (added alongside plot_ecology.py).

Outputs (written to --outdir):
  viral_pressure_prometheus.{png,svg}  proviral load + R-loop-prone fraction split by
                                       Prometheus presence, plus Prometheus prevalence
                                       binned by R-loop pressure with Wilson CIs (the
                                       anti-phage / R-loop-targeting hypothesis).
  niche_effect_forest.{png,svg}        phylogenetic-logistic effect of each environment/
                                       niche predictor on Prometheus presence, with CIs
                                       and BH significance (phylogeny-controlled).
  prevalence_niche_phylum.{png,svg}    prevalence bubble grid (niche x phylum): where the
                                       system concentrates, separating habitat from clade.

Effects are the phylogeny-controlled ones; raw splits are shown only as distributions,
never as the inferential claim.
"""
import os
import re
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

apply_style()


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


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0, c - h), min(1, c + h)


# --------------------------------------------------------------------------- #
def plot_viral_pressure(master, viral, assoc, outdir):
    prefix = os.path.join(outdir, "viral_pressure_prometheus")
    if master is None or "prometheus_present" not in master.columns:
        placeholder(prefix, "no data: master table missing/empty")
        return
    if viral is None:
        placeholder(prefix, "no data: viral_summary.tsv missing/empty")
        return
    m = master.copy()
    m["genome_id"] = m["genome_id"].astype(str)
    viral = viral.copy(); viral["genome_id"] = viral["genome_id"].astype(str)
    d = m.merge(viral, on="genome_id", how="left")
    d["present"] = _bool(d["prometheus_present"]).astype(int)
    for c in ("n_proviruses", "frac_rloop_proviruses"):
        if c not in d.columns:
            d[c] = np.nan
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if d["present"].nunique() < 2 or d["n_proviruses"].notna().sum() < 4:
        placeholder(prefix, "no data: need presence variation and viral values")
        return

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.7))
    # A: proviral load by presence
    for ax, col, lab in ((axes[0], "n_proviruses", "proviruses / genome"),
                         (axes[1], "frac_rloop_proviruses", "fraction R-loop-prone")):
        groups = [d.loc[d["present"] == g, col].dropna().values for g in (0, 1)]
        parts = ax.violinplot(groups, showmeans=False, showmedians=True, widths=0.8)
        for b, c in zip(parts["bodies"], (OKABE_ITO["grey"], OKABE_ITO["blue"])):
            b.set_facecolor(c); b.set_alpha(0.6)
        ax.set_xticks([1, 2]); ax.set_xticklabels(["absent", "present"], fontsize=6)
        ax.set_ylabel(lab); ax.set_title(lab, fontsize=6)
    # C: prevalence binned by R-loop pressure with Wilson CIs
    ax = axes[2]
    dd = d.dropna(subset=["frac_rloop_proviruses"])
    if len(dd) >= 8:
        try:
            dd = dd.assign(_bin=pd.qcut(dd["frac_rloop_proviruses"], q=min(4, dd["frac_rloop_proviruses"].nunique()),
                                        duplicates="drop"))
        except Exception:
            dd = dd.assign(_bin=pd.cut(dd["frac_rloop_proviruses"], bins=4))
        xs, ys, los, his = [], [], [], []
        for i, (b, g) in enumerate(dd.groupby("_bin", observed=True)):
            p, lo, hi = wilson(int(g["present"].sum()), len(g))
            xs.append(i); ys.append(p); los.append(p - lo); his.append(hi - p)
        ax.errorbar(xs, ys, yerr=[los, his], fmt="o-", color=OKABE_ITO["vermillion"],
                    markersize=4, capsize=2, linewidth=1)
        ax.set_xticks(xs); ax.set_xticklabels([f"Q{i+1}" for i in xs], fontsize=6)
        ax.set_ylim(0, 1)
    ax.set_xlabel("R-loop pressure quartile"); ax.set_ylabel("Prometheus prevalence")
    ax.set_title("prevalence vs R-loop pressure", fontsize=6)
    # annotate phyloglm effect if available
    note = ""
    if assoc is not None and {"predictor", "estimate", "p_value"}.issubset(assoc.columns):
        row = assoc[assoc["predictor"].astype(str) == "frac_rloop_proviruses"]
        if len(row):
            est = pd.to_numeric(row["estimate"], errors="coerce").iloc[0]
            p = pd.to_numeric(row.get("p_adj", row["p_value"]), errors="coerce").iloc[0]
            if np.isfinite(est):
                note = f"phyloglm β={est:.2f}" + (f", p={p:.3g}" if np.isfinite(p) else "")
    fig.suptitle("Viral pressure & R-loop targeting" + (f"   ({note})" if note else ""), fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_niche_forest(assoc, outdir):
    prefix = os.path.join(outdir, "niche_effect_forest")
    if assoc is None or not {"predictor", "estimate"}.issubset(assoc.columns):
        placeholder(prefix, "no data: ecology_associations.tsv missing/empty")
        return
    df = assoc.copy()
    df["estimate"] = pd.to_numeric(df["estimate"], errors="coerce")
    df = df[np.isfinite(df["estimate"])].copy()
    if df.empty:
        placeholder(prefix, "no data: no fitted ecology effects")
        return
    se = pd.to_numeric(df.get("std_error", np.nan), errors="coerce")
    df["lo"] = df["estimate"] - 1.96 * se
    df["hi"] = df["estimate"] + 1.96 * se
    padj = pd.to_numeric(df.get("p_adj", df.get("p_value", np.nan)), errors="coerce")
    df["sig"] = padj < 0.05
    df["label"] = df["predictor"].astype(str).str.replace("^niche_|^group_", "", regex=True)
    df = df.sort_values("estimate")
    y = np.arange(len(df))
    color = np.where(df["sig"].fillna(False),
                     np.where(df["estimate"] >= 0, OKABE_ITO["blue"], OKABE_ITO["vermillion"]),
                     OKABE_ITO["grey"])

    fig, ax = plt.subplots(figsize=(4.6, max(2.6, len(df) * 0.26)))
    ax.hlines(y, df["lo"].values, df["hi"].values, color=color, linewidth=0.9)
    ax.scatter(df["estimate"].values, y, c=color, s=16, zorder=4, edgecolor="white", linewidth=0.3)
    ax.axvline(0, color="black", linestyle="--", linewidth=0.6)
    ax.set_yticks(y); ax.set_yticklabels(df["label"].values, fontsize=5)
    if "n" in df.columns:
        for yi, nv in zip(y, df["n"].values):
            ax.text(ax.get_xlim()[1], yi, f" n={nv}", va="center", fontsize=4, color="#555")
    ax.set_xlabel("phyloglm log-odds (environment/niche → Prometheus)")
    ax.set_title("Ecological associations (phylogeny-controlled)")
    handles = [mpatches.Patch(color=OKABE_ITO["blue"], label="enriched (BH<0.05)"),
               mpatches.Patch(color=OKABE_ITO["vermillion"], label="depleted (BH<0.05)"),
               mpatches.Patch(color=OKABE_ITO["grey"], label="n.s.")]
    ax.legend(handles=handles, frameon=False, fontsize=5, loc="lower right")
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_prevalence_niche_phylum(master, outdir):
    prefix = os.path.join(outdir, "prevalence_niche_phylum")
    if master is None or "prometheus_present" not in master.columns:
        placeholder(prefix, "no data: master table missing/empty")
        return
    m = master.copy()
    niche_cols = [c for c in m.columns if c.startswith("niche_") or c.startswith("group_")]
    if not niche_cols or "gtdb_phylum" not in m.columns:
        placeholder(prefix, "no data: no niche/group columns or gtdb_phylum")
        return
    if "pass_primary" in m.columns:
        m = m[_bool(m["pass_primary"])]
    m["present"] = _bool(m["prometheus_present"]).astype(int)
    for c in niche_cols:
        m[c] = pd.to_numeric(m[c], errors="coerce").fillna(0).astype(int)
    # pick top niches (by occupancy) and top phyla (by n genomes)
    top_niches = sorted(niche_cols, key=lambda c: -m[c].sum())[:12]
    top_phyla = m["gtdb_phylum"].value_counts().head(12).index.tolist()
    xs, ys, sizes, colors = [], [], [], []
    for i, niche in enumerate(top_niches):
        for j, ph in enumerate(top_phyla):
            g = m[(m[niche] == 1) & (m["gtdb_phylum"] == ph)]
            if len(g) == 0:
                continue
            prev = g["present"].mean()
            xs.append(j); ys.append(i); sizes.append(8 + 30 * (len(g) / max(1, len(m)) ) ** 0.4)
            colors.append(prev)
    if not xs:
        placeholder(prefix, "no data: no niche x phylum cells with genomes")
        return
    fig, ax = plt.subplots(figsize=(max(4.5, 0.45 * len(top_phyla) + 2),
                                    max(3.0, 0.34 * len(top_niches) + 1.5)))
    sc = ax.scatter(xs, ys, s=sizes, c=colors, cmap="viridis", vmin=0, vmax=1,
                    edgecolor="#888", linewidth=0.3)
    ax.set_xticks(range(len(top_phyla))); ax.set_xticklabels(top_phyla, rotation=55, ha="right", fontsize=5)
    ax.set_yticks(range(len(top_niches)))
    ax.set_yticklabels([re.sub("^niche_|^group_", "", n) for n in top_niches], fontsize=5)
    ax.set_title("Prometheus prevalence by niche x phylum (size ~ n genomes)")
    fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, label="prevalence")
    save_fig(fig, prefix)


def main():
    ap = argparse.ArgumentParser(description="Extra ecology figures")
    ap.add_argument("--master", required=True)
    ap.add_argument("--viral", required=True)
    ap.add_argument("--associations", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    master = _load(args.master); viral = _load(args.viral); assoc = _load(args.associations)
    plot_viral_pressure(master, viral, assoc, args.outdir)
    plot_niche_forest(assoc, args.outdir)
    plot_prevalence_niche_phylum(master, args.outdir)


if __name__ == "__main__":
    main()
