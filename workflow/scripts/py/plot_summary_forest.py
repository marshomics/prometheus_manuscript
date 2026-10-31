"""plot_summary_forest.py -- Stage 7 summary forest of the CONFIRMATORY results.

Reads synthesis/results_registry.tsv and draws one panel per confirmatory row. The
confirmatory estimates live on different scales (log-odds of the phyloglm methanogen
effect, the RELAX intensity ratio K, transfer and gain counts), so each panel has its
OWN x-axis and its own null reference (registry `null_value`: 0 for log-odds, 1 for K
on a log axis, the Mk-simulated null mean for the gain enrichment; none for a count
without a test). A shared axis would place a log-odds of 1 next to a count of 30.
Confirmatory rows whose estimate could not be computed are shown as 'not estimable'
(with the registry note), never silently dropped.

Outputs (--outdir): summary_forest.{formats}
A header-only registry or one with no confirmatory row (genuinely empty input) gives a
labelled placeholder; a missing registry or one without the SCHEMA columns is an error.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO

import argparse
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

apply_style()

REQUIRED = ["family", "test", "term", "estimate", "p_value", "p_adj", "hypothesis"]
LOG_SCALES = {"ratio", "odds_ratio"}
FAMILY_ORDER = ["association", "ancestral", "reconciliation", "selection",
                "cooccurrence", "ecology", "disorder"]
LABELS = {  # (family, test, term) -> readable panel label
    ("association", "phyloglm_confounders", "methanogen"):
        "Methanogen association (phyloglm,\nadj. completeness, genome size, defense load)",
    ("association", "phyloglm_methanogen", "methanogen"):
        "Methanogen association (phyloglm, unadjusted)",
    ("association", "fitPagel", "correlated_evolution"):
        "Correlated evolution (Pagel 1994)",
    ("ancestral", "methanogen_gain_enrichment",
     "methanogen_independent_acquisition_enrichment"):
        "Gains on methanogen branches\n(vs Mk-structured null)",
    ("selection", "RELAX", "K"): "Selection intensity on methanogen\nbranches (RELAX K)",
    ("reconciliation", "generax", "transfer"): "Horizontal transfers (GeneRax DTL)",
}
SCALE_LABELS = {
    "log_odds": "log-odds (null 0)",
    "ratio": "ratio, log axis (null 1)",
    "odds_ratio": "odds ratio, log axis (null 1)",
    "count": "count",
    "coefficient": "coefficient (null 0)",
    "D_statistic": "D (0 Brownian, 1 random)",
}


def fnum(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(f) else f


def load_registry(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"registry not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{path} lacks registry column(s) {missing}")
    return df


def confirmatory_rows(df):
    conf = df[df["hypothesis"].str.strip() == "confirmatory"].copy()
    order = {f: i for i, f in enumerate(FAMILY_ORDER)}
    conf["_o"] = conf["family"].map(lambda f: order.get(f, len(order)))
    return conf.sort_values("_o", kind="stable").drop(columns="_o")


def panel_label(r):
    key = (r["family"], r["test"], r["term"])
    return LABELS.get(key, f"{r['family']}: {r['test']} ({r['term']})")


def p_text(r):
    p, padj = fnum(r.get("p_value")), fnum(r.get("p_adj"))
    if p is None:
        return "estimate only\n(no test)"
    m = r.get("fdr_m", "")
    meth = r.get("fdr_method", "") or "BH"
    tail = f" ({meth}, m={m})" if m not in ("", None) else ""
    return f"p = {p:.2g}\np_adj = {padj:.2g}{tail}" if padj is not None else f"p = {p:.2g}"


def draw_panel(ax, r, alpha):
    est = fnum(r.get("estimate"))
    lo, hi = fnum(r.get("ci_low")), fnum(r.get("ci_high"))
    null = fnum(r.get("null_value"))
    scale = r.get("estimate_scale", "") or ""
    padj = fnum(r.get("p_adj"))
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    if est is None or math.isinf(est):
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
        why = (r.get("note", "") or "no estimate")[:70]
        ax.text(0.5, 0.5, f"not estimable: {why}", transform=ax.transAxes,
                ha="center", va="center", fontsize=5.5, color=OKABE_ITO["vermillion"])
        return "not_estimable"
    vals = [v for v in (est, lo, hi, null) if v is not None and math.isfinite(v)]
    logx = scale in LOG_SCALES and all(v > 0 for v in vals)
    if logx:
        ax.set_xscale("log")
        lv = np.log10(vals)
        span = max(lv.max() - lv.min(), 0.3)
        ax.set_xlim(10 ** (lv.min() - 0.15 * span), 10 ** (lv.max() + 0.15 * span))
        ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
    else:
        if scale == "count" and min(vals) >= 0:
            vals = vals + [0.0]           # counts are read against zero
        vmin, vmax = min(vals), max(vals)
        span = (vmax - vmin) or max(abs(vmax), 1.0)
        ax.set_xlim(vmin - 0.15 * span if vmin < 0 or scale != "count" else 0.0,
                    vmax + 0.15 * span)
    if null is not None:
        ax.axvline(null, color="black", lw=0.6, ls="--", zorder=1)
    tested = fnum(r.get("p_value")) is not None
    sig = tested and padj is not None and padj < alpha
    col = OKABE_ITO["blue"] if sig else (OKABE_ITO["grey"] if tested else "black")
    if lo is not None and hi is not None and math.isfinite(lo) and math.isfinite(hi):
        ax.hlines(0, lo, hi, color=col, lw=1.1, zorder=2)
    ax.plot([est], [0], marker="o", ms=4.5, zorder=3, color=col,
            markerfacecolor=col if tested else "white", markeredgewidth=0.9)
    ax.set_ylim(-1, 1)
    ax.tick_params(axis="x", labelsize=5.5, length=2)
    ax.set_xlabel(SCALE_LABELS.get(scale, scale or "estimate"), fontsize=5.5, labelpad=1)
    return "significant" if sig else ("tested" if tested else "estimate_only")


def plot_summary_forest(registry, outdir, alpha=0.05, formats=("png", "svg")):
    prefix = os.path.join(outdir, "summary_forest")
    if registry is None or registry.empty:
        placeholder(prefix, "no data: results registry has no tests", formats)
        return []
    conf = confirmatory_rows(registry)
    if conf.empty:
        placeholder(prefix, "no confirmatory result in the results registry", formats)
        return []
    n = len(conf)
    fig, axes = plt.subplots(n, 1, figsize=(6.3, 0.35 + 0.62 * n), squeeze=False)
    states = []
    for ax, (_, r) in zip(axes[:, 0], conf.iterrows()):
        states.append(draw_panel(ax, r, alpha))
        ax.text(-0.04, 0.5, panel_label(r), transform=ax.transAxes, ha="right",
                va="center", fontsize=6)
        ax.text(1.03, 0.5, p_text(r), transform=ax.transAxes, ha="left",
                va="center", fontsize=5.5)
    handles = [
        Line2D([], [], marker="o", ls="", color=OKABE_ITO["blue"],
               label=f"p_adj < {alpha:g}"),
        Line2D([], [], marker="o", ls="", color=OKABE_ITO["grey"],
               label=f"p_adj \u2265 {alpha:g}"),
        Line2D([], [], marker="o", ls="", color="black", markerfacecolor="white",
               label="estimate only (no test)"),
        Line2D([], [], ls="--", color="black", lw=0.6, label="null value"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=5.5,
               bbox_to_anchor=(0.5, 1.0))
    methods = sorted({m for m in conf.get("fdr_method", pd.Series(dtype=str)) if m})
    regime = "/".join(methods) if methods else "BH"
    fig.suptitle(f"Confirmatory results (one FDR regime: {regime} within test family)",
                 fontsize=7, y=1.0 + 0.12 / n)
    fig.subplots_adjust(left=0.42, right=0.80, hspace=1.25, top=1 - 0.30 / (0.35 + 0.62 * n))
    save_fig(fig, prefix, formats)
    return states


def main(argv=None):
    ap = argparse.ArgumentParser(description="Summary forest of confirmatory results")
    ap.add_argument("--registry", required=True, help="synthesis/results_registry.tsv")
    ap.add_argument("--outdir", required=True, help="output directory")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="p_adj threshold for the 'significant' colour (default 0.05)")
    ap.add_argument("--formats", default="png,svg", help="comma list (default png,svg)")
    args = ap.parse_args(argv)
    if not 0 < args.alpha < 1:
        raise SystemExit(f"--alpha must be in (0, 1), got {args.alpha}")
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    os.makedirs(args.outdir, exist_ok=True)
    states = plot_summary_forest(load_registry(args.registry), args.outdir,
                                 args.alpha, formats)
    print(f"[plot_summary_forest] {len(states)} confirmatory panel(s): "
          f"{', '.join(states) or 'placeholder'}", file=sys.stderr)


if __name__ == "__main__":
    main()
