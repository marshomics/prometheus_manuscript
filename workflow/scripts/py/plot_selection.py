"""plot_selection.py -- selection by Prometheus domain.

Output (in --outdir): selection_by_domain.{png,svg}
  left : fraction of codon sites under pervasive positive selection (FEL, beta > alpha,
         p <= the run's threshold) and under episodic selection (MEME) per domain, with
         Wilson 95% CIs and the number of sites per domain (effect + uncertainty + n);
  right: per-site FEL beta/alpha by domain (finite values only; sites with alpha = 0 and
         beta > 0 are counted in the panel title instead of being plotted at an arbitrary y).
Domains are ordered by their first alignment column (N->C); 'NA' = linker / unassigned.
A site table with zero positive sites is still a result and is plotted (zeros), not
replaced by a placeholder. Header-only / missing input -> labelled placeholder.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402

apply_style()


def wilson(k, n, z=1.959964):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def _load(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])
    return None if df.empty else df


def domain_order(sites):
    first = sites.groupby("domain")["alignment_col"].min().sort_values()
    order = [d for d in first.index if d != "NA"]
    return order + (["NA"] if "NA" in first.index else [])


def per_domain(sites):
    rows = []
    for d in domain_order(sites):
        sub = sites[sites["domain"] == d]
        n = len(sub)
        k_fel = int((sub["classification"] == "positive").sum())
        k_meme = int((sub.get("meme_classification", pd.Series(dtype=str)) == "episodic").sum())
        rows.append({"domain": d, "n": n, "k_fel": k_fel, "k_meme": k_meme,
                     "has_meme": "meme_classification" in sub.columns and sub["meme_classification"].notna().any()})
    return pd.DataFrame(rows)


def plot_selection_by_domain(summary, sites, outdir):
    prefix = os.path.join(outdir, "selection_by_domain")
    if sites is None:
        placeholder(prefix, "no data: site_table.tsv missing or header-only (selection not run)")
        return
    for c in ("alignment_col", "domain", "classification"):
        if c not in sites.columns:
            placeholder(prefix, f"no data: column '{c}' missing from site_table.tsv")
            return
    sites = sites.copy()
    sites["alignment_col"] = pd.to_numeric(sites["alignment_col"], errors="coerce")
    sites["domain"] = sites["domain"].fillna("NA").astype(str)
    tab = per_domain(sites)
    thr = ""
    if summary is not None and {"method", "note"} <= set(summary.columns):
        fel = summary[summary["method"] == "FEL"]["note"].astype(str)
        if len(fel):
            thr = fel.iloc[0].split(";")[0]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.1, 1]})
    x = np.arange(len(tab))
    w = 0.38
    series = [("k_fel", "FEL pervasive (+)", OKABE_ITO["vermillion"])]
    if tab["has_meme"].any():
        series.append(("k_meme", "MEME episodic", OKABE_ITO["blue"]))
    for j, (col, lab, colr) in enumerate(series):
        frac = tab[col] / tab["n"].where(tab["n"] > 0)
        lo, hi = zip(*[wilson(k, n) for k, n in zip(tab[col], tab["n"])])
        off = (j - (len(series) - 1) / 2) * w
        ax.bar(x + off, frac, width=w, color=colr, label=lab, edgecolor="none")
        ax.errorbar(x + off, frac, yerr=[frac - np.array(lo), np.array(hi) - frac], fmt="none",
                    ecolor="black", elinewidth=0.6, capsize=1.5)
        for xi, k, h in zip(x + off, tab[col], hi):     # count k above the upper 95% CI
            ax.text(xi, h + 0.01, str(k), ha="center", va="bottom", fontsize=5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d}\n(n={n})" for d, n in zip(tab["domain"], tab["n"])], fontsize=6)
    ax.set_ylabel("fraction of codon sites")
    ax.set_title(f"Selected sites per domain ({thr})" if thr else "Selected sites per domain", fontsize=7)
    ax.legend(frameon=False, fontsize=6)

    ratio = pd.to_numeric(sites.get("beta_alpha_or_dnds", pd.Series(dtype=float)).replace("inf", np.inf),
                          errors="coerce")
    finite = sites.assign(r=ratio)[np.isfinite(ratio) & (ratio > 0)]
    n_inf = int(np.isinf(ratio).sum())
    data = [finite.loc[finite["domain"] == d, "r"].values for d in tab["domain"]]
    pos = [i for i, v in enumerate(data) if len(v) > 1]
    if pos:
        parts = ax2.violinplot([data[i] for i in pos], positions=pos, showmedians=True, widths=0.7)
        for b in parts["bodies"]:
            b.set_facecolor(OKABE_ITO["grey"])
            b.set_alpha(0.6)
    rng = np.random.default_rng(42)
    for i, v in enumerate(data):
        if len(v):
            ax2.scatter(i + rng.uniform(-0.12, 0.12, size=len(v)), v, s=3, alpha=0.5,
                        color=OKABE_ITO["black"], zorder=3)
    ax2.set_yscale("log")
    ax2.yaxis.set_major_locator(mticker.LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax2.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax2.axhline(1.0, color="black", linewidth=0.6, linestyle="--")
    ax2.set_xticks(x)
    ax2.set_xticklabels(tab["domain"], fontsize=6)
    ax2.set_ylabel("FEL beta/alpha (log)")
    ax2.set_title(f"Per-site beta/alpha ({n_inf} site(s) with alpha=0, beta>0 not shown)", fontsize=6)
    save_fig(fig, prefix)


def main():
    ap = argparse.ArgumentParser(description="Selection figures")
    ap.add_argument("--summary", required=True, help="evolution/selection/selection_summary.tsv")
    ap.add_argument("--sites", required=True, help="evolution/selection/site_table.tsv")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plot_selection_by_domain(_load(args.summary), _load(args.sites), args.outdir)


if __name__ == "__main__":
    main()
