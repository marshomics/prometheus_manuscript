"""Cross-domain R-loop / provirus mechanism figure (publication-ready).

The R-loop-targeting activity of Prometheus was demonstrated on Pseudomonas phages (a
bacterium). This figure asks whether the proviral-R-loop-vs-Prometheus signal the
archaeal side relies on also holds in bacteria, where the mechanism is established:

  viral_rloop_crossdomain.{formats}
    - fraction of R-loop-prone proviruses by Prometheus presence, in Archaea and in
      Bacteria (distributions);
    - Prometheus prevalence vs R-loop pressure (quantile bins) in both domains, with
      Wilson 95% CIs -- the dose-response, overlaid.
  --out-table (optional): the plotted numbers.

Genomes are the pass_primary set the phyloglm test uses. Presence follows the
Contract-v2 rule for the tier set T (--primary-tiers): present = best tier in T, absent
= no validated locus, other tiers excluded (counted in the panel titles). The R-loop-
prone FRACTION is undefined for a genome with no provirus, so those genomes are left out
of the fraction panels whether aggregate_viral encodes it as '' (current) or 0
(older versions; 0/0 is not 0). The inferential
claim is the phylogeny-controlled, completeness/size-adjusted phyloglm effect (from
ecology_stats), quoted in the title with its BH-adjusted p (p_adj) when available.
Empty inputs -> labelled placeholder.
"""
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO  # noqa: E402
import utils  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

TABLE_HEADER = ["domain_of_life", "bin", "frac_rloop_low", "frac_rloop_high", "n", "k",
                "prevalence", "ci_low", "ci_high"]


def _load(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return None if df.empty else df


def _bool(s):
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def wilson(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def presence(master, tiers):
    """Contract v2 item 6: 1 present (tier in T), 0 absent (tier ""), NaN excluded."""
    if "prometheus_tier" in master.columns:
        tier = master["prometheus_tier"].astype(str).str.strip()
        out = pd.Series(np.nan, index=master.index)
        out[tier == ""] = 0.0
        out[tier.isin(tiers)] = 1.0
        return out
    return _bool(master["prometheus_present"]).astype(float)


def merge(master, viral, tiers):
    """pass_primary genomes with presence (NaN = excluded) + viral columns; or None."""
    if master is None or viral is None or "genome_id" not in viral.columns:
        return None
    m = master.copy()
    if "pass_primary" in m.columns:
        m = m[_bool(m["pass_primary"])]
    m["genome_id"] = m["genome_id"].astype(str)
    m["present"] = presence(m, tiers)
    v = viral.copy()
    v["genome_id"] = v["genome_id"].astype(str)
    d = m.merge(v, on="genome_id", how="inner")
    for c in ("n_proviruses", "frac_rloop_proviruses"):
        d[c] = pd.to_numeric(d[c], errors="coerce") if c in d.columns else np.nan
    d.loc[~(d["n_proviruses"] > 0), "frac_rloop_proviruses"] = np.nan   # 0/0 is undefined
    return d if len(d) else None


def prevalence_curve(d, n_bins=4):
    """Prevalence by quantile bin of frac_rloop_proviruses (excluded genomes dropped)."""
    dd = d.dropna(subset=["frac_rloop_proviruses", "present"])
    if len(dd) < 8 or dd["frac_rloop_proviruses"].nunique() < 2:
        return None
    bins = pd.qcut(dd["frac_rloop_proviruses"], q=min(n_bins, dd["frac_rloop_proviruses"].nunique()),
                   duplicates="drop")
    out = []
    for i, (iv, g) in enumerate(dd.groupby(bins, observed=True)):
        p, lo, hi = wilson(int(g["present"].sum()), len(g))
        out.append((i, iv.left, iv.right, len(g), int(g["present"].sum()), p, lo, hi))
    return out


def bact_beta(assoc):
    if assoc is None or "predictor" not in assoc.columns:
        return ""
    row = assoc[assoc["predictor"].astype(str) == "frac_rloop_proviruses"]
    if not len(row):
        return ""
    est = pd.to_numeric(row["estimate"], errors="coerce").iloc[0]
    if not np.isfinite(est):
        return ""
    lab, p = "", float("nan")
    for col, name in (("p_adj", "p_adj"), ("p_value", "p")):
        if col in row.columns:
            p = pd.to_numeric(row[col], errors="coerce").iloc[0]
            if np.isfinite(p):
                lab = name
                break
    return f"bacteria phyloglm β={est:.2f}" + (f", {lab}={p:.3g}" if lab else "")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arch-master", required=True)
    ap.add_argument("--arch-viral", required=True)
    ap.add_argument("--bact-master", required=True)
    ap.add_argument("--bact-viral", required=True)
    ap.add_argument("--bact-assoc", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--primary-tiers", default="A")
    ap.add_argument("--formats", default="png,svg")
    ap.add_argument("--out-table", default="")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]
    tiers = [t.strip() for t in args.primary_tiers.split(",") if t.strip()]
    prefix = os.path.join(args.outdir, "viral_rloop_crossdomain")

    arch = merge(_load(args.arch_master), _load(args.arch_viral), tiers)
    bact = merge(_load(args.bact_master), _load(args.bact_viral), tiers)
    curves = {}
    for name, d in (("Archaea", arch), ("Bacteria", bact)):
        if d is not None:
            curves[name] = prevalence_curve(d)
    if args.out_table:
        rows = [[name, i, round(lo, 6), round(hi, 6), n, k, round(p, 6), round(l, 6), round(h, 6)]
                for name, cur in curves.items() if cur for (i, lo, hi, n, k, p, l, h) in cur]
        pd.DataFrame(rows, columns=TABLE_HEADER).to_csv(args.out_table, sep="\t", index=False)
    if arch is None and bact is None:
        placeholder(prefix, "no data: no archaeal or bacterial viral tables", formats)
        return

    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9))
    for ax, (name, d, col) in zip(axes[:2], (("Archaea", arch, OKABE_ITO["blue"]),
                                             ("Bacteria", bact, OKABE_ITO["vermillion"]))):
        drew = False
        n_exc = int(d["present"].isna().sum()) if d is not None else 0
        if d is not None:
            groups = [d.loc[d["present"] == g, "frac_rloop_proviruses"].dropna().values
                      for g in (0.0, 1.0)]
            if all(len(x) >= 2 for x in groups):
                parts = ax.violinplot(groups, showmedians=True, widths=0.8)
                for b, c in zip(parts["bodies"], (OKABE_ITO["grey"], col)):
                    b.set_facecolor(c)
                    b.set_alpha(0.6)
                ax.set_xticks([1, 2])
                ax.set_xticklabels([f"absent\n(n={len(groups[0])})",
                                    f"present\n(n={len(groups[1])})"], fontsize=6)
                drew = True
        if not drew:
            ax.text(0.5, 0.5, f"{name}: too few genomes with proviruses", ha="center",
                    va="center", fontsize=6, transform=ax.transAxes)
            ax.set_xticks([])
        ax.set_ylim(0, 1)
        ax.set_ylabel("fraction R-loop-prone proviruses", fontsize=6)
        ax.set_title(name + (f" ({n_exc} excl., tier outside {'+'.join(tiers)})" if n_exc else ""),
                     fontsize=6)
    ax = axes[2]
    drew = False
    for name, col in (("Archaea", OKABE_ITO["blue"]), ("Bacteria", OKABE_ITO["vermillion"])):
        cur = curves.get(name)
        if not cur:
            continue
        xs = [c[0] for c in cur]
        prev = [c[5] for c in cur]
        ax.errorbar(xs, prev, yerr=[[c[5] - c[6] for c in cur], [c[7] - c[5] for c in cur]],
                    fmt="o-", color=col, markersize=4, capsize=2, linewidth=1, label=name)
        drew = True
    if drew:
        ax.set_ylim(0, 1)
        ax.legend(frameon=False, fontsize=6)
    else:
        ax.text(0.5, 0.5, "too few genomes for a dose-response", ha="center", va="center",
                fontsize=6, transform=ax.transAxes)
    ax.set_xlabel("R-loop pressure (quantile bin of R-loop-prone fraction)", fontsize=6)
    ax.set_ylabel(f"Prometheus prevalence (tier {'+'.join(tiers)})", fontsize=6)
    ax.set_title("prevalence vs R-loop pressure", fontsize=7)
    beta = bact_beta(_load(args.bact_assoc))
    fig.suptitle("R-loop targeting across domains" + (f"   ({beta})" if beta else ""), fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, prefix, formats)


if __name__ == "__main__":
    main()
