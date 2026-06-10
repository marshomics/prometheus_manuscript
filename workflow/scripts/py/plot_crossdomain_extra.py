"""Extra cross-domain figure: bacteria vs archaea, prevalence + architecture.

Output (written to --outdir):
  crossdomain_comparison.{png,svg}   overall Prometheus prevalence (Wilson CI) in
                                     bacteria vs archaea, and the fraction of validated
                                     loci carrying each canonical domain in each domain
                                     of life (is archaeal Prometheus architecturally
                                     degenerate relative to the characterised bacterial
                                     form?).

Prevalence is shown side-by-side (not pooled) because bacterial and archaeal sampling
differ; the architecture panel anchors archaeal loci to the experimentally validated
bacterial reference architecture.
"""
import os
import math
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

apply_style()
CANON = ["DUF4011", "DUF3320", "DNA2", "nuclease", "REase"]


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


def _prev(master):
    if master is None or "prometheus_present" not in master.columns:
        return None
    m = master
    if "pass_primary" in m.columns:
        m = m[_bool(m["pass_primary"])]
    n = len(m)
    if n == 0:
        return None
    return wilson(int(_bool(m["prometheus_present"]).sum()), n) + (n,)


def _domain_fraction(typ):
    """Fraction of loci carrying each canonical domain (weighted by n_loci)."""
    if typ is None or "architecture" not in typ.columns:
        return None
    n = pd.to_numeric(typ.get("n_loci", 1), errors="coerce").fillna(1)
    total = n.sum()
    if total == 0:
        return None
    frac = {}
    for dom in CANON:
        col = f"has_{dom}"
        if col in typ.columns:
            has = typ[col].astype(str).str.lower().isin(("true", "1", "yes"))
        else:
            has = typ["architecture"].astype(str).apply(lambda a: dom in a.split("-"))
        frac[dom] = float((n[has].sum()) / total)
    return frac


def main():
    ap = argparse.ArgumentParser(description="Cross-domain comparison figure")
    ap.add_argument("--arch-master", required=True)
    ap.add_argument("--bact-master", required=True)
    ap.add_argument("--arch-typology", required=True)
    ap.add_argument("--bact-typology", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, "crossdomain_comparison")

    ap_ = _prev(_load(args.arch_master)); bp_ = _prev(_load(args.bact_master))
    af = _domain_fraction(_load(args.arch_typology)); bf = _domain_fraction(_load(args.bact_typology))
    if ap_ is None and bp_ is None and af is None and bf is None:
        placeholder(prefix, "no data: no bacterial/archaeal tables available"); return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.4, 2.9), gridspec_kw={"width_ratios": [1, 1.4]})
    # prevalence
    labels, prevs, los, his, ns, cols = [], [], [], [], [], []
    for name, pv, col in (("Archaea", ap_, OKABE_ITO["blue"]), ("Bacteria", bp_, OKABE_ITO["vermillion"])):
        if pv is not None:
            p, lo, hi, n = pv
            labels.append(f"{name}\n(n={n})"); prevs.append(p); los.append(p - lo); his.append(hi - p); cols.append(col)
    if prevs:
        x = np.arange(len(prevs))
        ax1.bar(x, prevs, color=cols, alpha=0.85, width=0.6)
        ax1.errorbar(x, prevs, yerr=[los, his], fmt="none", color="black", capsize=3, linewidth=0.8)
        ax1.set_xticks(x); ax1.set_xticklabels(labels, fontsize=6)
    ax1.set_ylim(0, 1); ax1.set_ylabel("Prometheus prevalence (Wilson 95% CI)")
    ax1.set_title("prevalence by domain", fontsize=7)
    # architecture fractions
    x = np.arange(len(CANON)); w = 0.38
    if af:
        ax2.bar(x - w / 2, [af.get(d, 0) for d in CANON], width=w, label="Archaea", color=OKABE_ITO["blue"], alpha=0.85)
    if bf:
        ax2.bar(x + w / 2, [bf.get(d, 0) for d in CANON], width=w, label="Bacteria", color=OKABE_ITO["vermillion"], alpha=0.85)
    ax2.set_xticks(x); ax2.set_xticklabels(CANON, rotation=40, ha="right", fontsize=6)
    ax2.set_ylim(0, 1); ax2.set_ylabel("fraction of loci with domain")
    ax2.set_title("domain architecture by domain of life", fontsize=7)
    ax2.legend(frameon=False, fontsize=6)
    fig.suptitle("Cross-domain comparison: prevalence & architecture", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, prefix)


if __name__ == "__main__":
    main()
