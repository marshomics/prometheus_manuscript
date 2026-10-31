"""Extra cross-domain figure: bacteria vs archaea, prevalence + architecture.

Output (written to --outdir):
  crossdomain_comparison.{formats}   Prometheus prevalence (Wilson 95% CI) in bacteria
                                     and archaea, and the fraction of validated loci
                                     carrying each configured domain category in each
                                     domain of life (is archaeal Prometheus
                                     architecturally reduced relative to the
                                     characterised bacterial form?).
  --out-table (optional)             the plotted numbers.

Prevalence follows the Contract-v2 presence rule for the tier set T (--primary-tiers,
config stats.confidence_tiers_primary): among pass_primary genomes a genome is PRESENT
when its best tier is in T, ABSENT when it has no validated locus (tier ""), and
EXCLUDED (not in the denominator; counted on the axis label) when its best tier is
outside T -- a B/C-only genome is not a confident absence in a tier-A prevalence.
Domains are shown side by side, never pooled, because sampling differs. Domain
categories come from the --domains spec (canonical N->C order) or, without it, from the
has_* columns of the typology tables. Placeholder only when every input is empty.
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

TABLE_HEADER = ["panel", "domain_of_life", "category", "tier_set", "n", "k", "n_excluded",
                "value", "ci_low", "ci_high"]


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
    utils.log("WARNING: master has no prometheus_tier; using prometheus_present (any tier)")
    return _bool(master["prometheus_present"]).astype(float)


def prevalence(master, tiers):
    """(p, lo, hi, n_in_denominator, k_present, n_excluded) or None."""
    if master is None or ("prometheus_tier" not in master.columns
                          and "prometheus_present" not in master.columns):
        return None
    m = master[_bool(master["pass_primary"])] if "pass_primary" in master.columns else master
    if m.empty:
        return None
    pr = presence(m, tiers)
    n_exc = int(pr.isna().sum())
    pr = pr.dropna()
    n, k = int(len(pr)), int(pr.sum())
    if n == 0:
        return None
    return wilson(k, n) + (n, k, n_exc)


def categories(spec, typologies):
    if spec:
        return utils.canonical_order(spec)
    cats = []
    for t in typologies:
        if t is not None:
            cats += [c[4:] for c in t.columns if c.startswith("has_") and c[4:] not in cats]
    return cats


def domain_fraction(typ, cats):
    """{category: (fraction of loci, n_loci)} weighted by n_loci, or None."""
    if typ is None or "architecture" not in typ.columns:
        return None
    n = pd.to_numeric(typ["n_loci"], errors="coerce").fillna(0) if "n_loci" in typ.columns \
        else pd.Series(1.0, index=typ.index)
    total = float(n.sum())
    if total == 0:
        return None
    toks = typ["architecture"].astype(str).apply(lambda a: set(a.split("-")) if a else set())
    frac = {}
    for cat in cats:
        col = f"has_{cat}"
        has = _bool(typ[col]) if col in typ.columns else toks.apply(lambda s, c=cat: c in s)
        frac[cat] = (float(n[has].sum()) / total, int(total))
    return frac


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arch-master", required=True)
    ap.add_argument("--bact-master", required=True)
    ap.add_argument("--arch-typology", required=True)
    ap.add_argument("--bact-typology", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--domains", default="", help="prometheus.domains JSON spec")
    ap.add_argument("--primary-tiers", default="A", help="tier set T, e.g. 'A' or 'A,B'")
    ap.add_argument("--formats", default="png,svg")
    ap.add_argument("--out-table", default="")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    apply_style()
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]
    tiers = [t.strip() for t in args.primary_tiers.split(",") if t.strip()]
    tier_lbl = "+".join(tiers)
    prefix = os.path.join(args.outdir, "crossdomain_comparison")
    spec = utils.load_domain_spec(args.domains) if args.domains else {}

    arch_t, bact_t = _load(args.arch_typology), _load(args.bact_typology)
    cats = categories(spec, (arch_t, bact_t))
    ap_ = prevalence(_load(args.arch_master), tiers)
    bp_ = prevalence(_load(args.bact_master), tiers)
    af, bf = domain_fraction(arch_t, cats), domain_fraction(bact_t, cats)
    rows = []
    for name, pv in (("Archaea", ap_), ("Bacteria", bp_)):
        if pv is not None:
            p, lo, hi, n, k, n_exc = pv
            rows.append(["prevalence", name, "Prometheus", tier_lbl, n, k, n_exc,
                         round(p, 6), round(lo, 6), round(hi, 6)])
    for name, fr in (("Archaea", af), ("Bacteria", bf)):
        for cat, (f, n) in (fr or {}).items():
            rows.append(["architecture", name, cat, "A+B+C", n, round(f * n), 0, round(f, 6), "", ""])
    if args.out_table:
        pd.DataFrame(rows, columns=TABLE_HEADER).to_csv(args.out_table, sep="\t", index=False)
    if ap_ is None and bp_ is None and af is None and bf is None:
        placeholder(prefix, "no data: no bacterial/archaeal tables available", formats)
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.4, 2.9),
                                   gridspec_kw={"width_ratios": [1, 1.4]})
    labels, prevs, los, his, cols = [], [], [], [], []
    for name, pv, col in (("Archaea", ap_, OKABE_ITO["blue"]),
                          ("Bacteria", bp_, OKABE_ITO["vermillion"])):
        if pv is not None:
            p, lo, hi, n, k, n_exc = pv
            labels.append(f"{name}\n{k}/{n}" + (f"\n({n_exc} excl.)" if n_exc else ""))
            prevs.append(p)
            los.append(p - lo)
            his.append(hi - p)
            cols.append(col)
    if prevs:
        x = np.arange(len(prevs))
        ax1.bar(x, prevs, color=cols, alpha=0.85, width=0.6)
        ax1.errorbar(x, prevs, yerr=[los, his], fmt="none", color="black", capsize=3,
                     linewidth=0.8)
        ax1.set_xticks(x)
        ax1.set_xticklabels(labels, fontsize=6)
    ax1.set_ylim(0, 1)
    ax1.set_ylabel(f"Prometheus prevalence, tier {tier_lbl} (Wilson 95% CI)")
    ax1.set_title("prevalence by domain of life", fontsize=7)
    x = np.arange(len(cats))
    w = 0.38
    if af:
        ax2.bar(x - w / 2, [af[c][0] for c in cats], width=w, color=OKABE_ITO["blue"],
                alpha=0.85, label=f"Archaea (n={next(iter(af.values()))[1]} loci)")
    if bf:
        ax2.bar(x + w / 2, [bf[c][0] for c in cats], width=w, color=OKABE_ITO["vermillion"],
                alpha=0.85, label=f"Bacteria (n={next(iter(bf.values()))[1]} loci)")
    ax2.set_xticks(x)
    ax2.set_xticklabels(cats, rotation=40, ha="right", fontsize=6)
    ax2.set_ylim(0, 1)
    ax2.set_ylabel("fraction of validated loci with domain")
    ax2.set_title("domain architecture by domain of life", fontsize=7)
    if af or bf:
        ax2.legend(frameon=False, fontsize=6)
    fig.suptitle("Cross-domain comparison: prevalence & architecture", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, prefix, formats)


if __name__ == "__main__":
    main()
