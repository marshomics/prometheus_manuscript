"""Cross-domain R-loop / provirus mechanism figure (publication-ready).

The R-loop-targeting activity of Prometheus was demonstrated on Pseudomonas phages (a
bacterium). This figure asks whether the proviral-R-loop-vs-Prometheus signal the
archaeal side relies on also holds in bacteria, where the mechanism is established:

  viral_rloop_crossdomain.{png,svg}
    - fraction of R-loop-prone proviruses by Prometheus presence, in Archaea and in
      Bacteria (distributions);
    - Prometheus prevalence vs R-loop pressure (quartiles) in both domains, with Wilson
      95% CIs -- the dose-response, overlaid.

The inferential claim is the phylogeny-controlled, completeness/size-adjusted phyloglm
effect (from ecology_stats), annotated in the title; the splits are shown only as
distributions. Empty inputs -> labelled placeholder.
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


def merge(master, viral):
    if master is None or viral is None or "prometheus_present" not in master.columns:
        return None
    m = master.copy(); m["genome_id"] = m["genome_id"].astype(str)
    v = viral.copy(); v["genome_id"] = v["genome_id"].astype(str)
    d = m.merge(v, on="genome_id", how="left")
    d["present"] = _bool(d["prometheus_present"]).astype(int)
    for c in ("n_proviruses", "frac_rloop_proviruses"):
        d[c] = pd.to_numeric(d.get(c, np.nan), errors="coerce")
    return d


def prevalence_curve(d):
    """(xs, prev, lo, hi) of Prometheus prevalence by R-loop-pressure quartile."""
    dd = d.dropna(subset=["frac_rloop_proviruses"])
    if len(dd) < 8:
        return None
    try:
        dd = dd.assign(_b=pd.qcut(dd["frac_rloop_proviruses"],
                                  q=min(4, dd["frac_rloop_proviruses"].nunique()), duplicates="drop"))
    except Exception:
        dd = dd.assign(_b=pd.cut(dd["frac_rloop_proviruses"], bins=4))
    xs, prev, lo, hi = [], [], [], []
    for i, (_, g) in enumerate(dd.groupby("_b", observed=True)):
        p, l, h = wilson(int(g["present"].sum()), len(g))
        xs.append(i); prev.append(p); lo.append(p - l); hi.append(h - p)
    return xs, prev, lo, hi


def _bact_beta(assoc):
    if assoc is None or "predictor" not in assoc.columns:
        return ""
    row = assoc[assoc["predictor"].astype(str) == "frac_rloop_proviruses"]
    if not len(row):
        return ""
    est = pd.to_numeric(row["estimate"], errors="coerce").iloc[0]
    p = pd.to_numeric(row.get("p_adj", row.get("p_value")), errors="coerce").iloc[0]
    if not np.isfinite(est):
        return ""
    return f"bacteria phyloglm β={est:.2f}" + (f", p={p:.3g}" if np.isfinite(p) else "")


def main():
    ap = argparse.ArgumentParser(description="Cross-domain R-loop mechanism figure")
    ap.add_argument("--arch-master", required=True)
    ap.add_argument("--arch-viral", required=True)
    ap.add_argument("--bact-master", required=True)
    ap.add_argument("--bact-viral", required=True)
    ap.add_argument("--bact-assoc", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, "viral_rloop_crossdomain")

    arch = merge(_load(args.arch_master), _load(args.arch_viral))
    bact = merge(_load(args.bact_master), _load(args.bact_viral))
    if arch is None and bact is None:
        placeholder(prefix, "no data: no archaeal or bacterial viral tables")
        return

    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9))
    # panels 1-2: R-loop-prone fraction by presence, per domain
    for ax, (name, d, col) in zip(axes[:2],
                                  (("Archaea", arch, OKABE_ITO["blue"]),
                                   ("Bacteria", bact, OKABE_ITO["vermillion"]))):
        if d is not None and d["present"].nunique() == 2 and d["frac_rloop_proviruses"].notna().sum() >= 4:
            groups = [d.loc[d["present"] == g, "frac_rloop_proviruses"].dropna().values for g in (0, 1)]
            if all(len(x) for x in groups):
                parts = ax.violinplot(groups, showmedians=True, widths=0.8)
                for b, c in zip(parts["bodies"], (OKABE_ITO["grey"], col)):
                    b.set_facecolor(c); b.set_alpha(0.6)
            ax.set_xticks([1, 2]); ax.set_xticklabels(["absent", "present"], fontsize=6)
        else:
            ax.text(0.5, 0.5, f"{name}: n/a", ha="center", va="center", fontsize=6, transform=ax.transAxes)
        ax.set_ylim(0, 1); ax.set_ylabel("fraction R-loop-prone proviruses", fontsize=6)
        ax.set_title(name, fontsize=7)
    # panel 3: prevalence vs R-loop pressure, both domains overlaid
    ax = axes[2]
    drew = False
    for name, d, col in (("Archaea", arch, OKABE_ITO["blue"]), ("Bacteria", bact, OKABE_ITO["vermillion"])):
        if d is None:
            continue
        cur = prevalence_curve(d)
        if cur is None:
            continue
        xs, prev, lo, hi = cur
        ax.errorbar(xs, prev, yerr=[lo, hi], fmt="o-", color=col, markersize=4,
                    capsize=2, linewidth=1, label=name)
        drew = True
    if drew:
        ax.set_ylim(0, 1); ax.legend(frameon=False, fontsize=6)
    ax.set_xlabel("R-loop pressure quartile", fontsize=6)
    ax.set_ylabel("Prometheus prevalence", fontsize=6)
    ax.set_title("prevalence vs R-loop pressure", fontsize=7)

    beta = _bact_beta(_load(args.bact_assoc))
    fig.suptitle("R-loop targeting across domains" + (f"   ({beta})" if beta else ""), fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_fig(fig, prefix)


if __name__ == "__main__":
    main()
