"""Extra QC figure: the covariates the headline association conditions on.

Output (written to --outdir):
  prometheus_covariates.{png,svg}   genome size, defense load (EXCLUDING Prometheus),
                                    and completeness split by Prometheus presence.

This is a diagnostic, not an inference: it shows visually why the headline phyloglm
adjusts for genome size, total defense load (minus the focal system), and completeness
-- if Prometheus-present genomes are systematically larger or defense-richer, an
unadjusted association would be confounded. Distributions only; no test is drawn here.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

apply_style()


def _bool(s):
    return s.astype(str).str.lower().isin(("true", "1", "yes"))


def main():
    ap = argparse.ArgumentParser(description="QC covariate diagnostic")
    ap.add_argument("--master", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, "prometheus_covariates")

    if not os.path.exists(args.master):
        placeholder(prefix, "no data: master table missing"); return
    try:
        m = pd.read_csv(args.master, sep="\t")
    except Exception:
        placeholder(prefix, "no data: master table unreadable"); return
    if m.empty or "prometheus_present" not in m.columns:
        placeholder(prefix, "no data: master table empty / no prometheus_present"); return

    m["present"] = _bool(m["prometheus_present"]).astype(int)
    defense_col = "n_defense_systems_excl_prometheus" if "n_defense_systems_excl_prometheus" in m.columns \
        else "n_defense_systems"
    panels = [("genome_size", "genome size (bp)"),
              (defense_col, "defense systems (excl. Prometheus)"),
              ("completeness", "completeness (%)")]
    panels = [(c, lab) for c, lab in panels if c in m.columns]
    if not panels or m["present"].nunique() < 2:
        placeholder(prefix, "no data: need presence variation and covariate columns"); return

    fig, axes = plt.subplots(1, len(panels), figsize=(2.3 * len(panels) + 0.5, 2.7))
    if len(panels) == 1:
        axes = [axes]
    for ax, (col, lab) in zip(axes, panels):
        vals = pd.to_numeric(m[col], errors="coerce")
        groups = [vals[(m["present"] == g) & vals.notna()].values for g in (0, 1)]
        if all(len(g) for g in groups):
            parts = ax.violinplot(groups, showmedians=True, widths=0.8)
            for b, c in zip(parts["bodies"], (OKABE_ITO["grey"], OKABE_ITO["blue"])):
                b.set_facecolor(c); b.set_alpha(0.6)
        ax.set_xticks([1, 2]); ax.set_xticklabels(["absent", "present"], fontsize=6)
        ax.set_ylabel(lab, fontsize=6)
    fig.suptitle("Covariates conditioned on in the headline association", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, prefix)


if __name__ == "__main__":
    main()
