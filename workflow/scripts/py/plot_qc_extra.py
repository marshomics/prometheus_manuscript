"""Extra QC figure: the covariates the headline association conditions on.

Output (written to --outdir):
  prometheus_covariates.{png,svg}   genome size (Mb), defense load (EXCLUDING Prometheus
                                    when n_defense_systems_excl_prometheus is present), and
                                    CheckM2 completeness, split by Prometheus presence.

This is a diagnostic, not an inference: it shows visually why the headline phyloglm
adjusts for genome size, total defense load (minus the focal system), and completeness
-- if Prometheus-present genomes are systematically larger or defense-richer, an
unadjusted association would be confounded. Distributions only; no test is drawn here.
A violin is drawn only for a group with >= 2 distinct values (a kernel density of a
constant or single value is undefined); otherwise the points and their median are shown.
The CLI (--master, --outdir) is fixed by rules/plots_extra.smk.
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
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def _draw_group(ax, pos, vals, colour):
    if len(vals) == 0:
        return
    if len(np.unique(vals)) >= 2:
        parts = ax.violinplot([vals], positions=[pos], showmedians=True, widths=0.8)
        for b in parts["bodies"]:
            b.set_facecolor(colour)
            b.set_alpha(0.6)
    else:
        jitter = np.linspace(-0.12, 0.12, len(vals)) if len(vals) > 1 else np.zeros(1)
        ax.scatter(pos + jitter, vals, s=6, color=colour, zorder=3, edgecolors="none")
        ax.hlines(float(np.median(vals)), pos - 0.3, pos + 0.3, color="black", lw=0.8)


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
    if "n_defense_systems_excl_prometheus" in m.columns:
        defense = ("n_defense_systems_excl_prometheus", "defense systems (excl. Prometheus)", 1.0)
    else:
        defense = ("n_defense_systems", "defense systems (incl. Prometheus)", 1.0)
    panels = [("genome_size", "genome size (Mb)", 1e-6), defense,
              ("completeness", "CheckM2 completeness (%)", 1.0)]
    panels = [p for p in panels if p[0] in m.columns]
    if not panels or m["present"].nunique() < 2:
        placeholder(prefix, "no data: need presence variation and covariate columns"); return

    fig, axes = plt.subplots(1, len(panels), figsize=(2.3 * len(panels) + 0.5, 2.7))
    axes = np.atleast_1d(axes)
    for ax, (col, lab, scale) in zip(axes, panels):
        vals = pd.to_numeric(m[col], errors="coerce") * scale
        groups = [vals[(m["present"] == g) & vals.notna()].to_numpy() for g in (0, 1)]
        for pos, g, colour in zip((1, 2), groups, (OKABE_ITO["grey"], OKABE_ITO["blue"])):
            _draw_group(ax, pos, g, colour)
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"absent\n(n={len(groups[0])})", f"present\n(n={len(groups[1])})"], fontsize=6)
        ax.set_xlim(0.4, 2.6)
        ax.set_ylabel(lab, fontsize=6)
    fig.suptitle("Covariates conditioned on in the headline association", fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    save_fig(fig, prefix)


if __name__ == "__main__":
    main()
