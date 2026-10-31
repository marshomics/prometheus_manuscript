"""Diagnostic figure for Phase-1 preprocessing.

Panel A: CheckM2 completeness vs contamination of every input genome that has both
estimates, coloured by the FULL Phase-1 quality filter (all five thresholds: the
`passing` column of filter_report.tsv), with the completeness / contamination threshold
lines. Genomes above the contamination axis limit are counted, not hidden.
Panel B: the filtering funnel input -> GTDB-Tk classified -> in the kept domain -> with a
GTDB species -> passing QC -> species representatives (counts from filter_report.tsv).
Without --report the scatter falls back to the two completeness/contamination thresholds
and says so in its legend.

Output (to --outdir): preprocess_quality.{formats}
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

apply_style()


def _read(path):
    if not path or not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _truthy(s):
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--report", default="", help="filter_report.tsv (filter_collapse.py --report)")
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--max-contamination", type=float, default=5.0)
    ap.add_argument("--keep-domain", default="Archaea")
    ap.add_argument("--formats", nargs="+", default=["png", "svg"])
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, "preprocess_quality")
    fmts = tuple(args.formats)

    ann = _read(args.annotated)
    if ann.empty or "genome_id" not in ann.columns:
        placeholder(prefix, "no data: annotated_metadata.tsv missing or empty", formats=fmts)
        return
    rep = _read(args.report)
    sp = _read(args.species)
    n_species = len(sp) if "genome_id" in sp.columns else 0

    comp = pd.to_numeric(ann.get("completeness"), errors="coerce")
    cont = pd.to_numeric(ann.get("contamination"), errors="coerce")
    if not rep.empty and "passing" in rep.columns:
        pmap = dict(zip(rep["genome_id"], _truthy(rep["passing"])))
        passing = ann["genome_id"].map(pmap).fillna(False).astype(bool)
        pass_label = "pass (all 5 QC thresholds)"
    else:
        passing = ((comp >= args.min_completeness) & (cont <= args.max_contamination)).fillna(False)
        pass_label = "pass (completeness/contamination only)"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 1.15]})
    m = comp.notna() & cont.notna()
    finite = cont[m]
    ytop = max(2.0 * args.max_contamination, min(100.0, 1.05 * float(finite.max()) if len(finite) else 0.0), 1.0)
    above = int((finite > ytop).sum())
    for mask, colour, lab, z in ((m & ~passing, OKABE_ITO["grey"], "fail", 2),
                                 (m & passing, OKABE_ITO["blue"], pass_label, 3)):
        ax1.scatter(comp[mask], cont[mask].clip(upper=ytop), s=3, alpha=0.5, color=colour,
                    label=f"{lab} (n={int(mask.sum()):,})", edgecolors="none", rasterized=True, zorder=z)
    ax1.axvline(args.min_completeness, ls="--", lw=0.6, color="k", zorder=1)
    ax1.axhline(args.max_contamination, ls="--", lw=0.6, color="k", zorder=1)
    ax1.set_xlim(0, 100.5)
    ax1.set_ylim(0, ytop)
    ax1.set_xlabel("CheckM2 completeness (%)")
    ax1.set_ylabel("CheckM2 contamination (%)")
    title = "Input genome quality"
    if above:
        title += f" ({above:,} above {ytop:g}% drawn at the top)"
    ax1.set_title(title)
    n_no_est = int((~m).sum())
    if n_no_est:
        ax1.text(0.02, 0.97, f"{n_no_est:,} genome(s) without CheckM2 estimates",
                 transform=ax1.transAxes, ha="left", va="top", fontsize=6)
    ax1.legend(frameon=False, fontsize=6, markerscale=2, loc="center left")

    if not rep.empty and "status" in rep.columns:
        st = rep["status"]
        gstat = rep["gtdbtk_status"] if "gtdbtk_status" in rep.columns else pd.Series("", index=rep.index)
        stages = [
            ("input", len(rep)),
            ("GTDB-Tk classified", int((gstat == "classified").sum())),
            (f"{args.keep_domain or 'any domain'}", int((st != "excluded_domain").sum())),
            ("GTDB species", int((~st.isin(["excluded_domain", "no_species"])).sum())),
            ("pass QC", int(st.isin(["representative", "species_member"]).sum())),
            ("species reps", int((st == "representative").sum())),
        ]
    else:
        stages = [("input", len(ann)), (pass_label.split(" (")[0], int(passing.sum())),
                  ("species reps", n_species)]
    labels = [s for s, _ in stages]
    vals = [v for _, v in stages]
    colours = [OKABE_ITO["grey"]] * (len(vals) - 1) + [OKABE_ITO["green"]]
    x = np.arange(len(vals))
    ax2.bar(x, vals, color=colours, width=0.65, edgecolor="none")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=30, ha="right")
    ax2.set_ylabel("Genomes")
    ax2.set_title("Filtering funnel")
    ax2.set_ylim(0, max(vals + [1]) * 1.15)
    for i, v in enumerate(vals):
        ax2.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=6)
    save_fig(fig, prefix, formats=fmts)


if __name__ == "__main__":
    main()
