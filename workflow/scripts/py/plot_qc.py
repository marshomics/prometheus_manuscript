"""plot_qc.py — QC figures for the analysis genomes (qc/genome_metadata.tsv).

Outputs (written to --outdir, one file per --formats entry):
  qc_completeness_contamination   CheckM2 completeness vs contamination coloured by
                                  quality_tier, with the configured tier/floor lines
                                  (qc.tier_high_*, qc.completeness_floor/contamination_ceiling)
  qc_assembly_stats               histograms: genome size (Mb, from bp), contig N50
                                  (log10 bp), number of contigs (log10), GC (%; BBMap
                                  gc_avg is a fraction and is converted)

MAG/isolate genome_type is deliberately not used (it is never produced by the pipeline).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO, TIER_COLORS

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

apply_style()
TIER_ORDER = ["A", "B", "C"]


def _load_meta(path):
    if not path or not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t", dtype={"genome_id": str})
    except pd.errors.EmptyDataError:
        return None
    return None if df.empty else df


def plot_completeness_contamination(df, outdir, args):
    prefix = os.path.join(outdir, "qc_completeness_contamination")
    fmts = tuple(args.formats)
    if df is None:
        placeholder(prefix, "no data: genome_metadata.tsv missing or empty", formats=fmts)
        return
    if not {"completeness", "contamination"}.issubset(df.columns):
        placeholder(prefix, "no data: completeness/contamination columns absent", formats=fmts)
        return
    d = df.assign(completeness=pd.to_numeric(df["completeness"], errors="coerce"),
                  contamination=pd.to_numeric(df["contamination"], errors="coerce"))
    d = d.dropna(subset=["completeness", "contamination"])
    if d.empty:
        placeholder(prefix, "no data: all rows have NaN completeness/contamination", formats=fmts)
        return

    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    tiers = (d["quality_tier"].fillna("NA").astype(str) if "quality_tier" in d.columns
             else pd.Series(["NA"] * len(d), index=d.index))
    seen = tiers.unique().tolist()
    plotted = [t for t in TIER_ORDER if t in seen] + sorted(t for t in seen if t not in TIER_ORDER)
    handles = []
    for tier in plotted:
        mask = tiers == tier
        colour = TIER_COLORS.get(tier, OKABE_ITO["grey"])
        ax.scatter(d.loc[mask, "completeness"], d.loc[mask, "contamination"], c=colour, s=8,
                   alpha=0.75, linewidths=0, zorder=3)
        handles.append(mlines.Line2D([], [], color=colour, marker="o", ls="", ms=3,
                                     label=f"tier {tier} (n={int(mask.sum())})"))
    ax.axvline(args.tier_high_completeness, color="black", lw=0.5, ls="--", zorder=2)
    ax.axhline(args.tier_high_contamination, color="black", lw=0.5, ls="--", zorder=2)
    ax.axvline(args.completeness_floor, color="black", lw=0.5, ls=":", zorder=2)
    ax.axhline(args.contamination_ceiling, color="black", lw=0.5, ls=":", zorder=2)
    handles += [mlines.Line2D([], [], color="black", lw=0.5, ls="--",
                              label=f"tier A: >= {args.tier_high_completeness:g}% / <= {args.tier_high_contamination:g}%"),
                mlines.Line2D([], [], color="black", lw=0.5, ls=":",
                              label=f"floor: >= {args.completeness_floor:g}% / <= {args.contamination_ceiling:g}%")]
    xmin = min(float(d["completeness"].min()), args.completeness_floor) - 2
    ax.set_xlim(max(0.0, xmin), 100.5)
    ax.set_ylim(0, max(float(d["contamination"].max()), args.contamination_ceiling) * 1.1)
    ax.set_xlabel("CheckM2 completeness (%)")
    ax.set_ylabel("CheckM2 contamination (%)")
    ax.set_title(f"Analysis genome quality (n={len(d)})")
    ax.legend(handles=handles, frameon=False, fontsize=5.5, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              borderaxespad=0.0)
    save_fig(fig, prefix, formats=fmts)


def plot_assembly_stats(df, outdir, args):
    prefix = os.path.join(outdir, "qc_assembly_stats")
    fmts = tuple(args.formats)
    if df is None:
        placeholder(prefix, "no data: genome_metadata.tsv missing or empty", formats=fmts)
        return
    if not {"genome_size", "n50", "n_contigs", "gc"}.intersection(df.columns):
        placeholder(prefix, "no data: no assembly stat columns found", formats=fmts)
        return

    def gc_pct(x):  # BBMap gc_avg is a fraction; tolerate a table that already holds %
        return x * 100.0 if x.max() <= 1.0 else x

    panels = [
        ("genome_size", "Genome size (Mb)", False, lambda x: x / 1e6),
        ("n50", "Contig N50 (log10 bp)", True, lambda x: x),
        ("n_contigs", "Contigs (log10 count)", True, lambda x: x),
        ("gc", "GC content (%)", False, gc_pct),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0))
    fig.subplots_adjust(hspace=0.5, wspace=0.35)
    for ax, (col, xlabel, log_x, transform) in zip(axes.flat, panels):
        vals = pd.to_numeric(df[col], errors="coerce").dropna() if col in df.columns else pd.Series(dtype=float)
        if vals.empty:
            ax.text(0.5, 0.5, f"{col}\nnot available", ha="center", va="center",
                    transform=ax.transAxes, fontsize=6)
            ax.axis("off")
            continue
        v = transform(vals)
        if log_x:
            v = np.log10(v.clip(lower=1))
        ax.hist(v, bins=min(30, max(5, int(np.sqrt(len(v))) + 1)), color=OKABE_ITO["blue"],
                alpha=0.8, edgecolor="none")
        ax.axvline(float(np.median(v)), color="black", lw=0.6, ls="--")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Genomes")
        med = float(np.median(transform(vals)))
        med_txt = f"{med:,.0f}" if abs(med) >= 1000 else f"{med:.3g}"
        ax.set_title(f"n={len(vals)}, median {med_txt}", fontsize=6)
    save_fig(fig, prefix, formats=fmts)


def main():
    parser = argparse.ArgumentParser(description="QC figures")
    parser.add_argument("--meta", required=True, help="qc/genome_metadata.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--completeness-floor", type=float, default=50.0)
    parser.add_argument("--contamination-ceiling", type=float, default=10.0)
    parser.add_argument("--tier-high-completeness", type=float, default=90.0)
    parser.add_argument("--tier-high-contamination", type=float, default=5.0)
    parser.add_argument("--formats", nargs="+", default=["png", "svg"])
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    df = _load_meta(args.meta)
    plot_completeness_contamination(df, args.outdir, args)
    plot_assembly_stats(df, args.outdir, args)


if __name__ == "__main__":
    main()
