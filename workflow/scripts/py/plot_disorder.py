"""Disorder figures for the ProA/Prometheus disorder module.

Outputs (written to --outdir, one file per --formats entry):
  disorder_along_protein    mean consensus disorder per alignment column (columns with
                            occupancy >= --min-occupancy) over a domain track drawn as
                            contiguous runs; disorder should pile up in linkers/tails
                            and drop inside folded domains (a built-in sanity check).
  disorder_by_group         per-genome disorder (mean over the genome's tier-T loci) by
                            methanogen status -- non-methanogen / ANME / methanogen,
                            ANME kept separate -- and, with a bacterial table, by domain
                            of life. Genomes outside the tier set are EXCLUDED and
                            counted (Contract v2 item 6).
  disorder_linker_vs_domain disorder in domain vs linker/tail columns (same occupancy
                            filter).
Domain categories and their order come from the prometheus.domains spec (--domains),
else from the domain values present; there is no hardcoded category list. The
inferential claims are the phylolm effects in disorder_associations.tsv, not these
distributions.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import utils  # noqa: E402

apply_style()
NA_COLOR = "#DDDDDD"
CATEGORY_PALETTE = [OKABE_ITO[k] for k in ("orange", "green", "vermillion", "skyblue",
                                           "purple", "blue", "yellow", "black")]
TRUE_SET = ("true", "1", "yes", "t")


def _load(path):
    """TSV -> DataFrame; None for no path / header-only. A missing named file is an error."""
    if not path:
        return None
    if not os.path.isfile(path):
        raise SystemExit(f"[plot_disorder] input {path} not found")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return None if df.empty else df


def _bool(s):
    return s.astype(str).str.strip().str.lower().isin(TRUE_SET)


def _resp(df):
    return "frac_disordered_nonlcr" if "frac_disordered_nonlcr" in df.columns else "frac_disordered"


def domain_categories(spec, bycol):
    """Ordered category list: the spec's canonical order, else order of appearance."""
    if spec:
        return utils.canonical_order(spec)
    if bycol is None or "domain" not in bycol.columns:
        return []
    seen = []
    for d in bycol["domain"].astype(str):
        if d not in ("NA", "", "nan") and d not in seen:
            seen.append(d)
    return seen


def domain_colors(categories):
    return {c: CATEGORY_PALETTE[i % len(CATEGORY_PALETTE)] for i, c in enumerate(categories)}


def domain_runs(cols, doms):
    """Contiguous (start, end, domain) runs over sorted alignment columns."""
    runs = []
    for c, d in zip(cols, doms):
        if runs and runs[-1][2] == d and c == runs[-1][1] + 1:
            runs[-1][1] = c
        else:
            runs.append([c, c, d])
    return [tuple(r) for r in runs]


def _occupancy_filter(df, min_occ):
    if "occupancy" in df.columns and min_occ > 0:
        occ = pd.to_numeric(df["occupancy"], errors="coerce")
        return df[occ >= min_occ]
    return df


def tier_filter_loci(per_locus, master, tiers):
    """Contract v2 item 6 for a per-locus disorder table and tier set T.
    -> (loci of genomes whose best tier is in T, restricted to tier-T loci;
        number of genomes EXCLUDED because their best tier is outside T)."""
    if not tiers:
        return per_locus, 0
    pl = per_locus
    if "confidence_tier" in pl.columns:
        pl = pl[pl["confidence_tier"].astype(str).str.strip().isin(tiers)]
    n_excl = 0
    if master is not None and "prometheus_tier" in master.columns:
        mt = master.set_index(master["genome_id"].astype(str))["prometheus_tier"].astype(str).str.strip()
        best = pl["genome_id"].astype(str).map(mt)
        pl = pl[best.isin(tiers)]
        n_excl = int(((mt != "") & ~mt.isin(tiers)).sum())
    return pl, n_excl


def _strip(ax, groups, labels, colors):
    """Violin where a group has >= 3 distinct values, plus every genome as a point."""
    rng = np.random.default_rng(12345)
    for i, (g, c) in enumerate(zip(groups, colors), 1):
        if len(g) >= 3 and np.unique(g).size >= 2:
            parts = ax.violinplot([g], positions=[i], showmedians=True, widths=0.8)
            for b in parts["bodies"]:
                b.set_facecolor(c); b.set_alpha(0.45)
        if len(g):
            ax.scatter(i + rng.uniform(-0.12, 0.12, len(g)), g, s=4, color=c, zorder=3,
                       linewidths=0)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels([f"{lab}\n(n={len(g)})" for lab, g in zip(labels, groups)], fontsize=6)


# --------------------------------------------------------------------------- #
def plot_along_protein(bycol, outdir, cats, formats, min_occ):
    prefix = os.path.join(outdir, "disorder_along_protein")
    if bycol is None or not {"alignment_col", "mean_disorder"}.issubset(bycol.columns):
        placeholder(prefix, "no data: disorder-by-column table missing/empty", formats)
        return
    df = bycol.copy()
    df["x"] = pd.to_numeric(df["alignment_col"], errors="coerce")
    df["y"] = pd.to_numeric(df["mean_disorder"], errors="coerce")
    df = _occupancy_filter(df.dropna(subset=["x", "y"]), min_occ).sort_values("x")
    if df.empty:
        placeholder(prefix, f"no data: no column with occupancy >= {min_occ}", formats)
        return
    colors = domain_colors(cats)
    fig, (ax, axd) = plt.subplots(2, 1, figsize=(7.0, 3.0), sharex=True,
                                  gridspec_kw={"height_ratios": [6, 1], "hspace": 0.06})
    ax.vlines(df["x"].values, 0, df["y"].values, color=OKABE_ITO["blue"], linewidth=0.4,
              alpha=0.6)
    if len(df) >= 15:
        roll = df["y"].rolling(15, center=True, min_periods=5).mean()
        ax.plot(df["x"].values, roll.values, color=OKABE_ITO["black"], linewidth=0.8,
                label="15-column running mean")
    ax.set_ylim(0, 1); ax.set_ylabel("mean consensus disorder")
    ax.set_title(f"Disorder along the Prometheus alignment (columns with >= "
                 f"{min_occ:.0%} of scored loci)")
    doms = df["domain"].astype(str).tolist() if "domain" in df.columns else ["NA"] * len(df)
    for s, e, d in domain_runs(df["x"].astype(int).tolist(), doms):
        axd.add_patch(plt.Rectangle((s - 0.5, 0), e - s + 1, 1,
                                    color=colors.get(d, NA_COLOR), linewidth=0))
    axd.set_xlim(df["x"].min() - 1, df["x"].max() + 1)
    axd.set_yticks([]); axd.set_ylim(0, 1); axd.set_xlabel("alignment column")
    handles = [Patch(color=colors[c], label=c) for c in cats]
    handles.append(Patch(color=NA_COLOR, label="linker/tail (NA)"))
    line_handles, _ = ax.get_legend_handles_labels()
    ax.legend(handles=handles + line_handles, frameon=False, fontsize=5,
              ncol=min(6, len(handles) + len(line_handles)), loc="upper right")
    save_fig(fig, prefix, formats)


# --------------------------------------------------------------------------- #
def plot_by_group(per_locus, master, bact_locus, outdir, formats, tiers):
    prefix = os.path.join(outdir, "disorder_by_group")
    if per_locus is None or "genome_id" not in per_locus.columns:
        placeholder(prefix, "no data: disorder_per_locus.tsv missing/empty", formats)
        return
    resp = _resp(per_locus)
    pl, n_excl = tier_filter_loci(per_locus, master, tiers)
    pl = pl.assign(_v=pd.to_numeric(pl[resp], errors="coerce")).dropna(subset=["_v"])
    pg = pl.groupby("genome_id")["_v"].mean()
    panels = []
    if master is not None and {"genome_id", "methanogen"}.issubset(master.columns):
        m = master.copy()
        m["genome_id"] = m["genome_id"].astype(str)
        if "pass_primary" in m.columns:
            m = m[_bool(m["pass_primary"])]
        m = m.set_index("genome_id")
        j = pd.DataFrame({"v": pg}).join(m, how="inner")      # genomes absent from master dropped
        anme = _bool(j["is_anme"]) if "is_anme" in j.columns else pd.Series(False, index=j.index)
        meth = _bool(j["methanogen"]) & ~anme
        groups = [j.loc[~meth & ~anme, "v"].values, j.loc[anme, "v"].values,
                  j.loc[meth, "v"].values]
        labels = ["non-methanogen", "ANME", "methanogen"]
        palette = [OKABE_ITO["grey"], OKABE_ITO["orange"], OKABE_ITO["blue"]]
        keep = [i for i, g in enumerate(groups) if len(g) or labels[i] != "ANME"]
        if sum(len(groups[i]) > 0 for i in keep) >= 1:
            panels.append(("methanogen status", [groups[i] for i in keep],
                           [labels[i] for i in keep], [palette[i] for i in keep]))
    if bact_locus is not None and "genome_id" in bact_locus.columns:
        bresp = _resp(bact_locus)
        bl = bact_locus.assign(_v=pd.to_numeric(bact_locus[bresp], errors="coerce"))
        if tiers and "confidence_tier" in bl.columns:
            bl = bl[bl["confidence_tier"].astype(str).str.strip().isin(tiers)]
        bpg = bl.dropna(subset=["_v"]).groupby("genome_id")["_v"].mean()
        panels.append(("domain of life", [pg.values, bpg.values], ["Archaea", "Bacteria"],
                       [OKABE_ITO["blue"], OKABE_ITO["vermillion"]]))
    if not panels:
        placeholder(prefix, "no data: no grouping (methanogen / bacterial) available", formats)
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(2.8 * len(panels) + 0.5, 2.8),
                             squeeze=False)
    for ax, (title, groups, labels, colors) in zip(axes[0], panels):
        _strip(ax, groups, labels, colors)
        ax.set_ylabel(f"per-genome disorder ({resp})", fontsize=6)
        ax.set_title(title, fontsize=7)
    tier_txt = f"tiers {','.join(tiers)}; {n_excl} genome(s) outside the tier set excluded" \
        if tiers else "all scored loci"
    fig.suptitle(f"ProA disorder by group ({tier_txt}; tests are phylogeny-controlled)",
                 fontsize=7)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    save_fig(fig, prefix, formats)


# --------------------------------------------------------------------------- #
def plot_linker_vs_domain(bycol, outdir, formats, min_occ):
    prefix = os.path.join(outdir, "disorder_linker_vs_domain")
    if bycol is None or not {"mean_disorder", "domain"}.issubset(bycol.columns):
        placeholder(prefix, "no data: disorder-by-column table missing domain/disorder", formats)
        return
    df = bycol.copy()
    df["y"] = pd.to_numeric(df["mean_disorder"], errors="coerce")
    df = _occupancy_filter(df.dropna(subset=["y"]), min_occ)
    df["grp"] = np.where(df["domain"].astype(str).isin(["NA", "nan", ""]), "linker/tail", "domain")
    groups = [df.loc[df["grp"] == g, "y"].values for g in ("domain", "linker/tail")]
    if not all(len(g) for g in groups):
        placeholder(prefix, "no data: need both domain and linker columns", formats)
        return
    fig, ax = plt.subplots(figsize=(3.2, 2.8))
    _strip(ax, groups, ["folded domain", "linker/tail"], [OKABE_ITO["green"], OKABE_ITO["orange"]])
    ax.set_ylabel("mean consensus disorder per column"); ax.set_ylim(0, 1)
    ax.set_title(f"Disorder: domain vs linker columns (occupancy >= {min_occ:.0%})", fontsize=7)
    save_fig(fig, prefix, formats)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Disorder figures")
    ap.add_argument("--by-column", required=True)
    ap.add_argument("--per-locus", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--bact-per-locus", default="")
    ap.add_argument("--outdir", required=True)
    # --- optional (new; backward compatible) ---
    ap.add_argument("--domains", default="", help="prometheus.domains spec as JSON")
    ap.add_argument("--formats", default="png,svg", help="comma-separated, e.g. png,svg")
    ap.add_argument("--tiers", default="",
                    help="tier set T for the group figure (e.g. 'A'); '' = all scored loci")
    ap.add_argument("--min-occupancy", type=float, default=0.5,
                    help="columns need >= this fraction of scored loci to be plotted")
    args = ap.parse_args(argv)
    formats = tuple(f.strip() for f in args.formats.split(",") if f.strip())
    if not formats:
        raise SystemExit("[plot_disorder] --formats is empty")
    tiers = [t.strip() for t in args.tiers.split(",") if t.strip()]
    os.makedirs(args.outdir, exist_ok=True)
    spec = utils.load_domain_spec(args.domains) if args.domains else {}
    bycol = _load(args.by_column)
    per_locus = _load(args.per_locus)
    master = _load(args.master)
    bact = _load(args.bact_per_locus)
    cats = domain_categories(spec, bycol)
    if bycol is not None and "domain" in bycol.columns and spec:
        unknown = sorted(set(bycol["domain"].astype(str)) - set(cats) - {"NA", "", "nan"})
        if unknown:
            raise SystemExit(f"[plot_disorder] domain values {unknown} are not categories "
                             "of the --domains spec")
    plot_along_protein(bycol, args.outdir, cats, formats, args.min_occupancy)
    plot_by_group(per_locus, master, bact, args.outdir, formats, tiers)
    plot_linker_vs_domain(bycol, args.outdir, formats, args.min_occupancy)
    return 0


if __name__ == "__main__":
    sys.exit(main())
