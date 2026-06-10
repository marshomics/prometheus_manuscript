"""Extra diversity / architecture figures (added alongside plot_diversity.py).

Outputs (written to --outdir):
  architecture_combinations.{png,svg}  UpSet-style: which of the 5 canonical domains
                                       co-occur, and how many loci per combination.
  architecture_by_phylum.{png,svg}     architecture x phylum count heatmap (which
                                       architectural variants live in which lineages).
  subfamily_identity_curve.{png,svg}   number of Prometheus subfamilies vs clustering
                                       identity + cluster-size rank-abundance.

Architecture variation (domain loss/insertion/fusion, reordering) is biologically
meaningful because the domains can have independent histories; these views make the
combinatorial structure and its taxonomic spread legible.
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
CANON = ["DUF4011", "DUF3320", "DNA2", "nuclease", "REase"]


def _load(path):
    if not path or not os.path.exists(path):
        return None
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return None
    return None if df.empty else df


# --------------------------------------------------------------------------- #
def plot_architecture_combinations(typ, outdir):
    prefix = os.path.join(outdir, "architecture_combinations")
    if typ is None or "architecture" not in typ.columns:
        placeholder(prefix, "no data: architecture_typology.tsv missing/empty")
        return
    df = typ.copy()
    if "n_loci" not in df.columns:
        df["n_loci"] = 1
    df = df.sort_values("n_loci", ascending=False).head(18)
    # domain-presence matrix from has_* columns (fall back to parsing the arch string)
    present = np.zeros((len(df), len(CANON)), dtype=bool)
    for i, (_, r) in enumerate(df.iterrows()):
        toks = set(str(r["architecture"]).split("-"))
        for j, dom in enumerate(CANON):
            col = f"has_{dom}"
            if col in df.columns:
                present[i, j] = str(r[col]).strip().lower() in ("true", "1", "yes")
            else:
                present[i, j] = dom in toks

    df = df.iloc[::-1].reset_index(drop=True)          # largest at top
    present = present[::-1]
    y = np.arange(len(df))
    fig, (axm, axb) = plt.subplots(
        1, 2, figsize=(5.2, max(3.0, len(df) * 0.26)), sharey=True,
        gridspec_kw={"width_ratios": [1.0, 1.6], "wspace": 0.05})
    # dot matrix
    for j in range(len(CANON)):
        axm.scatter(np.full(len(df), j), y,
                    c=[OKABE_ITO["blue"] if present[i, j] else "#E8E8E8" for i in range(len(df))],
                    s=26, edgecolor="none")
        for i in range(len(df)):
            pass
    # connect present domains within a row
    for i in range(len(df)):
        js = np.where(present[i])[0]
        if len(js) > 1:
            axm.plot([js.min(), js.max()], [i, i], color=OKABE_ITO["blue"], linewidth=0.8, zorder=1)
    axm.set_xticks(range(len(CANON)))
    axm.set_xticklabels(CANON, rotation=45, ha="right", fontsize=5)
    axm.set_yticks(y); axm.set_yticklabels([])
    axm.set_xlim(-0.5, len(CANON) - 0.5)
    axm.set_title("domains present", fontsize=6)
    axm.spines["left"].set_visible(False)
    # n_loci bars
    axb.barh(y, df["n_loci"].values, color=OKABE_ITO["blue"], alpha=0.85, height=0.6)
    for yi, v in zip(y, df["n_loci"].values):
        axb.text(v, yi, f" {int(v)}", va="center", fontsize=5)
    axb.set_xlabel("n loci")
    axb.set_title("architecture frequency", fontsize=6)
    fig.suptitle("Prometheus domain-architecture combinations", fontsize=8)
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_architecture_by_phylum(catalog, master, outdir):
    prefix = os.path.join(outdir, "architecture_by_phylum")
    if catalog is None or not {"genome_id", "architecture"}.issubset(catalog.columns):
        placeholder(prefix, "no data: prometheus_catalog.tsv missing/empty")
        return
    cat = catalog.copy()
    if "confidence_tier" in cat.columns:
        cat = cat[cat["confidence_tier"].astype(str) != "Rejected"]
    if cat.empty:
        placeholder(prefix, "no data: no retained loci in catalog")
        return
    if master is not None and {"genome_id", "gtdb_phylum"}.issubset(master.columns):
        ph = dict(zip(master["genome_id"].astype(str), master["gtdb_phylum"].astype(str)))
        cat["phylum"] = cat["genome_id"].astype(str).map(ph).fillna("")
    else:
        placeholder(prefix, "no data: master table lacks gtdb_phylum")
        return
    cat = cat[cat["phylum"].str.len() > 0]
    if cat.empty:
        placeholder(prefix, "no data: no phylum assignment for loci")
        return
    top_arch = cat["architecture"].value_counts().head(12).index
    top_phyla = cat["phylum"].value_counts().head(12).index
    sub = cat[cat["architecture"].isin(top_arch) & cat["phylum"].isin(top_phyla)]
    mat = sub.pivot_table(index="architecture", columns="phylum", values="genome_id",
                          aggfunc="count", fill_value=0).reindex(index=top_arch, columns=top_phyla, fill_value=0)

    fig, ax = plt.subplots(figsize=(max(4.0, 0.4 * mat.shape[1] + 2),
                                    max(3.0, 0.32 * mat.shape[0] + 1.5)))
    im = ax.imshow(mat.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=55, ha="right", fontsize=5)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=5)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            if v:
                ax.text(j, i, str(int(v)), ha="center", va="center", fontsize=4.5,
                        color="white" if v < mat.values.max() * 0.6 else "black")
    ax.set_title("Architecture x phylum (retained-locus counts)")
    fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="n loci")
    save_fig(fig, prefix)


# --------------------------------------------------------------------------- #
def plot_subfamily_identity_curve(sub, outdir):
    prefix = os.path.join(outdir, "subfamily_identity_curve")
    if sub is None:
        placeholder(prefix, "no data: subfamilies.tsv missing/empty")
        return
    clu_cols = [c for c in sub.columns if c.startswith("cluster_")]
    if not clu_cols:
        placeholder(prefix, "no data: no cluster_* columns in subfamilies.tsv")
        return
    ids, counts = [], []
    for c in sorted(clu_cols, key=lambda x: float(x.split("_")[1])):
        try:
            ident = float(c.split("_")[1])
        except ValueError:
            continue
        n = sub[c].replace("", np.nan).dropna().nunique()
        ids.append(ident); counts.append(n)
    if not ids:
        placeholder(prefix, "no data: cluster columns unparseable")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.0, 2.8))
    ax1.plot(ids, counts, "-o", color=OKABE_ITO["blue"], markersize=4)
    for x, yv in zip(ids, counts):
        ax1.text(x, yv, f" {yv}", fontsize=5, va="bottom")
    ax1.set_xlabel("clustering identity"); ax1.set_ylabel("n subfamilies")
    ax1.set_title("subfamily count vs identity")
    ax1.set_ylim(bottom=0)
    # rank-abundance at the lowest (broadest) identity
    finest = sorted(clu_cols, key=lambda x: float(x.split("_")[1]))[0]
    sizes = sub[finest].replace("", np.nan).dropna().value_counts().sort_values(ascending=False).values
    if len(sizes):
        ax2.plot(np.arange(1, len(sizes) + 1), sizes, "-o", color=OKABE_ITO["vermillion"], markersize=3)
        ax2.set_yscale("log")
    ax2.set_xlabel("subfamily rank"); ax2.set_ylabel("size (loci)")
    ax2.set_title(f"rank-abundance @ {finest.split('_')[1]} id")
    save_fig(fig, prefix)


def main():
    ap = argparse.ArgumentParser(description="Extra diversity/architecture figures")
    ap.add_argument("--typology", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--subfamilies", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plot_architecture_combinations(_load(args.typology), args.outdir)
    plot_architecture_by_phylum(_load(args.catalog), _load(args.master), args.outdir)
    plot_subfamily_identity_curve(_load(args.subfamilies), args.outdir)


if __name__ == "__main__":
    main()
