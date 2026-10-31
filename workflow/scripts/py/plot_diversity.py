"""plot_diversity.py — subfamily rarefaction and architecture-by-taxon figures.

Outputs (written to --outdir):
  subfamily_rarefaction.{png,svg}  sample-based (incidence) rarefaction of subfamilies
                                   against genomes sampled (design §2, "rarefaction of
                                   subfamilies vs genomes sampled")
  architecture_by_taxon.{png,svg}  per-locus domain architecture (catalog) x GTDB phylum
  --curve TSV (optional)           the numbers behind the rarefaction figure

Rarefaction. The sampling unit is the genome and the universe is EVERY analysis genome
(master_table genome_id), genomes without Prometheus included, so the x-axis is sampling
effort. With N genomes and subfamily j present in k_j of them, the expected number of
subfamilies in a random subset of n genomes is
    E[S_n] = sum_j 1 - C(N - k_j, n) / C(N, n)          (exact; Gotelli & Colwell 2001,
             Ecol Lett 4:379-391, doi:10.1046/j.1461-0248.2001.00230.x; Colwell et al. 2012,
             J Plant Ecol 5:3-21, doi:10.1093/jpe/rtr044)
computed with log-gamma. The shaded band is the 2.5-97.5% range over --permutations
random genome orders (seeded with --seed). Chao2 (incidence; Chao 1987, Biometrics 43:783, doi:10.2307/2531532; bias-corrected
form when Q2 = 0) is annotated as a lower-bound richness estimate. The level is the
lowest identity column (the subfamily level, cluster_0.30) unless --level is given.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from plotting_style import apply_style, save_fig, placeholder, OKABE_ITO
import utils  # noqa: F401

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.special import gammaln

apply_style()


def _load(path, name):
    """DataFrame or None for a missing / header-only table; malformed files raise."""
    if not path or not os.path.exists(path):
        return None, f"{name} missing"
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df.empty:
        return None, f"{name} empty"
    return df, None


def cluster_columns(df):
    cols = [c for c in df.columns if c.startswith("cluster_")]
    out = []
    for c in cols:
        try:
            out.append((float(c.split("_", 1)[1]), c))
        except ValueError:
            continue
    return [c for _, c in sorted(out)]


def expected_richness(incidence_k, n_total, sizes):
    """Exact sample-based rarefaction E[S_n] for each n in sizes."""
    k = np.asarray(incidence_k, dtype=float)
    out = []
    lgN = gammaln(n_total + 1)
    for n in sizes:
        if n <= 0:
            out.append(0.0)
            continue
        # log C(N-k, n) - log C(N, n); C(N-k, n) = 0 when N-k < n
        with np.errstate(invalid="ignore"):
            ok = (n_total - k) >= n
            lr = np.full(k.shape, -np.inf)
            nk = n_total - k[ok]
            lr[ok] = (gammaln(nk + 1) - gammaln(n + 1) - gammaln(nk - n + 1)
                      - (lgN - gammaln(n + 1) - gammaln(n_total - n + 1)))
        out.append(float(np.sum(1.0 - np.exp(lr))))
    return np.array(out)


def chao2(incidence_k, n_total):
    k = np.asarray(incidence_k)
    s_obs = int((k > 0).sum())
    q1, q2 = int((k == 1).sum()), int((k == 2).sum())
    m = float(n_total)
    if m <= 1:
        return float(s_obs)
    if q2 > 0:
        return s_obs + ((m - 1) / m) * q1 * q1 / (2.0 * q2)
    return s_obs + ((m - 1) / m) * q1 * (q1 - 1) / 2.0


def rarefaction_table(sub, universe, level, n_perm, seed):
    """DataFrame(n_genomes, expected, lo, hi) + summary dict."""
    sub = sub[sub[level].astype(str) != ""]
    genomes = list(dict.fromkeys(universe))
    gidx = {g: i for i, g in enumerate(genomes)}
    outside = sorted(set(sub["genome_id"]) - set(gidx))
    if outside:
        raise SystemExit(f"[plot_diversity] FATAL: {len(outside)} subfamily genome(s) missing "
                         f"from the master table (e.g. {outside[:3]})")
    inc = sub.groupby(level)["genome_id"].apply(lambda s: sorted(set(s)))
    k = np.array([len(v) for v in inc.values])
    N = len(genomes)
    step = max(1, N // 60)
    sizes = sorted(set(list(range(0, N + 1, step)) + [N]))
    exp = expected_richness(k, N, sizes)
    rng = np.random.default_rng(seed)
    members = [np.array([gidx[g] for g in v]) for v in inc.values]
    curves = np.zeros((n_perm, len(sizes)))
    for b in range(n_perm):
        rank = rng.permutation(N)              # rank[g] = sampling order of genome g
        first = np.array([rank[m].min() for m in members]) if members else np.array([])
        curves[b] = [(first < n).sum() for n in sizes]
    lo = np.percentile(curves, 2.5, axis=0) if n_perm else exp
    hi = np.percentile(curves, 97.5, axis=0) if n_perm else exp
    df = pd.DataFrame({"n_genomes": sizes, "expected_subfamilies": exp,
                       "perm_lo95": lo, "perm_hi95": hi})
    summ = {"N": N, "S_obs": int(len(k)), "Q1": int((k == 1).sum()),
            "Q2": int((k == 2).sum()), "chao2": chao2(k, N),
            "n_carrying": int(sub["genome_id"].nunique())}
    return df, summ


def plot_rarefaction(subfamilies, master, outdir, seed, n_perm, level, curve_out):
    prefix = os.path.join(outdir, "subfamily_rarefaction")
    if subfamilies is None:
        placeholder(prefix, "no data: subfamilies.tsv missing or empty")
        if curve_out:
            utils.write_empty(curve_out, "n_genomes\texpected_subfamilies\tperm_lo95\tperm_hi95")
        return
    cols = cluster_columns(subfamilies)
    if "genome_id" not in subfamilies.columns or not cols:
        raise SystemExit("[plot_diversity] FATAL: subfamilies.tsv needs genome_id and "
                         f"cluster_* columns (has {list(subfamilies.columns)})")
    level = level or cols[0]
    if level not in subfamilies.columns:
        raise SystemExit(f"[plot_diversity] FATAL: --level {level} not in {cols}")
    if master is not None and "genome_id" in master.columns:
        universe, xlab = master["genome_id"].tolist(), "Genomes sampled (all analysis genomes)"
    else:
        universe = subfamilies["genome_id"].tolist()
        xlab = "Genomes sampled (Prometheus-carrying genomes only; master table absent)"
    df, s = rarefaction_table(subfamilies, universe, level, n_perm, seed)
    if curve_out:
        df.to_csv(curve_out, sep="\t", index=False)

    ident = level.split("_", 1)[1]
    fig, ax = plt.subplots(figsize=(3.6, 2.8))
    ax.fill_between(df["n_genomes"], df["perm_lo95"], df["perm_hi95"], alpha=0.25,
                    color=OKABE_ITO["blue"], linewidth=0,
                    label=f"95% range ({n_perm} random genome orders)")
    ax.plot(df["n_genomes"], df["expected_subfamilies"], color=OKABE_ITO["blue"],
            linewidth=1.2, label="expected (exact rarefaction)")
    ax.axhline(s["chao2"], color=OKABE_ITO["vermillion"], linewidth=0.8, linestyle="--",
               label=f"Chao2 = {s['chao2']:.1f}")
    ax.set_xlabel(xlab)
    ax.set_ylabel(f"Subfamilies (MMseqs2 >= {float(ident):.0%} identity)")
    ax.set_title(f"Subfamily rarefaction: S_obs = {s['S_obs']} in {s['n_carrying']} of "
                 f"{s['N']} genomes (singletons Q1 = {s['Q1']}, Q2 = {s['Q2']})", fontsize=6)
    ax.set_xlim(0, max(1, s["N"]))
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=5, loc="lower right")
    save_fig(fig, prefix)


def plot_architecture_by_taxon(catalog, master, outdir):
    prefix = os.path.join(outdir, "architecture_by_taxon")
    if catalog is None:
        placeholder(prefix, "no data: prometheus_catalog.tsv not given or empty")
        return
    if master is None or not {"genome_id", "gtdb_phylum"}.issubset(master.columns):
        placeholder(prefix, "no data: master table missing or lacks gtdb_phylum")
        return
    need = {"genome_id", "architecture"}
    if not need.issubset(catalog.columns):
        raise SystemExit(f"[plot_diversity] FATAL: catalog lacks {sorted(need - set(catalog.columns))}")
    cat = catalog
    if "confidence_tier" in cat.columns:
        cat = cat[cat["confidence_tier"] != "Rejected"]
    cat = cat.assign(architecture=cat["architecture"].replace("", "none"))
    merged = cat.merge(master[["genome_id", "gtdb_phylum"]], on="genome_id", how="left")
    merged["gtdb_phylum"] = merged["gtdb_phylum"].replace("", "unclassified").fillna("unclassified")
    if merged.empty:
        placeholder(prefix, "no data: no retained loci")
        return
    top_phyla = merged["gtdb_phylum"].value_counts().head(10).index.tolist()
    df = merged[merged["gtdb_phylum"].isin(top_phyla)].copy()
    top_arch = df["architecture"].value_counts().head(8).index.tolist()
    df["architecture"] = df["architecture"].where(df["architecture"].isin(top_arch), "other")
    arch_list = top_arch + (["other"] if (df["architecture"] == "other").any() else [])
    palette = list(OKABE_ITO.values())
    color = {a: palette[i % len(palette)] for i, a in enumerate(arch_list)}
    counts = (df.groupby(["gtdb_phylum", "architecture"]).size().unstack(fill_value=0)
              .reindex(index=top_phyla, columns=arch_list, fill_value=0))
    totals = counts.sum(axis=1)
    fracs = counts.div(totals, axis=0) * 100

    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    bottom = np.zeros(len(fracs))
    for a in arch_list:
        ax.bar(range(len(fracs)), fracs[a].values, bottom=bottom, color=color[a],
               edgecolor="none", width=0.7)
        bottom += fracs[a].values
    for i, n in enumerate(totals.values):
        ax.text(i, 101, f"n={int(n)}", ha="center", va="bottom", fontsize=4.5)
    ax.set_xticks(range(len(fracs)))
    ax.set_xticklabels([p.replace("p__", "") for p in fracs.index], rotation=45, ha="right",
                       fontsize=5)
    ax.set_ylabel("% of retained loci")
    ax.set_ylim(0, 108)
    ax.set_title("Domain architecture by phylum (top 10 phyla, retained loci)", fontsize=6)
    ax.legend(handles=[mpatches.Patch(color=color[a], label=a) for a in arch_list],
              title="Architecture (N->C)", frameon=False, fontsize=4.5, title_fontsize=5,
              bbox_to_anchor=(1.01, 1), loc="upper left")
    save_fig(fig, prefix)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--subfamilies", required=True, help="diversity/subfamilies.tsv")
    parser.add_argument("--typology", required=True, help="diversity/architecture_typology.tsv")
    parser.add_argument("--master", required=True, help="scaffolds/master_table.tsv")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--catalog", default="", help="prometheus_catalog.tsv (per-locus architecture)")
    parser.add_argument("--permutations", type=int, default=200)
    parser.add_argument("--level", default="", help="cluster column (default: lowest identity)")
    parser.add_argument("--curve", default="", help="optional TSV of the rarefaction curve")
    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    subfamilies, _ = _load(args.subfamilies, "subfamilies")
    _load(args.typology, "typology")      # validated for readability (not plotted here)
    master, _ = _load(args.master, "master_table")
    catalog, _ = _load(args.catalog, "catalog") if args.catalog else (None, None)

    plot_rarefaction(subfamilies, master, args.outdir, args.seed, args.permutations,
                     args.level, args.curve)
    plot_architecture_by_taxon(catalog, master, args.outdir)


if __name__ == "__main__":
    main()
