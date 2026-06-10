"""Build a slim bacterial master table for the cross-domain comparative analyses.

One row per bacterial genome, with the columns the reused archaeal scripts need:
  genome_id, tree_tip, gtdb_domain..gtdb_species, completeness, contamination,
  pass_primary, prometheus_present, prometheus_tier, prometheus_n_loci,
  n_defense_systems, n_defense_systems_excl_prometheus

`compute_prevalence.py` and `cooccurrence.R` are reused unchanged on this table
(prevalence keys on pass_primary + prometheus_present/tier + gtdb_* ; co-occurrence
keys on genome_id + tree_tip). Defense-system counts come from the per-genome systems
files split out by `extract_bacterial_prometheus.py`, so the EXCL-Prometheus covariate
is consistent with the archaeal definition. Empty-input safe.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402
import build_master_table as bmt  # reuse parse_taxonomy/load_prometheus/count_defense_systems

HEADER = [
    "genome_id", "tree_tip",
    "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
    "gtdb_family", "gtdb_genus", "gtdb_species",
    "completeness", "contamination", "genome_size", "pass_primary",
    "prometheus_present", "prometheus_tier", "prometheus_n_loci",
    "n_defense_systems", "n_defense_systems_excl_prometheus",
]


def load_qc(path):
    """genome_id -> {completeness, contamination, genome_size} (strings)."""
    import pandas as pd
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    except Exception:
        return out
    gcol = "genome_id" if "genome_id" in df.columns else df.columns[0]
    comp = next((c for c in df.columns if c.lower() in ("completeness", "comp")), None)
    cont = next((c for c in df.columns if c.lower() in ("contamination", "contam")), None)
    size = next((c for c in df.columns if c.lower() in
                 ("genome_size", "length", "total_length", "genome_length", "size", "contig_bp")), None)
    for _, r in df.iterrows():
        gid = str(r[gcol]).strip()
        if gid:
            out[gid] = {"completeness": bmt._blank(r.get(comp, "")) if comp else "",
                        "contamination": bmt._blank(r.get(cont, "")) if cont else "",
                        "genome_size": bmt._blank(r.get(size, "")) if size else ""}
    return out


def assembly_size(genomes_dir, genome_id, ext):
    """Total residue count of an assembly (streamed; sequence chars only). '' if absent."""
    if not genomes_dir:
        return ""
    path = os.path.join(genomes_dir, f"{genome_id}.{ext}")
    if not os.path.exists(path):
        for e in ("fasta", "fna", "fa"):
            alt = os.path.join(genomes_dir, f"{genome_id}.{e}")
            if os.path.exists(alt):
                path = alt
                break
    if not os.path.exists(path):
        return ""
    try:
        n = 0
        with open(path) as fh:
            for line in fh:
                if not line.startswith(">"):
                    n += len(line.strip())
        return n if n else ""
    except Exception:
        return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--qc", required=True)
    ap.add_argument("--catalog", required=True, help="bacterial_prometheus_catalog.tsv")
    ap.add_argument("--defense-dir", required=True, help="split per-genome systems dir")
    ap.add_argument("--genomes-dir", default="", help="optional assemblies dir (for genome_size)")
    ap.add_argument("--genome-ext", default="fasta")
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pandas as pd
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    taxonomy = bmt.parse_taxonomy(args.taxonomy)   # accepts user_genome/classification
    # also accept a genome_id-keyed taxonomy table
    if not taxonomy:
        try:
            tdf = pd.read_csv(args.taxonomy, sep="\t", dtype=str, comment="#")
            gcol = "genome_id" if "genome_id" in tdf.columns else (
                "user_genome" if "user_genome" in tdf.columns else tdf.columns[0])
            ccol = next((c for c in tdf.columns if c.lower() == "classification"), None)
            for _, r in tdf.iterrows():
                gid = str(r[gcol]).strip()
                if not gid:
                    continue
                ranks = {col: "" for _, col in bmt.GTDB_RANKS}
                cstr = str(r.get(ccol, "")).strip() if ccol else ""
                if cstr and cstr.lower() not in ("nan", "unclassified", "na", "none"):
                    for token in [p.strip() for p in cstr.split(";")]:
                        for pre, col in bmt.GTDB_RANKS:
                            if token.startswith(pre):
                                ranks[col] = token[len(pre):].strip()
                taxonomy[gid] = ranks
        except Exception:
            pass

    qc = load_qc(args.qc)
    prometheus = bmt.load_prometheus(args.catalog)

    genomes = list(dict.fromkeys(list(qc.keys()) + list(taxonomy.keys())))
    if not genomes:
        utils.write_empty(args.out, "\t".join(HEADER))
        utils.log("No bacterial genomes (empty qc/taxonomy); header-only table.")
        return

    rows = []
    for gid in genomes:
        tax = taxonomy.get(gid, {col: "" for _, col in bmt.GTDB_RANKS})
        q = qc.get(gid, {"completeness": "", "contamination": ""})
        comp = q["completeness"]
        try:
            pass_primary = (float(comp) >= args.min_completeness) if comp != "" else False
        except (ValueError, TypeError):
            pass_primary = False
        prom = prometheus.get(gid, {"present": False, "n_loci": 0, "best": ""})
        n_def, n_def_excl = bmt.count_defense_systems(args.defense_dir, gid)
        gsize = q.get("genome_size", "") or assembly_size(args.genomes_dir, gid, args.genome_ext)
        rows.append({
            "genome_id": gid, "tree_tip": gid,
            "gtdb_domain": tax.get("gtdb_domain", ""), "gtdb_phylum": tax.get("gtdb_phylum", ""),
            "gtdb_class": tax.get("gtdb_class", ""), "gtdb_order": tax.get("gtdb_order", ""),
            "gtdb_family": tax.get("gtdb_family", ""), "gtdb_genus": tax.get("gtdb_genus", ""),
            "gtdb_species": tax.get("gtdb_species", ""),
            "completeness": comp, "contamination": q["contamination"], "genome_size": gsize,
            "pass_primary": pass_primary,
            "prometheus_present": prom["present"], "prometheus_tier": prom["best"],
            "prometheus_n_loci": prom["n_loci"],
            "n_defense_systems": n_def, "n_defense_systems_excl_prometheus": n_def_excl,
        })

    pd.DataFrame(rows, columns=HEADER).to_csv(args.out, sep="\t", index=False)
    n_prom = sum(1 for r in rows if r["prometheus_present"])
    utils.log(f"bacterial table: {len(rows)} genomes ({n_prom} Prometheus-present).")


if __name__ == "__main__":
    main()
