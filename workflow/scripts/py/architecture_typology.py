"""Enumerate distinct Prometheus domain architectures (design §2.2).

Groups validated (non-Rejected) loci by their N->C domain-architecture string and
reports, per architecture, the locus count, one example locus, and which of the
five canonical domains it contains. The architecture strings are taken from the
Stage-1 catalog (already computed and tier-filtered); the domtbl + domain JSON
are accepted so the typology can be recomputed if a catalog lacks the column.

Architecture variation (canonical four/five-domain form, domain losses,
insertions, fusions/fissions, reorderings) is biologically meaningful because the
domains may have independent evolutionary histories (tested in Stage 5).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import json
import pandas as pd

HEADER = (
    "architecture\tn_loci\texample_locus\t"
    "has_DUF4011\thas_DUF3320\thas_DNA2\thas_nuclease\thas_REase"
)
# Fixed schema columns (independent of how many domains the config defines).
CANONICAL = ["DUF4011", "DUF3320", "DNA2", "nuclease", "REase"]


def recompute_from_domtbl(domtbl_path, domains, loci):
    """Fallback: rebuild architecture per locus from the domtbl + domain JSON.

    Mirrors validate_prometheus.architecture_for so the typology is identical to
    the catalog when recomputed. Only `loci` (the retained locus_ids) are kept.
    """
    # import the sibling implementation to avoid divergence
    import validate_prometheus as vp
    hits = utils.parse_domtbl(domtbl_path)
    out = {}
    for locus in loci:
        arch, _, _ = vp.architecture_for(hits.get(locus, []), domains)
        out[locus] = arch
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--domains", required=True, help="JSON domain definitions")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    try:
        domains = json.loads(args.domains) if args.domains else {}
    except json.JSONDecodeError:
        domains = {}

    try:
        cat = pd.read_csv(args.catalog, sep="\t")
    except Exception:
        cat = pd.DataFrame()

    if cat.empty or "locus_id" not in cat.columns:
        utils.write_empty(args.out, HEADER)
        utils.log("Empty/invalid catalog; wrote header-only typology.")
        return

    # Only retained loci (Rejected loci are not part of the curated catalog).
    if "confidence_tier" in cat.columns:
        cat = cat[cat["confidence_tier"].astype(str) != "Rejected"]
    if cat.empty:
        utils.write_empty(args.out, HEADER)
        utils.log("No non-Rejected loci; wrote header-only typology.")
        return

    # Architecture per locus: prefer the catalog column, else recompute.
    if "architecture" in cat.columns:
        arch_by_locus = {
            str(r["locus_id"]): ("" if pd.isna(r["architecture"]) else str(r["architecture"]))
            for _, r in cat.iterrows()
        }
    else:
        loci = [str(x) for x in cat["locus_id"].tolist()]
        arch_by_locus = recompute_from_domtbl(args.domtbl, domains, loci)

    # Group identical architecture strings; keep first-seen locus as example.
    groups = {}  # arch -> {"n":int, "example":str}
    for locus, arch in arch_by_locus.items():
        g = groups.setdefault(arch, {"n": 0, "example": locus})
        g["n"] += 1

    rows = []
    for arch in sorted(groups.keys(), key=lambda a: (-groups[a]["n"], a)):
        present_tokens = set(arch.split("-")) if arch else set()
        rows.append({
            "architecture": arch,
            "n_loci": groups[arch]["n"],
            "example_locus": groups[arch]["example"],
            "has_DUF4011": "DUF4011" in present_tokens,
            "has_DUF3320": "DUF3320" in present_tokens,
            "has_DNA2": "DNA2" in present_tokens,
            "has_nuclease": "nuclease" in present_tokens,
            "has_REase": "REase" in present_tokens,
        })

    pd.DataFrame(rows, columns=HEADER.split("\t")).to_csv(
        args.out, sep="\t", index=False)
    utils.log(f"Typology: {len(rows)} distinct architectures over "
              f"{sum(g['n'] for g in groups.values())} retained loci.")


if __name__ == "__main__":
    main()
