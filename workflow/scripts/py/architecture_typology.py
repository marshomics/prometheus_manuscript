"""Enumerate distinct Prometheus domain architectures (design §2.2).

Retained (non-Rejected) catalog loci are grouped by their N->C architecture, recomputed
from the domtbl through the shared spec helpers (utils.architecture_tokens /
is_canonical_order / required_categories / other_domains), i.e. exactly the function
validate_prometheus used for the catalog's `architecture` column (a disagreement between
the two is an input error). Per architecture the table reports:

  architecture, n_loci, example_locus           (SCHEMA)
  has_DUF4011, has_DUF3320, has_DNA2, has_nuclease, has_REase
      (SCHEMA columns, kept for compatibility; `has_nuclease` is True only if the spec
      defines a 'nuclease' category -- the current spec has none: ProA's nuclease is the
      REase_MTES_1575 domain -> has_REase)
  has_<category> for any further spec category, then (additive):
  core_complete    every `required: true` category present
  canonical_order  each category once and in spec (N->C) order
  rearranged       core present but out of order / duplicated
  variant_class    canonical_full | canonical_core (an optional domain, e.g. DUF3320,
                   absent) | rearranged | core_loss | none
  missing_required, n_with_fusion, fusion_domains ("Pkinase:2;ResIII+Helicase_C:1": Pfam
                   hits outside the spec, the fusions/insertions of design §2.2),
  n_tier_A, n_tier_B, n_tier_C
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

SCHEMA_COLUMNS = ["architecture", "n_loci", "example_locus", "has_DUF4011", "has_DUF3320",
                  "has_DNA2", "has_nuclease", "has_REase"]
EXTRA_COLUMNS = ["core_complete", "canonical_order", "rearranged", "variant_class",
                 "missing_required", "n_with_fusion", "fusion_domains",
                 "n_tier_A", "n_tier_B", "n_tier_C"]
HEADER = "\t".join(SCHEMA_COLUMNS)   # backward-compatible constant


def columns_for(spec):
    extra_has = [f"has_{c}" for c in utils.canonical_order(spec)
                 if f"has_{c}" not in SCHEMA_COLUMNS]
    return SCHEMA_COLUMNS + extra_has + EXTRA_COLUMNS


def classify(tokens, spec):
    present = set(tokens)
    required = utils.required_categories(spec)
    missing = [c for c in required if c not in present]
    canonical = bool(tokens) and utils.is_canonical_order(tokens, spec)
    if not tokens:
        cls = "none"
    elif missing:
        cls = "core_loss"
    elif not canonical:
        cls = "rearranged"
    elif all(c in present for c in utils.canonical_order(spec)):
        cls = "canonical_full"
    else:
        cls = "canonical_core"
    return {"core_complete": not missing, "canonical_order": canonical,
            "rearranged": (not missing) and bool(tokens) and not canonical,
            "variant_class": cls, "missing_required": ",".join(missing)}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--domains", required=True, help="JSON domain definitions")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    try:
        spec = utils.load_domain_spec(args.domains)
    except ValueError as exc:
        raise SystemExit(f"[architecture_typology] FATAL: invalid --domains JSON: {exc}")
    if not spec:
        raise SystemExit("[architecture_typology] FATAL: empty --domains spec")
    cols = columns_for(spec)

    cat = pd.read_csv(args.catalog, sep="\t", dtype=str, keep_default_na=False)
    if "locus_id" not in cat.columns:
        raise SystemExit(f"[architecture_typology] FATAL: {args.catalog} has no locus_id")
    if "confidence_tier" in cat.columns:
        cat = cat[cat["confidence_tier"] != "Rejected"]
    if cat.empty:
        utils.write_empty(args.out, "\t".join(cols))
        utils.log("No retained loci; wrote header-only typology.")
        return

    hits = utils.parse_domtbl(args.domtbl)
    groups = {}
    mismatch = []
    for r in cat.to_dict("records"):
        locus = r["locus_id"]
        h = hits.get(locus, [])
        toks = utils.architecture_tokens(h, spec)
        arch = "-".join(toks)
        if "architecture" in r and r["architecture"] != arch:
            mismatch.append(f"{locus}: catalog '{r['architecture']}' vs domtbl '{arch}'")
        key = arch or "none"
        g = groups.setdefault(key, {"n": 0, "example": locus, "tokens": toks,
                                    "fusions": {}, "n_fusion": 0,
                                    "tiers": {"A": 0, "B": 0, "C": 0}})
        g["n"] += 1
        other = utils.other_domains(h, spec)
        if other:
            g["n_fusion"] += 1
            k = "+".join(other)
            g["fusions"][k] = g["fusions"].get(k, 0) + 1
        t = r.get("confidence_tier", "")
        if t in g["tiers"]:
            g["tiers"][t] += 1
    if mismatch:
        raise SystemExit(f"[architecture_typology] FATAL: {len(mismatch)} catalog "
                         f"architecture(s) disagree with the domtbl (e.g. {mismatch[:2]}); "
                         "the catalog and domtbl come from different runs or specs")

    rows = []
    for arch in sorted(groups, key=lambda a: (-groups[a]["n"], a)):
        g = groups[arch]
        present = set(g["tokens"])
        row = {"architecture": arch, "n_loci": g["n"], "example_locus": g["example"]}
        for c in ("DUF4011", "DUF3320", "DNA2", "nuclease", "REase"):
            row[f"has_{c}"] = c in present
        for c in utils.canonical_order(spec):
            row[f"has_{c}"] = c in present
        row.update(classify(g["tokens"], spec))
        row["n_with_fusion"] = g["n_fusion"]
        row["fusion_domains"] = ";".join(f"{k}:{v}" for k, v in
                                         sorted(g["fusions"].items(), key=lambda kv: (-kv[1], kv[0])))
        for t in ("A", "B", "C"):
            row[f"n_tier_{t}"] = g["tiers"][t]
        rows.append(row)
    pd.DataFrame(rows, columns=cols).to_csv(args.out, sep="\t", index=False)
    utils.log(f"Typology: {len(rows)} distinct architectures over "
              f"{sum(g['n'] for g in groups.values())} retained loci.")


if __name__ == "__main__":
    main()
