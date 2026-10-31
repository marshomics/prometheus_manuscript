"""Merge bacterial niche/oxygen/biome labels onto the bacterial master table.

Turns a niche metadata TSV (genome_id + columns such as oxygen, biome, host_associated)
into `niche_*` 0/1 columns that `ecology_stats.R` autodetects:
  * a binary column (values 0/1, true/false, yes/no in any case) -> `niche_<col>`;
  * a categorical column -> one-hot `niche_<col>_<level>` (level vs. the other levels);
  * a column already named niche_*/group_* keeps its name (no double prefix).
Missing is NOT absent: a genome with no label for a column gets NA (empty) in that
column's indicators, and genomes that are absent from the niche table, or have no label
at all, are DROPPED from the output (the niche test is run on annotated genomes only).
Coding them 0 would turn "unknown oxygen status" into "not anaerobic" and bias the
phylogenetic regression toward the null. Duplicated genome rows are unioned (a genome
labelled with two levels is 1 for both, the species-level presence semantics of the
archaeal group_/niche_ columns). Levels whose slugs collide (e.g. 'Host-associated' and
'host associated') are merged into one indicator. A missing --niche file is an error.
"""
from __future__ import annotations
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

TRUE_SET = {"1", "true", "yes", "t", "y"}
FALSE_SET = {"0", "false", "no", "f", "n"}
MISSING_SET = {"", "nan", "na", "n/a", "none", "null", "unknown", "?"}


def _slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s).strip()).strip("_").lower()


def _col_name(col):
    s = _slug(col)
    return s if s.startswith(("niche_", "group_")) else f"niche_{s}"


def encode(niche, gcol):
    """Return ({indicator_name: {genome: 0/1/None}}, annotated_genomes)."""
    out = {}
    annotated = set()
    genomes = niche[gcol].tolist()
    for col in niche.columns:
        if col == gcol:
            continue
        vals = [str(v).strip() for v in niche[col].fillna("").tolist()]
        low = [v.lower() for v in vals]
        known = [v for v in low if v not in MISSING_SET]
        if not known:
            utils.log(f"note: niche column '{col}' has no labels; skipped")
            continue
        is_binary = set(known) <= (TRUE_SET | FALSE_SET)
        if is_binary:
            name = _col_name(col)
            rec = out.setdefault(name, {})
            for g, v in zip(genomes, low):
                if v in MISSING_SET:
                    rec.setdefault(g, None)
                    continue
                annotated.add(g)
                x = 1 if v in TRUE_SET else 0
                rec[g] = max(x, rec.get(g) or 0)
            continue
        # categorical -> one-hot; levels merged by slug
        levels = {}
        for v, lv in zip(vals, low):
            if lv not in MISSING_SET:
                levels.setdefault(_slug(v), set()).add(v)
        merged = [k for k, s in levels.items() if len(s) > 1]
        if merged:
            utils.log(f"note: column '{col}': levels merged by name: "
                      f"{[sorted(levels[k]) for k in merged][:3]}")
        labelled = {}
        for g, v, lv in zip(genomes, vals, low):
            if lv in MISSING_SET:
                labelled.setdefault(g, set())
            else:
                labelled.setdefault(g, set()).add(_slug(v))
                annotated.add(g)
        for lev in sorted(levels):
            name = f"{_col_name(col)}_{lev}"
            rec = out.setdefault(name, {})
            for g, labs in labelled.items():
                rec[g] = (1 if lev in labs else 0) if labs else None
    return out, annotated


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", required=True)
    ap.add_argument("--niche", required=True, help="niche metadata TSV (genome_id + labels)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pandas as pd
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    if not os.path.exists(args.niche):
        utils.log(f"FATAL: niche metadata {args.niche} does not exist")
        sys.exit(1)

    master = pd.read_csv(args.master, sep="\t", dtype=str, keep_default_na=False)
    if "genome_id" not in master.columns:
        utils.log(f"FATAL: {args.master} has no genome_id column")
        sys.exit(1)
    master["genome_id"] = master["genome_id"].astype(str).str.strip()

    niche = pd.read_csv(args.niche, sep="\t", dtype=str, comment="#", keep_default_na=False)
    niche.columns = [c.strip() for c in niche.columns]
    if niche.shape[1] < 2:
        utils.log(f"FATAL: niche metadata {args.niche} needs genome_id + >= 1 label column")
        sys.exit(1)
    gcol = "genome_id" if "genome_id" in niche.columns else niche.columns[0]
    niche[gcol] = niche[gcol].astype(str).str.strip()
    niche = niche[niche[gcol] != ""]
    clash = [c for c in niche.columns if c != gcol and _col_name(c) in master.columns]
    if clash:
        utils.log(f"FATAL: niche column(s) {clash} would overwrite master columns")
        sys.exit(1)

    cols, annotated = encode(niche, gcol)
    if len(master) and annotated and not (annotated & set(master["genome_id"])):
        utils.log(f"FATAL: none of the {len(annotated)} labelled niche genomes is in the "
                  f"bacterial master (e.g. '{sorted(annotated)[0]}' vs "
                  f"'{master['genome_id'].iloc[0]}'): genome ids must be identical")
        sys.exit(1)
    unknown = sorted(set(niche[gcol]) - set(master["genome_id"]))
    if unknown:
        utils.log(f"WARNING: {len(unknown)} niche-table genome(s) are not in the bacterial "
                  f"master and are ignored (e.g. {unknown[:3]})")
    keep = master["genome_id"].isin(annotated)
    n_drop = int((~keep).sum())
    out = master[keep].copy()
    for name in sorted(cols):
        mapping = cols[name]
        out[name] = [("" if mapping.get(g) is None else str(int(mapping[g])))
                     for g in out["genome_id"]]
    out.to_csv(args.out, sep="\t", index=False)
    utils.log(f"bacterial niche master: {len(cols)} niche_* indicator(s) for {len(out)} "
              f"annotated genome(s); {n_drop} genome(s) without niche labels dropped -> {args.out}")


if __name__ == "__main__":
    main()
