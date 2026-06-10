"""Merge bacterial niche/oxygen/biome labels onto the bacterial master table.

Turns a free-form niche metadata TSV (genome_id + columns such as oxygen, biome,
host_associated) into binary `niche_*` presence columns that `ecology_stats.R`
autodetects, then left-joins them onto the bacterial master (missing -> 0). Binary
0/1 columns pass through as `niche_<col>`; categorical columns are one-hot encoded as
`niche_<col>_<value>`. If no niche file is given, the master is passed through
unchanged (the niche control then finds no predictors and writes a header-only result).
"""
from __future__ import annotations
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402


def _slug(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s).strip()).strip("_").lower()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--master", required=True)
    ap.add_argument("--niche", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pandas as pd
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    try:
        master = pd.read_csv(args.master, sep="\t", dtype=str)
    except Exception:
        utils.write_empty(args.out, "genome_id")
        return
    if master.empty or "genome_id" not in master.columns:
        master.to_csv(args.out, sep="\t", index=False)
        return
    master["genome_id"] = master["genome_id"].astype(str).str.strip()

    if not args.niche or not os.path.exists(args.niche):
        master.to_csv(args.out, sep="\t", index=False)
        utils.log("no niche metadata; bacterial master passed through unchanged.")
        return

    niche = pd.read_csv(args.niche, sep="\t", dtype=str, comment="#")
    gcol = "genome_id" if "genome_id" in niche.columns else niche.columns[0]
    niche[gcol] = niche[gcol].astype(str).str.strip()

    out_cols = {}
    for col in niche.columns:
        if col == gcol:
            continue
        vals = niche[col].fillna("").astype(str).str.strip()
        uniq = set(v for v in vals.unique() if v and v.lower() not in ("nan", "na"))
        is_binary = uniq.issubset({"0", "1", "true", "false", "yes", "no"})
        if is_binary:
            name = f"niche_{_slug(col)}"
            out_cols[name] = dict(zip(niche[gcol], vals.map(
                lambda v: 1 if str(v).lower() in ("1", "true", "yes") else 0)))
        else:
            for level in sorted(uniq):
                name = f"niche_{_slug(col)}_{_slug(level)}"
                out_cols[name] = dict(zip(niche[gcol], (vals == level).astype(int)))

    for name, mapping in out_cols.items():
        master[name] = master["genome_id"].map(mapping).fillna(0).astype(int)

    master.to_csv(args.out, sep="\t", index=False)
    utils.log(f"bacterial niche master: added {len(out_cols)} niche_* column(s) for "
              f"{len(master)} genomes -> {args.out}")


if __name__ == "__main__":
    main()
