"""Compute Prometheus prevalence by GTDB rank with Wilson 95% CIs."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import math
import pandas as pd

HEADER = (
    "rank\ttaxon\tn_genomes\tn_present\tprevalence\t"
    "ci_low\tci_high\tmedian_completeness\ttier_set"
)

GTDB_RANKS = [
    "gtdb_phylum",
    "gtdb_class",
    "gtdb_order",
    "gtdb_family",
    "gtdb_genus",
]


def wilson_ci(n: int, k: int, z: float = 1.96):
    """Wilson score interval for k successes in n trials."""
    if n == 0:
        return 0.0, 0.0
    p_hat = k / n
    denom = 1 + z * z / n
    centre = (p_hat + z * z / (2 * n)) / denom
    half_width = (z * math.sqrt(p_hat * (1 - p_hat) / n + z * z / (4 * n * n))) / denom
    ci_low = max(0.0, centre - half_width)
    ci_high = min(1.0, centre + half_width)
    return ci_low, ci_high


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--primary-tiers", required=True,
                    help="Comma-separated tier set, e.g. 'A' or 'A,B'")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    tier_set_str = args.primary_tiers
    primary_tiers = set(t.strip() for t in tier_set_str.split(",") if t.strip())

    try:
        master = pd.read_csv(args.master, sep="\t")
    except Exception:
        utils.write_empty(args.out, HEADER)
        return

    if master.empty:
        utils.write_empty(args.out, HEADER)
        return

    # Restrict denominator to pass_primary == True genomes
    if "pass_primary" in master.columns:
        # Handle string 'True'/'False' and actual booleans
        pp = master["pass_primary"].map(
            lambda x: str(x).strip().lower() in ("true", "1", "yes")
        )
        df = master[pp].copy()
    else:
        df = master.copy()

    if df.empty:
        utils.write_empty(args.out, HEADER)
        return

    # Determine prometheus_present column
    if "prometheus_present" not in df.columns:
        utils.write_empty(args.out, HEADER)
        return

    # present = prometheus_present True AND prometheus_tier in primary_tiers
    def is_present(row):
        pp_val = str(row.get("prometheus_present", "")).strip().lower()
        if pp_val not in ("true", "1", "yes"):
            return False
        tier = str(row.get("prometheus_tier", "")).strip()
        if not tier:
            return False
        return tier in primary_tiers

    df["_present"] = df.apply(is_present, axis=1)

    rows = []
    for rank in GTDB_RANKS:
        if rank not in df.columns:
            continue
        for taxon, grp in df.groupby(rank):
            taxon_str = str(taxon).strip()
            if not taxon_str or taxon_str.lower() in ("nan", "na", ""):
                continue
            n_genomes = len(grp)
            n_present = int(grp["_present"].sum())
            prevalence = n_present / n_genomes if n_genomes > 0 else 0.0
            ci_low, ci_high = wilson_ci(n_genomes, n_present)
            comp_col = "completeness"
            if comp_col in grp.columns:
                med_comp = pd.to_numeric(grp[comp_col], errors="coerce").median()
                med_comp = round(float(med_comp), 3) if not math.isnan(med_comp) else ""
            else:
                med_comp = ""
            rows.append({
                "rank": rank,
                "taxon": taxon_str,
                "n_genomes": n_genomes,
                "n_present": n_present,
                "prevalence": round(prevalence, 6),
                "ci_low": round(ci_low, 6),
                "ci_high": round(ci_high, 6),
                "median_completeness": med_comp,
                "tier_set": tier_set_str,
            })

    if not rows:
        utils.write_empty(args.out, HEADER)
        return

    out_df = pd.DataFrame(rows, columns=HEADER.split("\t"))
    out_df.to_csv(args.out, sep="\t", index=False)
    utils.log(f"Wrote {len(rows)} taxon-rank prevalence rows -> {args.out}")


if __name__ == "__main__":
    main()
