"""Prometheus prevalence by GTDB rank with Wilson 95% CIs (design Stage 4.1).

Denominator: genomes with pass_primary == True. Presence follows SCHEMA Contract v2
item 6 for a tier set T: PRESENT if prometheus_tier is in T, ABSENT if prometheus_tier
is "" (no non-rejected locus), EXCLUDED -- dropped from that tier set's denominator and
counted in n_excluded -- when the genome's best validated tier lies outside T. A B-only
genome is therefore never a confident absence in the tier-A prevalence.

Two tier sets are reported side by side (design 4.1 "at tier A and tier A+B"):
  primary      --primary-tiers      -> n_genomes, n_present, prevalence, ci_low, ci_high,
                                       tier_set, n_excluded
  sensitivity  --sensitivity-tiers  -> *_sensitivity columns (default A,B)

Rows: rank "overall" (taxon "all": every pass_primary genome) then gtdb_phylum ..
gtdb_genus. Genomes unassigned at a rank (blank / NA / "unclassified") are left out of
that rank's rows and counted in the log. median_completeness is taken over the
primary-set denominator genomes; n_genomes_total counts every pass_primary genome of
the taxon (present + absent + excluded).

Empty master (zero rows) -> header-only output (exit 0). Missing columns, unparseable
booleans, unknown tiers or prometheus_present/prometheus_tier disagreement -> exit 1.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import math
from statistics import NormalDist

import pandas as pd

BASE_COLS = ["rank", "taxon", "n_genomes", "n_present", "prevalence", "ci_low", "ci_high",
             "median_completeness", "tier_set"]
EXTRA_COLS = ["n_excluded", "n_genomes_total", "tier_set_sensitivity", "n_genomes_sensitivity",
              "n_present_sensitivity", "n_excluded_sensitivity", "prevalence_sensitivity",
              "ci_low_sensitivity", "ci_high_sensitivity"]
COLUMNS = BASE_COLS + EXTRA_COLS
HEADER = "\t".join(COLUMNS)

GTDB_RANKS = ["gtdb_phylum", "gtdb_class", "gtdb_order", "gtdb_family", "gtdb_genus"]
UNASSIGNED = {"", "nan", "na", "n/a", "none", "null", "unclassified"}
VALID_TIERS = {"A", "B", "C"}
Z975 = NormalDist().inv_cdf(0.975)          # 1.959963984540054
TRUE = {"true", "t", "1", "yes", "y", "1.0"}
FALSE = {"false", "f", "0", "no", "n", "0.0"}


def wilson_ci(n: int, k: int, z: float = Z975):
    """Wilson score 95% interval for k successes in n trials (NaN, NaN when n == 0)."""
    if n <= 0:
        return math.nan, math.nan
    p_hat = k / n
    denom = 1 + z * z / n
    centre = (p_hat + z * z / (2 * n)) / denom
    half = z * math.sqrt(p_hat * (1 - p_hat) / n + z * z / (4 * n * n)) / denom
    lo = 0.0 if k == 0 else max(0.0, centre - half)
    hi = 1.0 if k == n else min(1.0, centre + half)
    return lo, hi


def parse_tiers(spec: str, flag: str) -> list:
    tiers = []
    for t in spec.split(","):
        t = t.strip()
        if t and t not in tiers:
            tiers.append(t)
    if not tiers or not set(tiers) <= VALID_TIERS:
        sys.exit(f"[compute_prevalence] {flag} must be a comma-separated subset of A,B,C (got '{spec}')")
    return tiers


def parse_bool(series: pd.Series, col: str, blank=None) -> pd.Series:
    s = series.astype(str).str.strip().str.lower()
    out = pd.Series(pd.NA, index=s.index, dtype="object")
    out[s.isin(TRUE)] = True
    out[s.isin(FALSE)] = False
    if blank is not None:
        out[s == ""] = blank
    bad = out.isna()
    if bad.any():
        sys.exit(f"[compute_prevalence] column '{col}': {int(bad.sum())} non-boolean value(s), "
                 f"e.g. '{s[bad].iloc[0]}'")
    return out.astype(bool)


def tier_response(tier: pd.Series, tiers: list) -> pd.Series:
    """Contract v2 item 6: 1 present, 0 absent, NA excluded."""
    t = tier.astype(str).str.strip()
    bad = ~t.isin(VALID_TIERS | {""})
    if bad.any():
        sys.exit(f"[compute_prevalence] prometheus_tier: unexpected value(s) "
                 f"{sorted(set(t[bad]))} (allowed: '', A, B, C)")
    out = pd.Series(pd.NA, index=t.index, dtype="Int64")
    out[t == ""] = 0
    out[t.isin(tiers)] = 1
    return out


def summarise(grp: pd.DataFrame) -> dict:
    yp, ys = grp["_y_primary"], grp["_y_sens"]
    n_p, k_p = int(yp.notna().sum()), int(yp.fillna(0).sum())
    n_s, k_s = int(ys.notna().sum()), int(ys.fillna(0).sum())
    lo_p, hi_p = wilson_ci(n_p, k_p)
    lo_s, hi_s = wilson_ci(n_s, k_s)
    comp = grp.loc[yp.notna(), "_completeness"].dropna()
    return {
        "n_genomes": n_p, "n_present": k_p,
        "prevalence": round(k_p / n_p, 6) if n_p else math.nan,
        "ci_low": round(lo_p, 6) if n_p else math.nan,
        "ci_high": round(hi_p, 6) if n_p else math.nan,
        "median_completeness": round(float(comp.median()), 3) if len(comp) else math.nan,
        "n_excluded": int(yp.isna().sum()), "n_genomes_total": int(len(grp)),
        "n_genomes_sensitivity": n_s, "n_present_sensitivity": k_s,
        "n_excluded_sensitivity": int(ys.isna().sum()),
        "prevalence_sensitivity": round(k_s / n_s, 6) if n_s else math.nan,
        "ci_low_sensitivity": round(lo_s, 6) if n_s else math.nan,
        "ci_high_sensitivity": round(hi_s, 6) if n_s else math.nan,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--master", required=True)
    ap.add_argument("--primary-tiers", required=True,
                    help="Comma-separated primary tier set, e.g. 'A' or 'A,B'")
    ap.add_argument("--sensitivity-tiers", default="A,B",
                    help="Comma-separated sensitivity tier set (default A,B)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    primary = parse_tiers(args.primary_tiers, "--primary-tiers")
    sens = parse_tiers(args.sensitivity_tiers, "--sensitivity-tiers")
    lab_p, lab_s = ",".join(primary), ",".join(sens)

    master = pd.read_csv(args.master, sep="\t", dtype=str, keep_default_na=False)
    required = ["genome_id", "pass_primary", "prometheus_present", "prometheus_tier", "completeness"]
    missing = [c for c in required if c not in master.columns]
    if missing:
        sys.exit(f"[compute_prevalence] master lacks required column(s): {', '.join(missing)}")
    ranks = [r for r in GTDB_RANKS if r in master.columns]
    if not ranks:
        sys.exit("[compute_prevalence] master has none of the gtdb_phylum..gtdb_genus columns")
    if master.empty:
        utils.write_empty(args.out, HEADER)  # nonfatal: empty input (zero genomes)
        utils.log("master has zero genomes; wrote header-only output")
        return
    if master["genome_id"].duplicated().any():
        sys.exit("[compute_prevalence] duplicated genome_id in master")

    present_flag = parse_bool(master["prometheus_present"], "prometheus_present")
    has_tier = master["prometheus_tier"].astype(str).str.strip() != ""
    incons = present_flag != has_tier
    if incons.any():
        sys.exit(f"[compute_prevalence] {int(incons.sum())} genome(s) where prometheus_present "
                 f"disagrees with prometheus_tier (e.g. {master.loc[incons, 'genome_id'].iloc[0]})")

    pp = parse_bool(master["pass_primary"], "pass_primary", blank=False)
    df = master[pp].copy()
    df["_y_primary"] = tier_response(df["prometheus_tier"], primary)
    df["_y_sens"] = tier_response(df["prometheus_tier"], sens)
    comp = df["completeness"].astype(str).str.strip()
    df["_completeness"] = pd.to_numeric(comp.where(~comp.str.lower().isin(UNASSIGNED)), errors="coerce")
    bad_comp = df["_completeness"].isna() & ~comp.str.lower().isin(UNASSIGNED)
    if bad_comp.any():
        sys.exit(f"[compute_prevalence] completeness: non-numeric value '{comp[bad_comp].iloc[0]}'")
    utils.log(f"{len(master)} genomes; pass_primary={len(df)}; primary tiers {lab_p}: "
              f"present={int(df['_y_primary'].fillna(0).sum())}, "
              f"excluded={int(df['_y_primary'].isna().sum())}; sensitivity tiers {lab_s}: "
              f"present={int(df['_y_sens'].fillna(0).sum())}, excluded={int(df['_y_sens'].isna().sum())}")

    rows = []
    if len(df):
        rows.append({"rank": "overall", "taxon": "all", **summarise(df)})
    for rank in ranks:
        taxa = df[rank].astype(str).str.strip()
        assigned = ~taxa.str.lower().isin(UNASSIGNED)
        if (~assigned).any():
            utils.log(f"{rank}: {int((~assigned).sum())} pass_primary genome(s) unassigned at this "
                      f"rank; left out of its rows")
        sub = df[assigned].assign(_taxon=taxa[assigned])
        for taxon, grp in sub.groupby("_taxon", sort=True):
            rows.append({"rank": rank, "taxon": taxon, **summarise(grp)})

    if not rows:
        utils.write_empty(args.out, HEADER)  # nonfatal: no pass_primary genomes
        utils.log("no pass_primary genomes; wrote header-only output")
        return
    out = pd.DataFrame(rows)
    out["tier_set"] = lab_p
    out["tier_set_sensitivity"] = lab_s
    out[COLUMNS].to_csv(args.out, sep="\t", index=False, na_rep="NA")
    utils.log(f"wrote {len(out)} rank/taxon prevalence rows -> {args.out}")


if __name__ == "__main__":
    main()
