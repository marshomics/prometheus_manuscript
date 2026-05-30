"""Collect hypothesis tests from all analysis families, apply per-family BH FDR, tag hypothesis class."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import math
import pandas as pd

HEADER = (
    "family\ttest\tterm\testimate\tci_low\tci_high\t"
    "p_value\tp_adj\thypothesis"
)

CONFIRMATORY_TERMS = [
    "methanogen",
    "relax",
    "transfer",
    "independent_acquisition",
]


def _coerce_p(v):
    """Return a finite float p-value, or None if missing/non-numeric/NaN."""
    if v is None:
        return None
    try:
        p = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(p):
        return None
    return p


def bh_correction(pvals):
    """Benjamini-Hochberg FDR-adjusted p-values, preserving input order.

    Only the entries that carry a usable (finite) p-value are corrected; the
    number of tests m is the count of those entries, NOT len(pvals). Counting
    missing/NaN entries in m would inflate the BH denominator and make every
    adjusted p-value too large. Missing entries map back to None.
    """
    out = [None] * len(pvals)
    valid = [(i, p) for i, p in ((i, _coerce_p(v)) for i, v in enumerate(pvals))
             if p is not None]
    m = len(valid)
    if m == 0:
        return out
    # Step-up over the m valid tests, largest p first.
    valid.sort(key=lambda ip: ip[1])
    prev_adj = 1.0
    for rank, (orig_i, p) in enumerate(reversed(valid), start=1):
        adj = p * m / (m - rank + 1)
        adj = min(adj, prev_adj, 1.0)
        prev_adj = adj
        out[orig_i] = adj
    return out


def is_confirmatory(term: str) -> bool:
    t_lower = str(term).lower()
    return any(ct in t_lower for ct in CONFIRMATORY_TERMS)


def safe_float(val):
    try:
        f = float(val)
        if math.isnan(f):
            return None
        return f
    except (TypeError, ValueError):
        return None


def normalize_association(df: pd.DataFrame):
    """Normalize association_tests.tsv into canonical rows."""
    # SCHEMA: test, model, term, estimate, std_error, statistic, p_value, ci_low, ci_high, n, note
    rows = []
    for _, row in df.iterrows():
        rows.append({
            "family": "association",
            "test": str(row.get("test", row.get("model", ""))).strip(),
            "term": str(row.get("term", "")).strip(),
            "estimate": safe_float(row.get("estimate")),
            "ci_low": safe_float(row.get("ci_low")),
            "ci_high": safe_float(row.get("ci_high")),
            "p_value": safe_float(row.get("p_value")),
        })
    return rows


def normalize_cooccurrence(df: pd.DataFrame):
    """Normalize cooccurrence.tsv into canonical rows."""
    # SCHEMA: system_a, system_b, n_a, n_b, n_both, method, statistic, p_value, p_adj, direction
    rows = []
    for _, row in df.iterrows():
        term = f"{row.get('system_a','')}_vs_{row.get('system_b','')}"
        rows.append({
            "family": "cooccurrence",
            "test": str(row.get("method", "")).strip(),
            "term": term.strip(),
            "estimate": safe_float(row.get("statistic")),
            "ci_low": None,
            "ci_high": None,
            "p_value": safe_float(row.get("p_value")),
        })
    return rows


def normalize_ecology(df: pd.DataFrame):
    """Normalize ecology_associations.tsv into canonical rows."""
    # SCHEMA: predictor, level, model, estimate, std_error, p_value, p_adj, n, note
    rows = []
    for _, row in df.iterrows():
        term_parts = [str(row.get("predictor", "")), str(row.get("level", ""))]
        term = "_".join(p for p in term_parts if p and p.lower() not in ("nan", ""))
        rows.append({
            "family": "ecology",
            "test": str(row.get("model", "")).strip(),
            "term": term.strip(),
            "estimate": safe_float(row.get("estimate")),
            "ci_low": None,
            "ci_high": None,
            "p_value": safe_float(row.get("p_value")),
        })
    return rows


def normalize_selection(df: pd.DataFrame):
    """Normalize selection_summary.tsv into canonical rows."""
    # SCHEMA: method, statistic, value, p_value, n_sites_or_branches, note
    rows = []
    for _, row in df.iterrows():
        rows.append({
            "family": "selection",
            "test": str(row.get("method", "")).strip(),
            "term": str(row.get("statistic", "")).strip(),
            "estimate": safe_float(row.get("value")),
            "ci_low": None,
            "ci_high": None,
            "p_value": safe_float(row.get("p_value")),
        })
    return rows


def load_table(path: str, normalizer):
    """Load and normalize a table; return [] on any error."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return []
        return normalizer(df)
    except Exception as e:
        utils.log(f"WARNING: could not load {path}: {e}")
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--association", required=True)
    ap.add_argument("--cooccurrence", required=True)
    ap.add_argument("--ecology", required=True)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--prevalence", required=True)
    ap.add_argument("--fdr", default="BH")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    all_rows = []
    all_rows.extend(load_table(args.association, normalize_association))
    all_rows.extend(load_table(args.cooccurrence, normalize_cooccurrence))
    all_rows.extend(load_table(args.ecology, normalize_ecology))
    all_rows.extend(load_table(args.selection, normalize_selection))
    # prevalence is not a test table per se; include if needed but schema has no p_value rows

    if not all_rows:
        utils.write_empty(args.out, HEADER)
        return

    df = pd.DataFrame(all_rows)

    # Apply BH FDR within each family
    df["p_adj"] = None
    for family, grp in df.groupby("family"):
        pvals = grp["p_value"].tolist()
        padj = bh_correction(pvals)
        df.loc[grp.index, "p_adj"] = padj

    # hypothesis column
    df["hypothesis"] = df["term"].apply(
        lambda t: "confirmatory" if is_confirmatory(t) else "exploratory"
    )

    # Reorder columns to match schema
    out_cols = HEADER.split("\t")
    for c in out_cols:
        if c not in df.columns:
            df[c] = None
    df = df[out_cols]

    df.to_csv(args.out, sep="\t", index=False)
    utils.log(f"Wrote {len(df)} rows to results registry -> {args.out}")


if __name__ == "__main__":
    main()
