"""results_registry.py -- Stage 7.2: every hypothesis test in one table, one FDR regime.

Collects the tests written by the analysis stages, applies ONE multiple-testing
correction (Benjamini-Hochberg by default) WITHIN each test family, and tags every row
confirmatory or exploratory by an EXACT (family, test, term) key.

Families and their producers (columns per workflow/scripts/SCHEMA.md):
  association     distribution/association_tests.tsv        (association_tests.R)
  cooccurrence    distribution/cooccurrence.tsv             (cooccurrence.R)
  ecology         ecology/ecology_associations.tsv          (ecology_stats.R)
  methanogenesis  ecology/methanogenesis_associations.tsv   (methanogenesis_stats.R; optional
                  input, always exploratory, roles assigned by the producer)
  selection       evolution/selection/selection_summary.tsv (run_selection.py)
  reconciliation  evolution/reconciliation/reconciliation_events.tsv (run_reconciliation.py)
  ancestral       evolution/ancestral/gains_losses.tsv      (ancestral_states.R)
  disorder        disorder/disorder_associations.tsv        (disorder_stats.R; opt-in, always exploratory)

Confirmatory set (DEFAULT_CONFIRMATORY; override with --confirmatory 'k1;k2'). One pre-registered
primary test per hypothesis, matched exactly (case-insensitive equality, no substrings):
  association:phyloglm_confounders:methanogen   methanogen effect on Prometheus presence after
                                                phylogeny + completeness + genome size + defense
                                                load excl. Prometheus (association_tests.R HEADLINE)
  ancestral:methanogen_gain_enrichment:methanogen_independent_acquisition_enrichment
                                                gains on methanogen branches vs the Mk-structured null
  selection:RELAX:K                             selection-intensity shift on methanogen branches
  reconciliation:generax:transfer               horizontal transfer (GeneRax DTL; no p-value)
Everything else is exploratory; its `role` says why: `baseline` (naive_glm / naive_fisher:
assume independent genomes, reported for comparison only), `descriptive` (no p-value:
event/site counts), `secondary` (other tests, e.g. the unadjusted phyloglm, fitPagel, the D
statistic), `sensitivity` (a non-primary tier set: the same test on overlapping data).

The methanogenesis family (design 6.2: pathway capacity within methanogens) is always
exploratory and keeps the producer's role: `primary` = the habitat-adjusted partial effect
of each pathway indicator (test habitat_adjusted, term = predictor), `secondary` = the
marginal model, `descriptive` = the stochastic-mapping change count (test simmap; estimate =
mean number of changes, CI = its 2.5-97.5 % range, scale count) and the corHMM four-model
comparison (test corhmm; estimate = statistic = delta AIC, independent minus correlated,
scale statistic), `sensitivity` = the sens_* re-fits. Its FDR family is the primary rows
only; a p-value withheld by the producer's replication or VIF guard is NA and not counted.

FDR: per family, over the rows that carry a usable p-value and are neither `sensitivity`
nor `baseline` rows (methanogenesis: `primary` rows only); m (`fdr_m`) counts only those
rows, so missing p-values never inflate it. Baselines stay out because their p-values are invalid by design (non-independence) and
BH is rank-based: their very small p-values would LOWER the adjusted p of the real tests.
Methods: BH (default), BY, holm, bonferroni -- identical to R's p.adjust; anything else is
an error.

Failure policy (Contract v2): a missing or unreadable input, a non-empty table without
its SCHEMA columns, a p-value outside [0, 1], or a conflicting duplicate key is an error.
A header-only input is a genuinely empty family. With --strict-confirmatory a confirmatory
key that a non-empty family should contain but does not is an error, because a label
mismatch would otherwise silently tag the headline test exploratory.

Output: synthesis/results_registry.tsv -- the SCHEMA columns
  family, test, term, estimate, ci_low, ci_high, p_value, p_adj, hypothesis
followed by the additive columns
  role, estimate_scale, null_value, ci_method, statistic, std_error, n,
  fdr_method, fdr_m, tier_set, model, note
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import math
import re

import pandas as pd

CORE_COLS = ["family", "test", "term", "estimate", "ci_low", "ci_high",
             "p_value", "p_adj", "hypothesis"]
EXTRA_COLS = ["role", "estimate_scale", "null_value", "ci_method", "statistic",
              "std_error", "n", "fdr_method", "fdr_m", "tier_set", "model", "note"]
OUT_COLS = CORE_COLS + EXTRA_COLS
HEADER = "\t".join(OUT_COLS)

FAMILY_ORDER = ["association", "ancestral", "reconciliation", "selection",
                "cooccurrence", "ecology", "methanogenesis", "disorder"]
ALWAYS_EXPLORATORY = {"disorder",          # sequence-predicted; never confirmatory
                      "methanogenesis"}    # within-methanogen pathway tests (design 6.2)
PRODUCER_ROLE_FAMILIES = {"methanogenesis"}   # role comes from the producer's `role` column
FDR_ROLES_ONLY = {"methanogenesis": {"primary"}}   # FDR over these roles only
BASELINE_TESTS = {"naive_glm", "naive_fisher"}   # association_tests.R independence baselines
NO_FDR_ROLES = {"sensitivity", "baseline"}

DEFAULT_CONFIRMATORY = (
    "association:phyloglm_confounders:methanogen",
    "ancestral:methanogen_gain_enrichment:methanogen_independent_acquisition_enrichment",
    "selection:RELAX:K",
    "reconciliation:generax:transfer",
)
# A non-empty association / ancestral table always carries its confirmatory key (the
# producers write every test row, with a 'failed:' note when a fit fails), so a missing
# key there is a label mismatch. Selection and reconciliation keys are required only when
# the producer ran that method/tool (RELAX is configurable; ALE is optional).
ALWAYS_REQUIRED_FAMILIES = {"association", "ancestral"}

Z975 = 1.959963984540054                   # qnorm(0.975)
P_TOL = 1e-9                               # float slack on [0, 1] p-values
NA_STRINGS = {"", "na", "nan", "none", "null", "<na>"}


# --------------------------------------------------------------------------- #
# value helpers
# --------------------------------------------------------------------------- #
def safe_float(val):
    """Finite or infinite float, or None for missing / non-numeric / NaN."""
    if val is None:
        return None
    if isinstance(val, str):
        s = val.strip()
        if s.lower() in NA_STRINGS:
            return None
        val = s
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    if math.isnan(f):
        return None
    return f


def safe_str(val) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    return "" if s.lower() in NA_STRINGS else s


def check_p(p, where):
    """Validate a p-value: None stays None; values outside [0, 1] are an error."""
    if p is None:
        return None
    if p < -P_TOL or p > 1 + P_TOL or math.isinf(p):
        raise ValueError(f"p-value {p!r} outside [0, 1] in {where}")
    return min(max(p, 0.0), 1.0)


def wald_ci(est, se):
    if est is None or se is None or math.isinf(est) or math.isinf(se) or se < 0:
        return None, None
    return est - Z975 * se, est + Z975 * se


def canonical_tiers(s) -> str:
    """'B,A' / 'A+B' / 'A;B' -> 'A,B'. '' when no tier set is given."""
    toks = [t for t in re.split(r"[,+;/\s]+", safe_str(s)) if t]
    return ",".join(sorted(dict.fromkeys(toks)))


def parse_key(s: str):
    parts = [p.strip() for p in str(s).split(":", 2)]
    if len(parts) != 3 or not all(parts):
        raise ValueError(f"confirmatory key {s!r} is not 'family:test:term'")
    return tuple(parts)


def key_of(row):
    return (row["family"].casefold(), row["test"].casefold(), row["term"].casefold())


# --------------------------------------------------------------------------- #
# multiple-testing correction (matches R's stats::p.adjust)
# --------------------------------------------------------------------------- #
FDR_ALIASES = {"bh": "BH", "fdr": "BH", "fdr_bh": "BH", "benjamini-hochberg": "BH",
               "by": "BY", "fdr_by": "BY", "holm": "holm", "bonferroni": "bonferroni"}


def normalize_method(method: str) -> str:
    m = FDR_ALIASES.get(str(method).strip().lower())
    if m is None:
        raise ValueError(f"unsupported --fdr method {method!r} "
                         f"(supported: BH, BY, holm, bonferroni)")
    return m


def p_adjust(pvals, method="BH"):
    """Adjusted p-values in input order. Entries that are None/NaN stay None and are
    excluded from m (R's p.adjust: n defaults to the number of non-NA p-values)."""
    method = normalize_method(method)
    out = [None] * len(pvals)
    valid = [(i, p) for i, p in enumerate(safe_float(v) for v in pvals) if p is not None]
    m = len(valid)
    if m == 0:
        return out
    if m == 1:                            # p.adjust returns a single p unchanged
        out[valid[0][0]] = valid[0][1]
        return out
    if method == "bonferroni":
        for i, p in valid:
            out[i] = min(1.0, m * p)
        return out
    if method == "holm":
        asc = sorted(valid, key=lambda ip: ip[1])
        run = 0.0
        for rank, (i, p) in enumerate(asc, start=1):
            run = max(run, (m + 1 - rank) * p)
            out[i] = min(1.0, run)
        return out
    q = sum(1.0 / k for k in range(1, m + 1)) if method == "BY" else 1.0
    desc = sorted(valid, key=lambda ip: ip[1], reverse=True)
    run = math.inf
    for j, (i, p) in enumerate(desc):
        rank = m - j                      # ascending rank of this p
        run = min(run, q * m / rank * p)
        out[i] = min(1.0, run)
    return out


def bh_correction(pvals):
    """Backward-compatible alias: Benjamini-Hochberg with NA excluded from m."""
    return p_adjust(pvals, "BH")


# --------------------------------------------------------------------------- #
# input loading (strict)
# --------------------------------------------------------------------------- #
def load_table(path: str, name: str, required_cols, optional=False):
    """Read a producer TSV as strings. Returns a DataFrame (possibly empty) or None when an
    optional input is not configured. Missing/unreadable files and non-empty tables
    lacking their SCHEMA columns are errors."""
    if not path:
        if optional:
            return None
        raise ValueError(f"--{name}: no path given")
    if not os.path.exists(path):
        raise FileNotFoundError(f"--{name}: {path} does not exist")
    if os.path.getsize(path) == 0:
        raise ValueError(f"--{name}: {path} is a 0-byte file (producers always write a "
                         f"header; a 0-byte table means the producer crashed)")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df.empty:
        utils.log(f"{name}: {path} is header-only (no tests in this family)")
        return df
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"--{name}: {path} lacks SCHEMA column(s) {missing}; "
                         f"found {list(df.columns)}")
    return df


def base_row(family, test, term, **kw):
    row = {c: None for c in OUT_COLS}
    row.update(family=family, test=safe_str(test), term=safe_str(term),
               tier_set="", model="", note="", ci_method="none",
               estimate_scale="none")
    row.update(kw)
    return row


def _ci(row, est, lo, hi, se):
    """Producer CI when given; else a Wald (normal) CI from the standard error."""
    if lo is not None and hi is not None:
        row.update(ci_low=lo, ci_high=hi, ci_method="producer")
    elif est is not None and se is not None:
        wl, wh = wald_ci(est, se)
        if wl is not None:
            row.update(ci_low=wl, ci_high=wh, ci_method="wald_z")


# --------------------------------------------------------------------------- #
# per-family normalisers
# --------------------------------------------------------------------------- #
ASSOC_COLS = ["test", "term", "estimate", "p_value"]


def normalize_association(df):
    rows = []
    for i, r in df.iterrows():
        test = safe_str(r.get("test")) or safe_str(r.get("model"))
        est, se = safe_float(r.get("estimate")), safe_float(r.get("std_error"))
        if test.casefold() == "d_statistic":
            scale, null = "D_statistic", 1.0       # p tests D vs random (D=1); 0 = Brownian
        elif test.casefold() == "fitpagel":
            scale, null = "none", None             # LR test; statistic = likelihood ratio
        else:
            scale, null = "log_odds", 0.0          # glm / phyloglm logistic coefficients
        row = base_row("association", test, r.get("term"),
                       estimate=est, std_error=se,
                       statistic=safe_float(r.get("statistic")),
                       p_value=check_p(safe_float(r.get("p_value")),
                                       f"association_tests row {i + 1}"),
                       n=safe_float(r.get("n")), model=safe_str(r.get("model")),
                       note=safe_str(r.get("note")), estimate_scale=scale,
                       null_value=null, tier_set=canonical_tiers(r.get("tier_set")))
        _ci(row, est, safe_float(r.get("ci_low")), safe_float(r.get("ci_high")),
            se if scale == "log_odds" else None)
        rows.append(row)
    return rows


COOC_COLS = ["system_a", "system_b", "method", "statistic", "p_value"]
_OR_RE = re.compile(r"OR\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|Inf|-Inf|NaN|NA)")


def _odds_ratio(r):
    if "odds_ratio" in r.index:
        return safe_float(r.get("odds_ratio"))
    m = _OR_RE.search(safe_str(r.get("direction")))
    return safe_float(m.group(1)) if m else None


def normalize_cooccurrence(df):
    rows = []
    for i, r in df.iterrows():
        a, b = safe_str(r.get("system_a")), safe_str(r.get("system_b"))
        counts = ", ".join(f"{c}={safe_str(r.get(c))}" for c in ("n_a", "n_b", "n_both")
                           if c in r.index)
        note = "; ".join(x for x in (safe_str(r.get("direction")), counts) if x)
        rows.append(base_row(
            "cooccurrence", r.get("method"), f"{a}_vs_{b}",
            estimate=_odds_ratio(r), estimate_scale="odds_ratio", null_value=1.0,
            statistic=safe_float(r.get("statistic")),
            p_value=check_p(safe_float(r.get("p_value")), f"cooccurrence row {i + 1}"),
            note=note, tier_set=canonical_tiers(r.get("tier_set"))))
    return rows


PRED_COLS = ["predictor", "level", "model", "estimate", "p_value"]


def _pred_term(r):
    parts = [safe_str(r.get("predictor")), safe_str(r.get("level"))]
    return "_".join(p for p in parts if p)


def normalize_ecology(df):
    rows = []
    for i, r in df.iterrows():
        est, se = safe_float(r.get("estimate")), safe_float(r.get("std_error"))
        row = base_row("ecology", r.get("model"), _pred_term(r),
                       estimate=est, std_error=se, estimate_scale="log_odds",
                       null_value=0.0,
                       p_value=check_p(safe_float(r.get("p_value")), f"ecology row {i + 1}"),
                       n=safe_float(r.get("n")), model=safe_str(r.get("model")),
                       note=safe_str(r.get("note")),
                       tier_set=canonical_tiers(r.get("tier_set")))
        _ci(row, est, None, None, se)       # phyloglm: Wald z, consistent with its p
        rows.append(row)
    return rows


MGEN_COLS = ["predictor", "analysis", "role", "estimate", "p_value"]
MGEN_ROLES = {"primary", "secondary", "descriptive", "sensitivity"}
MGEN_ROLE_OF = {"marginal": "secondary", "habitat_adjusted": "primary", "simmap": "descriptive",
                "corhmm": "descriptive"}


def _mgen_bool(v):
    s = safe_str(v).casefold()
    return {"true": True, "false": False}.get(s)


def normalize_methanogenesis(df):
    """methanogenesis_stats.R: phyloglm rows (log-odds), the simmap change count (count)
    and the corHMM delta AIC (statistic). The producer's role is kept and checked."""
    rows = []
    for i, r in df.iterrows():
        where = f"methanogenesis row {i + 1}"
        analysis, pred = safe_str(r.get("analysis")), safe_str(r.get("predictor"))
        role = safe_str(r.get("role")).casefold()
        if role not in MGEN_ROLES:
            raise ValueError(f"{where}: role {r.get('role')!r} not in {sorted(MGEN_ROLES)}")
        # sens_* rows are sensitivity re-fits (new ones may be added by the producer);
        # any other unknown analysis is an error rather than a guessed normalisation
        want = MGEN_ROLE_OF.get(analysis, "sensitivity" if analysis.startswith("sens_") else None)
        if want is None:
            raise ValueError(f"{where}: unknown analysis {analysis!r} (known: "
                             f"{sorted(MGEN_ROLE_OF)} and sens_*)")
        if role != want:
            raise ValueError(f"{where}: analysis {analysis!r} must have role {want!r}, got {role!r}")
        p = check_p(safe_float(r.get("p_value")), where)
        note = safe_str(r.get("note"))
        tier = canonical_tiers(r.get("tier_set"))
        n = safe_float(r.get("n"))
        model = safe_str(r.get("model"))
        if analysis == "simmap":
            if p is not None:
                raise ValueError(f"{where}: the simmap row is descriptive and carries no p-value")
            est = safe_float(r.get("n_changes_mean"))
            rep = _mgen_bool(r.get("replicated"))
            extra = (f"replicated={'TRUE' if rep else 'FALSE'}" if rep is not None
                     else "replication not assessed")
            row = base_row("methanogenesis", analysis, pred, estimate=est,
                           estimate_scale="count", n=n, model=model, role=role,
                           note="; ".join(x for x in (extra, note) if x), tier_set=tier)
            lo, hi = safe_float(r.get("n_changes_lo")), safe_float(r.get("n_changes_hi"))
            if est is not None and lo is not None and hi is not None:
                row.update(ci_low=lo, ci_high=hi, ci_method="producer")
        elif analysis == "corhmm":
            if p is not None:
                raise ValueError(f"{where}: the corhmm row is descriptive and carries no p-value")
            d_aic = safe_float(r.get("delta_aic"))
            best = safe_str(r.get("best_model"))
            row = base_row("methanogenesis", analysis, pred, estimate=d_aic, statistic=d_aic,
                           estimate_scale="statistic", n=n, model=model, role=role,
                           note="; ".join(x for x in ((f"best_model={best}" if best else ""),
                                                      note) if x), tier_set=tier)
        else:
            est, se = safe_float(r.get("estimate")), safe_float(r.get("std_error"))
            row = base_row("methanogenesis", analysis, pred, estimate=est, std_error=se,
                           estimate_scale="log_odds", null_value=0.0, p_value=p, n=n,
                           model=model, role=role, note=note, tier_set=tier)
            _ci(row, est, safe_float(r.get("ci_low")), safe_float(r.get("ci_high")), se)
        rows.append(row)
    return rows


def normalize_disorder(df):
    """phylolm (Gaussian, t-based p): no normal-theory CI is synthesised here."""
    rows = []
    for i, r in df.iterrows():
        rows.append(base_row(
            "disorder", r.get("model"), _pred_term(r),
            estimate=safe_float(r.get("estimate")), std_error=safe_float(r.get("std_error")),
            estimate_scale="coefficient", null_value=0.0,
            p_value=check_p(safe_float(r.get("p_value")), f"disorder row {i + 1}"),
            n=safe_float(r.get("n")), model=safe_str(r.get("model")),
            note=safe_str(r.get("note")), tier_set=canonical_tiers(r.get("tier_set"))))
    return rows


SEL_COLS = ["method", "statistic", "value", "p_value"]
COUNT_STATS = {"sites_positive", "sites_episodic", "branches_selected"}


def _is_p_label(stat: str) -> bool:
    s = stat.casefold()
    return s in ("p", "p_value", "pvalue", "p-value") or s.endswith("_p")


def normalize_selection(df):
    rows = []
    for i, r in df.iterrows():
        method, stat = safe_str(r.get("method")), safe_str(r.get("statistic"))
        val = safe_float(r.get("value"))
        if method.casefold() == "relax" and stat.casefold() == "k":
            scale, null = "ratio", 1.0             # K>1 intensified, K<1 relaxed
        elif _is_p_label(stat):
            scale, null, val = "none", None, None  # 'value' repeats the p (BUSTED)
        elif stat in COUNT_STATS:
            scale, null = "count", None
        else:
            scale, null = "statistic", None
        rows.append(base_row(
            "selection", method, stat, estimate=val, estimate_scale=scale, null_value=null,
            p_value=check_p(safe_float(r.get("p_value")), f"selection row {i + 1}"),
            n=safe_float(r.get("n_sites_or_branches")), note=safe_str(r.get("note")),
            tier_set=canonical_tiers(r.get("tier_set"))))
    return rows


RECON_COLS = ["tool", "event", "count"]
# Notes run_reconciliation.py writes on count-0 rows when a tool did NOT produce a
# reconciliation; such a 0 is "not estimated", never "zero transfers".
_RECON_FAIL_RE = re.compile(r"fail|not on path|skip|disabled|exception|no parseable|"
                            r"unparsed|trivial|too few|only \d+ mapped|unavailable",
                            re.IGNORECASE)


def normalize_reconciliation(df):
    """One row per tool: the Transfer count (no p-value; DTL point estimate)."""
    rows = []
    for i, r in df.iterrows():
        if safe_str(r.get("event")).casefold() != "transfer":
            continue
        count, note = safe_float(r.get("count")), safe_str(r.get("note"))
        frac = safe_float(r.get("fraction"))
        if (count is None or count == 0) and _RECON_FAIL_RE.search(note):
            count, note = None, f"not estimated: {note}"
        elif frac is not None:
            note = "; ".join(x for x in (f"vertical_fraction=(S+SL)/(S+SL+T+TL)={frac:g}", note) if x)
        rows.append(base_row("reconciliation", r.get("tool"), "transfer",
                             estimate=count, estimate_scale="count", note=note,
                             tier_set=canonical_tiers(r.get("tier_set"))))
    return rows


ANC_COLS = ["metric", "value"]


def normalize_ancestral(df):
    """metric/value table -> the enrichment test + the acquisition count, per tier set."""
    groups = {}                           # tier_set -> (values, notes)
    for _, r in df.iterrows():
        tier = canonical_tiers(r.get("tier_set"))
        vals, notes = groups.setdefault(tier, ({}, {}))
        k = safe_str(r.get("metric"))
        if k in vals:
            raise ValueError(f"gains_losses.tsv: duplicated metric {k!r} "
                             f"(tier_set={tier or '-'})")
        vals[k], notes[k] = r.get("value"), safe_str(r.get("note"))
    rows = []
    for tier, (vals, notes) in groups.items():
        if "enrichment_p" in vals:
            obs = safe_float(vals.get("gains_on_methanogen_branches"))
            exp = safe_float(vals.get("expected_gains_methanogen"))
            fold = (f"; fold_enrichment=obs/null={obs / exp:.4g}"
                    if obs is not None and exp not in (None, 0.0) else "")
            rows.append(base_row(
                "ancestral", "methanogen_gain_enrichment",
                "methanogen_independent_acquisition_enrichment",
                estimate=obs, null_value=exp, estimate_scale="count",
                p_value=check_p(safe_float(vals.get("enrichment_p")),
                                "gains_losses enrichment_p"),
                note=(notes.get("enrichment_p", "") + fold).strip("; "), tier_set=tier))
        if "n_independent_acquisitions" in vals:
            rows.append(base_row(
                "ancestral", "independent_acquisitions", "independent_acquisition_count",
                estimate=safe_float(vals.get("n_independent_acquisitions")),
                estimate_scale="count", note=notes.get("n_independent_acquisitions", ""),
                tier_set=tier))
    return rows


FAMILIES = {
    # family: (cli name, required columns, normaliser, optional input)
    "association": ("association", ASSOC_COLS, normalize_association, False),
    "cooccurrence": ("cooccurrence", COOC_COLS, normalize_cooccurrence, False),
    "ecology": ("ecology", PRED_COLS, normalize_ecology, False),
    "methanogenesis": ("methanogenesis", MGEN_COLS, normalize_methanogenesis, True),
    "selection": ("selection", SEL_COLS, normalize_selection, False),
    "reconciliation": ("reconciliation", RECON_COLS, normalize_reconciliation, False),
    "ancestral": ("ancestral", ANC_COLS, normalize_ancestral, False),
    "disorder": ("disorder", PRED_COLS, normalize_disorder, True),
}


# --------------------------------------------------------------------------- #
# assembly: tier sets, duplicates, tagging, FDR, confirmatory check
# --------------------------------------------------------------------------- #
def mark_tier_sets(rows, family, primary):
    sets = sorted({r["tier_set"] for r in rows if r["tier_set"]})
    if not sets:
        return
    if not primary:
        if family in PRODUCER_ROLE_FAMILIES:
            return                        # roles (incl. the tier sensitivity) are explicit
        if len(sets) > 1:
            raise ValueError(f"{family}: rows for several tier sets {sets}; pass "
                             f"--primary-tiers so only the primary set enters the FDR")
        return
    demoted = 0
    for r in rows:
        if r["tier_set"] and r["tier_set"] != primary:
            demoted += r.get("role") not in (None, "sensitivity")
            r["role"] = "sensitivity"
    if demoted and family in PRODUCER_ROLE_FAMILIES:
        utils.log(f"WARNING: {family}: {demoted} row(s) computed for a tier set other than "
                  f"the registry primary ({primary}) were demoted to sensitivity")


def drop_duplicates(rows, family):
    seen, out = {}, []
    for r in rows:
        k = key_of(r) + (r["tier_set"],)
        if k in seen:
            prev = seen[k]
            same = all(prev[c] == r[c] for c in ("estimate", "p_value", "ci_low", "ci_high"))
            if not same:
                raise ValueError(f"{family}: conflicting duplicate rows for key "
                                 f"{r['family']}:{r['test']}:{r['term']} "
                                 f"(tier_set={r['tier_set'] or '-'})")
            utils.log(f"WARNING: {family}: identical duplicate row dropped "
                      f"({r['test']}:{r['term']})")
            continue
        seen[k] = r
        out.append(r)
    return out


def tag(rows, confirmatory):
    for r in rows:
        if r.get("role") == "sensitivity" or r["family"] in PRODUCER_ROLE_FAMILIES:
            r["hypothesis"] = "exploratory"   # producer-role families keep their role
            continue
        if key_of(r) in confirmatory and r["family"] not in ALWAYS_EXPLORATORY:
            r["hypothesis"], r["role"] = "confirmatory", "primary"
        else:
            r["hypothesis"] = "exploratory"
            t = r["test"].casefold()
            if r["family"] == "association" and (t in BASELINE_TESTS or t.startswith("naive_")):
                r["role"] = "baseline"
            elif r["p_value"] is None:
                r["role"] = "descriptive"
            else:
                r["role"] = "secondary"


def in_fdr_family(r):
    only = FDR_ROLES_ONLY.get(r["family"])
    if only is not None:
        return r["role"] in only
    return r["role"] not in NO_FDR_ROLES


def apply_fdr(rows, method):
    fam = [r for r in rows if in_fdr_family(r) and r["p_value"] is not None]
    adj = p_adjust([r["p_value"] for r in fam], method)
    for r, a in zip(fam, adj):
        r["p_adj"], r["fdr_m"], r["fdr_method"] = a, len(fam), method


def key_required(key, fam_rows, raw_nonempty):
    """Should a non-empty family contain this confirmatory key?"""
    fam, test, _term = key
    if fam in ALWAYS_REQUIRED_FAMILIES:
        return raw_nonempty               # judged on the producer table, not on the rows
    if not fam_rows:                      # e.g. reconciliation without a Transfer row
        return False
    return any(r["test"].casefold() == test for r in fam_rows)


def check_confirmatory(by_family, confirmatory, strict, raw_nonempty=None):
    raw_nonempty = raw_nonempty or {}
    present = {key_of(r) for rows in by_family.values() for r in rows
               if r.get("hypothesis") == "confirmatory"}
    problems = []
    for key in confirmatory:
        if key[0] in ALWAYS_EXPLORATORY:
            utils.log(f"WARNING: {':'.join(key)} ignored: family '{key[0]}' is always "
                      f"exploratory")
            continue
        if key in present:
            continue
        fam_rows = by_family.get(key[0], [])
        labels = sorted({f"{r['test']}:{r['term']}" for r in fam_rows})
        tiers = sorted({r["tier_set"] for r in fam_rows if r["tier_set"]})
        msg = (f"confirmatory key {':'.join(key)} not found; family '{key[0]}' has "
               f"{len(fam_rows)} row(s) with test:term {labels[:12]}" +
               (f"; tier sets {tiers} (only the primary tier set is confirmatory)"
                if tiers else ""))
        if key_required(key, fam_rows, raw_nonempty.get(key[0], bool(fam_rows))):
            problems.append(msg)
            utils.log(("ERROR: " if strict else "WARNING: ") + msg)
        else:
            utils.log(f"note: {msg} (family empty or method/tool not run)")
    n_conf = sum(1 for rows in by_family.values() for r in rows
                 if r.get("hypothesis") == "confirmatory")
    n_all = sum(len(v) for v in by_family.values())
    if n_all and not n_conf:
        utils.log("WARNING: the registry has tests but NO confirmatory row; check the "
                  "--confirmatory keys against the producers' test/term labels")
    if strict and problems:
        raise SystemExit("results_registry: confirmatory label mismatch:\n  " +
                         "\n  ".join(problems))
    return n_conf


def build_registry(paths, method="BH", primary_tiers="", confirmatory=None,
                   strict=False):
    """paths: {family: path}. Returns the registry as a DataFrame (OUT_COLS order)."""
    method = normalize_method(method)
    primary = canonical_tiers(primary_tiers)
    conf = {tuple(p.casefold() for p in parse_key(k))
            for k in (confirmatory or DEFAULT_CONFIRMATORY)}
    by_family, raw_nonempty = {}, {}
    for fam in FAMILY_ORDER:
        cli, cols, norm, optional = FAMILIES[fam]
        df = load_table(paths.get(fam, ""), cli, cols, optional=optional)
        raw_nonempty[fam] = df is not None and not df.empty
        rows = norm(df) if raw_nonempty[fam] else []
        mark_tier_sets(rows, fam, primary)
        rows = drop_duplicates(rows, fam)
        tag(rows, conf)
        apply_fdr(rows, method)
        by_family[fam] = rows
    n_conf = check_confirmatory(by_family, conf, strict, raw_nonempty)
    all_rows = [r for fam in FAMILY_ORDER for r in by_family[fam]]
    utils.log(f"{len(all_rows)} test row(s), {n_conf} confirmatory; FDR={method} within "
              f"family; per-family m: " +
              ", ".join(f"{f}={sum(1 for r in by_family[f] if r['fdr_m'])}"
                        for f in FAMILY_ORDER))
    reg = pd.DataFrame(all_rows, columns=OUT_COLS)
    for c in ("n", "fdr_m"):              # counts print as integers, not '25.0'
        v = pd.to_numeric(reg[c], errors="coerce")
        if v.notna().any() and (v.dropna() % 1 == 0).all():
            reg[c] = v.astype("Int64")
    return reg


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--association", required=True)
    ap.add_argument("--cooccurrence", required=True)
    ap.add_argument("--ecology", required=True)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--reconciliation", required=True)
    ap.add_argument("--ancestral", required=True)
    ap.add_argument("--disorder", default="", help="optional disorder_associations.tsv")
    ap.add_argument("--methanogenesis", default="",
                    help="optional methanogenesis_associations.tsv (methanogenesis_stats.R)")
    ap.add_argument("--fdr", default="BH", help="BH (default), BY, holm or bonferroni")
    ap.add_argument("--primary-tiers", default="",
                    help="primary tier set (e.g. 'A'); rows of other tier sets become "
                         "sensitivity rows outside the FDR family")
    ap.add_argument("--confirmatory", default="",
                    help="SEMICOLON-separated family:test:term keys replacing the built-in "
                         "set (test labels such as 'phylolm(disorder~pred, lambda)' contain "
                         "commas)")
    ap.add_argument("--strict-confirmatory", action="store_true",
                    help="error (not warn) when a non-empty family lacks its "
                         "confirmatory key")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    conf = [k for k in (s.strip() for s in args.confirmatory.split(";")) if k] or None
    paths = {f: getattr(args, FAMILIES[f][0]) for f in FAMILY_ORDER}
    reg = build_registry(paths, method=args.fdr, primary_tiers=args.primary_tiers,
                         confirmatory=conf, strict=args.strict_confirmatory)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    if reg.empty:
        utils.write_empty(args.out, HEADER)  # every family genuinely empty
    else:
        reg.to_csv(args.out, sep="\t", index=False, na_rep="")
    utils.log(f"wrote {len(reg)} row(s) -> {args.out}")


if __name__ == "__main__":
    main()
