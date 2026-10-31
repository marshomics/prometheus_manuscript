"""Predict intrinsic disorder on validated Prometheus/ProA proteins by CONSENSUS.

A residue is called disordered ONLY where two predictors agree:
  * metapredict (pinned network version, default V3; deep-learning consensus disorder).
    Disordered = score >= --meta-threshold (metapredict's own convention; 0.5 is the
    documented V2/V3 threshold, metapredict/parameters.py METAPREDICT_V3_THRESHOLD).
  * the SECOND predictor (--second-predictor):
      iupred3      IUPred3 (energy estimation; methodologically independent), staged by
                   path (config disorder.iupred3_path; free academic licence, not
                   installable). Disordered = score > --iupred3-threshold (IUPred's
                   "score above 0.5" convention).
      alphapredict predicted AlphaFold2 pLDDT from the alphaPredict network, run through
                   metapredict's bundled copy (metapredict.predict_pLDDT(pLDDT_version='V1')
                   loads alpha_fold_networkV7_hs100_nL2_200e_all_prot.pt, byte-identical
                   to alphaPredict 1.2's network). Disordered = pLDDT < --plddt-cutoff.
                   This is a CORRELATED fallback, not an independent consensus: the
                   metapredict V2/V3 networks were themselves trained on legacy
                   metapredict scores combined with (predicted / AF2) pLDDT
                   (metapredict/backend/network_parameters.py). It is labelled as such.
      none         metapredict only (explicit opt-in; labelled metapredict_only).
      auto         (default) iupred3 when --iupred3-path is set, else alphapredict.
A predictor that is requested but missing or crashing makes this script FAIL (no
silent downgrade to a weaker consensus, no silently dropped loci).

Defensibility controls:
  * intact-only (default) -- truncated/fragmented ORFs show artefactual disorder at the
    broken ends, so only loci with catalog catalytic_status == 'intact' AND
    fragmented == False are scored (anything else, including a missing value, is
    excluded). The catalog's confidence_tier is carried per locus so the downstream
    statistics/figures apply the Contract-v2 tier rule (item 6).
  * low-complexity awareness -- segmasker (BLAST+, SEG window 12 / locut 2.2 / hicut 2.5)
    masks LCRs; frac_disordered_nonlcr is the consensus disorder fraction of the
    NON-LCR residues (LCR residues are excluded from numerator and denominator).
  * host-proteome background (optional) -- a genome_id + background_frac[_disordered]
    table; disorder_minus_background = frac_disordered - background, and the background
    itself is carried (host_background) so the statistics can adjust for it.
  * residues: Prodigal's trailing '*' is removed (MAFFT drops it too, so positions stay
    identical to the family alignment); B/U/X/Z are converted with protfasta's
    STANDARD_CONVERSION (B->N, U->C, X->G, Z->Q -- what metapredict's own FASTA readers
    do; the metapredict encoder raises on any non-canonical letter). Any other
    character fails the run. n_nonstandard records the per-locus count.

Outputs:
  disorder_per_locus.tsv    one row per scored locus (the comparative-trait table)
  disorder_per_residue.tsv  per-residue consensus track (alignment mapping + figure)
Genuinely empty input (no proteins, or no locus passing the intact filter) gives
header-only outputs; every other problem is an error.
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

# SCHEMA columns first (order unchanged), then additive columns (Contract v2 item 5).
LOCUS_COLUMNS = [
    "locus_id", "genome_id", "length_aa", "intact", "predictors",
    "frac_disordered", "frac_disordered_metapredict", "frac_disordered_second",
    "n_idrs", "longest_idr", "frac_in_idrs", "frac_lcr", "frac_disordered_nonlcr",
    "disorder_minus_background",
    "confidence_tier", "n_nonstandard", "host_background", "predictor_detail",
]
LOCUS_HEADER = "\t".join(LOCUS_COLUMNS)
RESIDUE_COLUMNS = ["locus_id", "position", "meta_score", "second_score", "consensus",
                   "lcr", "aa"]
RESIDUE_HEADER = "\t".join(RESIDUE_COLUMNS)

STANDARD_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
# protfasta 0.1.26 _configs.STANDARD_CONVERSION (minus the deletions of '*', '-', ' ',
# which would shift residue positions): used by metapredict's FASTA readers.
CONVERSION = {"B": "N", "U": "C", "X": "G", "Z": "Q"}
SECOND_LABEL = {"iupred3": "metapredict+iupred3", "alphapredict": "metapredict+alphapredict",
                "none": "metapredict_only"}
PLDDT_NETWORK = "V1"      # metapredict's copy of the alphaPredict (V7) network
TRUE_SET = ("true", "1", "yes", "t")
FALSE_SET = ("false", "0", "no", "f")


class DisorderError(RuntimeError):
    """Missing tool, crashed predictor or inconsistent input: the rule must fail."""


def as_bool(v) -> bool:
    s = str(v).strip().lower()
    if s in TRUE_SET:
        return True
    if s in FALSE_SET:
        return False
    raise DisorderError(f"not a boolean flag value: {v!r}")


# =========================================================================
# Pure logic (unit-tested without the predictors installed)
# =========================================================================
def sanitize_sequence(locus: str, seq: str):
    """-> (clean_sequence, n_converted). Strips trailing '*', converts B/U/X/Z,
    raises DisorderError on any other non-canonical character."""
    s = str(seq).strip().upper().rstrip("*")
    bad = sorted({c for c in s if c not in STANDARD_AA and c not in CONVERSION})
    if bad:
        raise DisorderError(f"{locus}: unsupported residue character(s) {bad} "
                            "(only the 20 canonical amino acids, B/U/X/Z and a trailing "
                            "'*' are accepted)")
    if not s:
        raise DisorderError(f"{locus}: empty sequence")
    n_conv = sum(1 for c in s if c in CONVERSION)
    return "".join(CONVERSION.get(c, c) for c in s), n_conv


def _check_track(name, scores, n, lo, hi, locus):
    if scores is None:
        raise DisorderError(f"{locus}: {name} returned no scores")
    vals = [float(x) for x in scores]
    if len(vals) != n:
        raise DisorderError(f"{locus}: {name} returned {len(vals)} scores for a "
                            f"{n}-residue sequence")
    for v in vals:
        if not math.isfinite(v) or v < lo - 1e-6 or v > hi + 1e-6:
            raise DisorderError(f"{locus}: {name} score {v} outside [{lo}, {hi}]")
    return vals


def binarize(scores, kind: str, threshold: float):
    """Per-residue 0/1 call with each tool's documented convention:
    metapredict: score >= thr; iupred3: score > thr; plddt: pLDDT < thr."""
    if kind == "metapredict":
        return [1 if s >= threshold else 0 for s in scores]
    if kind == "iupred3":
        return [1 if s > threshold else 0 for s in scores]
    if kind == "plddt":
        return [1 if s < threshold else 0 for s in scores]
    raise ValueError(f"unknown predictor kind {kind!r}")


def consensus_binary(meta_scores, second_scores, meta_thr, second_thr,
                     second_is_plddt, plddt_cut, locus="sequence", n=None):
    """Per-residue consensus: 1 iff metapredict AND the second predictor call disorder.

    second_scores None = single-predictor mode (consensus == metapredict). Every track
    must have exactly n (= sequence length; default len(meta_scores)) finite values in
    range (else DisorderError). Returns (consensus, meta_bin, second_bin) with
    second_bin = [None]*n in single mode.
    """
    n = len(meta_scores) if n is None else n
    meta = _check_track("metapredict", meta_scores, n, 0.0, 1.0, locus)
    meta_bin = binarize(meta, "metapredict", meta_thr)
    if second_scores is None:
        return list(meta_bin), meta_bin, [None] * n
    if second_is_plddt:
        sec = _check_track("pLDDT", second_scores, n, 0.0, 100.0, locus)
        second_bin = binarize(sec, "plddt", plddt_cut)
    else:
        sec = _check_track("IUPred3", second_scores, n, 0.0, 1.0, locus)
        second_bin = binarize(sec, "iupred3", second_thr)
    cons = [1 if (a and b) else 0 for a, b in zip(meta_bin, second_bin)]
    return cons, meta_bin, second_bin


def call_idrs(binary, min_len):
    """Return [(start, end)] 1-based inclusive runs of 1s with length >= min_len."""
    idrs, start = [], None
    for i, b in enumerate(binary):
        if b and start is None:
            start = i
        elif not b and start is not None:
            if i - start >= min_len:
                idrs.append((start + 1, i))
            start = None
    if start is not None and len(binary) - start >= min_len:
        idrs.append((start + 1, len(binary)))
    return idrs


def summarize(consensus, meta_bin, second_bin, lcr_bin, min_idr_len):
    """Per-locus disorder metrics from the binary tracks (None for an empty track).

    frac_disordered_nonlcr = consensus-disordered non-LCR residues / non-LCR residues
    ('' when every residue is LCR)."""
    n = len(consensus)
    if n == 0:
        return None
    lcr_bin = list(lcr_bin) if lcr_bin is not None else [0] * n
    if len(lcr_bin) != n:
        raise DisorderError(f"LCR track length {len(lcr_bin)} != {n}")
    idrs = call_idrs(consensus, min_idr_len)
    in_idr = sum(e - s + 1 for s, e in idrs)
    longest = max((e - s + 1 for s, e in idrs), default=0)
    n_lcr = sum(lcr_bin)
    n_nonlcr = n - n_lcr
    dis_nonlcr = sum(1 for c, l in zip(consensus, lcr_bin) if c and not l)
    have_second = any(b is not None for b in second_bin)
    return {
        "length_aa": n,
        "frac_disordered": round(sum(consensus) / n, 4),
        "frac_disordered_metapredict": round(sum(meta_bin) / n, 4),
        "frac_disordered_second": (round(sum(b for b in second_bin if b is not None) / n, 4)
                                   if have_second else ""),
        "n_idrs": len(idrs),
        "longest_idr": longest,
        "frac_in_idrs": round(in_idr / n, 4),
        "frac_lcr": round(n_lcr / n, 4),
        "frac_disordered_nonlcr": (round(dis_nonlcr / n_nonlcr, 4) if n_nonlcr else ""),
    }


def first_profile(res):
    """IUPred library result -> per-residue list. Accepts (scores, glob_text) tuples
    (IUPred2A/IUPred3 style) and a bare score sequence; anything else is an error."""
    def is_seq(x):      # list / tuple / >=1-d numpy array (not str, bytes, dict, scalar)
        if isinstance(x, (list, tuple)):
            return True
        shape = getattr(x, "shape", None)
        return shape is not None and len(shape) >= 1 and hasattr(x, "tolist")

    def is_num(v):
        return isinstance(v, (int, float)) or (hasattr(v, "__float__") and not is_seq(v))
    if is_seq(res) and len(res) > 0 and is_seq(res[0]):
        return list(res[0])
    if is_seq(res) and all(is_num(v) for v in res):
        return list(res)
    raise DisorderError(f"unrecognised IUPred3 return value of type {type(res).__name__}")


def parse_segmasker_interval(text: str):
    """segmasker -outfmt interval -> {id: [(start0, end0), ...]} (0-based inclusive).
    Header lines echo the FASTA defline ('>id description'; 'lcl|' stripped)."""
    out, cur = {}, None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            cur = line[1:].split()[0]
            if cur.startswith("lcl|"):
                cur = cur[4:]
            if cur in out:
                raise DisorderError(f"segmasker output repeats sequence id {cur}")
            out[cur] = []
            continue
        if cur is None:
            raise DisorderError(f"segmasker interval line before any header: {line!r}")
        parts = [p.strip() for p in line.split("-")]
        if len(parts) != 2 or not all(p.isdigit() for p in parts):
            raise DisorderError(f"unparseable segmasker interval line: {line!r}")
        a, b = int(parts[0]), int(parts[1])
        if b < a:
            raise DisorderError(f"segmasker interval end < start: {line!r}")
        out[cur].append((a, b))
    return out


# =========================================================================
# Tool callers (each raises DisorderError on a missing tool or a crash)
# =========================================================================
def load_metapredict(threads: int):
    """Import metapredict (mandatory) and cap torch's intra-op threads to the rule's."""
    try:
        import metapredict as meta
        import torch
    except Exception as exc:  # noqa: BLE001
        raise DisorderError(f"metapredict could not be imported ({exc}); install the "
                            "workflow/envs/disorder.yaml environment") from exc
    torch.set_num_threads(max(1, int(threads)))
    return meta


def metapredict_batch(meta, seqs: dict, version: str, device: str):
    """{locus: seq} -> {locus: [scores]} with metapredict.predict_disorder (dict input
    returns {key: [sequence, np.ndarray]} in metapredict 3.x)."""
    try:
        res = meta.predict_disorder(dict(seqs), version=version, device=device,
                                    normalized=True, round_values=True,
                                    return_numpy=True, show_progress_bar=False)
    except Exception as exc:  # noqa: BLE001
        raise DisorderError(f"metapredict.predict_disorder failed: {exc}") from exc
    return _unpack_batch(res, seqs, "metapredict")


def plddt_batch(meta, seqs: dict, device: str):
    try:
        res = meta.predict_pLDDT(dict(seqs), pLDDT_version=PLDDT_NETWORK, device=device,
                                 normalized=True, round_values=True, return_numpy=True,
                                 show_progress_bar=False)
    except Exception as exc:  # noqa: BLE001
        raise DisorderError(f"metapredict.predict_pLDDT failed: {exc}") from exc
    return _unpack_batch(res, seqs, "pLDDT")


def _unpack_batch(res, seqs, name):
    if not isinstance(res, dict) or set(res) != set(seqs):
        raise DisorderError(f"{name}: batch output keys do not match the input loci")
    out = {}
    for k, s in seqs.items():
        v = res[k]
        if not (isinstance(v, (list, tuple)) and len(v) == 2 and v[0] == s):
            raise DisorderError(f"{name}: unexpected batch entry for {k}")
        out[k] = [float(x) for x in v[1]]
    return out


def load_iupred3(path: str):
    if not path or not os.path.isdir(path):
        raise DisorderError(f"IUPred3 requested but disorder.iupred3_path {path!r} is not "
                            "a directory")
    lib = os.path.join(path, "iupred3_lib.py")
    if not os.path.isfile(lib):
        raise DisorderError(f"{lib} not found (stage the IUPred3 package directory)")
    if path not in sys.path:
        sys.path.insert(0, path)
    try:
        import importlib
        return importlib.import_module("iupred3_lib")
    except Exception as exc:  # noqa: BLE001
        raise DisorderError(f"could not import iupred3_lib from {path}: {exc}") from exc


def iupred3_scores(lib, seq: str, mode: str, smoothing: str, locus: str):
    kwargs = {"smoothing": smoothing} if smoothing else {}
    try:
        res = lib.iupred(seq, mode, **kwargs)
    except Exception as exc:  # noqa: BLE001
        raise DisorderError(f"{locus}: iupred3_lib.iupred(seq, {mode!r}"
                            f"{', smoothing=' + repr(smoothing) if smoothing else ''}) "
                            f"failed: {exc}") from exc
    return first_profile(res)


def low_complexity_positions(seqs: dict, workdir: str):
    """segmasker over the (sanitised) sequences -> {locus: set(1-based LCR positions)}.
    Sequences are written with neutral ids (q1..qN) so no FASTA-id parsing quirk can
    break the mapping back to loci."""
    exe = shutil.which("segmasker")
    if exe is None:
        raise DisorderError("segmasker (BLAST+) is not on PATH but low-complexity "
                            "masking is on (disorder.mask_low_complexity)")
    keys = list(seqs)
    fa = os.path.join(workdir, "lcr_input.faa")
    with open(fa, "w") as fh:
        for i, k in enumerate(keys, 1):
            fh.write(f">q{i}\n{seqs[k]}\n")
    res = subprocess.run([exe, "-in", fa, "-infmt", "fasta", "-outfmt", "interval"],
                         capture_output=True, text=True)
    if res.returncode != 0:
        raise DisorderError(f"segmasker failed (exit {res.returncode}): "
                            f"{res.stderr.strip()[-2000:]}")
    parsed = parse_segmasker_interval(res.stdout)
    expected = {f"q{i}" for i in range(1, len(keys) + 1)}
    if set(parsed) != expected:
        raise DisorderError(f"segmasker reported {len(parsed)} sequences for "
                            f"{len(keys)} inputs")
    out = {}
    for i, k in enumerate(keys, 1):
        n = len(seqs[k])
        pos = set()
        for a, b in parsed[f"q{i}"]:
            if b >= n:
                raise DisorderError(f"segmasker interval {a}-{b} beyond {k} length {n}")
            pos.update(range(a + 1, b + 2))   # 0-based inclusive -> 1-based
        out[k] = pos
    return out


# =========================================================================
# Inputs
# =========================================================================
def load_catalog(path: str):
    """locus_id -> {catalytic_status, fragmented, confidence_tier} (strict columns)."""
    import pandas as pd
    if not path or not os.path.isfile(path):
        raise DisorderError(f"catalog {path!r} not found")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    need = {"locus_id", "catalytic_status", "fragmented", "confidence_tier"}
    miss = need - set(df.columns)
    if miss:
        raise DisorderError(f"catalog {path} lacks column(s) {sorted(miss)}")
    if df["locus_id"].duplicated().any():
        raise DisorderError(f"catalog {path} has duplicated locus_id values")
    return {r["locus_id"]: {"catalytic_status": r["catalytic_status"].strip().lower(),
                            "fragmented": r["fragmented"].strip().lower(),
                            "confidence_tier": r["confidence_tier"].strip()}
            for _, r in df.iterrows()}


def is_intact(entry) -> bool:
    """intact catalytic modules AND explicitly not fragmented (unknown -> excluded)."""
    return entry["catalytic_status"] == "intact" and entry["fragmented"] in FALSE_SET


def load_background(path: str):
    """genome_id -> background disorder fraction (strict; '' path -> {})."""
    import pandas as pd
    if not path:
        return {}
    if not os.path.isfile(path):
        raise DisorderError(f"host background table {path!r} not found")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    vcol = next((c for c in ("background_frac_disordered", "background_frac")
                 if c in df.columns), None)
    if "genome_id" not in df.columns or vcol is None:
        raise DisorderError(f"{path} needs columns genome_id and background_frac_disordered "
                            "(or background_frac)")
    out = {}
    for _, r in df.iterrows():
        g = r["genome_id"].strip()
        try:
            v = float(r[vcol])
        except ValueError:
            raise DisorderError(f"{path}: non-numeric background for {g}: {r[vcol]!r}")
        if not (0.0 <= v <= 1.0):
            raise DisorderError(f"{path}: background for {g} outside [0, 1]: {v}")
        if g in out:
            raise DisorderError(f"{path}: duplicated genome_id {g}")
        out[g] = v
    return out


def resolve_second(choice: str, iupred3_path: str) -> str:
    choice = (choice or "auto").strip().lower()
    if choice == "auto":
        return "iupred3" if iupred3_path else "alphapredict"
    if choice not in SECOND_LABEL:
        raise DisorderError(f"--second-predictor must be auto|iupred3|alphapredict|none, "
                            f"not {choice!r}")
    if choice == "iupred3" and not iupred3_path:
        raise DisorderError("--second-predictor iupred3 needs --iupred3-path")
    return choice


def write_locus_table(path, locus_rows):
    import pandas as pd
    pd.DataFrame(locus_rows, columns=LOCUS_COLUMNS).to_csv(path, sep="\t", index=False)


def _fmt(v):
    return "" if v is None else f"{round(float(v), 4):g}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--catalog", default="")
    ap.add_argument("--intact-only", default="True")
    ap.add_argument("--meta-threshold", type=float, default=0.5)
    ap.add_argument("--iupred3-path", default="")
    ap.add_argument("--iupred3-mode", default="long", choices=["long", "short"])
    ap.add_argument("--iupred3-threshold", type=float, default=0.5)
    ap.add_argument("--plddt-cutoff", type=float, default=50.0)
    ap.add_argument("--min-idr-len", type=int, default=12)
    ap.add_argument("--mask-lcr", default="True")
    ap.add_argument("--host-background", default="")
    ap.add_argument("--out-locus", required=True)
    ap.add_argument("--out-residue", required=True)
    # --- optional (new; backward compatible) ---
    ap.add_argument("--second-predictor", default="auto",
                    help="auto (iupred3 if --iupred3-path else alphapredict) | iupred3 | "
                         "alphapredict | none")
    ap.add_argument("--iupred3-smoothing", default="",
                    help="passed as iupred3_lib.iupred(..., smoothing=...) when set; "
                         "'' = the IUPred3 library default")
    ap.add_argument("--metapredict-version", default="V3",
                    help="metapredict disorder network (pinned; default V3)")
    ap.add_argument("--device", default="cpu",
                    help="torch device for metapredict (cpu = bit-reproducible)")
    ap.add_argument("--threads", type=int, default=1, help="torch intra-op threads")
    args = ap.parse_args(argv)

    if not (0.0 < args.meta_threshold < 1.0):
        raise DisorderError("--meta-threshold must be in (0, 1)")
    if not (0.0 < args.iupred3_threshold < 1.0):
        raise DisorderError("--iupred3-threshold must be in (0, 1)")
    if not (0.0 < args.plddt_cutoff < 100.0):
        raise DisorderError("--plddt-cutoff must be in (0, 100)")
    if args.min_idr_len < 1:
        raise DisorderError("--min-idr-len must be >= 1")
    intact_only = as_bool(args.intact_only)
    mask_lcr = as_bool(args.mask_lcr)
    for p in (args.out_locus, args.out_residue):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    if not os.path.isfile(args.proteins):
        raise DisorderError(f"proteins FASTA {args.proteins!r} not found")
    records = list(utils.read_fasta(args.proteins))
    ids = [r[0] for r in records]
    if len(ids) != len(set(ids)):
        dup = sorted({i for i in ids if ids.count(i) > 1})
        raise DisorderError(f"duplicated protein id(s) in {args.proteins}: {dup[:5]}")
    second = resolve_second(args.second_predictor, args.iupred3_path)
    label = SECOND_LABEL[second]
    if not records:
        utils.write_empty(args.out_locus, LOCUS_HEADER)
        utils.write_empty(args.out_residue, RESIDUE_HEADER)
        utils.log("no proteins in the input FASTA; header-only disorder outputs.")
        return 0

    catalog = load_catalog(args.catalog) if args.catalog else None
    if intact_only and catalog is None:
        raise DisorderError("--intact-only needs --catalog")
    if catalog is not None:
        missing = [i for i in ids if i not in catalog]
        if missing:
            raise DisorderError(f"{len(missing)} protein(s) absent from the catalog "
                                f"(e.g. {missing[:3]})")
    keep, n_skipped = [], 0
    for locus, seq in records:
        if intact_only and not is_intact(catalog[locus]):
            n_skipped += 1
            continue
        keep.append((locus, seq))
    utils.log(f"intact filter: {len(keep)} kept, {n_skipped} skipped "
              f"(intact_only={intact_only}).")
    if not keep:
        utils.write_empty(args.out_locus, LOCUS_HEADER)
        utils.write_empty(args.out_residue, RESIDUE_HEADER)
        utils.log("no locus passed the intact filter; header-only disorder outputs.")
        return 0

    clean, n_conv = {}, {}
    for locus, seq in keep:
        clean[locus], n_conv[locus] = sanitize_sequence(locus, seq)
    background = load_background(args.host_background)
    if background:
        lacking = sorted({utils.genome_of(l) for l in clean} - set(background))
        if lacking:
            raise DisorderError(f"host background table lacks {len(lacking)} scored "
                                f"genome(s) (e.g. {lacking[:3]})")

    # fail fast: every requested tool is located before the (long) metapredict batch
    lib = load_iupred3(args.iupred3_path) if second == "iupred3" else None
    if mask_lcr and shutil.which("segmasker") is None:
        raise DisorderError("segmasker (BLAST+) is not on PATH but low-complexity "
                            "masking is on (disorder.mask_low_complexity)")
    meta = load_metapredict(args.threads)
    detail = (f"metapredict {getattr(meta, '__version__', '?')} {args.metapredict_version} "
              f"(score>={args.meta_threshold})")
    if second == "alphapredict":
        detail += f"; pLDDT {PLDDT_NETWORK} = alphaPredict network (pLDDT<{args.plddt_cutoff:g})"
    elif second == "iupred3":
        detail += (f"; IUPred3 {args.iupred3_mode} smoothing={args.iupred3_smoothing or 'default'} "
                   f"(score>{args.iupred3_threshold})")
    detail += f"; min_idr_len={args.min_idr_len}; LCR={'segmasker' if mask_lcr else 'off'}"
    utils.log(f"{detail}; device={args.device}; consensus label: {label}")
    meta_scores = metapredict_batch(meta, clean, args.metapredict_version, args.device)
    second_scores = {}
    if second == "alphapredict":
        second_scores = plddt_batch(meta, clean, args.device)
    elif second == "iupred3":
        utils.log(f"IUPred3 library {getattr(lib, '__file__', '?')}")
        for locus, seq in clean.items():
            second_scores[locus] = iupred3_scores(lib, seq, args.iupred3_mode,
                                                  args.iupred3_smoothing, locus)
    lcr = {}
    if mask_lcr:
        with tempfile.TemporaryDirectory(
                prefix=".segmasker_", dir=os.path.dirname(os.path.abspath(args.out_locus))) as wd:
            lcr = low_complexity_positions(clean, wd)

    locus_rows = []
    with open(args.out_residue, "w") as res_fh:      # streamed: ~L rows per locus
        res_fh.write(RESIDUE_HEADER + "\n")
        for locus, seq in clean.items():
            meta_s = meta_scores[locus]
            sec_s = second_scores.get(locus) if second != "none" else None
            n = len(seq)
            cons, meta_bin, second_bin = consensus_binary(
                meta_s, sec_s, args.meta_threshold, args.iupred3_threshold,
                second == "alphapredict", args.plddt_cutoff, locus=locus, n=n)
            lcr_set = lcr.get(locus, set())
            lcr_bin = [1 if (i + 1) in lcr_set else 0 for i in range(n)]
            summ = summarize(cons, meta_bin, second_bin, lcr_bin, args.min_idr_len)
            gid = utils.genome_of(locus)
            bg = background.get(gid)
            locus_rows.append({
                "locus_id": locus, "genome_id": gid, "length_aa": summ["length_aa"],
                "intact": (is_intact(catalog[locus]) if catalog is not None else ""),
                "predictors": label,
                "frac_disordered": summ["frac_disordered"],
                "frac_disordered_metapredict": summ["frac_disordered_metapredict"],
                "frac_disordered_second": summ["frac_disordered_second"],
                "n_idrs": summ["n_idrs"], "longest_idr": summ["longest_idr"],
                "frac_in_idrs": summ["frac_in_idrs"], "frac_lcr": summ["frac_lcr"],
                "frac_disordered_nonlcr": summ["frac_disordered_nonlcr"],
                "disorder_minus_background": (round(summ["frac_disordered"] - bg, 4)
                                              if bg is not None else ""),
                "confidence_tier": catalog[locus]["confidence_tier"] if catalog else "",
                "n_nonstandard": n_conv[locus],
                "host_background": bg if bg is not None else "",
                "predictor_detail": detail,
            })
            res_fh.write("".join(
                f"{locus}\t{i + 1}\t{_fmt(meta_s[i])}\t"
                f"{_fmt(sec_s[i]) if sec_s is not None else ''}\t"
                f"{cons[i]}\t{lcr_bin[i]}\t{seq[i]}\n" for i in range(n)))
    write_locus_table(args.out_locus, locus_rows)
    utils.log(f"disorder: scored {len(locus_rows)} loci ({label}); intact filter "
              f"skipped {n_skipped}; LCR masking={'segmasker' if mask_lcr else 'off'}.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except DisorderError as exc:
        utils.log(f"ERROR: {exc}")
        sys.exit(1)
