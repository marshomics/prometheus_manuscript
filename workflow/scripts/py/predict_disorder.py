"""Predict intrinsic disorder on validated Prometheus/ProA proteins by CONSENSUS.

A residue is called disordered ONLY where two methods agree, so no single predictor's
bias drives the result:
  * metapredict (deep-learning consensus disorder; pip/conda) -- always the first call;
  * IUPred3 (biophysical energy-estimation; INDEPENDENT methodology) when its path is
    staged (config disorder.iupred3_path), else alphaPredict (predicted AF2 pLDDT) as
    an installable but CORRELATED fallback (clearly labelled).

Defensibility controls applied here:
  * intact-only -- a truncated/fragmented locus shows artefactual disorder at broken
    ends, so by default only catalytically-intact, non-fragmented loci are scored;
  * low-complexity awareness -- LCR is not disorder; segmasker (BLAST+) masks LCR and
    we additionally report disorder EXCLUDING LCR;
  * (optional) host-proteome normalisation -- subtract each genome's background disorder
    so the metric is 'unusually disordered for its host', not lineage baseline.

Outputs:
  disorder_per_locus.tsv    one row per scored locus (the comparative-trait table)
  disorder_per_residue.tsv  per-residue consensus track (for alignment mapping + figure)

The pure scoring/summary logic is separated from the tool calls so it is unit-tested
without the heavy ML dependencies; if metapredict is unavailable the module degrades to
a header-only output (the DAG never breaks), like every other tool wrapper here.
"""
from __future__ import annotations
import argparse
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

LOCUS_HEADER = (
    "locus_id\tgenome_id\tlength_aa\tintact\tpredictors\t"
    "frac_disordered\tfrac_disordered_metapredict\tfrac_disordered_second\t"
    "n_idrs\tlongest_idr\tfrac_in_idrs\tfrac_lcr\tfrac_disordered_nonlcr\t"
    "disorder_minus_background"
)
RESIDUE_HEADER = "locus_id\tposition\tmeta_score\tsecond_score\tconsensus\tlcr"


# =========================================================================
# Pure logic (unit-tested without metapredict/IUPred3 installed)
# =========================================================================
def consensus_binary(meta_scores, second_scores, meta_thr, second_thr,
                     second_is_plddt, plddt_cut):
    """Per-residue 1 iff metapredict AND the second predictor both call disorder.

    second_is_plddt: the second predictor returns pLDDT (alphaPredict), so 'disordered'
    means pLDDT < plddt_cut; otherwise it is a disorder score and 'disordered' means
    score >= second_thr. Returns (consensus, meta_bin, second_bin) lists of 0/1.
    """
    n = len(meta_scores)
    meta_bin = [1 if (s is not None and s >= meta_thr) else 0 for s in meta_scores]
    if second_scores is None:
        # single-predictor mode: consensus == metapredict
        return list(meta_bin), meta_bin, [None] * n
    second_bin = []
    for s in second_scores:
        if s is None:
            second_bin.append(0)
        elif second_is_plddt:
            second_bin.append(1 if s < plddt_cut else 0)
        else:
            second_bin.append(1 if s >= second_thr else 0)
    m = min(n, len(second_bin))
    cons = [1 if (meta_bin[i] and second_bin[i]) else 0 for i in range(m)]
    return cons, meta_bin[:m], second_bin[:m]


def call_idrs(binary, min_len):
    """Return [(start,end)] 1-based inclusive runs of 1s with length >= min_len."""
    idrs = []
    start = None
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
    """Per-locus disorder metrics from the binary tracks."""
    n = len(consensus)
    if n == 0:
        return None
    n_dis = sum(consensus)
    idrs = call_idrs(consensus, min_idr_len)
    in_idr = sum(e - s + 1 for s, e in idrs)
    longest = max((e - s + 1 for s, e in idrs), default=0)
    lcr_bin = lcr_bin or [0] * n
    frac_lcr = sum(lcr_bin[:n]) / n
    dis_nonlcr = sum(1 for i in range(n) if consensus[i] and not lcr_bin[i]) / n
    second_frac = (sum(b for b in second_bin if b is not None) / n
                   if any(b is not None for b in second_bin) else "")
    return {
        "length_aa": n,
        "frac_disordered": round(n_dis / n, 4),
        "frac_disordered_metapredict": round(sum(meta_bin) / n, 4),
        "frac_disordered_second": (round(second_frac, 4) if second_frac != "" else ""),
        "n_idrs": len(idrs),
        "longest_idr": longest,
        "frac_in_idrs": round(in_idr / n, 4),
        "frac_lcr": round(frac_lcr, 4),
        "frac_disordered_nonlcr": round(dis_nonlcr, 4),
    }


# =========================================================================
# Tool callers (return None when the tool is unavailable)
# =========================================================================
def metapredict_scores(seq):
    try:
        import metapredict as meta
    except Exception:
        return None
    try:
        out = meta.predict_disorder(seq)
        # v2 returns a list; some builds return an object with .disorder
        return list(getattr(out, "disorder", out))
    except Exception as exc:  # noqa: BLE001
        utils.log(f"metapredict failed on a sequence: {exc}")
        return None


def alphapredict_plddt(seq):
    try:
        import alphaPredict as alpha
    except Exception:
        return None
    try:
        return list(alpha.predict(seq))
    except Exception:
        return None


def iupred3_scores(seq, iupred3_path, mode):
    """IUPred3 per-residue disorder via its library (preferred) or CLI subprocess."""
    if not iupred3_path or not os.path.isdir(iupred3_path):
        return None
    # library import (fast)
    try:
        if iupred3_path not in sys.path:
            sys.path.insert(0, iupred3_path)
        import iupred3_lib  # noqa
        res = iupred3_lib.iupred(seq, mode)
        scores = res[0] if isinstance(res, (list, tuple)) else res
        return list(scores)
    except Exception:
        pass
    # CLI fallback
    script = os.path.join(iupred3_path, "iupred3.py")
    if not os.path.exists(script):
        return None
    fd, fa = tempfile.mkstemp(suffix=".fasta"); os.close(fd)
    try:
        with open(fa, "w") as fh:
            fh.write(">q\n" + seq + "\n")
        res = subprocess.run([sys.executable, script, fa, mode],
                             capture_output=True, text=True)
        if res.returncode != 0:
            return None
        scores = []
        for line in res.stdout.splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) >= 3:
                try:
                    scores.append(float(f[2]))
                except ValueError:
                    continue
        return scores or None
    finally:
        try:
            os.remove(fa)
        except OSError:
            pass


def low_complexity_positions(faa_path):
    """segmasker -outfmt interval -> {locus_id: set(1-based LCR positions)} or None."""
    from shutil import which
    if which("segmasker") is None:
        return None
    try:
        res = subprocess.run(["segmasker", "-in", faa_path, "-outfmt", "interval"],
                             capture_output=True, text=True)
    except Exception:
        return None
    if res.returncode != 0:
        return None
    out, cur = {}, None
    for line in res.stdout.splitlines():
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out.setdefault(cur, set())
        elif cur is not None and "-" in line:
            try:
                a, b = [int(x) for x in line.replace("-", " ").split()[:2]]
                out[cur].update(range(a + 1, b + 2))  # interval 0-based inclusive -> 1-based
            except (ValueError, IndexError):
                continue
    return out


# =========================================================================
def intact_loci(catalog_path):
    """locus_ids that are catalytically intact AND not fragmented (truncation guard)."""
    keep = None
    if not catalog_path or not os.path.exists(catalog_path):
        return keep
    try:
        import pandas as pd
        df = pd.read_csv(catalog_path, sep="\t", dtype=str)
        if "locus_id" not in df.columns:
            return keep
        ok = df
        if "catalytic_status" in df.columns:
            ok = ok[ok["catalytic_status"].astype(str).str.lower() == "intact"]
        if "fragmented" in ok.columns:
            ok = ok[~ok["fragmented"].astype(str).str.lower().isin(("true", "1", "yes"))]
        return set(ok["locus_id"].astype(str))
    except Exception:
        return keep


def load_background(path):
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        import pandas as pd
        df = pd.read_csv(path, sep="\t")
        gcol = "genome_id" if "genome_id" in df.columns else df.columns[0]
        vcol = next((c for c in df.columns if "background" in c.lower() or "frac" in c.lower()), None)
        if vcol:
            for _, r in df.iterrows():
                out[str(r[gcol]).strip()] = float(r[vcol])
    except Exception:
        pass
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--catalog", default="")
    ap.add_argument("--intact-only", default="True")
    ap.add_argument("--meta-threshold", type=float, default=0.5)
    ap.add_argument("--iupred3-path", default="")
    ap.add_argument("--iupred3-mode", default="long")
    ap.add_argument("--iupred3-threshold", type=float, default=0.5)
    ap.add_argument("--plddt-cutoff", type=float, default=50.0)
    ap.add_argument("--min-idr-len", type=int, default=12)
    ap.add_argument("--mask-lcr", default="True")
    ap.add_argument("--host-background", default="")
    ap.add_argument("--out-locus", required=True)
    ap.add_argument("--out-residue", required=True)
    args = ap.parse_args()

    for p in (args.out_locus, args.out_residue):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    records = list(utils.read_fasta(args.proteins))
    if not records:
        utils.write_empty(args.out_locus, LOCUS_HEADER)
        utils.write_empty(args.out_residue, RESIDUE_HEADER)
        utils.log("no proteins; header-only disorder outputs.")
        return

    # metapredict is mandatory (the first half of every consensus); if absent, degrade.
    if metapredict_scores("MAAAA") is None:
        utils.write_empty(args.out_locus, LOCUS_HEADER)
        utils.write_empty(args.out_residue, RESIDUE_HEADER)
        utils.log("metapredict not available; header-only disorder outputs.")
        return

    # choose the second predictor
    use_iupred = bool(args.iupred3_path) and os.path.isdir(args.iupred3_path)
    second_is_plddt = False
    label = "metapredict_only"
    if use_iupred:
        label = "metapredict+iupred3"
    elif alphapredict_plddt("MAAAA") is not None:
        second_is_plddt = True
        label = "metapredict+alphapredict"
    utils.log(f"disorder consensus: {label}")

    intact = None
    if str(args.intact_only).lower() in ("true", "1", "yes"):
        intact = intact_loci(args.catalog)
    lcr = low_complexity_positions(args.proteins) if str(args.mask_lcr).lower() in ("true", "1", "yes") else None
    background = load_background(args.host_background)

    locus_rows, residue_lines = [], []
    n_skipped_intact = 0
    for locus, seq in records:
        if intact is not None and locus not in intact:
            n_skipped_intact += 1
            continue
        meta = metapredict_scores(seq)
        if meta is None:
            continue
        if use_iupred:
            second = iupred3_scores(seq, args.iupred3_path, args.iupred3_mode)
        elif second_is_plddt:
            second = alphapredict_plddt(seq)
        else:
            second = None
        cons, meta_bin, second_bin = consensus_binary(
            meta, second, args.meta_threshold, args.iupred3_threshold,
            second_is_plddt, args.plddt_cutoff)
        n = len(cons)
        lcr_set = (lcr or {}).get(locus, set())
        lcr_bin = [1 if (i + 1) in lcr_set else 0 for i in range(n)]
        summ = summarize(cons, meta_bin, second_bin, lcr_bin, args.min_idr_len)
        if summ is None:
            continue
        gid = utils.genome_of(locus)
        bg = background.get(gid)
        dmb = round(summ["frac_disordered"] - bg, 4) if bg is not None else ""
        locus_rows.append({
            "locus_id": locus, "genome_id": gid, "length_aa": summ["length_aa"],
            "intact": (intact is None) or (locus in (intact or set())),
            "predictors": label,
            "frac_disordered": summ["frac_disordered"],
            "frac_disordered_metapredict": summ["frac_disordered_metapredict"],
            "frac_disordered_second": summ["frac_disordered_second"],
            "n_idrs": summ["n_idrs"], "longest_idr": summ["longest_idr"],
            "frac_in_idrs": summ["frac_in_idrs"], "frac_lcr": summ["frac_lcr"],
            "frac_disordered_nonlcr": summ["frac_disordered_nonlcr"],
            "disorder_minus_background": dmb,
        })
        for i in range(n):
            sc2 = second_bin[i] if i < len(second_bin) else None
            ms = round(meta[i], 4) if i < len(meta) and meta[i] is not None else ""
            ss = (round(second[i], 4) if (second is not None and i < len(second)
                                          and second[i] is not None) else "")
            residue_lines.append(f"{locus}\t{i+1}\t{ms}\t{ss}\t{cons[i]}\t{lcr_bin[i]}")

    if not locus_rows:
        utils.write_empty(args.out_locus, LOCUS_HEADER)
        utils.write_empty(args.out_residue, RESIDUE_HEADER)
        utils.log(f"no loci scored (intact filter skipped {n_skipped_intact}); header-only.")
        return

    import pandas as pd
    pd.DataFrame(locus_rows, columns=LOCUS_HEADER.split("\t")).to_csv(
        args.out_locus, sep="\t", index=False)
    with open(args.out_residue, "w") as fh:
        fh.write(RESIDUE_HEADER + "\n")
        fh.write("\n".join(residue_lines) + "\n")
    utils.log(f"disorder: scored {len(locus_rows)} loci ({label}); "
              f"intact filter skipped {n_skipped_intact}.")


if __name__ == "__main__":
    main()
