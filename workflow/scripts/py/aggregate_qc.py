"""Aggregate QC metrics into genome_metadata.tsv and genomeInfo.csv."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd


HEADER = (
    "genome_id\tsource\tassembly_level\tgenome_type\t"
    "completeness\tcontamination\tgenome_size\tn_contigs\tn50\tgc\t"
    "n_rrna\tn_trna\tquality_tier\tpass_primary"
)
GENOMEINFO_HEADER = "genome,completeness,contamination"


def count_rrna(path: str) -> int:
    """Count non-comment, non-empty lines in a GFF (each = one rRNA gene)."""
    if not path or not os.path.exists(path):
        return 0
    n = 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                n += 1
    return n


def count_trna(path: str) -> int:
    """Count data rows in a tRNAscan-SE table (skip the 3 header lines)."""
    if not path or not os.path.exists(path):
        return 0
    n = 0
    with open(path) as fh:
        lines = fh.readlines()
    # Skip first 3 header lines; also skip separator lines (dashes)
    for line in lines[3:]:
        line = line.strip()
        if line and not line.startswith("-"):
            n += 1
    return n


def assign_tier(completeness: float, contamination: float,
                comp_floor: float, cont_ceil: float,
                comp_high: float, cont_high: float):
    if completeness < comp_floor or contamination > cont_ceil:
        return "fail"
    if completeness >= comp_high and contamination <= cont_high:
        return "A"
    if completeness >= 70:
        return "B"
    return "C"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkm2", required=True)
    ap.add_argument("--stats", required=True)
    ap.add_argument("--rrna-dir", required=True)
    ap.add_argument("--trna-dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--ext", required=True)
    ap.add_argument("--comp-floor", type=float, required=True)
    ap.add_argument("--cont-ceil", type=float, required=True)
    ap.add_argument("--comp-high", type=float, required=True)
    ap.add_argument("--cont-high", type=float, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--genomeinfo", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.genomeinfo)), exist_ok=True)

    # --- samples ---
    try:
        samples = pd.read_csv(args.samples, sep="\t")
    except Exception:
        utils.write_empty(args.out, HEADER)
        utils.write_empty(args.genomeinfo, GENOMEINFO_HEADER)
        return

    if samples.empty or "genome_id" not in samples.columns:
        utils.write_empty(args.out, HEADER)
        utils.write_empty(args.genomeinfo, GENOMEINFO_HEADER)
        return

    # --- CheckM2 ---
    checkm2 = {}
    try:
        cm = pd.read_csv(args.checkm2, sep="\t")
        # Normalise column names
        cm.columns = [c.strip() for c in cm.columns]
        for _, row in cm.iterrows():
            name = str(row.get("Name", "")).strip()
            try:
                comp = float(row["Completeness"])
                cont = float(row["Contamination"])
            except (KeyError, ValueError):
                comp, cont = 0.0, 100.0
            checkm2[name] = (comp, cont)
    except Exception:
        utils.log("WARNING: could not read CheckM2 report; defaulting completeness/contamination to 0/100")

    # --- seqkit stats ---
    stats = {}
    try:
        sk = pd.read_csv(args.stats, sep="\t")
        sk.columns = [c.strip() for c in sk.columns]
        # Find GC column (seqkit -a -T can emit 'GC(%)' or 'gc_content' or 'GC')
        gc_col = None
        for candidate in ["GC(%)", "gc_content", "GC", "gc(%)"]:
            if candidate in sk.columns:
                gc_col = candidate
                break
        for _, row in sk.iterrows():
            fname = str(row.get("file", "")).strip()
            base = os.path.basename(fname)
            # Strip extension(s) — remove trailing .{ext}
            ext = args.ext.lstrip(".")
            if base.endswith("." + ext):
                gid = base[: -(len(ext) + 1)]
            else:
                gid = os.path.splitext(base)[0]
            def _num(col, cast, default):
                # Cast one seqkit cell independently so a single bad column
                # (e.g. a missing N50) does not zero out the others.
                try:
                    return cast(row[col])
                except (KeyError, ValueError, TypeError):
                    return default
            genome_size = _num("sum_len", lambda v: int(float(v)), 0)
            n_contigs = _num("num_seqs", lambda v: int(float(v)), 0)
            n50 = _num("N50", lambda v: int(float(v)), 0)
            gc = _num(gc_col, float, 0.0) if gc_col else 0.0
            stats[gid] = (genome_size, n_contigs, n50, gc)
    except Exception:
        utils.log("WARNING: could not read seqkit stats; size/contig/N50/GC will be 0")

    # --- Build output ---
    rows = []
    gi_rows = []

    for _, srow in samples.iterrows():
        gid = str(srow["genome_id"]).strip()
        source = srow.get("source", "")
        asm_level = srow.get("assembly_level", "")
        gtype = srow.get("genome_type", "")

        comp, cont = checkm2.get(gid, (0.0, 100.0))
        genome_size, n_contigs, n50, gc = stats.get(gid, (0, 0, 0, 0.0))

        rrna_path = os.path.join(args.rrna_dir, f"{gid}.rrna.gff")
        trna_path = os.path.join(args.trna_dir, f"{gid}.trna.txt")
        n_rrna = count_rrna(rrna_path)
        n_trna = count_trna(trna_path)

        tier = assign_tier(comp, cont, args.comp_floor, args.cont_ceil,
                           args.comp_high, args.cont_high)
        pass_primary = (comp >= args.comp_floor and cont <= args.cont_ceil)

        rows.append({
            "genome_id": gid,
            "source": source,
            "assembly_level": asm_level,
            "genome_type": gtype,
            "completeness": comp,
            "contamination": cont,
            "genome_size": genome_size,
            "n_contigs": n_contigs,
            "n50": n50,
            "gc": gc,
            "n_rrna": n_rrna,
            "n_trna": n_trna,
            "quality_tier": tier,
            "pass_primary": pass_primary,
        })
        gi_rows.append({
            "genome": f"{gid}.{args.ext.lstrip('.')}",
            "completeness": comp,
            "contamination": cont,
        })

    if not rows:
        utils.write_empty(args.out, HEADER)
        utils.write_empty(args.genomeinfo, GENOMEINFO_HEADER)
        return

    out_df = pd.DataFrame(rows, columns=HEADER.split("\t"))
    out_df.to_csv(args.out, sep="\t", index=False)

    gi_df = pd.DataFrame(gi_rows, columns=["genome", "completeness", "contamination"])
    gi_df.to_csv(args.genomeinfo, index=False)

    utils.log(f"Wrote {len(rows)} genomes -> {args.out}")


if __name__ == "__main__":
    main()
