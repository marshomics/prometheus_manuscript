"""Ingest pre-computed bacterial DefenseFinder output and extract Prometheus candidates.

The bacterial side is provided as COMBINED DefenseFinder tables (one systems TSV and
one genes TSV covering all genomes, each with a genome column added) plus per-genome
`{genome_id}.faa` proteomes (and optionally `{genome_id}.fna` CDS). This script:

  1. splits BOTH combined tables into the per-genome DefenseFinder layout the archaeal
     scripts read -- {out_defense_dir}/{g}/{g}_defense_finder_systems.tsv and
     {g}_defense_finder_genes.tsv (genome column dropped, all other columns kept) -- so
     `validate_prometheus.py` and `build_defense_matrix.py` run UNCHANGED on bacteria.
     Every genome of the analysis universe (--genomes) gets a file, header-only when
     DefenseFinder found nothing there. AntiDefenseFinder rows (`activity` ==
     'Antidefense', DefenseFinder >= 2.0 run with -a) are dropped: the archaeal run never
     calls anti-defense systems, so keeping them would inflate the bacterial defense load.
  2. selects Prometheus systems by EXACT, case-insensitive equality of `type` or
     `subtype` with a configured name (comma-separated list allowed, e.g.
     'Prometheus'); no substring matching.
  3. pulls the member proteins of those systems from the per-genome faa (and the CDS
     from cds_dir), writing candidates + sources + raw hits with the combined id
     `utils.make_protein_id(genome, hit_id)` = '{genome}__{hit_id}' -- the SAME id for
     the protein and its CDS, as codon_align.py requires.

DefenseFinder's `hit_id` is the first whitespace token of the protein FASTA header it
was run on, so the proteome lookup is by that token. CDS are matched by the same token
first and, for NCBI `cds_from_genomic` files, by their `[protein_id=...]` tag.

With --out-genes-dir, every genome carrying Prometheus hits whose proteome headers are
Prodigal's (`>{contig}_{n} # start # end # strand # ID={k}_{n};partial={xy};...`, the
format of GTDB / GTDB-Tk proteomes) also gets {out_genes_dir}/{g}/{g}.genes.tsv in the
SCHEMA item-2 layout of the archaeal `rename_prodigal.py` (contig_length/transl_table
blank: not in protein headers). `validate_prometheus.py --genes-dir` then calls
fragmentation from Prodigal's `partial` code -- the same evidence as for archaea -- instead
of its length proxy (length < min_len), which misses e.g. a 1,500-aa ProA cut at a contig
edge. A proteome with any non-Prodigal header gets no table (no partially covered table,
which validate_prometheus would reject) and is reported as genes_table=none.

No silent loss: a Prometheus genome without a proteome, or a hit_id absent from its
proteome, is a fatal input error unless --allow-missing (then it is reported). Detection
is taken at face value here; architecture/placement/tiering are done by the reused
`validate_prometheus.py`, exactly as for archaea.
"""
from __future__ import annotations
import argparse
import glob
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

SOURCES_HEADER = "protein_id\tgenome_id\tdetect_defensefinder\tdf_score\tdf_subtype"
RAW_HEADER = "protein_id\tgenome_id\tmethod\tscore\tsubtype"
REPORT_HEADER = ("genome_id\tin_universe\tn_systems\tn_antidefense_dropped\t"
                 "n_prometheus_systems\tn_prometheus_genes\tn_proteins_found\t"
                 "n_cds_found\tstatus\tnote\tgenes_table")
GENOME_FALLBACK_COLS = ["genome_id", "genome", "user_genome", "sample", "assembly", "replicon"]
FAA_EXTS = ("faa", "faa.gz", "fasta", "fasta.gz", "fa", "fa.gz")
CDS_EXTS = ("fna", "fna.gz", "ffn", "ffn.gz", "fasta", "fasta.gz", "fa", "fa.gz")
PROTEIN_ID_TAG = re.compile(r"\[protein_id=([^\]\s]+)\]")
# Prodigal protein header (rename_prodigal.py HDR_RE + its ID/partial attributes)
PRODIGAL_HDR = re.compile(r"^(\S+) # (\d+) # (\d+) # (-?1) # ID=(\d+)_(\d+);partial=([01]{2})(;|$)")
GENES_COLS = ["protein_id", "contig", "contig_index", "gene_index", "start", "end", "strand",
              "partial", "contig_length", "length_aa", "transl_table"]   # == rename_prodigal


class InputError(Exception):
    """A malformed or inconsistent input: reported and turned into exit status 1."""


def _pick(df, *candidates):
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand and cand.lower() in low:
            return low[cand.lower()]
    return None


def _open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_fasta_headers(path):
    """Yield (first_token, full_header, seq) records; gzip-aware."""
    name, header, seq = None, None, []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    yield name, header, "".join(seq)
                header = line[1:].strip()
                name = header.split()[0] if header else ""
                seq = []
            elif line and not line.startswith("#"):
                seq.append(line.strip())
    if name is not None:
        yield name, header, "".join(seq)


def find_file(directory, genome, exts):
    for e in exts:
        p = os.path.join(directory, f"{genome}.{e}")
        if os.path.exists(p):
            return p
    return None


def load_proteome(path, want_genes=False):
    """(first token -> sequence, genes rows | None). A duplicated token is an input error.
    genes rows (GENES_COLS) are returned only when want_genes and EVERY header is a
    Prodigal header; otherwise None (validate_prometheus then uses its length proxy)."""
    out, rows, prodigal = {}, [], want_genes
    for tok, hdr, seq in read_fasta_headers(path):
        if tok in out:
            raise InputError(f"duplicated protein id '{tok}' in {path}")
        out[tok] = seq
        if prodigal:
            m = PRODIGAL_HDR.match(hdr or "")
            if not m or m.group(1) != tok:
                prodigal = False
                continue
            _name, start, end, strand, ci, gi, partial, _ = m.groups()
            contig = tok[: -(len(gi) + 1)] if tok.endswith("_" + gi) else tok
            rows.append([utils.sanitize(tok), contig, int(ci), int(gi), int(start), int(end),
                         "+" if strand == "1" else "-", partial, "",
                         len(seq.rstrip("*")), ""])
    if not (prodigal and rows):
        return out, None
    ids = [r[0] for r in rows]
    if len(set(ids)) != len(ids):
        raise InputError(f"protein ids in {path} collide after sanitising (genes table)")
    return out, rows


def clear_genes_tables(out_genes_dir):
    """Remove genes tables from a previous run so no stale genome survives."""
    n = 0
    for f in glob.glob(os.path.join(out_genes_dir, "*", "*.genes.tsv")):
        os.remove(f)
        n += 1
    return n


def load_cds(path):
    """(by_token, by_protein_id_tag) maps of CDS sequences."""
    by_tok, by_tag = {}, {}
    for tok, hdr, seq in read_fasta_headers(path):
        by_tok.setdefault(tok, seq)
        m = PROTEIN_ID_TAG.search(hdr or "")
        if m:
            by_tag.setdefault(m.group(1), seq)
    return by_tok, by_tag


def load_universe(path):
    """Genome ids of the analysis universe (genome_id column, else first column)."""
    import pandas as pd
    if not path:
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if df.shape[1] == 0:
        return []
    col = "genome_id" if "genome_id" in df.columns else df.columns[0]
    ids = [str(g).strip() for g in df[col].dropna().tolist() if str(g).strip()]
    return list(dict.fromkeys(ids))


def clear_split_tables(out_defense_dir):
    """Remove per-genome tables from a previous run so no stale genome survives."""
    n = 0
    for pat in ("*_defense_finder_systems.tsv", "*_defense_finder_genes.tsv"):
        for p in glob.glob(os.path.join(out_defense_dir, "*", pat)):
            os.remove(p)
            n += 1
    return n


def write_split(df, gcol, genomes, out_defense_dir, kind):
    """Write {dir}/{g}/{g}_defense_finder_{kind}.tsv for every genome in `genomes`
    (header-only where the genome has no rows). Returns {genome: n_rows}."""
    cols = [c for c in df.columns if c != gcol]
    counts = {}
    groups = {str(g): sub for g, sub in df.groupby(gcol, sort=False)} if len(df) else {}
    for g in genomes:
        gdir = os.path.join(out_defense_dir, g)
        os.makedirs(gdir, exist_ok=True)
        sub = groups.get(g)
        path = os.path.join(gdir, f"{g}_defense_finder_{kind}.tsv")
        if sub is None:
            with open(path, "w") as fh:
                fh.write("\t".join(cols) + "\n")
            counts[g] = 0
        else:
            sub[cols].to_csv(path, sep="\t", index=False)
            counts[g] = len(sub)
    return counts


def read_table(path, label):
    import pandas as pd
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:  # noqa: BLE001
        raise InputError(f"could not read {label} {path} ({exc})")
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--df-systems", required=True)
    ap.add_argument("--df-genes", required=True)
    ap.add_argument("--faa-dir", required=True)
    ap.add_argument("--cds-dir", default="", help="optional per-genome {g}.fna CDS dir; "
                    "when set, the matching Prometheus CDS are written to --out-cds")
    ap.add_argument("--systems-genome-col", default="genome_id")
    ap.add_argument("--genes-genome-col", default="genome_id")
    ap.add_argument("--system", default="Prometheus",
                    help="DefenseFinder type/subtype name(s), comma-separated; exact, "
                         "case-insensitive match")
    ap.add_argument("--genomes", default="",
                    help="optional analysis-universe TSV (genome_id column): every genome "
                         "gets a per-genome table; id overlap with the DF tables is checked")
    ap.add_argument("--allow-missing", default="False",
                    help="True = report (not fail on) Prometheus genomes whose proteome or "
                         "protein is missing")
    ap.add_argument("--out-defense-dir", required=True)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--out-sources", required=True)
    ap.add_argument("--out-raw", required=True)
    ap.add_argument("--out-cds", default="", help="combined Prometheus CDS (genome__hit_id)")
    ap.add_argument("--out-report", default="", help="per-genome extraction report TSV")
    ap.add_argument("--out-genes-dir", default="",
                    help="optional: write {dir}/{g}/{g}.genes.tsv (SCHEMA item 2) from "
                         "Prodigal proteome headers for genomes with Prometheus hits")
    args = ap.parse_args()
    try:
        run(args)
    except InputError as exc:
        utils.log(f"FATAL: {exc}")
        sys.exit(1)


def run(args) -> None:
    import pandas as pd
    allow_missing = str(args.allow_missing).strip().lower() in ("true", "1", "yes")
    for p in (args.out_candidates, args.out_sources, args.out_raw, args.out_cds, args.out_report):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    os.makedirs(args.out_defense_dir, exist_ok=True)

    sysdf = read_table(args.df_systems, "--df-systems")
    gdf = read_table(args.df_genes, "--df-genes")
    sgcol = _pick(sysdf, args.systems_genome_col, *GENOME_FALLBACK_COLS)
    ggcol = _pick(gdf, args.genes_genome_col, *GENOME_FALLBACK_COLS)
    if sgcol is None:
        raise InputError(f"no genome column in --df-systems (tried '{args.systems_genome_col}', "
                         f"{GENOME_FALLBACK_COLS})")
    if ggcol is None:
        raise InputError(f"no genome column in --df-genes (tried '{args.genes_genome_col}', "
                         f"{GENOME_FALLBACK_COLS})")
    sidcol, g_sidcol = _pick(sysdf, "sys_id"), _pick(gdf, "sys_id")
    tcol, stcol = _pick(sysdf, "type"), _pick(sysdf, "subtype")
    hitcol = _pick(gdf, "hit_id")
    scorecol = _pick(gdf, "hit_score", "score")
    if sidcol is None:
        raise InputError("--df-systems has no 'sys_id' column")
    if tcol is None and stcol is None:
        raise InputError("--df-systems needs a 'type' and/or 'subtype' column")
    if g_sidcol is None or hitcol is None:
        raise InputError("--df-genes needs 'sys_id' and 'hit_id' columns")

    for df, col, label in ((sysdf, sgcol, "--df-systems"), (gdf, ggcol, "--df-genes")):
        blank = df[col].isin(["", "nan", "NaN", "None"])
        if blank.any():
            raise InputError(f"{int(blank.sum())} row(s) of {label} have no genome id "
                             f"(column '{col}'); every row must name its genome")
        bad = sorted({g for g in df[col] if utils.ID_SEP in g})
        if bad:
            raise InputError(f"genome id(s) contain '{utils.ID_SEP}' (breaks the "
                             f"'genome__protein' id scheme): {bad[:3]}")
        bad = sorted({g for g in df[col] if utils.sanitize(g) != g})
        if bad:
            raise InputError(f"genome id(s) with characters outside [A-Za-z0-9_.-] in {label}: "
                             f"{bad[:3]} (ids must match file names and the species-tree tips)")

    # ---- anti-defense rows are not defense systems ------------------------------ #
    n_adf_sys = {}
    actcol = _pick(sysdf, "activity")
    if actcol is not None:
        adf = sysdf[actcol].str.lower() == "antidefense"
        for g, n in sysdf[adf].groupby(sgcol).size().items():
            n_adf_sys[str(g)] = int(n)
        if adf.any():
            utils.log(f"dropping {int(adf.sum())} AntiDefenseFinder system row(s) "
                      f"(activity == Antidefense)")
        sysdf = sysdf[~adf].copy()
    gactcol = _pick(gdf, "activity")
    if gactcol is not None:
        gdf = gdf[gdf[gactcol].str.lower() != "antidefense"].copy()

    df_genomes = list(dict.fromkeys(list(sysdf[sgcol]) + list(gdf[ggcol])))
    universe = load_universe(args.genomes) if args.genomes else None
    if universe is not None:
        uset = set(universe)
        if df_genomes and universe and not (set(df_genomes) & uset):
            raise InputError(
                f"none of the {len(df_genomes)} genome ids in the DefenseFinder tables occur in "
                f"the bacterial genome list ({args.genomes}; e.g. DF '{df_genomes[0]}' vs list "
                f"'{universe[0]}'): the ids must be identical everywhere")
        outside = [g for g in df_genomes if g not in uset]
        if outside:
            utils.log(f"WARNING: {len(outside)} genome(s) in the DefenseFinder tables are not in "
                      f"the bacterial genome list and are ignored (e.g. {outside[:3]})")
        split_genomes = list(dict.fromkeys(universe + df_genomes))
    else:
        uset = None
        split_genomes = df_genomes

    # ---- per-genome DefenseFinder layout (systems + genes) ---------------------- #
    n_stale = clear_split_tables(args.out_defense_dir)
    if n_stale:
        utils.log(f"removed {n_stale} per-genome table(s) from a previous run")
    if args.out_genes_dir:
        os.makedirs(args.out_genes_dir, exist_ok=True)
        n_stale_g = clear_genes_tables(args.out_genes_dir)
        if n_stale_g:
            utils.log(f"removed {n_stale_g} genes table(s) from a previous run")
    n_sys = write_split(sysdf, sgcol, split_genomes, args.out_defense_dir, "systems")
    write_split(gdf, ggcol, split_genomes, args.out_defense_dir, "genes")

    # ---- Prometheus systems: exact type/subtype match --------------------------- #
    names = {s.strip().lower() for s in str(args.system).split(",") if s.strip()}
    if not names:
        raise InputError("--system is empty")
    is_prom = pd.Series(False, index=sysdf.index)
    for c in (tcol, stcol):
        if c is not None:
            is_prom |= sysdf[c].str.lower().isin(names)
    prom = sysdf[is_prom]
    if uset is not None:
        prom = prom[prom[sgcol].isin(uset)]
    subtype_by_key = {}
    for _, r in prom.iterrows():
        subtype_by_key[(r[sgcol], r[sidcol])] = r[stcol] if stcol else r[tcol]
    n_prom_sys = prom.groupby(sgcol).size().to_dict() if len(prom) else {}

    # ---- Prometheus member genes -> (genome, hit_id, score, subtype) ----------- #
    per_genome_hits = {}
    raw_rows = []
    if subtype_by_key:
        mask = pd.Series([k in subtype_by_key for k in zip(gdf[ggcol], gdf[g_sidcol])],
                         index=gdf.index, dtype=bool)
        pg = gdf[mask]
        for _, r in pg.iterrows():
            genome, sid, hit_id = r[ggcol], r[g_sidcol], r[hitcol]
            if not hit_id or hit_id.lower() == "nan":
                raise InputError(f"Prometheus system {sid} of {genome} has a gene row "
                                 f"without hit_id")
            score = None
            if scorecol and r.get(scorecol, "") not in ("", "nan"):
                try:
                    score = float(r[scorecol])
                except ValueError:
                    score = None
            subtype = subtype_by_key[(genome, sid)]
            per_genome_hits.setdefault(genome, []).append((hit_id, score, subtype))
            raw_rows.append((utils.make_protein_id(genome, hit_id), genome, "defensefinder",
                             score, subtype))
        sys_without_genes = set(subtype_by_key) - set(zip(pg[ggcol], pg[g_sidcol]))
        if sys_without_genes:
            raise InputError(f"{len(sys_without_genes)} Prometheus system(s) have no rows in "
                             f"--df-genes (e.g. {sorted(sys_without_genes)[:2]}): the systems "
                             f"and genes tables do not match")

    # ---- proteins (and CDS) ------------------------------------------------------ #
    want_cds = bool(args.cds_dir)
    cand_records, cds_records, source_rows = [], [], []
    report = {}
    problems = []
    seen_pid = {}
    for genome, hits in per_genome_hits.items():
        rep = report.setdefault(genome, {"n_genes": 0, "n_prot": 0, "n_cds": 0,
                                         "status": "ok", "note": [], "genes_table": ""})
        best = {}
        for hit_id, score, subtype in hits:
            best.setdefault(hit_id, (score, subtype))
        rep["n_genes"] = len(best)
        faa = find_file(args.faa_dir, genome, FAA_EXTS)
        if faa is None:
            rep["status"] = "no_proteome"
            problems.append(f"{genome}: no proteome {args.faa_dir}/{genome}.faa")
            continue
        seqs, genes_rows = load_proteome(faa, want_genes=bool(args.out_genes_dir))
        if args.out_genes_dir:
            if genes_rows is None:
                rep["genes_table"] = "none"
                rep["note"].append("proteome headers not Prodigal: fragmentation by length proxy")
            else:
                gdir = os.path.join(args.out_genes_dir, utils.sanitize(genome))
                os.makedirs(gdir, exist_ok=True)
                with open(os.path.join(gdir, f"{utils.sanitize(genome)}.genes.tsv"), "w") as fh:
                    fh.write("\t".join(GENES_COLS) + "\n")
                    for r in genes_rows:
                        fh.write("\t".join(map(str, r)) + "\n")
                rep["genes_table"] = "prodigal"
        cds_tok, cds_tag = {}, {}
        if want_cds:
            cfna = find_file(args.cds_dir, genome, CDS_EXTS)
            if cfna is None:
                rep["note"].append("no CDS file")
            else:
                cds_tok, cds_tag = load_cds(cfna)
        for hit_id, (score, subtype) in best.items():
            if hit_id not in seqs:
                rep["status"] = "protein_missing"
                problems.append(f"{genome}: hit_id '{hit_id}' not in {faa}")
                continue
            pid = utils.make_protein_id(genome, hit_id)
            if pid in seen_pid and seen_pid[pid] != (genome, hit_id):
                raise InputError(f"protein ids collide after sanitising: {seen_pid[pid]} and "
                                 f"{(genome, hit_id)} -> '{pid}'")
            seen_pid[pid] = (genome, hit_id)
            seq = seqs[hit_id].rstrip("*")
            cand_records.append((pid, seq))
            rep["n_prot"] += 1
            if want_cds:
                cds = cds_tok.get(hit_id) or cds_tag.get(hit_id)
                if cds is None:
                    rep["note"].append(f"no CDS for {hit_id}")
                else:
                    cds_records.append((pid, cds.upper()))
                    rep["n_cds"] += 1
                    if len(cds) not in (3 * len(seq), 3 * len(seq) + 3):
                        rep["note"].append(f"CDS length {len(cds)} != 3x{len(seq)} aa for {hit_id}")
            source_rows.append({
                "protein_id": pid, "genome_id": genome, "detect_defensefinder": 1,
                "df_score": "" if score is None else score, "df_subtype": subtype})

    if problems:
        msg = (f"{len(problems)} Prometheus protein(s)/proteome(s) missing from --faa-dir "
               f"(e.g. {problems[:3]})")
        if not allow_missing:
            raise InputError(msg + "; fix the inputs or set bacteria.allow_missing_proteins")
        utils.log("WARNING: " + msg + " -- reported in --out-report and skipped")
    if want_cds and cand_records and not cds_records:
        raise InputError(f"bacteria.cds_dir is set but no CDS matched any of the "
                         f"{len(cand_records)} candidate proteins (headers must equal the "
                         f"protein ids or carry [protein_id=...])")
    n_no_cds = len(cand_records) - len(cds_records)
    if want_cds and n_no_cds:
        utils.log(f"WARNING: {n_no_cds} candidate(s) have no CDS; they are absent from the "
                  f"cross-domain codon alignment")

    # ---- write outputs ----------------------------------------------------------- #
    utils.write_fasta(sorted(set(cand_records)), args.out_candidates)
    if args.out_cds:
        utils.write_fasta(sorted(set(cds_records)), args.out_cds)
    if source_rows:
        pd.DataFrame(source_rows, columns=SOURCES_HEADER.split("\t")).to_csv(
            args.out_sources, sep="\t", index=False)
    else:
        utils.write_empty(args.out_sources, SOURCES_HEADER)
    if raw_rows:
        pd.DataFrame([{"protein_id": p, "genome_id": g, "method": m,
                       "score": "" if s is None else s, "subtype": st}
                      for p, g, m, s, st in raw_rows],
                     columns=RAW_HEADER.split("\t")).to_csv(args.out_raw, sep="\t", index=False)
    else:
        utils.write_empty(args.out_raw, RAW_HEADER)
    if args.out_report:
        with open(args.out_report, "w") as fh:
            fh.write(REPORT_HEADER + "\n")
            for g in split_genomes:
                rep = report.get(g, {"n_genes": 0, "n_prot": 0, "n_cds": 0, "status": "ok",
                                     "note": [], "genes_table": ""})
                in_u = "" if uset is None else str(g in uset)
                status = rep["status"] if (uset is None or g in uset) else "not_in_universe"
                fh.write("\t".join(map(str, [
                    g, in_u, n_sys.get(g, 0), n_adf_sys.get(g, 0), n_prom_sys.get(g, 0),
                    rep["n_genes"], rep["n_prot"], rep["n_cds"], status,
                    "; ".join(rep["note"])[:500], rep["genes_table"]])) + "\n")
    utils.log(f"Bacterial {args.system}: {len(cand_records)} candidate protein(s) in "
              f"{len(per_genome_hits)} genome(s); CDS for {len(cds_records)}; "
              f"{len(split_genomes)} per-genome defense tables written.")


if __name__ == "__main__":
    main()
