# Review: stage `preprocess_qc` (leaf-1.1.1)

Phase-1 preprocessing (`workflow/rules/preprocess.smk`: batch_genomes, gtdbtk_batch,
checkm2_batch, bbmap_batch, merge_annotations, filter_collapse, plot_preprocess) and Stage-0
gene calls / QC (`workflow/rules/qc.smk`: prodigal_genes, genome_metadata, plot_qc), their
scripts, `plot_qc_extra.py` (rule in plots_extra.smk, CLI frozen) and the envs
`qc.yaml`, `bbmap.yaml`, `checkm2.yaml` (new). Reviewed against the pristine copy at
commit 56f5437; `file:line` refers to the **pristine** file unless marked "(new)".

Suite: `bash tests/run_leaf.sh preprocess_qc` -> 27 tests, prints `PREPROCESS_QC TESTS PASSED`
(about 50 s; at most 2 threads). Real tools exercised: Prodigal 2.6.3 (planted-gene genome,
20-kb training boundary, rule-level runs), Snakemake 8.30 (the real rules of preprocess.smk +
qc.smk end to end), awk (the rule's own merge program). GTDB-Tk 2.4.0 / CheckM2 1.0.2 need
10-100 GB databases: they are verified from source at the pinned tags and replaced in the
rule-level test by stand-ins that write their exact file layout and columns. BBMap 39.06 was
run for real: the bioconda linux-64 build `bbmap-39.06-h92535d8_1` (Java classes + shell
scripts; its optional JNI library `jni/libbbtoolsjni.*` is not used by statswrapper.sh)
unpacked and run with conda-forge OpenJDK 17 (env
`prom-preprocess_qc-bbmap`), because no macOS build of 39.06 exists. The suite compares real
BBMap with the stand-in whenever `statswrapper.sh` is on the toolchain PATH or
`PROM_BBMAP_STATSWRAPPER` points at it (JVM capped at 2 processors; 27/27 passed that way here),
and always parses the captured real 39.06 output.

## Primary-source evidence used

| Tool (pinned) | Source | What was verified |
|---|---|---|
| GTDB-Tk 2.4.0 | github Ecogenomics/GTDBTk tag 2.4.0: `gtdbtk/config/output.py`, `main.py`, `classify.py`, `files/classify_summary.py`, `files/batchfile.py`, `biolib_lite/common.py`, `tools.py`, `cli.py`, `__main__.py` | summaries are `classify/{prefix}.{ar53,bac120}.summary.tsv`; classify.py `symlink_f` links both into the out_dir root; `identify/{prefix}.*.markers_summary.tsv` and `classify/ani_screen/*_skani.summary.tsv` also match `*summary.tsv`; the 20 summary columns (`user_genome` ... `warnings`, `N/A` for None); filtered genomes -> `Unclassified Archaea/Bacteria`, Prodigal/align failures -> `Unclassified` (bac120 summary); `_genomes_to_process` raises `GTDBTkExit` when a user id equals a reference id, and `get_reference_ids` adds the bare `GCA_/GCF_` accession of every `RS_/GB_` reference; ids may not contain `()[],;= `, paths may not contain spaces; batchfile = `path<TAB>id[<TAB>table]`; `--skip_ani_screen` still required unless `--mash_db` (main.py:1094); a completed ANI screen in `out_dir/gtdbtk.json` is reused; `GTDBTkExit` -> exit 1; `--version` performs a network version check |
| CheckM2 1.0.2 | github chklovski/CheckM2 tag 1.0.2: `bin/checkm2`, `checkm2/predictQuality.py`, `checkm2/prodigal.py`, `checkm2/fileManager.py`, `checkm2/versionControl.py`, `checkm2/version/version_hashes_1.0.2.json` (and `_1.1.0.json` at tag 1.1.0) | DIAMOND DB versions 0-2 accepted (v2 = zenodo 5571251; 1.1.0 adds v3 = zenodo 14897628, `incompatible_below_checkm2ver` 1.1.0); `predict --database_path` sha256-checks the DB and `sys.exit(1)`s on an unknown hash; `database --download` fetches the highest compatible version into `CheckM2_database/uniref100.KO.1.dmnd`; flags `--input/-i --output-directory -x --threads --force --database_path --remove_intermediates --version`; bins = files ending in `-x`; `Name` = `os.path.splitext(basename)[0]` (prodigal.py:25); 0-byte bins and bins with an empty protein file are skipped with only a warning; report columns `Name Completeness Contamination Completeness_Model_Used Translation_Table_Used Coding_Density Contig_N50 Average_Gene_Length Genome_Size GC_Content Total_Coding_Sequences Total_Contigs Max_Contig_Length Additional_Notes`; completeness/contamination `np.round(..., 2)`; missing DB / no bins -> `sys.exit(1)` |
| BBMap 39.06 | bioconda `bbmap-39.06-h92535d8_1` (ships `current/jgi/AssemblyStats2.java`, `AssemblyStatsWrapper.java`, `statswrapper.sh`) + a real run | format=6 prints ONE `#`-prefixed header per JVM (wrapper sets `header=f` after the first file), `N_UNDERSCORE=true` by default, wrapper adds `addname=t` -> `filename` column; `scaf_L50`/`ctg_L50` are N50 LENGTHS, `*_N50` the counts; `formatKB()` prints plain integers when `FORMAT>=2`; `gc_avg` = (G+C)/(A+C+G+T), a fraction; contigs split at `maxNs=10` N; an unreadable input throws `RuntimeException`, exit 1, after printing the earlier rows. Real run on a 50/30/20 kb + (10 kb+N20+5 kb) assembly: `n_contigs 5, contig_bp 115000, scaf_bp 115020, ctg_N50 2, ctg_L50 30000, gc_avg 0.50070` |
| Prodigal 2.6.3 | github hyattpd/Prodigal tag v2.6.3 `main.c`, `sequence.c` + real runs (prom-tools) | `MIN_SINGLE_GENOME 20000` -> exit 10; training length = every char `A`..`z` + 12 bp per record when >1 record (verified: 19,975 letters in 2 records -> exit 10, 19,976 -> ok); empty/header-only input -> exit 18; header `>{seqid}_{n} # start # end # 1\|-1 # ID={k}_{n};partial={xy};...`; GFF `# Sequence Data: seqnum=;seqlen=;seqhdr="..."` + `# Model Data: ...transl_table=` per record (also records without genes); in `-p meta` the user `-g` is never applied (all 50 meta models, incl. table-4 ones, are tried); `-v` exits 0; `-m` drops a gene spanning an N-gap entirely (tested; not adopted) |
| Snakemake 8.30 | `snakemake/shell.py:86` | every bash shell block runs under `set -euo pipefail;` |
| conda | `micromamba create --platform linux-64 --dry-run` | pristine qc.yaml does not solve; new qc/bbmap/checkm2 envs solve (G3) |

## Findings

### PQ-01 · critical · preprocess.smk:70-77 (gtdbtk_batch) — a crashed GTDB-Tk was swallowed
- **Issue**: `gtdbtk classify_wf ... || true` followed by `printf "user_genome\tclassification\n"` turned any GTDB-Tk failure (crash, missing reference data, id collision, OOM) into a header-only summary. The whole batch (~20k genomes) then had no taxonomy and was silently removed by the domain filter.
- **Evidence**: GTDB-Tk 2.4.0 `__main__.py` exits 1 on `GTDBTkExit`; the rule discarded that exit code.
- **Fix**: no `|| true`; fail fast when `db.gtdbtk_data` is not a directory or `GTDBTK_DATA_PATH` is unset; fail when GTDB-Tk exits 0 but writes no classify summary; record the GTDB-Tk version in the log.
- **Test**: `test_snakemake_phase1_phase2_end_to_end`, scenarios `gt` (simulated crash -> rule fails, no summary) and `gtdb` (absent data dir).

### PQ-02 · critical · preprocess.smk:72-76 (gtdbtk_batch) — summary merge read the wrong files
- **Issue**: `find <outdir> -name '*summary.tsv'` also matched GTDB-Tk's root-level symlinks to the two classify summaries (every genome twice), `identify/*.markers_summary.tsv` (a different table whose header could become the merged header) and `classify/ani_screen/*_skani.summary.tsv`. `tail -n +2` of a table with a different header silently mixed columns.
- **Evidence**: gtdbtk/config/output.py (`PATH_AR53_SUMMARY_OUT`, `PATH_*_MARKER_SUMMARY`, `PATH_*_ANI_SKANI`); classify.py `symlink_f(PATH_AR53_SUMMARY_OUT ..., os.path.join(out_dir, basename))`.
- **Fix**: merge exactly `classify/gtdbtk.ar53.summary.tsv` and `classify/gtdbtk.bac120.summary.tsv` (whichever exist) with an awk program that requires `user_genome` as the first header field, requires identical headers, appends a `marker_set` column and fails on any row without the batch id prefix. merge_annotations then rejects duplicated or unknown ids.
- **Test**: `test_gtdbtk_summary_merge_awk` (runs the rule's own awk; positive and negative controls) and the rule-level test (the stand-in writes the root symlinks and markers summaries; asserts exact header, no duplicate, all genomes).

### PQ-03 · critical · preprocess.smk:68-70 + batch_genomes.py — GTDB-Tk aborts on genomes named like GTDB references
- **Issue**: genomes were passed to GTDB-Tk by file stem = genome_id. GTDB-Tk 2.4.0 aborts the entire run if any user id equals a GTDB reference id, and that set includes the bare NCBI accession of every representative (`GCA_000008085.1`). An NCBI-derived input set named by accession makes every batch containing a GTDB reference genome fail (before PQ-01: the batch was silently lost).
- **Evidence**: gtdbtk/main.py `_genomes_to_process` (`invalid_genomes = set(genomic_files) & set(get_reference_ids())` -> `GTDBTkExit`); gtdbtk/tools.py `get_reference_ids` adds `raw_id[3:]` for `RS_`/`GB_` ids.
- **Fix**: batch_genomes.py writes `batch_XXXX/gtdbtk_batchfile.tsv` (`<staged fasta>\tuser_<genome_id>`; format from gtdbtk/files/batchfile.py) and the rule runs `classify_wf --batchfile`; the merge strips the prefix (rows without it are an error). CheckM2/BBMap keep reading the `{genome_id}.fna` symlinks.
- **Test**: `test_batch_genomes_layout_matches_common_smk` (batchfile content), `test_gtdbtk_summary_merge_awk` (prefix stripped / missing prefix rejected), rule-level test (the stand-in refuses ids without the prefix).

### PQ-04 · critical · preprocess.smk:102-108 (checkm2_batch) — a crashed CheckM2 was swallowed
- **Issue**: `checkm2 predict ... || true` plus a placeholder `Name Completeness Contamination` header: a crash or missing database left every genome of the batch without completeness, so all of them silently failed the quality filter.
- **Evidence**: CheckM2 1.0.2 predictQuality.py `sys.exit(1)` on missing DB / no bins.
- **Fix**: no `|| true`; fail when `db.checkm2` is missing or empty; fresh output directory; `cp` of `quality_report.tsv` fails if CheckM2 wrote none; version logged; own env (PQ-06).
- **Test**: rule-level scenarios `cm` (simulated crash) and `cmdb` (absent database).

### PQ-05 · critical · qc.smk:33-36 (prodigal_genes) — a crashed Prodigal produced an empty proteome
- **Issue**: `prodigal ... || true` then `touch` the outputs: a Prodigal failure (exit 10 single mode < 20 kb, exit 18 no sequence, ...) became an empty `.faa`, read downstream as "Prometheus absent".
- **Evidence**: Prodigal 2.6.3 main.c exit codes 10/18 (reproduced with the real binary).
- **Fix**: no `|| true`/`touch`; Prodigal writes `*.raw` files; rename_prodigal.py fails on zero predicted genes or inconsistent outputs; raw files removed only after success.
- **Test**: rule-level `BOGUS` (header-only genome: Prodigal exit 18 -> rule fails, no `.faa`) and `NOGENE` (25 kb without ORFs -> "0 genes" -> rule fails); `test_rename_prodigal_failure_modes`.

### PQ-06 · critical · envs/qc.yaml:7,16 — environment unsolvable
- **Issue**: `checkm2=1.0.2` (pins numpy 1.19.2 and a tensorflow 2.4 stack) cannot be solved with `python=3.11` / `numpy=1.26.4`; checkm2_batch, prodigal_genes and the provenance rule could not build their env. drep, barrnap, trnascan-se, assembly-stats, seqkit and pyrodigal were never used.
- **Evidence**: micromamba linux-64 dry-run solve (driver); the new checkm2-only env resolves python 3.8.15, numpy 1.19.2, tensorflow 2.4.0, diamond 2.0.4, lightgbm 3.2.1, scikit-learn 0.23.2. Database compatibility from CheckM2 source: `checkm2/version/version_hashes_1.0.2.json` lists DIAMOND DB versions 0-2 (v2 = zenodo 10.5281/zenodo.5571251); `version_hashes_1.1.0.json` adds v3 (zenodo 10.5281/zenodo.14897628, `incompatible_below_checkm2ver` 1.1.0); both releases read `CheckM2_database/uniref100.KO.1.dmnd`.
- **Fix**: `qc.yaml` = python 3.11 + prodigal 2.6.3 + pandas 2.2.2 + numpy 1.26.4; new `checkm2.yaml` = checkm2 1.0.2 + python 3.8 + diamond 2.0.4; its header records that the configured `db.checkm2` must be DIAMOND DB v<=2 (a v3 download needs CheckM2 >=1.1.0); checkm2_batch uses it; `bbmap.yaml` drops the unused, numpy-unpinned pandas.
- **Test**: G3 `tests/check_env_solve.sh` -> `ENV SOLVE OK` for all three; `test_rule_hygiene_and_static_contract` checks pins and env wiring.

### PQ-07 · major · run_bbmap_stats.py:69-77 — BBMap failures were warnings
- **Issue**: statswrapper.sh missing from PATH -> header-only output (every genome silently loses n50/n_contigs and fails QC); a non-zero exit only logged a WARNING, although BBMap aborts the JVM at the first unreadable file, so the remaining genomes of that 2000-genome chunk were silently missing. Short data lines were silently skipped.
- **Evidence**: real BBMap 39.06 run with a missing input: exit 1 after printing the earlier rows (`RuntimeException: Input file does not appear to be valid`).
- **Fix**: fail when the tool is absent, on any non-zero exit, on a header without `filename`, on a malformed line, and unless every genome of the chunk has exactly one row; run in the genome directory with bare file names.
- **Test**: `test_run_bbmap_stats_chunks_and_failures`, rule-level scenario `bb`; `test_run_bbmap_stats_real_bbmap_if_available` (real BBMap == stand-in on all compared columns).

### PQ-08 · major · merge_annotations.py:106-113 — missing BBMap column only warned
- **Issue**: a configured BBMap column absent from the statswrapper output (e.g. a mistyped `n50_col`) gave blank n50 for every genome and all of them failed the N50 filter after a WARNING. Non-numeric values (e.g. `45.1 KB` from a non-machine format) were coerced to NaN silently.
- **Evidence**: BBMap format=6 column names from AssemblyStats2.java and the real run.
- **Fix**: hard error for a missing column and for any non-empty value that does not parse as a number; error when `n50 > genome_size`; the CheckM2 `Contig_N50` (scaffold-level) cross-check is logged, with a warning when most genomes disagree by more than 2-fold (a count/length column mix-up).
- **Test**: `test_merge_annotations_hard_errors` (missing `ctg_L50`, `30 KB`, swapped columns).

### PQ-09 · major · merge_annotations.py:36-46,77-102 — inputs located by glob, read errors skipped
- **Issue**: `glob(batch_dir/batch_*/...)` also read batch directories left by an earlier, larger run (stale genomes re-entered the table); an unreadable file was skipped with a warning.
- **Evidence**: pristine `_read_concat` (36-46) globbed `batch_dir/batch_*/<file>` for all three tools (77, 92, 102) and logged `WARNING: could not read` on a read error (44-45); common.smk derives BATCHES from the current metadata, so a `batch_*` directory left by an earlier, larger run is outside the DAG but still matched the glob.
- **Fix**: the rule passes the exact per-batch file lists of the current DAG (`--gtdbtk/--checkm2/--bbmap`); `--batch-dir` remains as a legacy fallback; read errors are fatal.
- **Test**: rule-level test (explicit inputs; 4 batches); `test_merge_annotations_empty_input` (zero files).

### PQ-10 · major · merge_annotations.py:83, 97, 115 — duplicates and unknown ids silently resolved
- **Issue**: GTDB-Tk rows were stored in a dict (a duplicate silently overwrote), CheckM2/BBMap duplicates were dropped with `drop_duplicates`, and rows whose id was not an input genome (a broken file-stem/genome_id contract, e.g. `M1.fna`) were silently ignored, leaving the genome unannotated.
- **Evidence**: pristine line 83 stores GTDB-Tk rows in a dict keyed by `user_genome` (a later duplicate overwrites), lines 97/115 `drop_duplicates(subset=[idc])` keep the first CheckM2/BBMap row, and line 120 left-merges onto the metadata, which discards tool rows whose id is not an input genome; CheckM2 1.0.2 names a bin `os.path.splitext(basename)[0]` (prodigal.py:25), so a mis-staged `M1.fna.fna` is reported as `M1.fna`.
- **Fix**: duplicated ids and ids that are not input genomes are errors for all three tools; the batch manifest must match the metadata exactly.
- **Test**: `test_merge_annotations_hard_errors` (duplicate across batches, `M1.fna` stem, manifest mismatch).

### PQ-11 · major · preprocess.smk:147 — merge rule default `length_col: scaf_bp`
- **Issue**: the rule default (used when the config key is absent) was `scaf_bp`, contradicting config/SCHEMA/script (`contig_bp`): genome_size would include scaffold gap Ns.
- **Evidence**: real BBMap run: `scaf_bp 115020` vs `contig_bp 115000` for an assembly with one 20-N gap.
- **Fix**: rule default `contig_bp`.
- **Test**: `test_rule_hygiene_and_static_contract`; `test_merge_annotations_values` (genome_size = contig_bp); rule-level (`ARC_B1` genome_size excludes its 15 gap Ns, n_contigs = records + 1).

### PQ-12 · major · batch_genomes.py:33-38,54,59-67 — genome set and ids not the ones the DAG counts
- **Issue**: (a) the metadata was read without `comment='#'`, so comment lines became genomes, while common.smk now excludes them: the `.ready` files written could differ from BATCHES; (b) duplicated ids were silently dropped; (c) ids were `utils.sanitize`d for the symlink names, but merge_annotations joined on the raw ids, so a genome whose id had any character outside `[A-Za-z0-9_.-]` silently lost all annotations and was filtered out; (d) a missing FASTA was only counted as "dangling" and the genome silently vanished downstream; (e) symlinks of genomes no longer in a batch stayed in `batch_XXXX/genomes`, and GTDB-Tk/CheckM2/BBMap process every file there.
- **Evidence**: pristine line 33 reads the metadata without `comment`, while common.smk `_count_input_genomes` uses `comment="#", skip_blank_lines=True` (common.smk:67); line 38 `drop_duplicates` drops repeated ids silently; line 54 names links `utils.sanitize(id)` while merge_annotations.py:72-74 joins on the raw id; lines 64-65 only count dangling links; lines 60-62 replace same-name links but never remove links of genomes that left the batch, and CheckM2 1.0.2 treats every file ending in `-x` in the input directory as a bin (bin/checkm2).
- **Fix**: parse exactly as `common.smk._count_input_genomes` and check `--expected-n` (= N_INPUT) so the batches are exactly BATCHES; duplicates are errors; ids must satisfy `utils.sanitize(id) == id`, start alphanumeric and contain no `__`; every FASTA must exist, be non-empty, uncompressed and start with `>` (all problems listed, then exit 1); deterministic id order; stale pipeline symlinks removed, non-symlink entries refused; manifest written.
- **Test**: `test_batch_genomes_layout_matches_common_smk` (executes common.smk's own counting code), `test_batch_genomes_fails_on_bad_inputs` (missing, gzip, empty, not FASTA, 4 invalid ids, duplicate, count mismatch, inline `#`), `test_batch_genomes_rerun_refreshes_only_its_links`, `test_batch_genomes_empty_metadata`.

### PQ-13 · major · filter_collapse.py:184-191 — staging failures ignored, real files deleted
- **Issue**: a failed symlink was only logged (Phase 2 then lacks the representative's genome); `os.remove(dst)` deleted whatever file was at the destination, including a real FASTA.
- **Evidence**: pristine lines 186-187 call `os.remove(dst)` whenever `dst` exists (also a regular file), lines 189-190 turn an `OSError` into `WARNING: could not stage`, and the source FASTA is never checked before `os.symlink` (188).
- **Fix**: a missing source FASTA is an error; a non-symlink file at the destination is never overwritten (error); only links recorded in the script's own `.staged_links` ledger are removed when a genome stops being a representative; links are staged where Phase 2 reads them (`genomes_dir` / `genome_extension`).
- **Test**: `test_filter_collapse_staging_is_safe`; rule-level (staged links resolve to the input FASTAs).

### PQ-14 · major · filter_collapse.py:132 — n_genomes_in_species depended on the aggregation mode
- **Issue**: computed over the niche-aggregation set, so with `niche_aggregation: passing` it equalled n_passing_in_species.
- **Evidence**: pristine line 128 `agg_set = keyed if niche_aggregation == "all" else keyed[keyed["passing"]]` and line 132 `n_in_species = agg_set.groupby("species_key").size()`; line 133 computes `n_passing` from the same passing set, so the two columns coincide under `passing`.
- **Fix**: n_genomes_in_species counts every genome of the species, n_passing_in_species the passing ones, in both modes.
- **Test**: `test_filter_collapse_niche_union_modes` (2 / 1 in both modes).

### PQ-15 · major · filter_collapse.py:86 — invalid niche values silently converted
- **Issue**: `to_numeric(errors='coerce').fillna(0).clip(0,1)` turned `yes`/`TRUE` into 0 and `2` into 1 without notice: presence data corrupted.
- **Evidence**: pristine line 86 `pd.to_numeric(df[c], errors="coerce").fillna(0).clip(0, 1)`: `errors="coerce"` turns a non-numeric string into NaN, `fillna(0)` into 0, and `clip(0, 1)` maps 2 to 1; nothing is logged.
- **Fix**: accept 0/1, 0.0/1.0, true/false, yes/no (case-insensitive); blank/NA -> 0 (count logged); anything else is an error.
- **Test**: `test_filter_collapse_niche_union_modes` (blank -> 0; `2` -> error).

### PQ-16 · major · filter_collapse.py:94-107 — keep_domain fragile, filter could be skipped
- **Issue**: `keep_domain: "d__Archaea"` (the GTDB form) matched nothing and silently emptied the analysis; a missing `gtdb_domain` column only warned and skipped the filter (bacteria would enter the archaeal analysis).
- **Evidence**: merge_annotations stores the domain without its prefix (`RANKS` `("d__", "gtdb_domain")`, lines 18 and 30-32 strip it), and pristine filter_collapse line 98 compares it with `keep_domain.strip().lower()`, so `d__Archaea` matches no genome; lines 105-107 only warn and skip the filter when `gtdb_domain` is absent.
- **Fix**: `Archaea` and `d__Archaea` are equivalent; the column always exists (merge_annotations); genomes outside the domain are reported with the reason.
- **Test**: `test_filter_collapse_domain_unclassified_tiers_and_files` (`d__Archaea` == `Archaea`; `domain=Bacteria`).

### PQ-17 · major · merge_annotations.py:128 / filter_collapse.py:119-125 — threshold decided by float noise
- **Issue**: thresholds must be applied exactly as stated (`quality_score > 50`), but `completeness - 5*contamination` in binary floats is off by 1e-14: completeness 64.15, contamination 2.83 gives `50.00000000000001`, which passed although the exact value is 50 (must fail).
- **Evidence**: CheckM2 1.0.2 reports two decimals (`np.round(..., 2)`), so the exact score has at most two decimals.
- **Fix**: quality_score rounded to 1e-6 in merge_annotations and recomputed (with a warning on disagreement) in filter_collapse; comparisons unchanged (`>=`, `<=`, `<`, `>`, `>`).
- **Test**: `test_filter_collapse_thresholds_are_exact` (each threshold at and just past its boundary; 64.15/2.83 fails; 75.01/5.00 passes); `test_merge_annotations_values`.

### PQ-18 · major · rename_prodigal.py:30-48 — gene ids not Contract v2, no coordinate table, silent CDS drops
- **Issue**: ids were `{genome}_{N}` (genome-wide counter): no contig, so DefenseFinder gembase replicons and gene adjacency were impossible; no genes.tsv; CDS records without a matching protein were silently dropped.
- **Evidence**: pristine line 34 `nid = f"{g}_{i}"` (one counter per genome), lines 42-46 skip a CDS without a matching protein (`missing += 1; continue`) with only a warning (50-52), and no coordinate table is written; Contract v2 requires `{genome}_{contig#}_{gene#}` via `utils.make_gene_id`; Prodigal 2.6.3 headers carry `ID={k}_{n}` with k = record index and n = gene index on that record (main.c/sequence.c and the real run).
- **Fix**: ids `utils.make_gene_id(genome, k, n)` from Prodigal's `ID=k_n` (k = record index, n = gene index on the record), identical in `.faa`, `.fna`, GFF `ID=` and `genes.tsv` (`protein_id, contig, contig_index, gene_index, start, end, strand, partial, contig_length, length_aa` + `transl_table`). Proteins, CDS and GFF features must agree one-to-one, with equal coordinates and `len(CDS) = 3 x len(protein)`; otherwise exit 1. `genes.tsv` is a declared output, exactly `results/genes/{genome}/{genome}.genes.tsv`.
- **Test**: `test_prodigal_real_run_contract_v2_ids_and_genes_tsv` (real Prodigal on a genome with 21 planted genes on 4 records incl. a gene-less record and a right-edge partial: every row equals the design; CDS translate to the proteins; split/make_gene_id round-trip); rule-level test (same through the rule); `test_rename_prodigal_failure_modes`.

### PQ-19 · major · filter_collapse.py / merge_annotations.py — genomes GTDB-Tk/CheckM2 did not report were invisible
- **Issue**: genomes GTDB-Tk reported as `Unclassified*`, genomes absent from its output and genomes CheckM2 skipped simply disappeared at the domain or QC filter with no record.
- **Evidence**: GTDB-Tk 2.4.0 classify.py `add_filtered_genomes_to_summary` / `add_failed_genomes_to_summary`; CheckM2 1.0.2 skips empty bins / empty protein files.
- **Fix**: annotated_metadata.tsv gains `gtdbtk_status` (classified / unclassified / missing), `checkm2_status`, `bbmap_status` and GTDB-Tk/CheckM2 provenance columns; new `preprocess/filter_report.tsv` gives every input genome a status (`representative`, `species_member`, `failed_quality`, `no_species`, `excluded_domain`), the reasons, and its species representative; counts and warnings logged; the funnel plot uses it.
- **Test**: `test_merge_annotations_values`, `test_filter_collapse_domain_unclassified_tiers_and_files`, rule-level test (exact status of all 7 genomes, incl. one absent from GTDB-Tk output).

### PQ-20 · minor · qc.smk:30-32 — sub-20-kb fallback test did not match Prodigal's count
- **Issue**: `grep -v '^>' | tr -d '\n' | wc -c` counted `\r`, spaces and other characters and ignored the 12-bp spacers Prodigal adds between records, so CRLF files near 20 kb could choose `single` and hit Prodigal's exit 10 (then swallowed, PQ-05). Under `set -euo pipefail` a header-only file aborted the rule at the `grep` with no message.
- **Evidence**: Prodigal 2.6.3 sequence.c `read_seq_training`; boundary reproduced with the real binary (19,975 letters in 2 records -> slen 19,999 -> exit 10; 19,976 -> slen 20,000 -> ok).
- **Fix**: an awk count of sequence letters + 12 bp per record when more than one record (exactly Prodigal's slen for valid FASTA; it can only err toward `meta`), invalid `qc.prodigal_mode` rejected, decision logged.
- **Test**: `test_prodigal_training_length_matches_prodigal` (runs the rule's own awk; real Prodigal at both sides of the boundary); rule-level (`ARC_E1`, 7.7 kb -> `mode=meta`).

### PQ-21 · minor · qc.smk:33 — translation table implicit; meta mode ignores it
- **Issue**: no `-g`; archaea use table 11 (Prodigal's default, but unpinned in the rule). In `-p meta` Prodigal never applies `-g` and may pick a table-4 model for a record.
- **Evidence**: main.c meta loop over `meta[i].tinf` (`user_tt` only used for training files).
- **Fix**: explicit `-g {qc.prodigal_translation_table, default 11}`; the table Prodigal actually used per record (GFF `# Model Data ... transl_table=`) is written to `genes.tsv` (`transl_table`) and a non-11 table is logged as a warning. `-m` (mask N runs) was tested and not adopted: it drops a gap-spanning gene entirely instead of truncating it.
- **Test**: `test_rule_hygiene_and_static_contract`; genes.tsv assertions (`transl_table` 11).

### PQ-22 · minor · filter_collapse.py:40 — empty representatives.txt began with "genome_id"
- **Issue**: `utils.write_empty(out_reps_list, "genome_id")` made the header-less list (consumers: build_master_table `load_reps`, validate_prometheus `--reps`) contain a genome called `genome_id`.
- **Evidence**: pristine line 40 `utils.write_empty(args.out_reps_list, "genome_id")`; `utils.write_empty` writes the header line (utils.py:63-67); pristine build_master_table.py `load_reps` (51-59, "one genome_id per line") adds every non-blank, non-`#` line, so the empty case produced a representative named `genome_id`. The current consumers (edited by other leaves) skip a first-line `genome_id`; an empty file is correct for both versions.
- **Fix**: an empty file when there is no representative.
- **Test**: `test_filter_collapse_empty_input`, `test_snakemake_empty_input_set`.

### PQ-23 · minor · filter_collapse.py:181 — quality_tier hard-coded to "A"
- **Issue**: `quality_tier: "A"`, `pass_primary: True` for every representative, wrong as soon as the Phase-1 thresholds are relaxed.
- **Evidence**: pristine line 181 writes the constants `"quality_tier": "A", "pass_primary": True`; config.yaml:75-78 defines `qc.completeness_floor` 50, `qc.contamination_ceiling` 10, `qc.tier_high_completeness` 90 ("MIMAG high-quality threshold for tier A genomes") and `qc.tier_high_contamination` 5, which filter_collapse ignored.
- **Fix**: MIMAG completeness/contamination tiers from `qc.tier_high_*` (A), `qc.completeness_floor/contamination_ceiling` (B), else C; `pass_primary` from the floor/ceiling. Identical to the old values under the default thresholds (all A/True).
- **Test**: `test_filter_collapse_domain_unclassified_tiers_and_files` (A/B/C and True/True/False).

### PQ-24 · minor · filter_collapse.py:140-147 — representative tie-break
- **Issue**: ties were broken by completeness then genome_id only. Deterministic, but not the most defensible choice among equal-score genomes.
- **Evidence**: pristine lines 140-141 sort by `["species_key", "quality_score", "completeness", "genome_id"]`; with equal quality_score and completeness the contamination is equal too (score = completeness - 5 x contamination), so assembly contiguity (N50, contig count) was ignored and the smallest id won.
- **Fix**: quality_score desc, completeness desc, contamination asc, contig N50 desc, contigs asc, genome_id asc, stable sort.
- **Test**: `test_filter_collapse_representative_tie_breaks` (same answer for three input orders).

### PQ-25 · minor · merge_annotations.py:130-149 — "KB" heuristic instead of a real check
- **Issue**: a warning heuristic (`max(n50) < 1000 and median size > 100 kb`) guessed at unit problems.
- **Evidence**: AssemblyStats2.java `formatKB`: plain integer bp whenever FORMAT >= 2, so format=6 can never emit KB/MB.
- **Fix**: replaced by the strict numeric parse, the `n50 <= genome_size` check and the CheckM2 cross-check (PQ-08).
- **Test**: `test_merge_annotations_hard_errors`.

### PQ-26 · minor · batch_genomes.py / merge_annotations.py:72 — a mid-line '#' truncates a metadata row
- **Issue**: the table is read with `comment='#'` (as in common.smk), which drops everything after a `#` anywhere in a line, silently truncating filenames or the niche columns of that row. merge_annotations also read the metadata without the comment handling, so comment lines became genomes.
- **Evidence**: pandas 2.2.2 `read_csv(..., comment="#")` on the data row `G1<TAB>g1#v2.fna<TAB>1` returns `filename = "g1"` and `niche_gut = NaN` (run here); common.smk reads the metadata with `comment="#"` (common.smk:67, 90, 127); pristine merge_annotations.py:72 reads it without `comment`, so a comment line became a genome there.
- **Fix**: both scripts parse like common.smk and reject a `#` inside a data row (whole-line comments stay allowed).
- **Test**: `test_batch_genomes_fails_on_bad_inputs` (case `hash`), `test_merge_annotations_values` (comment and blank lines excluded).

### PQ-27 · minor · plot_preprocess.py:40,55,59-61 — "pass" did not mean pass
- **Issue**: the scatter coloured genomes by completeness/contamination only (not the 5-threshold filter); the funnel had no domain/species stages; `ylim(0, max)` was degenerate when every contamination was 0.
- **Evidence**: pristine line 40 defines `passing` from completeness and contamination only (filter_collapse also applies n_contigs, N50 and quality_score); line 55 `set_ylim(0, min(50, nanmax(cont)))` is `(0, 0)` when every contamination is 0, which matplotlib 3.8.4 only warns about and auto-expands (axis.py:1236, "Attempting to set identical low and high ylims"); lines 59-61 draw three bars (input, pass, species reps).
- **Fix**: colour by `filter_report.tsv` `passing` (fallback labelled as completeness/contamination only); funnel input -> GTDB-Tk classified -> domain -> species -> pass QC -> representatives; robust axis limits with the count of points above the limit and of genomes without estimates; configurable formats.
- **Test**: `test_plot_preprocess_values_formats_and_empty` (legend counts, funnel, formats, placeholder).

### PQ-28 · minor · plot_qc.py:70-73,97-114 — hard-coded thresholds, forbidden label, unit guesses
- **Issue**: MIMAG lines hard-coded at 50/90 and 5/10 regardless of config; colouring by `genome_type` (README/SCHEMA exclude MAG/isolate labels); genome size divided by 1e6 only if `max > 1e4` (axis always said Mb); BBMap `gc_avg` (a fraction) labelled "GC content" without unit.
- **Evidence**: pristine lines 70-73 draw literal 50/90 and 5/10 lines, while config.yaml:75-78 makes these thresholds configurable; lines 97-101 colour by `genome_type`, which README.md:208/267, SCHEMA.md:64 and samples.schema.yaml:5-7 exclude (the design note section 0.1 proposes recording isolate vs MAG/SAG, but README/SCHEMA settle on not using it); line 111 `x / 1e6 if x.max() > 1e4 else x` under a fixed "Mb" label; line 114 labels BBMap `gc_avg`, a fraction (AssemblyStats2.java; real run 0.50070), "GC content".
- **Fix**: lines and legend from `qc.*`; tier colours only; size always bp -> Mb; GC converted to %; n and median per panel; configurable formats; legend outside the axes (it overlapped the floor line).
- **Test**: `test_plot_qc_units_thresholds_and_empty` (legend text, `GC content (%)`, medians 41 % and 1.8 Mb, placeholders).

### PQ-29 · minor · plot_qc_extra.py:47-63 — degenerate groups and a mislabelled fallback
- **Issue**: `violinplot` was called for groups with one value or a constant value (a kernel density is undefined there); the defense panel said "excl. Prometheus" even when it fell back to `n_defense_systems`; genome size in bp.
- **Evidence**: pristine lines 62-63 call `violinplot` whenever both groups are non-empty; matplotlib 3.8.4 `Axes.violinplot` `_kde_method` (axes/_axes.py:8252-8257) returns the indicator `(X[0] == coords)` for a constant vector instead of a density, so a single-genome or constant group is drawn as a zero-width spike (no error); lines 47-50 keep the label "(excl. Prometheus)" when falling back to `n_defense_systems`; line 49 labels genome size in bp.
- **Fix**: violins only for groups with >= 2 distinct values, otherwise points + median; label states incl./excl.; Mb; n per group. CLI unchanged.
- **Test**: `test_plot_qc_extra_degenerate_groups` (single-genome and constant groups; placeholder without presence variation).

### PQ-30 · minor · preprocess.smk:220, qc.smk:44,62 — rules without log:
- **Issue**: plot_preprocess, genome_metadata and plot_qc had no `log:` (Contract hygiene).
- **Evidence**: pristine plot_preprocess (preprocess.smk:220-242), genome_metadata (qc.smk:44-59) and plot_qc (qc.smk:62-78) have no `log:` directive; the driver's tests/check_rules.py requires one in every rule.
- **Fix**: logs added; every shell writes to its log.
- **Test**: G2 `tests/check_rules.py` -> `RULES HYGIENE OK`.

### PQ-31 · minor · preprocess.smk:220-241, qc.smk:62-78 — plot outputs vs PLOT_FORMATS
- **Issue**: the rules declare `plotting.formats` outputs, but the scripts always wrote png + svg; any other config produced missing outputs.
- **Evidence**: outputs are `expand(..., ext=PLOT_FORMATS)` (preprocess.smk:226, qc.smk:66-67), but the shell passes no formats and `plotting_style.save_fig` defaults to `formats=("png", "svg")` (plotting_style.py:79).
- **Fix**: `--formats` passed from PLOT_FORMATS.
- **Test**: `test_plot_preprocess_values_formats_and_empty` (`--formats png` -> only png).

### PQ-32 · minor · aggregate_qc.py — dead code removed
- **Issue**: no rule called it; it implemented the removed seqkit/barrnap/tRNAscan/dRep-genomeInfo path, used the forbidden `genome_type`, and defaulted a missing CheckM2 result to 0/100 (a silent failure).
- **Evidence**: no rule, Snakefile or `.smk` file of the pristine tree references `aggregate_qc` (grep); aggregate_qc.py:12 writes a `genome_type` column, and lines 103, 106 and 153 default completeness/contamination to 0.0/100.0 when the CheckM2 report is unreadable, unparsable or lacks the genome.
- **Fix**: deleted (moved to Trash; still in the pristine copy). Wiring it in would duplicate genome_metadata.tsv, which filter_collapse already produces.
- **Test**: `test_py_compile_owned_scripts` asserts its absence.

### PQ-33 · minor · preprocess.smk gtdbtk_batch — stale GTDB-Tk state reused
- **Issue**: GTDB-Tk 2.4 reloads a completed ANI screen from `out_dir/gtdbtk.json`; after a batch layout change (sorted-id batching shifts genomes between batches) that state would be stale when the ANI screen is enabled.
- **Evidence**: gtdbtk/main.py classify_wf ("The ani_screen step has already been completed, we load existing results").
- **Fix**: the rule removes its own GTDB-Tk and CheckM2 output directories before running.
- **Test**: rule-level test (fresh directories; reruns in the failure scenarios).

### PQ-34 · minor · preprocess.smk:96-108 — CheckM2 intermediates
- **Issue**: CheckM2 keeps per-genome protein files and DIAMOND output (both grow with the batch size); only `quality_report.tsv` is used.
- **Evidence**: CheckM2 1.0.2 bin/checkm2:112 defines `--remove_intermediates` ("Remove all intermediate files (protein files, diamond output)"); the pristine rule (preprocess.smk:102-103) does not pass it and copies only `quality_report.tsv`.
- **Fix**: `--remove_intermediates` by default via `cfg('preprocess/checkm2_extra', '--remove_intermediates')` (flag verified in bin/checkm2 1.0.2).
- **Test**: rule-level (stand-in accepts exactly the 1.0.2 flags).

### PQ-35 · minor · preprocess.smk:125 — bbmap_batch reserved 8 threads
- **Issue**: statswrapper runs one single-threaded JVM per chunk; `T_MED` reserved 8 slots for nothing.
- **Evidence**: BBMap 39.06 statswrapper.sh:69 runs one `java ... jgi.AssemblyStatsWrapper` process; AssemblyStatsWrapper.java and AssemblyStats2.java start no worker threads (the only `Thread` call is `Thread.yield()`, AssemblyStatsWrapper.java:61), and run_bbmap_stats.py runs the chunks one after another (one `subprocess.run` per chunk).
- **Fix**: `T_SMALL`.
- **Test**: `test_rule_hygiene_and_static_contract` asserts `threads: T_SMALL` in bbmap_batch.

### PQ-36 · minor · batch_genomes.py — batch paths with spaces
- **Issue**: GTDB-Tk 2.4.0 rejects genome paths containing a space (main.py `_verify_file_path`), which would surface only after the job started.
- **Evidence**: GTDB-Tk 2.4.0 main.py:137-139 `_verify_file_path` raises `GTDBTkExit` when a genome path contains a space; the batchfile lists the staged links under the batch directory.
- **Fix**: batch_genomes fails immediately with that explanation.
- **Test**: `test_batch_genomes_fails_on_bad_inputs` (outdir `with space/b`: fails with the GTDB-Tk explanation and stages nothing).

## Verified correct, no change needed
- The five Phase-1 inequalities in filter_collapse were already `>=`, `<=`, `<`, `>`, `>` as the config comments state; kept (now with exact scores, PQ-17).
- Representative selection already took the whole top row (`drop_duplicates` after sorting, not `GroupBy.first`); kept.
- Niche union semantics (OR = max of 0/1 over the species' genomes; `all` by default) match the design note; kept (validated, PQ-15).
- CheckM2 `Name` = file stem, so `{genome_id}.fna` symlinks give `Name == genome_id` (checkm2/prodigal.py:25), also for ids containing dots; GTDB-Tk `remove_extension` likewise; BBMap prints the input name, reduced to its basename.
- BBMap mapping `gc=gc_avg`, `genome_size=contig_bp`, `n_contigs=n_contigs`, `n50=ctg_L50` (config, SCHEMA) is right: ctg_L50 is the contig N50 length in bp (real run).
- GTDB-Tk `--skip_ani_screen` (config default) is still valid and required without `--mash_db` in 2.4.0.
- Output formats consumed by other stages are unchanged: `taxonomy.tsv` (`user_genome, classification`; root_tree substring match on the classification), `representatives.txt` (one id per line, no header), `genome_metadata.tsv` (GENOME_META_HEADER, same columns and order), `species_metadata.tsv` columns (additive only).
- `validate_prometheus.py` (another leaf) reads genes.tsv with `dtype=str, keep_default_na=False` and uses `protein_id, contig_index, gene_index, partial`, so the 2-character partial codes survive; `extract_bacterial_prometheus.py` writes the bacterial genes.tsv with the same columns and `+`/`-` strands, so the formats agree across domains.
- `.test` dry runs: `bash tests/check_dag.sh` on a private copy of the repo with the final files -> `DAG CHECK OK` (`preprocess` 11 jobs, `all` before Phase 1 11 with no Phase-2 jobs, `all` 103, `crossdomain` 86, `disorder` 43; unsupported `plotting.formats` rejected at parse time).

## Not verifiable here (see `unresolved` in the leaf report)
- GTDB-Tk and CheckM2 were not executed (reference data 10-100 GB); their interfaces were verified from source at the pinned tags and through format-faithful stand-ins.
- BBMap 39.06 has no macOS build; it was run from the unpacked bioconda linux-64 build (statswrapper.sh is pure Java; the JNI library is not loaded) with conda-forge OpenJDK 17. BBMap is not in the shared toolchain, so the real-BBMap comparison test runs only when `statswrapper.sh` is on the toolchain PATH or `PROM_BBMAP_STATSWRAPPER` is set; the leaf sandbox cannot write the package into its conda env, so the gate run (G1) parses the captured real output instead.
- A direct matplotlib experiment on degenerate violins was not executed (the platform safety check could not run and asked for user approval); the behaviour was read from the matplotlib 3.8.4 source instead (PQ-29 evidence), and the suite exercises those inputs.
