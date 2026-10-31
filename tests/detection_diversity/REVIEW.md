# Review: Stage 1 detection/validation + Stage 2 diversity (leaf-1.1.2, `detection_diversity`)

Reviewed commit: `56f5437` (github.com/marshomics/prometheus_manuscript). Line numbers in
**Location** refer to that commit. Fixes are in the working tree; every fix is covered by
`tests/detection_diversity/run_tests.py` (26 tests; `bash tests/run_leaf.sh detection_diversity`:
about 136 s on the first run, which downloads the pinned models, and 103 s with the model cache
warm; at most 2 threads per tool).

Scope: `workflow/rules/{prometheus,diversity}.smk`, `workflow/envs/defensefinder.yaml`, and the
scripts `check_defensefinder.py`, `combine_proteomes.py`, `parse_defensefinder.py`,
`validate_prometheus.py`, `plot_prometheus_validation.py`, `cluster_diversity.py`,
`architecture_typology.py`, `plot_diversity.py`, `plot_diversity_extra.py`,
`R/plot_validation_tree.R`.

Severity scale (as in the review brief): **critical** = silently wrong biological result for
many or all genomes; **major** = wrong, lost or unreproducible result under realistic
conditions, or a crash converted into a plausible-looking output; **minor** = wrong only in
edge cases, a misleading output, or hygiene and provenance.

Real tools used as evidence and in the tests: DefenseFinder 3.0.0 + macsyfinder 2.1.4
(env `prom-detection_diversity-defensefinder`), defense-finder-models 3.1.0 + CasFinder 3.1.0
(downloaded from codeload.github.com by the test and installed with `macsydata install --target`),
HMMER 3.4, MAFFT 7.526, ClipKIT 2.14.0, IQ-TREE 2.3.6, MMseqs2 18.8cc5c, R 4.3.3 with ape 5.8.1 and
ggtree 3.10.0.

**Findings: 33** (critical 3, major 20, minor 10).

---

## A. DefenseFinder detection (rules `defensefinder_check`, `defensefinder`; `check_defensefinder.py`; `envs/defensefinder.yaml`)

### DD-01 [critical] A crashed DefenseFinder run became a "no Prometheus" genome
- **Location:** `workflow/rules/prometheus.smk:62` (`defense-finder run ... > {log} 2>&1 || true`), `:68-69` (header-only placeholders written whenever the output is missing).
- **Issue:** any DefenseFinder failure (missing or incompatible models, offline version check, out-of-memory) produced header-only tables. The genome then counted as Prometheus-negative in prevalence and every downstream association, with no error.
- **Evidence:** with the model set that `defense-finder update` installs today (CasFinder 3.1.1, 44 definitions with grammar `vers="2.1"`), macsyfinder 2.1.4 aborts: its `definition_parser` accepts only `"2.0"` ("version supported is '2.0'"; reproduced here). Under `|| true` every genome would have been scored negative.
- **Fix:** removed `|| true`. Placeholders are written only for a genuinely empty proteome (0-byte `.faa`), using DefenseFinder's own headers, and each placeholder line carries a `# nonfatal:` justification (Contract v2 item 5). Every table, run or ingested, is then verified by `check_defensefinder.py --check-outputs`.
- **Test:** `defensefinder_rule_ingest_mode`, `check_outputs_and_ingest_copy`, `real_defensefinder_preflight_run_parse_context`.

### DD-02 [major] Multi-contig genomes were searched as one replicon
- **Location:** `workflow/rules/prometheus.smk:62` (no `--db-type`).
- **Issue:** DefenseFinder `run` defaults to `--db-type ordered_replicon`, which treats the whole proteome as one ordered replicon. The last gene of contig *i* and the first gene of contig *i+1* become neighbours, so systems and gene-distance context can span contig junctions in MAGs and drafts.
- **Evidence:** DefenseFinder v3.0.0 `cli.py` (`--db-type` default `ordered_replicon`). In macsyfinder 2.1.4 gembase mode, the replicon is the sequence id up to its last `_`. Contract v2 ids `{genome}_{contig#}_{gene#}` therefore give replicon `{genome}_{contig#}`. A real run on a two-contig genome produced replicons `GX_1` and `GX_2`.
- **Fix:** `--db-type gembase`. `--check-outputs --expect-gembase` requires every `replicon` to equal its `hit_id` minus the gene field.
- **Test:** `check_outputs_and_ingest_copy` (replicon `G1` for hit `G1_1_1` is rejected), `real_defensefinder_preflight_run_parse_context`.

### DD-03 [major] Every DefenseFinder run depended on the GitHub API
- **Location:** `workflow/rules/prometheus.smk:62` (no `--skip-model-version-check`).
- **Issue:** without the flag, DefenseFinder looks up the latest models version on each run. The lookup uses `check_last_version_models()` → macsypy `remote.list_package_vers`, with a 30-day cache file in `$HOME/.defensefinder_model_lastversion`. That is a network dependency on compute nodes. It also does not pin anything: the run continues with whatever models are installed and only logs a warning.
- **Evidence:** DefenseFinder v3.0.0 `cli.py` lines 29-49 and 225-243 (archived source). GitHub API requests from this sandbox were rate-limited (HTTP 403) during the review.
- **Fix:** `--skip-model-version-check` added. The models version is pinned and verified once by the preflight (DD-05) instead.
- **Test:** `real_defensefinder_preflight_run_parse_context` (runs the rule's command line, which includes the flag).

### DD-04 [major] Floating detector pins; documented update path breaks the run
- **Location:** `workflow/envs/defensefinder.yaml:10-11` (`defense-finder>=2.0.0`, `macsyfinder>=2.1.4`), `:17` (recommends `defense-finder update`). The models were not pinned anywhere.
- **Issue:** the detector version, and therefore every Prometheus call, depended on install date. `defense-finder update` installs the latest models and CasFinder; today that includes CasFinder 3.1.1, which macsyfinder 2.1.4 cannot parse (DD-01).
- **Evidence:** bioconda provides defense-finder 2.0.0, 2.0.1 and 3.0.0 (noarch); the DefenseFinder setup requires `macsyfinder==2.1.4`. Tags: DF v3.0.0 2026-05-04; models 3.1.0 2026-05-14. `Prometheus.xml` (md5 `5eee0d48274738866021ff589b12ddaf`) and `Prometheus__ProA.hmm` (md5 `83b327e1edff1df8601e0dd98bca55a7`; LENG 1668, GA 600, NSEQ 505) are byte-identical in models 2.0.0+breakingchanges, 2.0.2 and 3.1.0. The pin therefore does not change the Prometheus model relative to the version the project was designed on.
- **Fix:** exact pins `defense-finder=3.0.0`, `macsyfinder=2.1.4`, `hmmer=3.4`, `pandas=2.2.3`. Models are pinned through rule config defaults: `prometheus/defensefinder_models_pin` "3.1.0", `defensefinder_casfinder_pin` "3.1.0", optional `defensefinder_models_dir`. The preflight enforces them. Install commands are documented in the env file (`macsydata install --target ... 'defense-finder-models==3.1.0'` and `'CasFinder==3.1.0'`). The linux-64 solve resolves defense-finder 3.0.0 `pyhdfd78af_0`, macsyfinder 2.1.4, hmmer 3.4, pandas 2.2.3.
- **Test:** `preflight_version_and_model_checks` (software, models and CasFinder pin mismatches all fail); ledger gate G3 (`check_env_solve.sh`).

### DD-05 [major] The preflight could pass without any usable model set
- **Location:** `workflow/scripts/py/check_defensefinder.py:62` (tries `defense-finder --version`, which is not an option; only the fallback `version` subcommand works), `:94-123` (looks for models under `sys.prefix/share`, not macsyfinder's search path), `:182-191` (only a WARNING when models are not found, then passes). The models version was never compared (`--models-version` appeared only in the message). CasFinder was not checked.
- **Issue:** the preflight never established that the installed DefenseFinder could call Prometheus.
- **Evidence:** DefenseFinder `cli.py` (`version` subcommand prints "Using DefenseFinder version X"; `run` exits with "Uncomplete defense-finder models" unless both defense-finder-models and CasFinder are installed). macsypy 2.1.4 searches `--models-dir`, else `$VIRTUAL_ENV/share/macsyfinder/models` or the last existing of `/usr/share` and `/usr/local/share`, plus `~/.macsyfinder/models`.
- **Fix:**
  - Software version via `defense-finder version`, against the exact pin and the 2.0.0 minimum.
  - Models resolved through macsypy's own registry (`_find_all_installed_packages`) and the same directory list.
  - Exact models and CasFinder pins; CasFinder required.
  - `Prometheus.xml` must declare `Prometheus__ProA` as mandatory, and its profile must exist.
  - All definitions must use grammar 2.0.
  - **Positive control:** DefenseFinder must call Prometheus on the `hmmemit -c` consensus of `Prometheus__ProA.hmm` (1668 aa) before any genome is searched.
  - Key/value provenance sentinel: versions, paths, md5s, profile header.
  - Ingest mode records the provided outputs as UNVERIFIED.
  - New optional `--threads` sets the positive control's DefenseFinder workers; the rule passes `{threads}` (30 s with 1 worker, 16 s with 2).
- **Test:** `preflight_version_and_model_checks` (pins, models < 2.0.0, CasFinder 3.1.1, wrong gene, accessory presence, missing profile, missing CasFinder, missing binary), `preflight_ingest_mode`, `real_defensefinder_preflight_run_parse_context` (real preflight about 45 s: "called Prometheus"; grammar `{'2.0': 666}` / `{'2.0': 44}`).

### DD-06 [major] Ingest mode silently replaced a missing provided table with an empty one
- **Location:** `workflow/rules/prometheus.smk:64-65` (`cp {params.provided}/{genome}/*_defense_finder_*.tsv ... 2>/dev/null || true`, then the placeholder at `:68-69`).
- **Issue:** a misnamed or missing pre-computed DefenseFinder output became a Prometheus-negative genome. Stale outputs from a different gene calling were accepted unchecked.
- **Evidence:** in bash, `2>/dev/null || true` discards both the `cp` error and its exit status. The next line (`[ -s {output} ] || printf ...`) then writes the header-only placeholder, so the rule exits 0 with an empty table.
- **Fix:** `check_defensefinder.py --ingest-copy` takes the exact file name or exactly one glob match, and errors on zero or several. The copied tables are then verified against this genome's proteome (`hit_id` ⊆ `.faa` ids, `sys_id` consistency, v2+ model paths).
- **Test:** `check_outputs_and_ingest_copy` ("found 0", "found 2"), `defensefinder_rule_ingest_mode` (a proteome that no longer matches → "absent from").

### DD-07 [minor] Placeholder headers did not match DefenseFinder's tables
- **Location:** `workflow/rules/prometheus.smk:68-69`.
- **Issue:** the systems placeholder lacked `activity` and `genes_count`. The genes placeholder had 9 of the 25 columns. Consumers keyed on DefenseFinder columns would break, or behave differently, only for empty genomes.
- **Evidence:** DefenseFinder posttreat (`df_systems.py`, `df_genes.py`) and the macsypy 2.1.4 `TsvSystemSerializer` header plus `type`/`subtype`/`activity`.
- **Fix:** `DF_SYSTEMS_HEADER` / `DF_GENES_HEADER` constants with the exact headers. `--check-outputs` requires the v2 columns.
- **Test:** `check_outputs_and_ingest_copy` ("lacks DefenseFinder v2 columns").

### DD-08 [minor] DefenseFinder index written into another stage's directory
- **Location:** `workflow/rules/prometheus.smk:62` (no `--index-dir`).
- **Issue:** the sequence index was written next to the input `.faa` under `results/genes/{g}/`, an output directory owned by the gene-calling stage.
- **Evidence:** macsyfinder 2.1.4 `--index-dir` help (`macsypy/scripts/macsyfinder.py:402-405`) describes it as the place to store the index when the sequence-db directory is not writable. The default location is therefore the sequence-db directory (`macsypy/database.py` `_index_dir`).
- **Fix:** `--index-dir {outdir}/index`.
- **Test:** `real_defensefinder_preflight_run_parse_context` (the rule command line with `--index-dir`).

## B. Candidate assembly (`combined_proteome`, `assemble_candidates`; `combine_proteomes.py`, `parse_defensefinder.py`)

### DD-09 [minor] Substring system matching, and every gene of a system taken as ProA
- **Location:** `workflow/scripts/py/parse_defensefinder.py:83` (`if sys_lower in t or sys_lower in st`), `:104-124` (all genes of a matching system become candidates).
- **Issue:** any type or subtype containing the configured string matched, and accessory genes of a multi-gene model would have been treated as ProA. This is latent with the current single-gene Prometheus model, but the configured system name is a free parameter.
- **Evidence:** the pristine code path. In models 3.1.0, Prometheus is a single-gene model (`Prometheus.xml`: one mandatory gene, `Prometheus__ProA`), so the defect is latent today. The test fixture shows that a `NotPrometheusLike`/`Prometheus_relative` system and a `Prometheus__Accessory` gene would both have been taken.
- **Fix:** exact equality on `type` or `subtype` (case-insensitive). Only hits whose `gene_name` equals `--gene` (default `{system}__ProA`; empty string = all genes) are used. Additive columns: `df_sys_id`, `df_gene_name`, `df_i_eval`, `df_profile_cov`, `df_seq_cov`.
- **Test:** `parse_defensefinder_exact_matching` (`Prometheus_relative` and `Prometheus__Accessory` excluded; `--gene ""` keeps both Prometheus genes).

### DD-10 [major] Genome universe from the filesystem; silent id and read fallbacks
- **Location:** `workflow/scripts/py/parse_defensefinder.py:58` (`os.listdir(df_dir)`), `:39-41` (`pmap.get(key, make_protein_id(...))`), `:70-75`, `:92-97`, `:121` (unreadable tables warned and skipped), `:183` (candidates without a sequence silently dropped from `candidates.faa` but kept in `candidate_sources.tsv`).
- **Issue:**
  - Stale genome directories from earlier runs were counted.
  - Hits whose protein is not in the pipeline proteome produced source rows without sequences.
  - An unreadable table made a genome silently negative.
- **Evidence:** the pristine code path. Snakemake does not delete the output directories of genomes later removed from the sample table, so `os.listdir` still sees them.
- **Fix:** genomes come from `--samples`. A missing table, an id absent from the proteome and an unreadable file are all errors. `--strict` keeps the "ZERO systems in every genome" sanity check.
- **Test:** `parse_defensefinder_failures_and_empty`, `parse_defensefinder_exact_matching` (the `STALE` directory is ignored).

### DD-11 [major] Missing per-genome proteomes were skipped silently
- **Location:** `workflow/scripts/py/combine_proteomes.py:25-35` (unreadable or empty sample table → empty outputs), `:43` (`utils.read_fasta` returns nothing for a missing path, so the genome vanished).
- **Issue:** a missing per-genome proteome, or an unreadable sample table, produced an incomplete or empty combined proteome without an error.
- **Evidence:** `utils.read_fasta` returns immediately, without error, for a missing or empty path (`workflow/scripts/py/utils.py:20-22` at 56f5437).
- **Fix:** a missing `.faa`/`.fna`, an unreadable sample table, a duplicated `genome_id` and a duplicated protein id are errors. An empty per-genome file (no genes) contributes 0 records. CDS mode (`--suffix fna`, used by `evolution.smk combined_cds`) is unchanged. New **reference mode** builds the validation-tree input (DD-16).
- **Test:** `combine_proteomes_proteome_mode`, `combine_proteomes_reference_mode`.

## C. Domain architecture (`prometheus_domain_scan`)

### DD-12 [critical] An unpressed Pfam silently removed every domain
- **Location:** `workflow/rules/prometheus.smk:142-144` (`hmmscan --cut_ga ... || true`).
- **Issue:** with an unpressed `Pfam-A.hmm`, hmmscan exits 1 and writes no domtbl. Every candidate then had no domains, so architecture, catalytic status and tiers were degraded for the whole catalogue without an error.
- **Evidence:** HMMER 3.4, run here: an unpressed database gives exit 1 with "use hmmpress first" and no `--domtblout`. A pressed database reproduces the fixture domtbl exactly (62 domain hits).
- **Fix:** the rule requires `.h3m/.h3i/.h3f/.h3p` and has no `|| true`. `--notextw -o domains.hmmscan.txt` keeps the per-domain alignments needed for residue-level catalytic evidence (DD-17). Zero candidates gives a comment-only domtbl.
- **Test:** `domain_scan_rule_fails_on_unpressed_pfam`, `real_hmmscan_matches_fixture_and_motifs`.

## D. Validation tree (`prometheus_validation_tree`)

### DD-13 [major] Sequence count broke on zero matches
- **Location:** `workflow/rules/prometheus.smk:175` (`n=$(grep -c '>' all.faa || echo 0)`).
- **Issue:** with no match, `grep -c` prints `0` and exits 1, so `n` becomes `"0\n0"` and the integer test errors.
- **Evidence:** grep exits 1 when no line is selected, while `-c` still prints `0` (POSIX grep). `$(grep -c ... || echo 0)` therefore yields two lines.
- **Fix:** `awk` counts for candidates and for all sequences. A `();` placeholder only for 0 candidates or fewer than 4 sequences (no testable internal edge).
- **Test:** `validation_tree_rule_real_tools` (zero-candidate branch).

### DD-14 [major] ClipKIT and IQ-TREE failures produced "no tree", and every candidate became Unplaced
- **Location:** `workflow/rules/prometheus.smk:178` (`clipkit ... || cp all.aln all.trim.aln`), `:179-180` (`iqtree2 ... || true`).
- **Issue:** ClipKIT exits 0 on its own errors, so the fallback never ran; IQ-TREE then failed, and `|| true` hid it.
- **Evidence:** ClipKIT v2.3.0 `args_processing.py` prints "Input file does not exist" and calls `sys.exit()` (status 0); confirmed with the local ClipKIT.
- **Fix:** explicit non-empty output checks after MAFFT, ClipKIT and IQ-TREE; no fallbacks.
- **Test:** `validation_tree_rule_real_tools`.

### DD-15 [major] Placement tree had no branch supports and no seed
- **Location:** `workflow/rules/prometheus.smk:179` (`iqtree2 -m LG+G -fast`).
- **Issue:** `-fast` gives no supports, so placement could not be judged statistically, and runs were not reproducible. UFBoot cannot simply be added: IQ-TREE 2.3.6 rejects `-B` with `-fast` ("Ultrafast bootstrap does not work with -fast, -te or -n option", exit 2).
- **Evidence:** IQ-TREE 2.3.6, run here: `-fast -B 1000` exits 2 with the message quoted above. With no `--seed`, repeated runs can return different trees.
- **Fix:** `-m {prometheus/validation_tree_model, default MFP} -B {tree/ufboot} --alrt {tree/alrt} --bnni --seed {seed}`; `runtime=1440`, `mem_mb=16000`.
- **Test:** `validation_tree_rule_real_tools` (real MAFFT/ClipKIT/IQ-TREE with `LG+G4` for speed; supports parsed; GA → Prometheus, GB → Paralog).

### DD-16 [minor] Reference labels did not survive IQ-TREE tip renaming
- **Location:** `workflow/rules/prometheus.smk:172-173` (labels written from raw FASTA headers).
- **Issue:** IQ-TREE rewrites `: = , ; ( ) '` in names to `_` (it keeps `|` and `/`). A reference whose id contains one of these no longer matched its label and was treated as a candidate.
- **Evidence:** IQ-TREE 2.3.6, run here: tip names containing `: = , ; ( ) '` are rewritten to `_`, while `|` and `/` are kept.
- **Fix:** `combine_proteomes.py` reference mode (`--candidates`, `--ref GROUP=FASTA`, `--out-labels`, `--out-map`) writes `utils.sanitize`d ids to both `all.faa` and `reference_labels.tsv`. Sanitising collisions, duplicates, non-sanitised candidate ids, and an empty reference set with candidates present are errors.
- **Test:** `combine_proteomes_reference_mode`, `validation_tree_rule_real_tools` (labels ⊆ tips, e.g. `tr_A0A0A7PKR5_X_1`), `placement_outgroup_rooting_and_support` (an IQ-TREE-renamed tip still matches).

## E. Integration and tiers (`validate_prometheus.py`)

### DD-17 [critical] Tier A silently required the optional DUF3320
- **Location:** `workflow/scripts/py/validate_prometheus.py:315` (`detect_df and n_domains >= 4 and ...`).
- **Issue:** with the Contract-v2 four-category spec (DUF4011, DNA2, REase, DUF3320 with `required: false`), `n_domains >= 4` makes DUF3320 mandatory for Tier A. The same was already true under the five-category spec at 56f5437: there, the `nuclease` category matched no ProA domain, so the maximum was 4. That contradicts the spec and the biology. Roughly half of complete, genuine Prometheus loci would be capped at Tier B.
- **Evidence:** van den Berg et al. 2024 (*Cell Host Microbe* 32:1427): Prometheus is defined by DUF4011 (PF13195), a DNA2-family SF1 helicase (COG1112) and REase_MTES_1575 (PF18741), with DUF3320 (PF11784) absent from 54 % of systems. `config/config.yaml` marks DUF3320 `required: false`. Running the pristine `architecture_for` on the 11 fixture proteins gives n = 4 for the 7 that carry DUF3320 and n = 3 for the 4 that lack it, so the latter could never reach Tier A.
- **Fix:** Tier A = DefenseFinder call + every **required** category present (`utils.required_categories`) in canonical N→C order + supported Prometheus placement (DD-20) + catalytic `intact` (DD-19) + not fragmented (DD-23) + `min_len ≤ length ≤ max_len`. DUF3320 is recorded (`has_DUF3320`) but optional. Tier B = supported Prometheus placement, or required core complete + defense context; the reason lists the basis and each shortfall. Tier C otherwise. Rejected = supported paralog placement, or a complete (non-fragmented) ORF < 0.5 × min_len or > 2 × max_len.
- **Test:** `tier_rules` (canonical core without DUF3320 → A; every single shortfall → B; paralog → Rejected; fragment never rejected by length), `validate_end_to_end_known_answers`.

### DD-18 [major] Substring domain matching, an obsolete category, and no order check
- **Location:** `workflow/scripts/py/validate_prometheus.py:48-126` (`domain_matches` / `best_domain_for_hit` by substring of Pfam names and accessions; hard-coded `nuclease` category no longer in the spec; architecture = categories ordered by first start, so repeats and rearrangements were invisible).
- **Issue:** Pfam hits were assigned to categories by permissive substring rules, including a category (`nuclease`) that the spec no longer has. The architecture listed each category once, at its first start, so repeats and rearrangements were invisible, and the canonical N→C order was never checked.
- **Evidence:** at 56f5437, the spec was matched by case-insensitive substrings (`DNA2: names [DNA2, AAA, PIF1, helicase]`; `REase: [restriction, REase, MTES, endonuclease]`; a separate `nuclease` category). So any `AAA*` family, Helicase_C (PF00271) or Dna2 (PF08696) counted as the DNA2 helicase module. Running the pristine `architecture_for` on the fixture domtbl assigns the N-terminal Helicase_C to DNA2. Both haloarchaeal ProA homologues (A0A1H8VVA9, A0A2A2F1Z1) therefore come out as the spurious rearrangement `DNA2-DUF4011-REase`. With the exact spec they are the canonical core `DUF4011-DNA2-REase`, with `other_domains` ResIII and Helicase_C.
- **Fix:** the shared spec helpers: `utils.hit_category` (exact names and accessions), `utils.architecture_tokens` (N→C, repeats kept), `utils.is_canonical_order`, `utils.other_domains`. Additive catalog columns: `missing_required`, `canonical_order`, `other_domains`, `has_<category>`.
- **Test:** `real_hmmscan_matches_fixture_and_motifs` (A0A060RCK4 → `DUF4011-DNA2-REase-DUF3320`; A0AAV5N545 → `other_domains = Pkinase`), `tier_rules`.

### DD-19 [major] "Catalytically intact" was inferred from domain presence alone
- **Location:** `workflow/scripts/py/validate_prometheus.py:129-150` (`catalytic_status`: "intact" = DNA2 + nuclease/REase present + ≥ 3 domains).
- **Issue:** no residue was examined, so active-site mutants and pseudogenised copies were reported intact.
- **Evidence:**
  - SF1 helicase Walker A `[AG]x4GK[ST]` and motif II (DE) as in the DNA2/AAA_11/AAA_12/AAA_19/AAA_30 families.
  - REase_MTES_1575 (PF18741) is a "Vsr REase Fold" (InterPro; PMID 23768067), not a canonical PD-(D/E)xK motif. Its HMM consensus has no `(D/E)xK`, so a PD-(D/E)xK regex would be wrong.
  - The conserved acidic match states D28 and D45 are uppercase in the consensus. HMMER `p7_hmm.c` upper-cases a consensus residue iff its match emission is ≥ 0.5.
- **Fix:** residue-level assessment from the hmmscan alignments (`domains.hmmscan.txt`, `--hmmscan-out`):
  - Helicase module: Walker A and motif II at HMM-anchored positions, with a window fallback (`G..G.GK[ST]`; `[ACFGILMVWY]{3}DE.{18,32}GD..Q`) that can only confirm presence.
  - Nuclease module: residues aligned to the PF18741 acidic states.
  - Overall status: `intact` (both modules intact), `partial` (one module absent or disrupted), `degenerate` (both), `unknown` (no hits, or a module unassessable). `catalytic_evidence` records the residue and position.
  - Without alignments, a labelled `domain_proxy` fallback is used.
- **Result on the 11 fixture proteins:** 9 intact; A0A011V4U2 partial (`D45:C@1642`); A0A0A6X990 partial (`D28:G@1461`); the ProA consensus intact via the window fallback (`MotifII:DE@993(window)`).
- **Test:** `real_hmmscan_matches_fixture_and_motifs`, `catalytic_status_edge_cases` (K→A Walker A → disrupted; both modules lost → degenerate).

### DD-20 [major] Placement by nearest reference on an unrooted, unsupported tree
- **Location:** `workflow/scripts/py/validate_prometheus.py:156-228` (`assign_by_tree`: Prometheus if the patristically nearest reference tip is a Prometheus reference).
- **Issue:** uneven reference sampling and long branches decided the call, and no support was required for either "Prometheus" or "Paralog", including the Rejected tier.
- **Evidence:** the IQ-TREE FAQ (iqtree.github.io/doc/Frequently-Asked-Questions, last updated June 2024) says one can be more confident in a clade with SH-aLRT ≥ 80 % and UFBoot ≥ 95 % (UFBoot: Minh et al. 2013, Hoang et al. 2018; SH-aLRT: Guindon et al. 2010). IQ-TREE developers also stress that these are conventions, since bootstrap-type supports measure sampling variance. The thresholds are therefore configurable, and every locus records its supports (`placement_support`).
- **Fix:**
  - Root on the paralog outgroup, taking the minimal paralog side of the bipartition.
  - Assign Prometheus or Paralog only when the candidate lies in a clade with SH-aLRT ≥ `prometheus/placement_min_alrt` (80) and UFBoot ≥ `placement_min_ufboot` (95) that contains at least one reference of that group and none of the other. Crown or stem position is recorded.
  - Otherwise `Unplaced`, with `placement_note`: `unsupported_ingroup`, `unsupported_outgroup`, `references_not_separable`, `not_in_tree` or `no_tree`.
  - New columns: `placement_support`, `placement_clade`, `placement_note`.
- **Test:** `placement_outgroup_rooting_and_support`, `validation_tree_rule_real_tools`. On a real tree, labelling the Pkinase-fusion ProA homologues as paralogs gives `references_not_separable`, which is correctly detected.

### DD-21 [minor] Missing dendropy or an unparsable tree silently made everything Unplaced
- **Location:** `workflow/scripts/py/validate_prometheus.py:172-184`.
- **Issue:** an environment without dendropy, or a corrupt tree file, turned every candidate into Unplaced (capping every tier) with only a log warning.
- **Evidence:** the pristine code: `except Exception: ... return assignments` (empty), and downstream `tree_assign.get(pid, "Unplaced")`.
- **Fix:** dendropy is a hard import; it is pinned in `envs/annotation.yaml` (dendropy 4.6.1). A malformed tree is an error. The `();` placeholder is recognised explicitly (`no_tree`).
- **Test:** `validate_empty_and_failures` (bad Newick fails), `placement_outgroup_rooting_and_support` (`();` → `{}`).

### DD-22 [major] Defense context was genome-level, and the window parameter was ignored
- **Location:** `workflow/scripts/py/validate_prometheus.py:256-280` (`defense_context_for_genomes`: ≥ 2 systems of any type, or ≥ 1 non-Prometheus system anywhere in the genome), `:356` (`--context-window` parsed but unused).
- **Issue:** two Prometheus copies counted as "defense context", and so did any system anywhere in a multi-megabase genome. Design §1.5 asks for a ±N-gene defense-island neighbourhood.
- **Evidence:** the pristine code, and design §1.5 (defense-island neighbourhood of ±`prometheus/context_window` = 10 ORFs).
- **Fix:** a ±N-gene neighbourhood on the same contig, taken from `genes.tsv` (`contig`, `gene_index`) with the gene id as fallback. The locus's own system and every system of the focal type are excluded; contigs are treated as linear. New columns: `context_method`, `context_n_systems`, `context_systems`.
- **Test:** `validate_end_to_end_known_answers` (RM 11 genes away → False; Gabija on another contig → False; CBASS 3 genes away → True; the other Prometheus copy is not context), `real_defensefinder_preflight_run_parse_context` (real Septu system 2 genes from a ProA homologue → True; the lone ProA on contig 2 → False).

### DD-23 [major] "Fragmented" meant "short", and contig-edge fragments were Rejected
- **Location:** `workflow/scripts/py/validate_prometheus.py:455` (`fragmented = length_aa < min_len`), `:306-310` (length < 0.5 × min_len → Rejected for every locus).
- **Issue:** a truncated 2,000-aa ORF at a contig edge was not "fragmented", while a complete short homologue was. Short contig-edge fragments were rejected, although SCHEMA says truncated loci are retained and labelled.
- **Evidence:** Prodigal reports truncation explicitly (`partial=10/01/11`; `00` = complete), and this flag is carried as `partial` in the Contract-v2 `genes.tsv`. SCHEMA.md defines `fragmented` as truncated, retained and labelled.
- **Fix:** `fragmented` comes from Prodigal's `partial` flag in `genes.tsv` (`fragment_evidence`), with a labelled length proxy only when the locus is missing from `genes.tsv`. Fragments are never rejected by length. New columns: `partial`, `fragment_evidence`, `length_in_bounds`.
- **Test:** `tier_rules` (300-aa fragment → B; complete 300-aa ORF → Rejected), `validate_end_to_end_known_answers` (G4: 350-aa ORF running off a contig end is fragmented and kept).

### DD-24 [minor] Unused `--reps`, stop codons counted, invalid spec accepted
- **Location:** `workflow/scripts/py/validate_prometheus.py:354` (`--reps` required but unused), `:441` (length includes a terminal `*`), `:368-372` (invalid `--domains` JSON → empty spec).
- **Issue:** `detection_sensitivity.tsv` covered non-representative genomes; ORFs with a terminal stop were 1 aa too long; a malformed domain spec ran as an empty spec.
- **Evidence:** the pristine code: the `--reps` value is never read; `len(seq)` is taken on sequences that may end in `*`; a `json.loads` failure falls back to an empty spec.
- **Fix:** the genome universe of `detection_sensitivity.tsv` is metadata ∩ representatives. When `--reps` is the metadata table itself, as in `crossdomain.smk bact_validate`, this is a no-op. A terminal `*` is stripped. Invalid JSON and candidate ids without a sequence are errors.
- **Test:** `validate_end_to_end_known_answers` (G6, not a representative, is excluded), `validate_empty_and_failures`.

## F. Stage 1 figures

### DD-25 [minor] Missing tiers were plotted as Tier C
- **Location:** `workflow/scripts/py/plot_prometheus_validation.py:54`, `:96` (`fillna("C")`), `:30-35` (read errors swallowed).
- **Issue:** a catalog row with a missing or malformed tier was drawn as Tier C in the detection summary.
- **Evidence:** SCHEMA.md defines `confidence_tier` ∈ {A, B, C, Rejected}; `fillna("C")` maps a missing value to C.
- **Fix:** an unknown or missing tier is an error, and read errors propagate.
- **Test:** `diversity_extra_and_prometheus_plots` ("unexpected confidence_tier").

### DD-26 [minor] Validation-tree figure always exited 0
- **Location:** `workflow/scripts/R/plot_validation_tree.R:197-204` (top-level `tryCatch` → placeholder + `quit(status = 0)`), `:30-32`, `:74-76` (missing packages silently downgraded).
- **Issue:** any error (missing package, malformed tree, bad labels) produced a placeholder figure and exit status 0, so the figure rule always succeeded. The tree was also drawn unrooted.
- **Evidence:** the pristine code: the error handler draws a placeholder, and the script ends with `quit(status = 0)` on every path.
- **Fix:** errors exit 1. A placeholder is drawn only for a trivial tree (`();` or < 3 tips). The display is rooted on the paralog MRCA clade (`edgelabel = TRUE`), with support dots from the IQ-TREE `SH-aLRT/UFBoot` labels. The CLI is unchanged (rule in `synthesis.smk`).
- **Test:** `plot_validation_tree_r` (real ape/ggtree: rooted figure, placeholder for `();`, non-zero exit on a malformed Newick).

## G. Stage 2 diversity (`diversity.smk`; `cluster_diversity.py`, `architecture_typology.py`, `plot_diversity.py`, `plot_diversity_extra.py`)

### DD-27 [major] Subfamily column names broke the rarefaction figure
- **Location:** `workflow/rules/diversity.smk:12` (`str(x)` turns 0.30 into "0.3"), `workflow/scripts/py/cluster_diversity.py:40` (`cluster_{i}`), `workflow/scripts/py/plot_diversity.py:39` (looks for `cluster_0.30`).
- **Issue:** `subfamilies.tsv` had `cluster_0.3 … cluster_0.9` instead of SCHEMA's `cluster_0.30 …`. `plot_diversity.py` never found its column, so the rarefaction figure was always a placeholder.
- **Evidence:** Python `str(0.30) == '0.3'`; SCHEMA.md names the columns `cluster_0.30 ...`.
- **Fix:** `%.2f` formatting in the rule and the script. Identities that collide after formatting are rejected.
- **Test:** `mmseqs_cluster_real_and_rule` (columns `cluster_0.30 … cluster_0.90`), `mmseqs_failure_single_and_empty` (`0.3,0.30` rejected).

### DD-28 [major] MMseqs2 failures produced blank or header-only subfamily tables
- **Location:** `workflow/scripts/py/cluster_diversity.py:59-72` (a failed run returns None → blank columns), `:113-116` (mmseqs not on PATH → header-only), `:109-112` (a single sequence → header-only, so the protein vanished).
- **Issue:** a clustering failure yielded empty subfamily assignments, which downstream rarefaction and typology read as data. A single validated protein disappeared from `subfamilies.tsv`.
- **Evidence:** the pristine code paths: `return None` on a non-zero exit or a missing output, and `shutil.which("mmseqs") is None` leads to a header-only table.
- **Fix:** a failed run or a missing binary is an error. Every input protein must be assigned to exactly one cluster per identity. A single sequence is its own cluster. Empty input gives a header-only table. The temporary directory is removed. New optional `--summary` (`diversity/subfamily_summary.tsv`).
- **Test:** `mmseqs_failure_single_and_empty`.

### DD-29 [major] Clustering relied on implicit MMseqs2 defaults
- **Location:** `workflow/scripts/py/cluster_diversity.py:56`.
- **Issue:** the default sensitivity `-s 4` is low for 30 %-identity subfamilies. Without `--cluster-reassign`, cascaded clustering can leave members below the identity threshold to their final representative. Coverage and identity modes were implicit.
- **Evidence:** MMseqs2 15-6f452 (the version pinned in `annotation.yaml`) source and help: `-s` default 4; `--cluster-reassign` default 0; `--cov-mode` 0; with `-c > 0` and `--min-seq-id > 0` the automatic alignment mode computes coverage and identity.
- **Fix:** explicit `--cov-mode 0` (bidirectional; `diversity/mmseqs_cov_mode`), `--alignment-mode 3`, `--seq-id-mode 0`, `--cluster-mode 0`, `--cluster-reassign 1`, `-s 7.5` (`diversity/mmseqs_sensitivity`).
- **Test:** `mmseqs_cluster_real_and_rule` (11 fixture ProA homologues: the 2.8–2.9-kaa haloarchaeal pair co-cluster at 0.50, the Pkinase fusions at 0.30, 11 singletons at 0.90; cluster counts are non-decreasing with identity).

### DD-30 [major] Typology used a hard-coded, outdated domain list
- **Location:** `workflow/scripts/py/architecture_typology.py:27` (`CANONICAL` including the removed `nuclease` category), `:37-41` (re-derived through the old substring matcher), `:56-59` (invalid JSON → `{}`), `:61-64` (unreadable catalog → header-only).
- **Issue:** the typology counted a category the spec no longer has, re-derived architectures with the substring matcher (DD-18), and turned unreadable inputs into empty outputs. Canonical, optional-loss, core-loss and fusion variants were not distinguished.
- **Evidence:** the pristine code. The Contract-v2 spec has no `nuclease` category, so `has_nuclease` was always False, and the matcher it imported was the substring matcher of DD-18.
- **Fix:** the spec helpers. The SCHEMA columns are kept; `has_nuclease` is retained for compatibility and is always False with the current spec. Additive columns: `has_<other categories>`, `core_complete`, `canonical_order`, `rearranged`, `variant_class` (`canonical_full`, `canonical_core` = an optional domain such as DUF3320 absent, `rearranged`, `core_loss`, `none`), `missing_required`, `n_with_fusion`, `fusion_domains`, `n_tier_A/B/C`. A catalog/domtbl id mismatch is an error.
- **Test:** `architecture_typology_known_answers`.

### DD-31 [major] Rarefaction over carriers only, and subfamily ids labelled as architectures
- **Location:** `workflow/scripts/py/plot_diversity.py:49-68` (curve over Prometheus-carrying genomes only, Monte-Carlo only), `:94-109` (architecture_by_taxon filled "architecture" with `cluster_0.30` subfamily ids).
- **Issue:** genomes without Prometheus were not sampling units, so the curve was a carrier accumulation curve and overstated saturation per genome sampled. The architecture-by-taxon figure plotted the wrong variable.
- **Evidence:** incidence-based rarefaction: E[S_n] = Σ_k [1 − C(N − f_k, n) / C(N, n)] (Gotelli & Colwell 2001, *Ecol. Lett.* 4:379–391, doi:10.1046/j.1461-0248.2001.00230.x; Colwell et al. 2012, *J. Plant Ecol.* 5:3–21, doi:10.1093/jpe/rtr044). Chao2 (Chao 1987, *Biometrics* 43:783, doi:10.2307/2531532; references checked in Crossref): S_obs + ((m − 1)/m) · Q1² / (2 Q2), with the bias-corrected form S_obs + ((m − 1)/m) · Q1(Q1 − 1) / 2 when Q2 = 0 (m genomes; Q1 and Q2 = subfamilies found in exactly 1 and 2 genomes).
- **Fix:** exact rarefaction over all master-table genomes, with a permutation band (`--permutations`, default 200; `--seed`) and Chao2 annotated as a lower bound. Per-locus architectures come from `--catalog`. Curve written to `diversity/subfamily_rarefaction.tsv` (`--curve`); `--level` selects the taxonomic rank.
- **Test:** `rarefaction_exact_values` (4 genomes, f = [2, 1]: E[S] = 0, 0.75, 4/3, 1.75, 2.0; Chao2 = 2.375).

### DD-32 [minor] Extra diversity figure had a hard-coded domain list
- **Location:** `workflow/scripts/py/plot_diversity_extra.py:26` (`CANON` with `nuclease`, an always-empty column, and DUF3320 before DNA2), `:29-35` (read errors swallowed).
- **Issue:** the domain-presence panel had an always-empty `nuclease` column and the wrong N→C order, and unreadable inputs became placeholders.
- **Evidence:** the pristine code. The spec's N→C order is DUF4011, DNA2, REase, DUF3320 (`config/config.yaml` `prometheus.domains`).
- **Fix:** `domain_categories()` derives the categories from the typology `has_*` columns and architectures in N→C order. Read errors propagate. The CLI is unchanged (rule in `plots_extra.smk`).
- **Test:** `diversity_extra_and_prometheus_plots` (`[DUF4011, DNA2, REase, DUF3320]`).

## H. Rule hygiene

### DD-33 [minor] Rules without `log:`
- **Location:** `prometheus.smk`: `combined_proteome`, `assemble_candidates`, `validate_prometheus`, `plot_prometheus`. `diversity.smk`: `architecture_typology`, `plot_diversity`.
- **Issue:** six rules wrote no log, so their failures left no record.
- **Evidence:** `tests/check_rules.py` (the project's rule-hygiene contract) requires `log:` on every rule.
- **Fix:** every rule has `log:` and redirects stdout/stderr to it. No `|| true` remains in the owned rules.
- **Test:** `rules_hygiene_and_no_hardcoded_domain_lists` (`tests/check_rules.py` → RULES HYGIENE OK; no hard-coded `CANON` lists or `"nuclease" in present` logic in owned scripts).

---

## Verified correct, no change needed
- `utils.parse_domtbl` column indices (i-Evalue, score, ali from/to). A real `hmmscan --cut_ga --domtblout` run reproduces `tests/fixtures/proA_homologs.pfam_subset.domtbl` hit for hit (`real_hmmscan_matches_fixture_and_motifs`).
- DefenseFinder output names `{input basename}_defense_finder_{systems,genes}.tsv` match the rule outputs (DefenseFinder posttreat; real run).
- DefenseFinder `type`/`subtype` = `model_fqn.split('/')[-2]` / `[-1]`; `hit_pos` is a global rank in the input file (macsypy). The context code uses `genes.tsv`, or gene ids, rather than `hit_pos` across contigs.
- MAFFT 7.526 silently strips `*`, so terminal stops need no pre-processing for the alignment.
- IQ-TREE support labels are `SH-aLRT/UFBoot`. The output tree is unrooted with a trifurcating root, which is why placement re-roots on the paralog outgroup. Identical sequences are kept.
- Haloarchaeal ProA homologues (2,796 and 2,893 aa, N-terminal ResIII/Helicase_C fusion) fall inside `max_len_aa: 3000`, close to the bound.
- The ProA `hmmemit -c` consensus (1668 aa) is called as Prometheus by DefenseFinder 3.0.0 with models 3.1.0. This is the positive control.

## Design decisions (scientific)
- **REase evidence.** The PF18741 acidic core (D28/D45 match states) is a model-derived hypothesis, not a mutagenesis-mapped active site; van den Berg et al. 2024 report no ProA nuclease mutagenesis. Loci with a substitution at these states are labelled `partial` (Tier B), not rejected.
- **Defense context is not required for Tier A.** It can raise a locus to Tier B, and it never rejects.
- **Stem placement counts for Tier A.** A supported stem position next to the Prometheus crown is accepted (crown/stem recorded), because archaeal Prometheus may be divergent from the mostly bacterial references.

## Open items (not fixed here)
- Design §2.4 (genomic-context diversity via clinker/cblaster) and the methanogen-only versus mixed subfamily composition (§2.1) are not implemented; they were not in the pipeline at `56f5437`.
- The real-DefenseFinder test needs network access to codeload.github.com to download the pinned models (≈ 67 MB). Offline it prints `skip`; set `PROM_DF_MODELS` to a pre-installed models root.
- Local tool versions differ from some env pins: MMseqs2 18-8cc5c (pin 15.6f452), ClipKIT 2.14 (pin 2.3.0), seqkit 2.14 (pin 2.8.2). The behaviours relied on (MMseqs2 parameters, ClipKIT exit 0 on errors) were checked in the pinned versions' source.
