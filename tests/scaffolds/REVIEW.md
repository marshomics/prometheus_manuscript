# Stage 3 (scaffolds) review — leaf-1.1.3

Scope: `workflow/rules/scaffolds.smk` (gtdbtk_align, species_tree, mcr_search,
methanogen_classify, build_master_table) and `workflow/scripts/py/{root_tree,mcr_search,
classify_methanogen,build_master_table}.py`. Baseline: commit 56f5437 (pristine copy);
line numbers below refer to that baseline. Suite: `bash tests/run_leaf.sh scaffolds`
(28 tests, ~45 s, real IQ-TREE 2.3.6 / HMMER 3.4 / MAFFT 7.526 / R ape / Snakemake 8.30;
GTDB-Tk replaced by a stub that writes GTDB-Tk 2.4.0's own file names).

## Primary sources used

| Tool / format | Version | Evidence |
|---|---|---|
| GTDB-Tk | 2.4.0 (env pin) | source at tag `2.4.0`: `gtdbtk/config/output.py` (paths), `markers.py` (`align`, `_apply_mask`, `genome_domain`, `_write_msa`), `main.py` (`_genomes_to_process`, `_verify_genome_id`), `cli.py` (flags), `files/batchfile.py` |
| IQ-TREE 2 | 2.3.6 (env pin; installed) | real runs: PMSF two-step, exit codes, name rewriting, seeds; IQ-TREE docs (Complex Models: PMSF guide tree "under the simpler LG+F+G model", `-ft`) |
| HMMER | 3.4 (env pin; installed) | real `hmmsearch` runs incl. `--cut_ga`; tblout cross-checked with Biopython 1.83 `SearchIO` (`hmmer3-tab`) |
| KofamScan | takaram/kofam_scan master | `lib/kofam_scan/result/parser.rb` (score columns), `result/hit.rb` (`score >= threshold`), `result/with_threshold_scale/hit.rb`, `executor.rb` (`hmmsearch -T 0`), `ko.rb` (ko_list format) |
| KOfam ko_list | genome.jp, downloaded 2026-10-30 | K00399 741.47 full (F = 1.000), K00401 506.80 full, K00402 359.37 full; 8,684 of 28,504 KOs use `score_type = domain` (fixture: `fixtures/kofam_ko_list_mcr_subset.tsv`) |
| TIGRFAM/NCBIfam | InterPro API | TIGR03256 McrA (met_CoM_red_alp), TIGR03257 McrB, TIGR03259 McrG, TIGR03260 McrD |
| DefenseFinder | v2.0.0 | `defense_finder_posttreat/df_systems.py` (systems columns incl. `activity`), `best_solution.py` (`type = model_fqn.split('/')[-2]`; `activity = 'Antidefense'` for ADF models) |
| defense-finder-models | 2.0.0 tag | `definitions/DefenseFinder/Prometheus/Prometheus.xml` -> type `Prometheus` |
| MAFFT | 7.526 (installed) | real `--add --keeplength` runs |
| DendroPy | 4.6.1 (env pin; installed) | probes of `reroot_at_edge`, `reroot_at_midpoint`, `collapse_basal_bifurcation`, `prune_taxa_with_labels`, Newick writer |
| R ape | prom-r (installed) | `read.tree` on quoted / unquoted labels, `is.rooted` on bifurcating vs trifurcating roots |
| Snakemake | 8.30 (installed) | dry runs of all branches; real execution of gtdbtk_align + species_tree |

## Findings

Severity: critical = silently wrong or empty results downstream; major = wrong biology,
statistics or data loss in realistic inputs; minor = robustness / reproducibility /
documentation.

### scaffolds.smk

**SC-01 — critical — scaffolds.smk:35,37,40**
- Issue: `gtdbtk identify ... || true` and `gtdbtk align ... || true`, then an empty gzip
  placeholder when no MSA was found. A crashed GTDB-Tk or a missing database gave an
  empty MSA, a `();` species tree and silently empty phylogenetic analyses.
- Evidence: GTDB-Tk raises `GTDBTkExit` (non-zero exit) on bad input and failures
  (`main.py` `_genomes_to_process`; `markers.py` "Unable to assign a domain to any genomes").
- Fix: no `|| true`; `db.gtdbtk_data` must be an existing directory; a missing MSA after a
  successful GTDB-Tk run is an error pointing at `gtdbtk_dropped.tsv`.
- Test: `test_rule_gtdbtk_align_shell_fails_loudly` (identify crash, no MSA, missing DB;
  the rendered rule shell run under `set -euo pipefail`).

**SC-02 — major — scaffolds.smk:39**
- Issue: `find -name '*user_msa.fasta.gz' | head -n1` can pick
  `gtdbtk.bac120.user_msa.fasta.gz`.
- Evidence: `markers.py` writes one `{prefix}.{bac120|ar53}.user_msa.fasta.gz` per domain
  that has genomes (`output.py` `PATH_*_USER_MSA` under `align/`; `_write_msa(zip_output=True)`
  appends `.gz`); `genome_domain` puts a genome in bac120 when its bac120 marker % >= its ar53 %.
- Fix: explicit `align/gtdbtk.{tree.marker_set}.user_msa.fasta.gz` (marker_set validated:
  ar53 | bac120).
- Test: `test_rules_gtdbtk_align_and_species_tree_end_to_end` (the stub writes both MSAs;
  the tree holds only the ar53 genomes; the bac120 genome is reported).

**SC-03 — major — scaffolds.smk:19,33**
- Issue: `tree.scope: domain` was a silent no-op: dropping `--skip_gtdb_refs` does not put
  references in the user MSA.
- Evidence: `markers.py`: `trimmed_user_msa = {k: v for k, v in trimmed_seqs.items() if k in
  user_msa}`; references only go to `{prefix}.{ms}.msa.fasta.gz`, headers `>{id} {taxonomy}`.
- Fix: scope=domain aligns with references and uses `gtdbtk.{ms}.msa.fasta.gz`. The rule
  strips header descriptions, extracts the reference taxonomy from the headers so an
  outgroup can be a GTDB reference taxon, roots, then prunes the references, so
  `species_tree.treefile` tips are always analysed genomes. The unpruned rooted tree is
  kept as `species_tree/species_tree.full.treefile`. An invalid scope fails at DAG build.
  Memory and runtime are scaled for domain scope (IQ-TREE docs: the first PMSF phase
  needs mixture-model RAM).
- Test: `test_domain_scope_reference_outgroup_then_prune`; Snakemake dry run of the domain
  branch (no `--skip_gtdb_refs`, `gtdbtk.ar53.msa.fasta.gz`).

**SC-04 — major — scaffolds.smk:34-40**
- Issue: genomes that GTDB-Tk leaves out of the MSA disappeared from every phylogenetic
  analysis without a trace.
- Evidence (markers.py@2.4.0):
  - `_apply_mask` prunes user genomes below `--min_perc_aa` and writes
    `align/{prefix}.ar53.filtered.tsv` ("Insufficient number of amino acids in MSA (x%)").
  - Genomes assigned to bac120 go to the bac120 files.
  - Genomes without markers go to `align/intermediate_results/{prefix}.align.failed.tsv`.
  - identify failures go to `identify/{prefix}.failed_genomes.tsv`.
- Fix: `gtdbtk_dropped.tsv` lists each dropped genome with GTDB-Tk's own reason;
  `species_tree_membership.tsv` gives `in_msa`, `in_tree` and `status` for every analysed
  genome; the master table gains a `tree_status` column. `--min_perc_aa` is explicit
  (`tree.min_perc_aa`, GTDB-Tk default 10).
- Tests: end-to-end rule test (all four reasons); `test_trivial_tree_and_membership_reasons`;
  master test (`tree_status`).

**SC-05 — major — scaffolds.smk:34**
- Issue: `--genome_dir` processes every `*.fna` in the staging directory, including stale
  genomes from an earlier Phase-1 run. A re-used identify directory can also feed earlier
  genomes into standalone `align`, which skips the consistency check
  (`genomes_to_process` is None there).
- Evidence: `main.py` `_genomes_to_process` (`os.listdir`); `markers.py` `align` consistency
  check runs only when `genomes_to_process` is set.
- Fix: `--batchfile` (absolute path `<TAB>` genome_id, `files/batchfile.py` format) built
  from `representatives.txt`; identify/align directories are recreated on every run.
- Test: end-to-end rule test (batchfile content and stub call line asserted).

**SC-06 — critical — scaffolds.smk:80,83,86,88**
- Issue: `|| true` on both IQ-TREE runs, plus `cp species.treefile ... || echo "();"`. An
  IQ-TREE crash (e.g. out of memory) silently produced `();`, which root_tree then passed
  through.
- Evidence: IQ-TREE 2.3.6 exits 2 on errors (unknown model, missing `-ft` file, `-B` with
  < 4 sequences), verified locally.
- Fix: removed. The only nonfatal branch is an MSA with fewer than 4 genomes (`# nonfatal:`),
  and the membership and rooting files record it.
- Test: `test_rule_species_tree_shell_iqtree_crash_and_small_msa`.

**SC-07 — major — scaffolds.smk:75**
- Issue: `n=$(grep -c '>' $msa || echo 0)` yields `0\n0` for an empty MSA (grep -c prints 0
  and exits 1), which breaks `[ "$n" -ge 4 ]`; it also counted `>` anywhere on a line.
- Fix: `awk` count of header lines.
- Test: small-MSA branch of `test_rule_species_tree_shell_iqtree_crash_and_small_msa`; G2
  hygiene check.

**SC-08 — minor — scaffolds.smk:79**
- Issue: the PMSF guide tree used `-m LG+G -fast`; the IQ-TREE docs obtain the guide tree
  under the simpler LG+F+G model.
- Fix: `tree.pmsf_guide_model` (default `LG+F+G`, full search).
- Verified with real IQ-TREE 2.3.6: `-m LG+C20+F+G -ft guide.treefile -B 1000 -alrt 1000`
  prints "COMPUTING SITE FREQUENCY MODEL" and writes `.sitefreq`.
- Test: `test_root_outgroup_on_real_iqtree_pmsf_tree`; end-to-end rule test (`.sitefreq`
  exists).

**SC-09 — minor — scaffolds.smk:49-50,98**
- Issue: `rooting_method.txt` (cited in the README) was written but not declared as an
  output.
- Fix: declared outputs `rooting`, `membership` and `full_tree`.
- Test: end-to-end rule test.

**SC-10 — minor — scaffolds.smk:161-192**
- Issue: `build_master_table` had no `log:`.
- Fix: every rule logs.
- Test: `test_rules_hygiene` (G2).

**SC-11 — minor — scaffolds.smk:110-112,135-141**
- Issue: `methanogen.mcr_kos` was unused, and `mcr_required` was passed but ignored by the
  script. With `anme_via_mcrA_tree: false` the reference files were still mandatory inputs.
- Fix: `--kos` (each configured subunit must have a profile); `--required` implemented;
  reference inputs are conditional.
- Tests: `test_classify_without_placement_and_required_subunits`; dry run of the ANME-off
  branch.

### root_tree.py

**RT-01 — critical — root_tree.py:25-27,104-110,163-166**
- Issue: every failure (dendropy missing, parse error, rooting exception) wrote the
  UNROOTED input as the canonical rooted tree. The R consumers then midpoint-rooted it
  silently, because they root only when `!ape::is.rooted`.
- Fix: these are errors (exit 1); only a genuinely trivial tree (`();`, fewer than 3 tips)
  is passed through, and that is recorded.
- Test: `test_unknown_tip_and_missing_dendropy_are_errors` (`python -S` hides dendropy).

**RT-02 — major — root_tree.py:141,147,151,154**
- Issue: SH-aLRT/UFBoot supports ended up on the wrong bipartitions after rerooting.
- Evidence: DendroPy 4.6.1 keeps internal labels on nodes when rerooting. Probe on an
  IQ-TREE 2.3.6 PMSF tree: 3 of its 5 labelled bipartitions carried the wrong label after
  rerooting; for example, the {G1,G2} split lost its 100/100. This is the general problem described by Czech, Huerta-Cepas & Stamatakis
  2017 (MBE 34:1535, doi:10.1093/molbev/msx055).
- Fix: labels are keyed by bipartition before rooting and re-assigned afterwards; both
  root children carry the root split's support.
- Tests: `test_support_labels_follow_bipartitions_after_rerooting`;
  `test_root_outgroup_on_real_iqtree_pmsf_tree` (every IQ-TREE split label preserved).

**RT-03 — major — root_tree.py:124-145**
- Issue: the outgroup MRCA was computed on IQ-TREE's arbitrary root. A monophyletic
  outgroup straddling that root was declared unrootable (midpoint fallback). A
  non-monophyletic outgroup was rooted at its MRCA, with ingroup tips inside the "outgroup",
  noted only as `non-monophyletic@MRCA` in a log.
- Fix: monophyly is tested as an unrooted bipartition. A configured outgroup that is absent,
  equal to all tips, or not monophyletic FAILS, with the intruding tips listed. A labelled
  midpoint fallback happens only with `--on-outgroup-failure midpoint`
  (`tree.outgroup_on_failure`). No outgroup configured gives
  `midpoint(no_outgroup_configured)` plus a warning. An already-rooted input is unrooted first.
- Test: `test_outgroup_failures_and_explicit_fallback`.

**RT-04 — major — root_tree.py:72**
- Issue: `outgroup_taxon` matched by substring (`tx in classification.lower()`), so
  `p__Altiarchaeota` also matched `p__Altiarchaeota_A`, which GTDB treats as a different
  lineage.
- Fix: exact rank-token matching (case-insensitive); a bare name matches that name at any
  rank; comma lists are allowed; the reference taxonomy is used for scope=domain.
- Tests: `test_outgroup_taxon_is_exact_rank_token_match`; the PMSF test (the
  `p__Altiarchaeota_A` tip is excluded from the outgroup).

**RT-05 — major — root_tree.py:158**
- Issue: the dendropy writer emitted a `[&R] ` prefix and single-quoted labels containing
  `/` (`'100/100'`). ape keeps the quote characters in node labels (verified), so the
  supports were non-numeric in R.
- Fix: own iterative Newick writer: no `[&R]`, no quotes unless Newick requires them, no
  root branch length, 10-decimal lengths.
- Test: PMSF test (`ape`: ROOTED TRUE, BINARY TRUE, numeric supports).

**RT-06 — major — root_tree.py:133,147,151,154**
- Issue: `reroot_at_midpoint` mishandles a midpoint that falls exactly on a node:
  - it turns a tip into the root node: `(A:1,B:1,C:1)` became `(A:1.0,C:1.0)B;`, losing tip B;
  - or it leaves a trifurcating root: `((A,B),(C,D),E)` gave root degree 3, which ape
    reports as unrooted.
- Fix: own exact midpoint (diameter path; the root goes on the edge holding D/2, with a
  zero-length side when D/2 is a node); the root is always bifurcating.
- Test: `test_midpoint_rooting_is_exact_and_preserves_distances` (25 random trees plus both
  pathological ones; all patristic distances are unchanged and both sides reach exactly D/2).

**RT-07 — major — (not handled)**
- Issue: IQ-TREE rewrites characters outside `[A-Za-z0-9_-.|/]` to `_` (25/25 characters
  predicted on a real run). The tip then differs from the genome_id, `tree_tip` is blank,
  and the genome is silently lost from the phylogenetic analyses.
- Fix: `--keep-tips` restores the original id when it maps back uniquely. An id that would
  need Newick quoting (`'`, `:`) is refused, since GTDB-Tk already forbids `()[],;=` and
  spaces. Unknown tips are an error.
- Test: `test_iqtree_renamed_tips_are_restored_to_genome_ids`.

**RT-08 — minor — root_tree.py:155**
- Issue: an unknown `--method` was silently midpoint-rooted.
- Fix: error.
- Test: failures test.

**RT-09 — minor — root_tree.py:35-42,94**
- Issue: the rooting record was a single free-text line, and log write errors were
  swallowed.
- Fix: `rooting_method.txt` keeps the first-line `method=...; tips=N`, followed by key=value
  lines (fallback reason, outgroup tips and monophyly, taxon match counts, pruned
  reference count, restored labels, root clade sizes, membership counts).
- Tests: assertions on these keys in the rooting tests.

### mcr_search.py

**MS-01 — major — mcr_search.py:15-24,214; config `methanogen.kofam_score_fraction: 0.5`**
- Issue: hits were accepted at half the KOfam adaptive threshold.
- Evidence:
  - KofamScan's significance is `score >= threshold` (`result/hit.rb`; the threshold scale
    is unset by default, and `-T` scales it).
  - K00399 has threshold 741.47 with F-measure 1.000 (ko_list), so a 0.5 fraction accepts
    hits >= 370.7 bits.
  - A complete (partial `00`) 250-aa McrA piece scores 497.8 against a fixture-built McrA
    HMM: accepted at 0.5, rejected at 1.0.
  - Divergent McrA homologues (e.g. alkyl-CoM reductases of alkane oxidisers) sit in that
    window.
- Fix:
  - The script default and the rule's `cfg()` default are 1.0 (KofamScan's definition).
  - Truncation is handled separately: a hit below the threshold is a credible FRAGMENT
    only when its gene runs off a contig end (Prodigal partial != `00` from
    `genes.tsv`, Contract v2) and it scores >= `partial_score_fraction` × threshold
    (default 0.5).
  - The config value itself must change (request R1).
- Tests: `test_mcr_search_half_kofam_threshold_is_permissive`;
  `test_mcr_search_real_hmmsearch_thresholds_fragments_copies` (fragment with flag `10`
  accepted as `partial_rescue`; the same fragment with `00` rejected).

**MS-02 — major — mcr_search.py:81-107,214**
- Issue: the ko_list `score_type` was ignored (the full score was always used).
- Evidence: KofamScan `parser.rb` uses tblout column 9 (best-1-domain score) for
  `score_type == domain` and column 6 (full score) for `full`; 8,684 KOs are `domain`.
- Fix: implemented as in KofamScan, including the 5-decimal rounding of the scaled
  threshold.
- Test: `test_mcr_search_kofam_domain_score_type_uses_best_domain` (tandem protein: full
  988 bits passes under `full`, best domain 496 fails under `domain`).

**MS-03 — critical — mcr_search.py:136-138,179-191,216**
- Issue:
  - An hmmsearch non-zero exit was only a WARNING, so its hits were lost.
  - No HMMs, no hmmsearch, or a missing ko_list produced a header-only table (zero
    methanogens) or fell back to an E-value floor of 1e-10.
- Fix: all of these are errors; there is no E-value fallback (a model without a curated
  threshold is an error); only an empty proteome gives a header-only table.
- Test: `test_mcr_search_failures_are_loud_and_empty_input_is_header_only`.

**MS-04 — major — mcr_search.py:45-56**
- Issue: the subunit came from the file stem only. A TIGRFAM profile named `TIGR03256.hmm`
  became gene `TIGR03256`, which classify ignores, so no methanogens were called
  (resources/README offers TIGRFAM as an option).
- Fix: map via KO, TIGRFAM accession (InterPro) or name. It is an error if no McrA profile
  exists or if a configured KO has no profile.
- Tests: `test_mcr_search_ga_mode_matches_hmmer_cut_ga` (TIGR03256); failures test.

**MS-05 — minor — mcr_search.py:59-78**
- Issue: GA mode compared only GA1 with the full score. HMMER 3.4 `--cut_ga` lists in
  `--tblout` a sequence that passes GA1 but has no domain >= GA2 (N = 0, verified).
- Fix: full score >= GA1 and best domain >= GA2 (the Pfam/pfam_scan membership rule; both
  scaled).
- Test: the GA test asserts the accepted set equals `--cut_ga --domtblout` targets, and the
  tandem protein is in `--tblout` but not accepted.

**MS-06 — minor — mcr_search.py:196-222**
- Issue: one protein could count as two subunits, inflating `n_mcr`.
- Fix: one subunit per protein (highest score/threshold); the loser is recorded as
  `superseded_by_<gene>`.
- Test: `test_mcr_search_one_subunit_per_protein`.

**MS-07 — minor — mcr_search.py:218**
- Issue: `--pmap` was ignored and the genome was taken from `genome_of`.
- Fix: the protein map is authoritative; a target absent from it is an error.
- Test: failures test.

**MS-08 — minor — (missing)**
- Issue: copy number, fragment status and the thresholds used were not reported.
- Fix: additive `mcr_hits.tsv` columns plus the audit table `mcr_hits_all.tsv` (every hit
  with E <= 1e-3 and its decision).
- Test: fragments test (G5 `n_copies = 2`, representative = higher-scoring copy; scores
  equal the independent Biopython parse within 0.05 bits).

### classify_methanogen.py

**CM-01 — critical — classify_methanogen.py:24-27,110-173,295**
- Issue: every ANME-check failure (missing or empty references, mafft or IQ-TREE failure,
  dendropy missing, McrA missing from the proteome) became `anme_check_skipped` with
  is_anme False, so ANME were silently counted as methanogens.
- Fix: these are errors; zero candidates is the only skip, and it is logged.
- Test: `test_classify_failures_are_loud_and_empty_input_is_header_only`.

**CM-02 — major — classify_methanogen.py:220,284-292**
- Issue: only the exact label `ANME` counted. An `other` nearest reference (e.g. an
  alkyl-CoM reductase) and subgroup labels such as `ANME-1` / `ANME-2a` both led to
  methanogen = True.
- Fix: methanogen requires `mcr_group == methanogen`; strict label vocabulary
  (`methanogen`, `other`, `ANME*`), with anything else an error; `other`, `ambiguous` and
  `unplaced` genomes are not methanogens, are counted, and are flagged.
- Tests: `test_classify_placement_with_real_mafft_iqtree` (labels `ANME-1`; the `other`
  candidate is not a methanogen); unit tests.

**CM-03 — major — classify_methanogen.py:176-221**
- Issue: assignment was by the patristically nearest reference alone, which is sensitive to
  long branches.
- Fix: nearest-labelled-clade rule on the midpoint-rooted McrA tree: walk up to the first
  node with references; one label is assigned, mixed labels give `ambiguous`. Candidates
  with fewer than 50 residues in reference columns are `unplaced`. The patristic nearest
  reference is kept as secondary evidence.
- Test: `test_nearest_labelled_clade_rule_units` (long-branch case: the clade says ANME
  although the nearest reference is a methanogen; mixed sister clade gives `ambiguous`;
  nested candidates).

**CM-04 — minor — classify_methanogen.py:149**
- Issue: IQ-TREE ran without `--seed`.
- Fix: `--seed` (config `seed`).
- Test: placement test.

**CM-05 — minor — classify_methanogen.py:118**
- Issue: the whole combined proteome was loaded into a dict in an 8 GB rule.
- Fix: stream it and keep only the candidate ids.
- Test: placement test.

**CM-06 — minor — classify_methanogen.py:230-231**
- Issue: `--required` was parsed but McrA was hard-coded.
- Fix: KO-to-subunit mapping; unknown tokens are an error.
- Test: required-subunits test.

**CM-07 — minor — classify_methanogen.py:245-251**
- Issue: only genomes with mcr hits were written.
- Fix: one row per genome of the protein map, with explicit `no_mcr`.
- Test: placement test (GX).

### build_master_table.py

**BM-01 — critical — build_master_table.py:66-81,91-96,120-125,145-150,170-175,204-225,240-250,353-365**
- Issue: silent empty joins:
  - a missing or unreadable side table became `{}` (e.g. a broken catalog meant no
    Prometheus anywhere);
  - a missing DefenseFinder file counted as 0 systems, and an unreadable one fell back to
    a line count;
  - missing dendropy or a tree parse error left `tree_tip` blank for every genome;
  - an unreadable meta table produced header-only outputs.
- Fix: all are errors. Side tables naming genomes outside the universe, duplicated
  genome_ids, and universe genomes missing from the species or taxonomy table are errors.
- Test: `test_master_table_inconsistent_inputs_fail`.

**BM-02 — major — build_master_table.py:212-213**
- Issue: the focal system was excluded by the substring `prometheus` in `type`.
- Evidence: in defense-finder-models 2.0.0 the type is exactly `Prometheus`
  (`definitions/DefenseFinder/Prometheus/Prometheus.xml`;
  `type = model_fqn.split('/')[-2]`).
- Fix: exact (case-insensitive) equality with `--system` (= `prometheus.defensefinder_system`).
- Test: master test (a `PrometheusLike` type stays in the excl-Prometheus count).

**BM-03 — major — (not handled)**
- Issue: DefenseFinder 2.0 with `-a` reports AntiDefenseFinder systems in the same table
  (`activity = Antidefense`), and they were counted as defense systems.
- Fix: rows with activity `Antidefense` are excluded from both counts.
- Test: master test.

**BM-04 — major — build_master_table.py:178-192**
- Issue: catalog tiers were not validated. A blank or unknown tier counted as present with
  tier `''`, which tier-based consumers read as absent.
- Fix: tiers must be A/B/C/Rejected (anything else is an error).
- Test: inconsistent-inputs test (tier `D`).

**BM-05 — major — Contract v2 item 6 (rev 2)**
- Issue: tier consumers coded B/C-only genomes as confident absences in the tier-A
  analysis.
- Fix: `prometheus_tier` = best non-rejected tier, `prometheus_present` = any non-rejected
  locus, `prometheus_n_loci` = count of non-rejected loci (confirmed). Added
  `prometheus_status_primary` (via `--primary-tiers`, fed by
  `stats.confidence_tiers_primary`) and `prometheus_status_sensitivity` (`--sensitivity-tiers`,
  default A,B), each present | absent | excluded.
- Test: `test_master_table_values_and_presence_rule` (regression: a B-only genome is
  EXCLUDED under A and PRESENT under A+B; C-only is excluded under both; Rejected-only is
  absent; A + Rejected is present with `n_loci = 1`).

**BM-06 — major — build_master_table.py:256,305-321**
- Issue: blank or NA ecology cells, and genomes absent from the species table, were
  imputed as 0 (absence). This contradicts the script's own "never silently imputed".
- Fix: a blank cell stays NA (empty); a universe genome absent from the species table is an
  error; non-binary numbers are coerced to 1 and counted in the log; a missingness table
  (`ecology_missingness.tsv`: n_1, n_0, n_NA per column) is written.
- Test: master test.

**BM-07 — minor — build_master_table.py:19-20,376**
- Issue: `--envo-map` was silently ignored, and config `ecology.envo_map` / `predictors`
  refer to free-text columns that the species table does not have.
- Fix: the ignore is logged and documented (design §3.4: the curated binary group/niche
  columns replace ENVO mapping).
- Test: master test (stderr).

**BM-08 — minor — (missing)**
- Issue: the README lists contig count and N50 as quality covariates, but they were absent
  from the master table.
- Fix: additive `n_contigs`, `n50` and `gc` columns.
- Test: master test.

## Verified correct (no change needed)

- GTDB-Tk 2.4.0 CLI:
  - `identify --batchfile/--out_dir/--cpus` and `align --identify_dir/--out_dir/--cpus/--min_perc_aa/--skip_gtdb_refs` (cli.py);
  - batchfile = 2 or 3 tab-separated columns;
  - genome ids may not contain `()[],;=` or spaces, nor clash with reference ids.
- IQ-TREE 2.3.6:
  - The FASTA sequence name is the first whitespace-delimited token.
  - `-B` needs >= 4 sequences.
  - The same `--seed` with `-T 1` is bit-identical; with `-T 2` branch lengths differ at the
    1e-6 level only (topology and supports identical). Seeds are set everywhere;
    exact reproducibility therefore also needs the same thread count.
- hmmsearch `--tblout` columns: target = 1, query = 3, full E/score = 5/6, best-domain
  E/score = 8/9 (KofamScan positions; Biopython parse agrees).
- MAFFT `--add --keeplength` keeps the reference width; a trailing `*` in a query is
  tolerated (it is sanitised anyway).
- DendroPy `reroot_at_edge(length1, length2)`: length1 goes to the tail (old parent) side
  and length2 to the head side. Only halves are used.
- `methanogenesis_subtype` is `unresolved` for methanogens and `NA` otherwise; no subtype is
  asserted from mcr genes (SCHEMA enum allows both). Kept.
- `representative`: the universe is the species representatives (genome_metadata), so it is
  True for all; kept as an explicit join.

## Not verified here

- GTDB-Tk itself was not executed: it is not installed locally and its reference data are
  about 100 GB. Its interface was checked against the 2.4.0 source, and the rule shell runs
  against a stub that writes the same file names.
- GeneRax parsing of the new species-tree format could not be tried: GeneRax's MPI runtime
  fails to start in this sandbox ("PMIx server's listener thread failed").
- The KOfam K00399 profile was not available (the profile set is 1.4 GB), so the
  0.5-vs-1.0 evidence combines the real ko_list threshold with a fixture-built McrA HMM.
- scope=domain resources (128 GB base, 7 days) are an estimate.
