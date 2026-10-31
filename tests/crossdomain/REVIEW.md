# Crossdomain review — leaf-1.1.8 (stage `crossdomain`)

Scope: the opt-in bacterial / cross-domain branch (`bacteria.enabled: true`).
Owned files: `workflow/rules/crossdomain.smk`; `workflow/scripts/py/{extract_bacterial_prometheus,
build_bacterial_table, build_bacterial_niche_master, crossdomain_analysis, plot_crossdomain_tree,
plot_crossdomain_extra, plot_viral_crossdomain, prune_tree_to_taxa}.py`; `tests/crossdomain/**`.
No owned R files (the R scripts this branch calls belong to other leaves; they are exercised
end-to-end here, not edited). Baseline: commit 56f5437 (pristine copy); line numbers in
**Location** refer to the pristine files.

Method: four passes (correctness, robustness/no-silent-failure, scientific defensibility,
tests), with tool behaviour verified by running the installed tools or reading their source:
IQ-TREE 2.3.6, MAFFT 7.526, ClipKIT 2.14.0, HMMER 3.4, HyPhy (FEL run on the e2e output),
Snakemake 8.30, dendropy 4.6.1, Biopython 1.83; DefenseFinder 2.0.0, geNomad 1.8.0 and
CheckV 1.0.3 from their source code (not installed).

Run: `bash tests/run_leaf.sh crossdomain` → `CROSSDOMAIN TESTS PASSED` (36 tests, ~100 s).

**Findings: 45** (3 critical, 29 major, 13 minor), CD-01 … CD-45.

| ID | Sev. | File | Issue (short) | Test(s) |
|---|---|---|---|---|
| CD-01 | critical | crossdomain.smk | cross-domain tree replaced by `();` (cp onto itself) | e2e_crossdomain_tree_and_analysis, rules_no_swallowed_failures_or_self_copy |
| CD-02 | critical | crossdomain.smk | tool failures swallowed (`\|\| true`, `\|\| cp`, Pfam guard) | e2e_crashed_tools_fail_their_rules, rules_no_swallowed_failures_or_self_copy |
| CD-03 | critical | crossdomain.smk | bacterial validation tree has no supports → tier A unreachable in bacteria | e2e_bacterial_branch, rules_parity_with_archaeal_calls |
| CD-04 | major | crossdomain.smk | default config (no bacterial CDS) aborts / archaea-only "cross-domain" dN/dS | e2e_no_bacterial_cds_selection_not_applicable |
| CD-05 | major | crossdomain.smk | `python` run inside the R-only conda env (niche) | rules_no_swallowed_failures_or_self_copy, e2e_niche_and_viral_modules |
| CD-06 | major | crossdomain.smk | validation universe taken from the QC table | e2e_bacterial_branch |
| CD-07 | major | extract_bacterial_prometheus.py | substring system matching | extract_exact_matching_split_and_ids |
| CD-08 | major | extract_bacterial_prometheus.py | genes table never split; no tables for zero-system genomes; stale tables | extract_exact_matching_split_and_ids, extract_stale_tables_removed_on_rerun |
| CD-09 | major | extract_bacterial_prometheus.py | AntiDefenseFinder rows counted as defense | extract_exact_matching_split_and_ids, e2e_bacterial_branch |
| CD-10 | major | extract_bacterial_prometheus.py | missing proteome/protein silently skipped; gzip read as text | extract_missing_proteome_fails_unless_allowed, extract_exact_matching_split_and_ids |
| CD-11 | major | extract_bacterial_prometheus.py | NCBI `cds_from_genomic` CDS never matched | extract_exact_matching_split_and_ids, e2e_codon_alignment_keeps_both_domains |
| CD-12 | major | extract / build_bacterial_table | genome ids not validated (blank, duplicate, `__`, zero overlap) | extract_empty_and_mismatch, table_duplicates_and_orphans_fail |
| CD-13 | major | prune_tree_to_taxa.py | GTDB tip matching / no relabel / silent trivial tree | prune_normalisation_collapse_and_rooting, prune_low_match_fails_and_empty_is_trivial, prune_gene_tree_exact_mode |
| CD-14 | major | crossdomain_analysis.py | rooting on a non-monophyletic outgroup; supports misassigned on reroot | analysis_nonmonophyletic_outgroup_unknown_and_trivial, analysis_fitch_rule_and_parser |
| CD-15 | major | crossdomain_analysis.py | direction from the ML Fitch root state alone | analysis_nested_clade_is_robust_bacterial_origin, analysis_sister_clade_is_ambiguous, analysis_archaeal_root_and_imbalance_control, analysis_bootstrap_direction_support |
| CD-16 | major | crossdomain_analysis.py | archaeal acquisitions counted on unsupported splits | analysis_support_collapse_counts |
| CD-17 | major | crossdomain_analysis.py | dendropy missing / unclassified tips handled silently | analysis_nonmonophyletic_outgroup_unknown_and_trivial |
| CD-18 | major | crossdomain.smk | reference ids not sanitised → label lookups fail after IQ-TREE | e2e_bacterial_branch |
| CD-19 | major | crossdomain_analysis.py | no synonymous-saturation control for cross-domain dN/dS | analysis_saturation_ng86, e2e_crossdomain_tree_and_analysis |
| CD-20 | major | crossdomain.smk | HyPhy gene tree ≠ codon-alignment taxa | e2e_codon_alignment_keeps_both_domains, e2e_hyphy_accepts_selection_inputs |
| CD-21 | major | build_bacterial_table.py | `pass_primary` ignores contamination | table_known_values, e2e_bacterial_branch |
| CD-22 | major | build_bacterial_table.py | taxonomy formats, duplicates, universe; coupling to build_master_table internals | table_known_values, table_duplicates_and_orphans_fail |
| CD-23 | major | build_bacterial_niche_master.py | unknown niche coded 0; missing niche file passed through | niche_master_missing_is_not_absent, niche_master_failures |
| CD-24 | major | plot_crossdomain_extra.py, plot_viral_crossdomain.py | Contract-v2 presence rule ignored (tier-agnostic prevalence) | plot_comparison_presence_rule_regression, plot_viral_presence_and_undefined_fraction |
| CD-25 | major | plot_viral_crossdomain.py | 0/0 R-loop fraction plotted as 0; unlabelled p | plot_viral_presence_and_undefined_fraction |
| CD-26 | major | crossdomain.smk | geNomad/CheckV placeholders with wrong headers mask crashes | e2e_niche_and_viral_modules, e2e_crashed_tools_fail_their_rules |
| CD-27 | major | extract / crossdomain.smk | bacterial fragmentation by length proxy, archaeal by Prodigal flags | extract_prodigal_genes_tables, e2e_bacterial_branch |
| CD-28 | major | crossdomain.smk | bacterial validation lacks the archaeal evidence/thresholds (motif alignments, focal system, placement cut-offs) | e2e_bacterial_branch, rules_parity_with_archaeal_calls |
| CD-29 | major | crossdomain.smk | tier sets/thresholds not passed to tier-based responses; no A+B sensitivity fits | rules_parity_with_archaeal_calls, e2e_niche_and_viral_modules, plot tests |
| CD-41 | major | build_bacterial_table.py | genome_size blank for compressed assemblies (decode error swallowed) | table_genome_size_from_assembly |
| CD-43 | major | crossdomain.smk | geNomad resumed stale intermediates (no `--restart`); no `--splits`; DB check = any directory | rules_parity_with_archaeal_calls, e2e_crashed_tools_fail_their_rules |
| CD-44 | major | crossdomain.smk | CheckV resumed a stale tmp/ (no `--restart`); DB check = any directory; no post-run table/header check | e2e_crashed_tools_fail_their_rules, rules_parity_with_archaeal_calls |
| CD-30 | minor | crossdomain.smk | `grep -c … \|\| echo 0` | rules_no_swallowed_failures_or_self_copy, e2e_all_empty_scenario |
| CD-31 | minor | crossdomain.smk | 11 rules without `log:` | rules_hygiene, e2e_niche_and_viral_modules |
| CD-32 | minor | crossdomain.smk | GTDB tree passed as `params` (no dependency/missing-file check) | e2e_bacterial_branch |
| CD-33 | minor | plot_crossdomain_extra.py | hard-coded domain list incl. `nuclease` | plot_comparison_presence_rule_regression |
| CD-34 | minor | plot scripts | `plotting.formats` ignored (png+svg hard-coded) | plot tests, e2e_crossdomain_tree_and_analysis |
| CD-35 | minor | plot_crossdomain_tree.py | any error → placeholder; unrooted tree drawn | plot_tree_figure_and_placeholder |
| CD-36 | minor | crossdomain.smk | viral summary `--samples` = QC table; no log | e2e_niche_and_viral_modules |
| CD-37 | minor | crossdomain.smk | hmmscan alignments dumped into the log | e2e_bacterial_branch |
| CD-38 | minor | crossdomain.smk | no empty-output checks after MAFFT/ClipKIT/IQ-TREE (ClipKIT exits 0 on errors) | e2e_crashed_tools_fail_their_rules |
| CD-39 | minor | crossdomain.smk | missing per-genome systems table counted as zero systems | e2e_bacterial_branch |
| CD-40 | minor | crossdomain.smk | codon-alignment interface: bacterial-only `--proteins`, whole archaeal CDS copied | e2e_codon_alignment_keeps_both_domains, rules_no_swallowed_failures_or_self_copy |
| CD-42 | minor | crossdomain.smk | cross-domain selection without `--seed` (HyPhy RNG seed) | rules_parity_with_archaeal_calls |
| CD-45 | minor | tests/crossdomain/run_tests.py | parity canary did not cover the geNomad/CheckV invocations | rules_parity_with_archaeal_calls |

---

## Critical

### CD-01 — cross-domain tree always replaced by `();`
- **Location:** crossdomain.smk:312, 319, 328-330 (`crossdomain_tree`).
- **Issue:** `output.tree` = `{CD}/tree/crossdomain.treefile` = the file IQ-TREE writes
  (`--prefix {params.workdir}/crossdomain`). The rule then ran
  `cp {params.workdir}/crossdomain.treefile {output.tree} 2>/dev/null || echo "();" > {output.tree}`.
- **Evidence:** `cp x x` exits 1 ("are identical (not copied)"; GNU: "are the same file"),
  reproduced under `bash -c 'set -euo pipefail; …'`, so the `||` branch overwrote the real
  tree on every run. Every cross-domain analysis therefore ran on `();`
  (`crossdomain_summary` status `incomplete`).
- **Fix:** no copy: IQ-TREE writes the declared outputs directly (+ `--wbtl` → `.ufboot`, all
  UFBoot trees for the robustness analysis), followed by explicit non-empty checks.
- **Test:** e2e_crossdomain_tree_and_analysis (10-tip tree, 1000 UFBoot trees, status
  complete); rules_no_swallowed_failures_or_self_copy (pattern can't return).

### CD-02 — tool failures swallowed
- **Location:** crossdomain.smk:74-78 (`[ -s pfam ]` guard, `hmmscan … || true`, `# no domains`
  placeholder), 105/301 (`clipkit … || cp` untrimmed fallback), 107/109 and 329 (`iqtree2 … ||
  true` + `();`), 526-529 (`genomad … || true`), 550-552 (`checkv … || true`).
- **Issue:** a crashed or misconfigured tool produced a plausible-looking empty result. A failed
  hmmscan gave every bacterial candidate "no domains", so all were Rejected and bacterial
  prevalence came out 0. A missing Pfam file skipped the scan silently. A ClipKIT error sent
  an untrimmed alignment to IQ-TREE.
- **Evidence:** e2e with failing `iqtree2`/`hmmscan` stubs on PATH: the pristine rules complete
  and write placeholders; the fixed rules exit non-zero and leave no output.
- **Fix:** no `|| true`/`|| cp`/`2>/dev/null`; placeholders only for genuinely empty input
  (zero candidates, < 4 sequences, empty viral FASTA), each marked `# nonfatal:`; a missing
  Pfam index or geNomad/CheckV database is a clear error.
- **Test:** e2e_crashed_tools_fail_their_rules; rules_no_swallowed_failures_or_self_copy;
  e2e_all_empty_scenario (empty inputs still complete).

### CD-03 — bacterial validation tree without supports (tier A unreachable in bacteria)
- **Location:** crossdomain.smk:106 (`iqtree2 -m LG+G -fast`, no `-B`/`--alrt`, no `--seed`).
- **Issue:** `validate_prometheus.py` (Contract v2) places a candidate in the Prometheus clade
  only when SH-aLRT ≥ `placement_min_alrt` (80) and UFBoot ≥ `placement_min_ufboot` (95). The
  archaeal `prometheus_validation_tree` now uses `-m MFP -B 1000 --alrt 1000 --bnni --seed`.
  A support-less bacterial tree makes every bacterial locus "unsupported", so tier A was
  unreachable for bacteria only. The cross-domain tier-A prevalence comparison would then
  read "absent in bacteria" by construction.
- **Evidence:** same e2e inputs, before vs after: with the `-fast` tree all 4 bacterial loci
  were `tree_assignment=Unplaced`, `reason=placement:unsupported_ingroup(<80…)`, tier C. With
  the supported tree they are `tree_assignment=Prometheus`, tier B (A is blocked only by the
  missing REase in the 650-aa test windows).
- **Fix:** `bact_validation_tree` mirrors the archaeal rule exactly: `combine_proteomes.py`
  reference mode (same id sanitisation and collision checks), MAFFT, ClipKIT kpic-smart-gap
  with output checks, IQ-TREE `-m {prometheus.validation_tree_model} -B {tree.ufboot} --alrt
  {tree.alrt} --bnni --seed`.
- **Test:** e2e_bacterial_branch (SH-aLRT/UFBoot labels present in the validation tree);
  rules_parity_with_archaeal_calls (the iqtree2/mafft/clipkit/combine_proteomes options are a
  superset of the archaeal call; `-fast` absent).

## Major

### CD-04 — no bacterial CDS (the config default) breaks or mislabels cross-domain dN/dS
- **Location:** crossdomain.smk:54 (`[ -f {output.cds} ] || : > {output.cds}`), 380-402;
  config.yaml `bacteria.cds_dir: ""`; common.smk `crossdomain_targets` always requests
  `crossdomain_selection_summary.tsv`.
- **Issue:** without bacterial CDS the codon alignment held archaea only, yet it was reported as
  the "cross-domain" selection test. With the revised `codon_align.py` (fails when more than
  `--max-drop-frac` of the expected loci lack a CDS, codon_align.py:255-263) the default
  configuration aborts the `crossdomain` target instead.
- **Fix:** `BACT_CDS = bool(bacteria.cds_dir)`. Without bacterial CDS, `crossdomain_codon_align`
  writes an empty codon alignment and a header-only report (`# nonfatal:`, reason logged).
  `run_selection.py` is then called with `--run False`, and the `note` column states
  `not run: bacteria.cds_dir not set (cross-domain dN/dS needs bacterial CDS)`. The column is
  located by header name, not by position.
- **Test:** e2e_no_bacterial_cds_selection_not_applicable.

### CD-05 — `python` run inside the R-only environment
- **Location:** crossdomain.smk:473-506 (`bact_niche`, `conda: rstats.yaml`).
- **Issue:** `python build_bacterial_niche_master.py` ran in `envs/rstats.yaml`, which has r-base
  and R packages only, with no python and no pandas. Under `--use-conda` the rule cannot run. The
  niche master was an undeclared output (a `params` path) and its stdout/stderr were not logged.
- **Fix:** split into `bact_niche_master` (python env, declared output, log) and `bact_niche`
  (R only).
- **Test:** rules_no_swallowed_failures_or_self_copy (no `python` in `bact_niche`);
  e2e_niche_and_viral_modules.

### CD-06 — validation universe taken from the QC table
- **Location:** crossdomain.smk:121, 141 (`meta=_bact("qc")`, `--meta {input.meta} --reps {input.meta}`).
- **Issue:** the genome universe (detection sensitivity, which genomes count as "absent") came
  from whatever QC table was configured. CheckM2 output, for instance, names the column `Name`,
  not `genome_id`. It also differed from the genome set the per-genome rules expand over.
- **Fix:** `bact_samples` writes `bacteria/bacterial_samples.tsv` (genome_id + QC columns) =
  `BACT_GENOMES`, and `bacterial_representatives.txt` in the archaeal representatives.txt
  format (one id per line). Both feed `--meta`/`--reps`, the bacterial table, the defense
  split and the viral summary.
- **Test:** e2e_bacterial_branch.

### CD-07 — substring system matching
- **Location:** extract_bacterial_prometheus.py:122-125 (`if sl in t or sl in st`).
- **Issue:** any type/subtype containing the configured name was counted as Prometheus (e.g.
  `Prometheus_II`).
- **Fix:** exact, case-insensitive equality with one of the configured names (comma list allowed).
- **Test:** extract_exact_matching_split_and_ids (Prometheus_II excluded, PROMETHEUS included).

### CD-08 — DefenseFinder genes table never split; tables missing or stale
- **Location:** extract_bacterial_prometheus.py:43-60, 111.
- **Issue:** only the systems table was split. `validate_prometheus.py` reads
  `{g}_defense_finder_genes.tsv` for the defense-island context and to exclude the locus' own
  system, so bacterial loci never had context (a tiering input). Genomes with zero systems got
  no table. Tables from a previous run survived.
- **Fix:** both tables split for every genome of the universe (header-only when empty; genome
  column dropped, all others kept); previous per-genome tables removed first.
- **Test:** extract_exact_matching_split_and_ids; extract_stale_tables_removed_on_rerun.

### CD-09 — AntiDefenseFinder rows counted as defense systems
- **Location:** extract_bacterial_prometheus.py (no `activity` filter).
- **Issue:** DefenseFinder ≥ 2.0 with `-a` writes anti-defense systems with `activity ==
  Antidefense` into the same table (best_solution.py:49-57; main.py:107). The archaeal run
  never calls them, so bacterial defense load and co-occurrence were inflated.
- **Fix:** Antidefense rows dropped from the split tables (counted per genome in the
  extraction report); the bacterial table counts `activity != Antidefense` only.
- **Test:** extract_exact_matching_split_and_ids; e2e_bacterial_branch (Anti_CBASS absent from
  the defense matrix).

### CD-10 — missing proteomes/proteins skipped silently; gzip read as text
- **Location:** extract_bacterial_prometheus.py:150-155, 171-182, 196-198; utils.read_fasta (plain `open`).
- **Issue:** a Prometheus genome without a proteome, or a hit_id absent from it, was dropped with
  only a warning. That shrinks bacterial prevalence silently. `faa.gz` was opened as text,
  so a compressed proteome looked empty.
- **Fix:** fatal unless `bacteria.allow_missing_proteins: true` (then reported per genome in
  `prometheus/extraction_report.tsv`); gzip-aware readers; duplicated ids are errors.
- **Test:** extract_missing_proteome_fails_unless_allowed; extract_exact_matching_split_and_ids (g1.faa.gz).

### CD-11 — NCBI CDS never matched
- **Location:** extract_bacterial_prometheus.py:188-198 (CDS looked up by first header token only).
- **Issue:** NCBI `cds_from_genomic` headers are `lcl|<acc>_cds_<protein_id>_<n> [protein_id=…]`.
  The first token never equals the protein accession DefenseFinder reports, so NCBI-annotated
  genomes silently lost all their CDS and dropped out of the dN/dS.
- **Fix:** lookup by token, then by the `[protein_id=…]` tag; CDS length checked against 3×aa.
- **Test:** extract_exact_matching_split_and_ids (CDS via tag); e2e_codon_alignment_keeps_both_domains
  (the WP_-annotated genome is codon-aligned).

### CD-12 — genome ids not validated
- **Location:** extract_bacterial_prometheus.py:53-54; build_bacterial_table.py (whole file).
- **Issue:** blank ids were skipped. Ids containing `__` (utils.ID_SEP) or unsafe characters
  break the `genome__locus` round trip (`make_protein_id`/`local_of`). Duplicated rows and a
  zero overlap between the DefenseFinder tables and the universe went unnoticed.
- **Fix:** fatal checks with examples (blank/duplicate/unsafe ids, `__`, zero overlap, catalog
  genomes outside the universe).
- **Test:** extract_empty_and_mismatch; table_duplicates_and_orphans_fail.

### CD-13 — GTDB tree pruning: matching, relabelling, silent trivial trees
- **Location:** prune_tree_to_taxa.py:46, 54, 63, 93 (trivial tree on empty input, missing
  dendropy, parse or retain failure); 67-71 (`replace("GB_", "").replace("RS_", "")`, anywhere
  in the string); 77-81 (zero/low overlap only logged); 85 (kept tips not relabelled).
- **Issue:** pruned tips kept GTDB labels (`RS_GCF_…`) while the master's `tree_tip` is the
  genome id, so the phylogenetic R tests could not match tips. Version, paired GCA↔GCF and
  non-representative genomes were unmatched. Any error gave `();`.
- **Fix:** tiered matching (exact → prefix_stripped → version_stripped → paired_accession →
  species_representative), ambiguous keys never used, one genome per tip (best tier, then
  quality), tips relabelled to genome_id, per-genome `--out-report`, fatal below
  `--min-match-fraction` (default 0.9). `--match exact` mode is used for gene trees.
- **Test:** prune_normalisation_collapse_and_rooting; prune_low_match_fails_and_empty_is_trivial;
  prune_gene_tree_exact_mode; e2e_bacterial_branch (paired_accession, prefix_stripped tips).

### CD-14 — rooting
- **Location:** crossdomain_analysis.py:178-195.
- **Issue:** the tree was rerooted at the MRCA of the paralog tips. If the outgroup is not
  monophyletic, that MRCA contains ingroup taxa and the root lands inside the ingroup. The
  fallback to midpoint was silent. dendropy keeps IQ-TREE support labels on nodes, so after
  rerooting the labels along the old root path describe different bipartitions.
- **Fix:** own newick reader with edge-attached supports. The tree is rooted on the edge
  separating all Paralog_ref tips from the ingroup; if no such edge exists,
  `outgroup_monophyletic=False` and no direction is called. A merged root bifurcation keeps
  the more conservative support.
- **Test:** analysis_nonmonophyletic_outgroup_unknown_and_trivial; analysis_fitch_rule_and_parser
  (reader == dendropy tip set and label count on the e2e tree).

### CD-15 — transfer direction from a single parsimony state
- **Location:** crossdomain_analysis.py:200-216.
- **Issue:** "bacteria→archaea" was called from the Fitch state at the ingroup root of the ML
  tree. With bacteria far outnumbering archaea that is the expected outcome under sampling
  alone, and no support or uncertainty entered the call.
- **Fix:** a direction is called only when every criterion holds: clean outgroup split; ML root
  state; the same state on the support-collapsed tree (SH-aLRT < 80 or UFBoot < 95 collapsed);
  ≥ 95% of the UFBoot replicate trees; ≥ 95% of balanced subsamples (larger domain thinned to
  the smaller, seeded). Otherwise the call is `ambiguous`, with the failed criteria listed.
  Midpoint rooting and excluding the references are reported as sensitivity analyses;
  `direction_method=topology-based`.
- **Test:** analysis_nested_clade_is_robust_bacterial_origin; analysis_sister_clade_is_ambiguous;
  analysis_archaeal_root_and_imbalance_control; analysis_bootstrap_direction_support.

### CD-16 — archaeal acquisitions counted on unsupported splits
- **Location:** crossdomain_analysis.py:82-104.
- **Issue:** maximal archaeal-only clades on the fully resolved ML tree; unsupported splits
  inflate the "independent acquisitions" count.
- **Fix:** adds `n_archaeal_clades_supported` (minimum archaeal groups compatible with the
  supported topology), `…_nested_supported`, `n_crossdomain_transitions_supported`, plus a
  per-clade table (`crossdomain_archaeal_clades.tsv`).
- **Test:** analysis_support_collapse_counts.

### CD-17 — silent degradation in the analysis
- **Location:** crossdomain_analysis.py:134-135 (no dendropy → placeholder summary), 155-163
  (unlabelled tips → `Unknown`, silently out of the character).
- **Fix:** no dendropy dependency; every tip must be classified (archaeal/bacterial protein ids
  or reference labels), otherwise error.
- **Test:** analysis_nonmonophyletic_outgroup_unknown_and_trivial.

### CD-18 — reference ids not sanitised
- **Location:** crossdomain.smk:99-101 (validation tree), 295-297 (cross-domain alignment): raw FASTA ids into MAFFT/IQ-TREE and raw ids in the label files.
- **Issue:** IQ-TREE 2.3.6 rewrites `:` `,` `(` `)` etc. in names (verified: `NAME:x,y` →
  `NAME_x_y`) but keeps `|`. Label lookups keyed on raw ids then miss the renamed tips, and the
  references stop defining the Prometheus/paralog clades.
- **Fix:** one sanitisation (`[^A-Za-z0-9_.-]` → `_`, = utils.sanitize) for FASTA and labels
  (validation tree via `combine_proteomes.py`), duplicate ids after sanitising are errors.
- **Test:** e2e_bacterial_branch (`sp_P51530_DNA2_HUMAN` in labels and tree).

### CD-19 — no synonymous-saturation control for cross-domain dN/dS
- **Location:** crossdomain.smk:380-441 (selection run without any divergence check).
- **Issue:** archaeal and bacterial ProA are separated by long branches; saturated dS makes
  dN/dS uninterpretable.
- **Evidence:** Biopython 1.83 `codonseq.cal_dn_ds(method="NG86")` returns -1 when dS is not
  estimable (saturated), verified.
- **Fix:** `crossdomain_analysis.py --codon` computes NG86 pairwise dS (within Archaea, within
  Bacteria, cross-domain; seeded pair sample) into `selection/crossdomain_saturation.tsv`. The
  summary gets `dnds_saturation_warning` when > 50% of cross-domain pairs have dS undefined or
  > `saturation_ds` (3.0).
- **Test:** analysis_saturation_ng86; e2e_crossdomain_tree_and_analysis.

### CD-20 — HyPhy gene tree does not match the codon alignment
- **Location:** crossdomain.smk:403-441 (`gene_tree = crossdomain_tree` incl. reference/paralog
  tips without CDS); run_selection.py:464-487 prunes only when `intact_only` is True.
- **Fix:** `crossdomain_selection_tree` prunes the cross-domain tree to exactly the
  codon-aligned taxa (`prune_tree_to_taxa.py --match exact --min-match-fraction 1.0`).
- **Test:** e2e_codon_alignment_keeps_both_domains (tips == codon ids);
  e2e_hyphy_accepts_selection_inputs (real HyPhy FEL run accepts the pair).

### CD-21 — `pass_primary` ignores contamination
- **Location:** build_bacterial_table.py:131-136.
- **Issue:** bacteria passed on completeness ≥ 90 alone, while the archaeal set also requires
  contamination ≤ 5. Prevalence was therefore compared at unmatched quality, and a contaminated
  MAG can carry a foreign ProA.
- **Fix:** `--max-contamination` (rule: `bacteria.max_contamination`, default 5.0); unknown
  completeness, or unknown contamination when the column exists, does not pass.
- **Test:** table_known_values; e2e_bacterial_branch (contamination 7.5 → False).

### CD-22 — taxonomy, duplicates, universe; coupling to another leaf's internals
- **Location:** build_bacterial_table.py:21 (`import build_master_table as bmt`), 33-56 (`load_qc`), 81-159 (`main`).
- **Issue:** the table relied on private helpers of build_master_table.py, whose API changed
  during this review (parse_taxonomy now raises without user_genome/classification). Tables
  with `gtdb_*` columns were not parsed, duplicated ids were not detected, and the rows were
  qc ∪ taxonomy instead of the analysis universe.
- **Fix:** self-contained, same definitions as the archaeal master (documented in the module
  docstring); accepts GTDB-Tk summaries, genome_id+classification or gtdb_* tables; duplicates
  and missing per-genome defense tables are errors; rows = `--samples`.
- **Test:** table_known_values; table_duplicates_and_orphans_fail.

### CD-23 — unknown niche coded as "not in niche"
- **Location:** build_bacterial_niche_master.py:44-48 (missing file → master passed through), 70 (`fillna(0)`).
- **Issue:** genomes without a label became 0 ("not anaerobic"), which biases the phylogenetic
  regression toward the null. A mistyped niche path silently produced a test with no niche
  columns.
- **Fix:** NA for missing labels, unannotated genomes dropped, duplicate rows unioned, slug
  collisions merged, missing `--niche` file or zero overlap fatal.
- **Test:** niche_master_missing_is_not_absent; niche_master_failures.

### CD-24 — figures ignore the Contract-v2 presence rule
- **Location:** plot_crossdomain_extra.py:54-63; plot_viral_crossdomain.py:55-61.
- **Issue:** prevalence counted `prometheus_present` (any tier, incl. C), while the archaeal
  numbers are tier-A. A B/C-only genome is not a confident absence at tier A.
- **Fix:** `--primary-tiers` (config `stats.confidence_tiers_primary`): present if the best tier
  is in T, absent if the tier is empty, excluded otherwise (left out of the denominator;
  counts shown on the axis/panel labels); `--out-table` writes the plotted numbers.
- **Test:** plot_comparison_presence_rule_regression (B-only genome excluded at tier A, present
  at A+B); plot_viral_presence_and_undefined_fraction.

### CD-25 — undefined R-loop fraction plotted as 0
- **Location:** plot_viral_crossdomain.py:55-66 (`merge`: fraction used as given), 84-94 (p label).
- **Issue:** `frac_rloop_proviruses` for a genome with no provirus is 0/0. The pristine
  aggregate_viral wrote 0, which piles provirus-free genomes into the lowest bin. The title
  printed `p_adj` or `p_value` without saying which.
- **Fix:** genomes with `n_proviruses == 0` excluded from the fraction panels (works with the
  current blank encoding and the old 0); the label names `p_adj`/`p`.
- **Test:** plot_viral_presence_and_undefined_fraction.

### CD-26 — geNomad/CheckV placeholders mask crashes
- **Location:** crossdomain.smk:526-529, 549-552.
- **Issue:** `|| true` plus placeholder tables whose headers do not match the tools. The
  geNomad 1.8.0 virus summary has 11 columns (summary.py:669-672; it is written even with zero
  viruses, lines 501-503); the placeholder had 5. CheckV 1.0.3 `quality_summary.tsv` has 14
  columns (quality_summary.py:168-183); the placeholder had 4. A crash looked like "no viruses".
- **Fix:** tool failures propagate; missing databases are explicit errors. CheckV on an empty
  FASTA exits 0 without output (utility.check_fasta), so only that genuinely empty case writes
  the full 14-column header (`# nonfatal:`).
- **Test:** e2e_niche_and_viral_modules (14 columns for a virus-free genome);
  e2e_crashed_tools_fail_their_rules (missing geNomad DB fails with a clear message).

### CD-27 — bacterial fragmentation by length proxy only
- **Location:** extract_bacterial_prometheus.py (no gene coordinates); crossdomain.smk (`bact_validate` without `--genes-dir`).
- **Issue:** archaeal loci are called fragmented from Prodigal's `partial` code (SCHEMA item 2,
  validate_prometheus.py:909-916). Bacterial loci fell back to `length < min_len`, which misses
  e.g. a 1,500-aa ProA cut at a contig edge. The result is a domain-specific bias in tier A.
- **Fix:** `--out-genes-dir` writes `bacteria/genes/{g}/{g}.genes.tsv` (rename_prodigal.py
  layout) from Prodigal headers (GTDB/GTDB-Tk proteomes) for genomes with hits.
  `bact_validate` passes `--genes-dir`. A proteome with any non-Prodigal header gets no table:
  validate_prometheus rejects partially covered tables, so it uses the length proxy instead,
  reported as `genes_table=none`.
- **Test:** extract_prodigal_genes_tables (layout == rename_prodigal.GENES_COLS, read back by
  validate_prometheus.load_genes_table); e2e_bacterial_branch (`prodigal_partial=10` →
  fragmented).

### CD-28 — bacterial validation with less evidence than archaeal
- **Location:** crossdomain.smk:113-151 (`bact_validate`), 61-81 (`bact_domain_scan`).
- **Issue:** the revised archaeal `validate_prometheus` call passes `--hmmscan-out` (residue-level
  Walker A / motif II / nuclease-core checks), `--system`, `--min-alrt/--min-ufboot` and
  `--genes-dir`. Without these, bacterial catalytic status and tiers come from different
  evidence and thresholds.
- **Fix:** `bact_domain_scan` writes `--notextw -o domains.hmmscan.txt` (checks the pressed Pfam
  index first) and `bact_validate` passes the same options (focal system = first configured
  bacterial name).
- **Test:** e2e_bacterial_branch (`== domain` alignment blocks); rules_parity_with_archaeal_calls.

### CD-29 — tier sets and thresholds not passed to tier-based responses
- **Location:** crossdomain.smk:178-194 (`bact_prevalence`), 253-277 (`bact_cooccurrence`),
  473-506 and 599-625 (ecology fits), 556-577 (`bact_rloop`), 380-441 (codon/selection).
- **Issue:** the bacterial calls relied on script defaults. A non-default
  `stats.confidence_tiers_primary`, `cooccurrence_min_count`, `phyloglm_btol`,
  `ecology.covariates`, R-loop model/score, codon table or GARD setting reached archaea but not
  bacteria, making the comparison inconsistent (Contract v2 item 6). `",".join(cfg(...))`
  would also turn a string `"A,B"` into `"A,,,B"`. The archaeal ecology rule re-fits at A+B;
  the bacterial rules did not.
- **Fix:** every option of the archaeal call is passed with the same config key (helper
  `_tiers`/`_csv` accepts lists or strings); `bact_niche`/`bact_ecology_stats` add the A+B
  re-fit (`*_sensitivity.tsv`); `bact_prevalence` passes `--sensitivity-tiers`.
- **Test:** rules_parity_with_archaeal_calls (reads the archaeal rules and fails on any option
  the bacterial call lacks, reasoned exemptions only: `--predictors`, `--viral` for the niche
  fit, `--relax-labelling`); e2e_niche_and_viral_modules (sensitivity outputs, tier log lines).

### CD-41 — `genome_size` silently blank for compressed assemblies
- **Location:** build_bacterial_table.py:57-78 (`assembly_size`: `open(path)` in text mode inside
  `try … except Exception: return ""`).
- **Issue:** when the QC table has no size column, genome size comes from the assembly. A
  gzip-compressed assembly (`bacteria.genome_extension: fna.gz`, which geNomad reads) raised a
  decode error that was swallowed. Every such genome got `genome_size = ""`, silently removing
  the genome-size covariate from the size-adjusted bacterial models (co-occurrence, ecology).
- **Fix:** gzip-aware streaming count (`.gz` variants also searched); no blanket `except`.
- **Test:** table_genome_size_from_assembly (plain 150 bp, gzip 131 bp, missing → '').

### CD-43 — geNomad: stale intermediates resumed, no memory control, weak DB check
- **Location:** crossdomain.smk:507-532 (`bact_genomad`; unchanged by my first revision).
- **Issue:** without `--restart`, geNomad 1.8.0 skips every step whose output already exists
  (annotate.py:138-196, find_proviruses.py:521, aggregated_classification.py:149-262). A re-run
  after a failed or interrupted run therefore silently resumes stale intermediates. There was no
  `--splits` for the MMseqs2 marker search (memory), and any existing directory passed as a
  database.
- **Evidence:** geNomad 1.8.0 source, cli.py:350-418 (`end-to-end INPUT OUTPUT DATABASE`;
  `--cleanup`, `--restart`, `--threads`, `--splits` IntRange(min=0) default 0). The archaeal
  `ecology.smk:genomad` (final) uses `--restart --splits {ecology.genomad.splits}` and checks
  `{db}/version.txt`.
- **Fix:** mirrors the archaeal rule: `version.txt` check with the same message, `--cleanup
  --restart --threads --splits {ecology.genomad.splits, default 0}`.
- **Test:** rules_parity_with_archaeal_calls (geNomad options); e2e (the fake geNomad parses
  options like geNomad 1.8.0 and rejects unknown ones); e2e_crashed_tools_fail_their_rules
  (missing database → explicit error).

### CD-44 — CheckV: stale tmp/ resumed, weak DB check, silent exit 0 not caught
- **Location:** crossdomain.smk:533-555 (`bact_checkv`).
- **Issue:** CheckV 1.0.3 deletes its `tmp/` only with `--restart` (end_to_end.py:53-57), so a
  re-run resumed stale intermediates. Without `--remove_tmp` (end_to_end.py:64-65) the
  intermediates stay. Any existing directory passed as the database. CheckV also exits 0
  without results on some fatal errors (empty or duplicate-id FASTA, missing executables);
  Snakemake only notices a *missing* table, not a malformed one.
- **Fix:** mirrors the archaeal rule: `genome_db`/`hmm_db` check with the same message,
  `--restart --remove_tmp`, then a required `quality_summary.tsv` whose first line equals the
  14-column CheckV 1.0.3 header (`envs/viral.yaml` pins checkv=1.0.3). The header-only table is
  still written only for an empty viral FASTA.
- **Test:** e2e_crashed_tools_fail_their_rules (CheckV exiting 0 without a table, a 4-column
  header, an incomplete database → rule fails with an explicit message); rules_parity_with_archaeal_calls
  (options, and the bacterial header == ecology.smk `CHECKV_QS_HEADER`).

## Minor

### CD-30 — `grep -c '>' … || echo 0`
- **Location:** crossdomain.smk:102, 298, 326.
- **Issue:** `grep -c` prints `0` and exits 1 when nothing matches, so `n` became `"0\n0"`. The
  following `[ "$n" -ge 4 ]` then fails with "integer expression expected" and silently takes
  the placeholder branch. The command also counts `>` anywhere in a line.
- **Fix:** awk record counter (`_AWK_COUNT`, counts header lines only).
- **Test:** rules_no_swallowed_failures_or_self_copy (`|| echo` banned); e2e_all_empty_scenario.

### CD-31 — rules without a log
- **Location:** bact_validate, bact_table, bact_prevalence, bact_architecture,
  bact_defense_matrix, crossdomain_analysis, plot_crossdomain, plot_crossdomain_extra,
  bact_rloop, bact_viral_summary, plot_viral_crossdomain (11 rules).
- **Fix:** every rule has `log:` and redirects its tools into it.
- **Test:** rules_hygiene (`tests/check_rules.py` → RULES HYGIENE OK); e2e_niche_and_viral_modules (log files exist).

### CD-32 — GTDB reference tree passed as `params`
- **Location:** crossdomain.smk:223 (`ref_tree=_bact("gtdb_bac_tree")`).
- **Issue:** no missing-file check by Snakemake and no rerun when the tree changes.
- **Fix:** declared input; per-genome tip-matching report as an output.
- **Test:** e2e_bacterial_branch.

### CD-33 — hard-coded domain categories
- **Location:** plot_crossdomain_extra.py:27 (`CANON = [... "nuclease", "REase"]`).
- **Issue:** SCHEMA item 3: categories come only from the domain spec, and there is no
  `nuclease` category.
- **Fix:** categories from `--domains` (utils canonical order) or the `has_*` columns present.
- **Test:** plot_comparison_presence_rule_regression.

### CD-34 — `plotting.formats` ignored
- **Location:** plot_crossdomain_tree.py:76, 120 (`for ext in ("png", "svg")`); same in the
  other two plot scripts; the rules declare `expand(..., ext=PLOT_FORMATS)` (crossdomain.smk:363,
  450, 636).
- **Issue:** with e.g. `formats: [pdf]` the rules fail with MissingOutputException.
- **Fix:** `--formats` from `PLOT_FORMATS`.
- **Test:** plot tests; e2e_crossdomain_tree_and_analysis (png + svg as configured).

### CD-35 — tree figure: silent placeholder, unrooted drawing
- **Location:** plot_crossdomain_tree.py:75-90, 125-127.
- **Issue:** any exception produced a placeholder figure (logged only). The ML tree was drawn
  as stored (unrooted), not with the rooting the analysis used.
- **Fix:** draws the outgroup-rooted tree written by crossdomain_analysis; placeholder only for a
  trivial tree; other errors propagate; supported nodes marked; labels only up to `--max-labels`.
- **Test:** plot_tree_figure_and_placeholder.

### CD-36 — viral summary universe and log
- **Location:** crossdomain.smk:578-598 (`samples=_bact("qc")`, no `log:`).
- **Fix:** `--samples bacterial_samples.tsv` (= the genomes the per-genome viral rules run on); logged.
- **Test:** e2e_niche_and_viral_modules (one row per universe genome).

### CD-37 — hmmscan alignments written into the log
- **Location:** crossdomain.smk:75-76 (`hmmscan … > {log} 2>&1`, no `-o`).
- **Issue:** per-target alignments (large) went to the log and were not available to validation.
- **Fix:** `-o domains.hmmscan.txt --notextw` (consumed by `validate_prometheus --hmmscan-out`, CD-28).
- **Test:** e2e_bacterial_branch.

### CD-38 — no output checks after MAFFT / ClipKIT / IQ-TREE
- **Location:** crossdomain.smk:103-109, 299-305, 327-333.
- **Evidence:** ClipKIT 2.14.0 exits 0 and writes nothing for a missing input ("Input file does
  not exist") and for an unreadable one ("Format type could not be read"), verified.
- **Fix:** `[ -s out ] || exit 1` after each tool (tree and UFBoot file for the cross-domain tree).
- **Test:** e2e_crashed_tools_fail_their_rules.

### CD-39 — missing per-genome systems table counted as zero systems
- **Location:** crossdomain.smk:235-252 (`build_defense_matrix.py` without `--require-all`).
- **Issue/fix:** the extractor now writes a (possibly header-only) table for every genome of
  the universe, so `--require-all` is the defensible setting: a missing table is an error, not
  "no defense systems".
- **Test:** e2e_bacterial_branch.

### CD-40 — codon-alignment interface
- **Location:** crossdomain.smk:380-402 (`--proteins {input.proteins}` = bacterial validated
  proteins only; `cat {input.arch_cds} {input.bact_cds}` copies every archaeal CDS).
- **Issue:** the pristine codon_align.py parsed `--proteins` but never used it (codon_align.py:57),
  so no data were lost then. Under the Contract-v2 codon_align, `--proteins` is the set of
  expected loci whose loss is fatal, so archaeal CDS failures would only have warned. The full
  archaeal CDS set (all genes of all genomes) was also copied on every run.
- **Fix:** combined archaeal + bacterial validated protein list; CDS restricted to those ids;
  per-record `--report`.
- **Test:** e2e_codon_alignment_keeps_both_domains (every expected locus `kept`, references
  `no_cds`, ARC3 not copied); rules_no_swallowed_failures_or_self_copy.

### CD-42 — cross-domain selection without the HyPhy seed
- **Location:** crossdomain.smk:403-441 (`crossdomain_selection`).
- **Issue:** `run_selection.py --seed` sets HyPhy's RNG seed (SetParameter(RANDOM_SEED), default
  12345; run_selection.py:465). The archaeal `evolution.smk:selection` passes `--seed
  {cfg seed}`; the bacterial call did not, so a non-default `seed` reached archaea only.
- **Fix:** `--seed {params.seed}` with `seed=cfg("seed", 12345)`, as the archaeal rule.
- **Test:** rules_parity_with_archaeal_calls (it reported exactly this gap when the evolution leaf
  added the option).

### CD-45 — parity canary did not cover the viral tool calls
- **Location:** tests/crossdomain/run_tests.py (`PARITY`).
- **Issue:** the parity test compared the shared scripts but not the geNomad/CheckV invocations,
  so the drift behind CD-43/CD-44 went unnoticed until the ecology leaf reported it.
- **Fix:** geNomad and CheckV entries added, plus an equality check of the CheckV header against
  ecology.smk.
- **Test:** rules_parity_with_archaeal_calls.

---

## Verified correct (no change needed)
- DefenseFinder 2.0.0 systems columns `sys_id, type, subtype, activity, sys_beg, sys_end,
  protein_in_syst, genes_count, name_of_profiles_in_sys` (df_systems.py:28); type/subtype from
  the model FQN (best_solution.py:49-57). The extraction keys on these.
- geNomad 1.8.0 output prefix = input file stem (summary.py:125-127; one extra split for
  compressed input), so `bact_genomad`'s `{genome}_summary/{genome}_virus*.{fna,tsv}` outputs
  are correct for `{genome}.{ext}` and `{genome}.{ext}.gz`.
- CheckV 1.0.3 CLI `checkv end_to_end <in> <out> -d <db> -t <n> [--restart] [--remove_tmp] [--quiet]`
  (end_to_end.py); geNomad 1.8.0 `end-to-end` options `--cleanup --restart --threads --splits`
  (cli.py:350-418).
- HMMER 3.4 `hmmscan --cut_ga --domtblout` output parsed by `utils.parse_domtbl` (real run in e2e).
- IQ-TREE 2.3.6 node labels are `SH-aLRT/UFBoot`; `--wbtl` writes all 1000 UFBoot trees.
- dendropy 4.6.1 reads GTDB quoted internal labels as node labels; `retain_taxa` collapses
  unifurcations; a force-rooted tree keeps its root (ape `read.tree` → `is.rooted` TRUE).
- HyPhy FEL accepts the cross-domain codon alignment with the pruned tree (e2e).
- Snakemake 8.30 dry run of the real `crossdomain` target on `.test` (Phase 2 simulated):
  86 jobs, exit 0.
- Output keys used outside this file are unchanged (`bact_validate.output.proteins/catalog` in
  disorder.smk; all paths requested by `crossdomain_targets` in common.smk).

## Contract v2 compliance
1. Gene ids: bacterial proteomes keep their DefenseFinder hit ids; combined ids are
   `{genome}__{hit_id}` via `utils.make_protein_id` (as validate/codon_align expect).
2. genes.tsv: bacterial tables from Prodigal headers when present (CD-27); otherwise the
   documented length proxy, reported per genome.
3. Domain spec: consumed only through the utils helpers; no `nuclease` (CD-33).
4. No silent failure: CD-02, CD-26, CD-38; every placeholder is for genuinely empty input and
   marked `# nonfatal:` (checked by tests/check_rules.py).
5. Additive interfaces only (see below); existing paths, rule names, output keys and column
   orders are kept.
6. Tier presence rule (rev 2): implemented in build_bacterial_table
   (`prometheus_status_primary/_sensitivity`), plot_crossdomain_extra and plot_viral_crossdomain;
   the configured tier sets are passed to compute_prevalence (primary + sensitivity),
   cooccurrence.R and ecology_stats.R (primary + A+B re-fit). Regression test with a B-only
   genome (excluded at tier A, present at A+B): plot_comparison_presence_rule_regression;
   status columns: table_known_values. crossdomain_analysis has no tier-based response (it
   uses validated loci only); confirmed.

## Interface changes (additive)
- New outputs: `bacteria/bacterial_samples.tsv`, `bacteria/bacterial_representatives.txt`,
  `bacteria/prometheus/extraction_report.tsv`, `bacteria/prometheus/domains.hmmscan.txt`,
  `bacteria/genes/{g}/{g}.genes.tsv` (side outputs of `bact_extract`, Prodigal-header genomes
  with hits), `bacteria/scaffolds/bacterial_species_tree_tips.tsv`,
  `bacteria/ecology/bacterial_niche_master.tsv` (now a declared output),
  `bacteria/ecology/{bacterial_niche,bacterial_viral}_associations_sensitivity.tsv`,
  `bacteria/ecology/viral/{g}/rloop_rlfs_calls.tsv`, `crossdomain/align/crossdomain.in.faa`,
  `crossdomain/tree/crossdomain.ufboot`, `crossdomain/tree/crossdomain_rooted.treefile`,
  `crossdomain/crossdomain_archaeal_clades.tsv`, `crossdomain/selection/{crossdomain_saturation.tsv,
  crossdomain_proteins.faa, crossdomain_codon_report.tsv, crossdomain_selection.nwk}`,
  `crossdomain/figures/{crossdomain_comparison.tsv, viral_rloop_crossdomain.tsv}`.
- New columns: bacterial master `detect_defensefinder`, `prometheus_status_primary`,
  `prometheus_status_sensitivity`; extraction report `genes_table`; tip table `archaeal_clade`;
  crossdomain_summary metrics `status, rooting, outgroup_monophyletic, outgroup_split_support,
  n_crossdomain_transitions_supported, n_archaeal_clades_supported,
  n_archaeal_clades_nested_supported, archaeal_monophyletic, ingroup_root_state,
  collapsed_root_state, n_bootstrap_trees, bootstrap_root_B/A_fraction, n_subsample,
  subsample_root_B/A_fraction, midpoint_root_state, root_state_excluding_references,
  inferred_direction_ml, direction_method, dnds_saturation_warning`.
- Semantics: `inferred_direction` is a robust call or `ambiguous` (ML-only call kept as
  `inferred_direction_ml`); bacterial `pass_primary` includes contamination ≤ 5; the niche
  master holds annotated genomes only, NA for missing labels.
- New optional flags (owned scripts): extract `--genomes --allow-missing --out-report
  --out-genes-dir`; build_bacterial_table `--samples --max-contamination --system
  --primary-tiers --sensitivity-tiers`; prune_tree_to_taxa `--keep-col --keep-fasta
  --out-report --match --taxonomy --ref-taxonomy --quality --min-match-fraction --no-relabel
  --keep-internal-labels`; crossdomain_analysis `--out-clades --out-rooted-tree --ufboot
  --max-boot-trees --min-ufboot --min-alrt --n-subsample --subsample-min --robust-fraction
  --ref-domain --seed --codon --out-saturation --n-pairs --saturation-ds`; plots `--domains
  --primary-tiers --formats --out-table`. build_bacterial_niche_master `--niche` is now required.
- Optional config keys (rule defaults shown): `bacteria.max_contamination` 5.0,
  `bacteria.allow_missing_proteins` false, `bacteria.min_tree_match_fraction` 0.9,
  `bacteria.gtdb_bac_taxonomy` "", `bacteria.covariates` (default `ecology.covariates`),
  `ecology.genomad.splits` (0, shared with archaea),
  `bacteria.crossdomain.{min_ufboot 95, min_alrt 80, n_subsample 200, subsample_min 5,
  robust_fraction 0.95, ref_domain Bacteria, saturation_pairs 200, saturation_ds 3.0,
  max_codon_drop_frac (default evolution.max_cds_drop_frac)}`.

## Requests for the driver (files outside this leaf)
1. **Critical, evolution.smk `family_tree`:** `cp {params.workdir}/family.treefile {output.tree}` and the
   `.ufboot` copy still copy IQ-TREE's own outputs onto themselves (output = workdir/family.*).
   With the `|| echo "();"` fallback removed, `cp` exits 1 under `set -euo pipefail`
   (verified: "are identical (not copied)"), so the rule fails after every IQ-TREE run. The
   pristine version (evolution.smk:60-61) instead silently replaced tree and UFBoot trees with
   `();`. Fix: delete both `cp` lines.
2. Add the optional config keys above to config.yaml/README (defaults live in the rules).
3. SCHEMA/README: document the new outputs, columns and semantics listed above.
4. common.smk `crossdomain_targets`: request the cross-domain selection outputs only when
   `bacteria.cds_dir` is set (currently handled by an explicit "not run" summary, CD-04).
5. validate_prometheus: keep `--genes-dir` per-genome optional (bacteria have tables only for
   Prodigal-header proteomes).
6. Environment: MAFFT needs an absolute, writable TMPDIR (relative TMPDIR silently changes the
   alignment); the toolchain HOME is not writable for leaves.
7. Snakemake 8.30 runs `conda info`/`conda --version` for rules with `conda:` even when
   environments are not used; the e2e harness strips those lines in its private copy.
8. `rules_parity_with_archaeal_calls` is a deliberate canary: when an archaeal call to a shared
   script gains an option, it fails until crossdomain.smk passes it too (or a reasoned
   exemption is added to `PARITY` in run_tests.py).
9. Manuscript/registry: cross-domain dN/dS is exploratory; report `dnds_saturation_warning`.
10. run_selection.py: prune the gene tree whenever the codon set differs from the tree tips, not
    only with `intact_only` (this branch prunes for itself, CD-20).
11. Genomes with non-11 genetic codes (e.g. Mycoplasmatota, table 4) will be dropped by
    codon_align as translation mismatches; per-genome codon tables would be needed to keep them.

## Unresolved / limitations
- geNomad and CheckV are not installed here; the e2e uses fakes that follow the verified output
  formats.
- Bacterial proteomes without Prodigal headers still use the length proxy for fragmentation
  (reported per genome as `genes_table=none`).
- Cross-domain dN/dS remains exploratory even with the saturation diagnostic; the direction
  inference is topology-based (no DTL reconciliation), as documented in crossdomain_analysis.
- Scratch from kept test runs (`/private/tmp/cdtest_*`, one `/tmp/ck*` ClipKIT probe) could not be
  deleted: removal needed user approval. Normal runs delete their scratch.

## Tests (tests/crossdomain/run_tests.py, 36 tests, ~100 s)
- Unit tests for every owned script, the static rule checks, the archaeal-parity check, and an
  end-to-end Snakemake run of the branch: a private copy of `workflow/` plus a test Snakefile
  with stub archaeal outputs, running real MAFFT, ClipKIT, IQ-TREE, HMMER, HyPhy and R.
  Scenarios: full, all-empty, crashed tools (failing stubs on PATH), and no bacterial CDS (a copy
  of the e2e workdir).
- Deviations: the `conda:` directives are stripped in the private copy (see request 7);
  geNomad/CheckV are fakes. The fixture proteins are 650-aa DUF4011-anchored windows of real ProA
  homologues, with DNA2 C-termini as outgroup, to keep IQ-TREE fast; the missing REase in these
  windows is why the e2e loci reach tier B, not A. The validation tree uses LG+G4 in the test
  config (rule default MFP).
- Scratch lives in `/tmp/cdtest_*` (absolute TMPDIR, HOME, MPLCONFIGDIR set) and is deleted
  unless `CROSSDOMAIN_KEEP=1`.
