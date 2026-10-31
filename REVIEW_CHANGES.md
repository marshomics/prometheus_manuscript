# Review changes

A full robustness, accuracy and scientific-defensibility review of this workflow. The
baseline is commit `56f543735352193e66994eea740fbe722e2a208c` of
`github.com/marshomics/prometheus_manuscript`. Every rule and script of every stage was
checked against primary sources: tool source code and documentation, the pinned
packages, and the Prometheus literature. Each defect was fixed in place and covered by a
test. This file is the index. The full record of each finding (issue, primary-source
evidence, fix, covering test) is in the stage's `tests/<stage>/REVIEW.md`. The finding
IDs below are those files' own IDs.

## Summary

| Stage (`tests/<stage>/REVIEW.md`) | Scope | Critical | Major | Minor | Total |
|---|---|---|---|---|---|
| [preprocess_qc](tests/preprocess_qc/REVIEW.md) | Phase 1 preprocessing + gene calls | 6 | 13 | 17 | 36 |
| [detection_diversity](tests/detection_diversity/REVIEW.md) | Detection, validation + diversity | 3 | 20 | 10 | 33 |
| [scaffolds](tests/scaffolds/REVIEW.md) | Taxonomy, species tree, methanogens, master table | 6 | 21 | 16 | 43 |
| [distribution](tests/distribution/REVIEW.md) | Prevalence, association, co-occurrence | 3 | 17 | 15 | 35 |
| [evolution](tests/evolution/REVIEW.md) | Gene trees, reconciliation, ancestral states, selection | 6 | 19 | 10 | 35 |
| [ecology](tests/ecology/REVIEW.md) | Viral pressure, R-loops, niche tests | 5 | 10 | 7 | 22 |
| [synthesis](tests/synthesis/REVIEW.md) | Registry, master tree, summary figures | 1 | 11 | 10 | 22 |
| [crossdomain](tests/crossdomain/REVIEW.md) | Bacterial / cross-domain branch (opt-in) | 3 | 29 | 13 | 45 |
| [disorder](tests/disorder/REVIEW.md) | Intrinsic disorder (opt-in) | 4 | 16 | 10 | 30 |
| cross-stage (`DRV-*`, below) | integration | 1 | 7 | 12 | 20 |
| **all** | | **38** | **163** | **120** | **321** |

The 301 stage findings come from nine stage reviews. The 20 cross-stage findings
(`DRV-*`) come from integrating them, mainly through the end-to-end integration test
(`tests/integration/`), which runs the real Snakemake DAG through both phases.

### The findings that changed results

- **Every gene tree was empty.** `family_tree` and `crossdomain_tree` copied IQ-TREE's
  outputs onto themselves. `cp` exits 1, so the `|| echo "();"` fallback replaced the ML
  gene tree and the UFBoot trees with `();` on every run. Reconciliation, selection, the
  tanglegram, the per-domain tree comparison and the whole cross-domain analysis then ran
  on an empty tree and reported "not estimated" or skipped (EVO-34, CD-01).
- **The codon alignment was empty whenever any locus started with GTG.** PAL2NAL ran with
  the standard code (table 1), which rejects GTG starts. It then aborts with its error on
  stderr, nothing on stdout and exit 0, so dN/dS (including the confirmatory RELAX test)
  was silently skipped. Fixed with code 11, plus id and order verification (EVO-01, DRV-17).
- **GTDB-Tk 2.4.0 aborts on genome ids equal to GTDB reference ids.** That includes bare
  NCBI accessions, the normal id form. Phase 1 (PQ-03) and the Phase-2 marker alignment
  (DRV-01) now pass prefixed ids and strip the prefix afterwards.
- **The headline association never ran** (D1). The confounder-adjusted phyloglm, the
  confirmatory methanogen test, failed on every input, and the failure was silent (D2).
- **Absent and excluded were conflated.** Under the primary tier set {A}, genomes whose
  only locus was tier B/C were counted as confident absences in prevalence, association,
  co-occurrence, ancestral states, ecology and disorder. The rule is now `present` /
  `absent` / `excluded` everywhere (Contract v2 item 6; BM-05, P1, D4, C2, EVO-11, ECO-12,
  D-24).
- **Tier A silently required the optional DUF3320**, which is absent from ~54% of
  Prometheus, so roughly half of all genuine loci could not reach tier A (DD-17).
- **Unknown niches were coded as absent** before the species union, so a species with no
  environment record became a negative in every niche test, and the missingness report
  always said zero (DRV-04).
- **Systemic silent failures.** `|| true`, `tryCatch(… quit(status = 0))`, placeholder
  tables and E-value fallbacks turned crashed tools, missing databases and unreadable
  inputs into "no Prometheus", "no methanogens" or "not estimated" across every stage
  (e.g. PQ-01, PQ-04, PQ-05, DD-01, DD-12, SC-01, SC-06, RT-01, MS-03, CM-01, BM-01,
  EVO-17, ECO-01–03, ECO-10, S1, CD-02, D-05, D-22). Now every rule fails loudly, and
  `tests/check_rules.py` enforces it.
- **Identical inputs gave different results.** Every MAFFT alignment ran its iterative refinement multithreaded,
  which MAFFT documents as non-reproducible, so alignments, trees, placement tiers and everything downstream
  could change between two runs of the same data. Fixed with `--threadit 0` (DRV-19).
- **Environments that could not be built:** the gene-calling env (PQ-06), and the disorder
  module's non-existent and unsatisfiable pins (D-01, D-02).

## Cross-stage (integration) findings

| ID | Severity | Location | Issue | Fix | Covering test |
|---|---|---|---|---|---|
| DRV-01 | critical | `rules/scaffolds.smk` `gtdbtk_align` | GTDB-Tk 2.4.0 (`tools.get_reference_ids`, `main._genomes_to_process`) raises `GTDBTkExit` when a user genome id equals a GTDB reference id, including bare `GCA_`/`GCF_` accessions. The Phase-2 marker alignment therefore failed for accession-named genomes (and in the baseline the failure was swallowed into an empty MSA). | Batchfile ids prefixed `user_`; the prefix is stripped from the MSA and from the dropped-genome table. | `tests/scaffolds` (a GTDB-Tk stub that refuses reference ids), `tests/integration` (NCBI-accession ids) |
| DRV-02 | major | `plotting.formats`; all figure rules and writers | The key was only half wired. `all` requested one file per configured format, but the Python and R writers always wrote png + svg (the value never reached them), the R writer had no pdf support, and some rules (e.g. `ancestral_states`) declared fixed png/svg outputs. Any non-default value therefore broke the run. | Validated `PLOT_FORMATS` in `common.smk`, exported to the scripts (`PROM_PLOT_FORMATS`, DPI, fonts); every figure rule declares one output per format; Python `default_formats()`, R `plot_formats()`. | `tests/check_dag.sh` (invalid format rejected), `tests/integration` (`[png, pdf]`), `tests/synthesis` |
| DRV-03 | major | `scripts/R/plotting_theme.R` `save_plot`; `scripts/py/plotting_style.py` `save_fig` | A figure that failed to render produced a warning and a missing or partial file, and the rule succeeded. | R: `stop()` on any render error and remove the partial file. Python: raise if a file is missing or empty. | `tests/synthesis` `test_save_plot_render_error_fails` |
| DRV-04 | major | `scripts/py/filter_collapse.py` (niche union) | Blank `group_*` / `niche_*` cells were coded 0 before the species union. This contradicts SCHEMA Contract v2 (species-level {1, 0, NA}) and the master table's NA handling. A species with no environment record became a negative in every niche test, and `ecology_missingness.tsv` always reported n_NA = 0. Observed in the integration run: 6 of 13 species. | The union keeps NA: 1 if any genome is 1, else 0 if any is 0, else blank; the count is logged. | `tests/preprocess_qc` `test_filter_collapse_niche_union_modes`, `tests/integration` (species table, master table, missingness counts, ecology n) |
| DRV-05 | minor | `scripts/R/plot_tree_master.R` (niche richness) | Once unknowns reach the master table, `rowSums(nm != 0)` made a partly unknown species "no data". | Count observed niches; NA only when every niche is unknown. | `tests/synthesis` master-tree test (G5, G7) |
| DRV-06 | major | `workflow/schemas/config.schema.yaml` | The schema left most keys untyped and unbounded, so invalid values (an unknown FDR method, an unsupported figure format, UFBoot < 1000, which IQ-TREE rejects, malformed tier lists, unknown enums) passed DAG build and failed late or silently. | Types, enums and bounds for every key; validated on the main and `.test` configs. | 7 negative cases rejected; `tests/check_dag.sh` |
| DRV-07 | major | `rules/common.smk` `all_targets` | `all` never requested `evolution/domains/rf_distances.tsv`, so the per-domain trees (design 5.2) never ran, even with `evolution.per_domain_trees: true` (leaf finding EVO-32, fixed in the shared file). | Target added when the key is true. | `tests/check_dag.sh`, `tests/integration` |
| DRV-08 | minor | `rules/common.smk` | A genome id that breaks Newick or the `genome__protein` separator, or a data row truncated by an inline `#`, failed hours later inside a tool (leaf findings in `batch_genomes.py`, applied to the shared DAG-build code). | `_ID_RE` (`[A-Za-z0-9][A-Za-z0-9_.-]*`, no `__`) and `_check_no_inline_hash` at DAG build. | `tests/preprocess_qc` (drift guard executes `common.smk`'s own helpers) |
| DRV-09 | minor | `rules/evolution.smk` `ancestral_states` | The ancestral figure was declared in a single fixed format. | One output per `plotting.formats` entry. | `tests/check_dag.sh`, `tests/integration` |
| DRV-10 | minor | `rules/plots_extra.smk` `plot_ecology_extra` | The rule did not pass `--primary-tiers` / `--formats`, so the figure used script defaults rather than the configured tier set. | Both passed from config. | `tests/check_dag.sh`, `tests/integration` |
| DRV-11 | minor | `scripts/py/results_registry.py` | The reconciliation row reported `vertical_fraction` without its definition. | The note now gives `vertical_fraction=(S+SL)/(S+SL+T+TL)`. | `tests/synthesis` |
| DRV-12 | minor | `config/config.yaml` | Dead keys (`db.interproscan`, `db.kofam_profiles`); `db.checkm2` / `db.checkv` defaults that did not match the directories the documented download commands create; many keys used by rules but undocumented. | Dead keys removed; defaults aligned (`…/CheckM2_database/uniref100.KO.1.dmnd`, `…/checkv-db-v1.5`); every key documented. | config validated against the schema |
| DRV-13 | major | `README.md` | False or stale claims: a "detection-sensitivity" covariate that does not exist; tier A requiring a four-domain architecture; midpoint rooting as a silent fallback; "every figure as PNG and SVG"; a dangling `AUDIT.md` link; `defense-finder update`, which installs an unpinned CasFinder that macsyfinder 2.1.4 rejects; geiger; wrong counts. | Rewritten against the code; every referenced config key exists. | key-existence check (62 keys) |
| DRV-14 | minor | `resources/README.md` | Genomes under `config/samples.tsv`, a `defensefinder_merged` mode, the `interproscan` / `kofam_profiles` keys and name-based domain matching, none of which exist. | Rewritten: pinned `macsydata install`, CheckM2 database ≤ v2, GTDB R220, KOfam `ko_list` required, exact domain matching, ingest-mode id requirement. | — |
| DRV-15 | minor | `Prometheus_archaea_pipeline_design.md` | The Stage-1 implementation note said tier A needs the four-domain architecture, and the design described many analyses that are implemented differently or not at all. | Note corrected; §19 "As implemented" lists every departure. | — |
| DRV-16 | minor | `workflow/scripts/SCHEMA.md` | Several tables were undocumented or stale. | Rewritten as Contract v2 with every table, column and missing-value convention. | per-stage contract tests, `tests/integration` |
| DRV-17 | minor | `scripts/py/codon_align.py` (docstring) | The PAL2NAL failure mechanism was described inaccurately. | Corrected: PAL2NAL 14.1 prints the inconsistency error to stderr (`pal2nal.pl:564`), writes nothing to stdout and exits 0. | `tests/evolution` |
| DRV-18 | minor | `tests/` | Test scratch used relative paths and was never cleaned. The `common.smk` drift guard broke when the shared file gained an import. There was no test of the assembled workflow. | Absolute scratch with cleanup (`PROM_KEEP_TMP`, `PROM_TEST_TMP`); the drift guard executes `common.smk`'s imports; new `tests/integration/` harness. | all suites |
| DRV-19 | major | every MAFFT call: `prometheus_validation_tree`, `bact_validation_tree`, `crossdomain_align`, `family_align`, `domain_trees.py`, `classify_methanogen.py` | `mafft --thread N` with iterative refinement (L-INS-i, or `--auto` on small inputs) is not reproducible run to run. MAFFT 7.526 gave 2 different alignments in 3 runs on identical input; in the integration world this flipped the validation-tree placement of all 7 candidates between Prometheus and Unplaced across reruns, because the clade support sits at the UFBoot threshold. The MAFFT multithreading documentation states that the iterative refinement can return different results for different runs. | `--threadit 0` on every call: refinement single-threaded, other stages still threaded. Three runs now give identical output. | replicate test in the review; `tests/integration` |
| DRV-20 | minor | `workflow/envs/reconciliation.yaml` | The env pinned `ete3=3.1.3`, which no script imports (checked in the baseline and the current tree). The evolution review dropped it but did not record the change. | Dependency removed; the env (python 3.11, generax 2.1.3, dendropy 4.6.1, pandas 2.2.2) solves for linux-64. | `tests/check_env_solve.sh` |

## Methanogenesis pathway capacity (added after the review, design 3.3 / 6.2)

The review left design 3.3 partly unimplemented: `methanogenesis_subtype` was `unresolved`
for every methanogen, so the master-tree subtype ring was uninformative and design 6.2 ("does
the methanogenesis subtype predict presence?") could not be asked. That gap is now closed as an
opt-in module (`methanogenesis.enabled`, default false; the default DAG is unchanged at 103
Phase-2 jobs).

**What it does.** 85 pathway KOs are searched in the methanogens with the same KofamScan-exact
acceptance as `mcr_search` (`kofam.py`, now shared), structured by the verbatim KEGG module
definitions (M00567, M00357, M00356, M00563, M00880) into five multi-label pathway calls plus
substrate flags. Calls are three-valued (present / absent / unresolved) and `not_applicable`
for non-methanogens. Rules, evidence and primary citations: `workflow/resources/methanogenesis/
rules/RULES.md`; statistics: `methanogenesis_stats.R`.

**Defects found and fixed while building it** (each is now covered by a test):

| ID | Severity | Issue | Fix | Covering test |
|---|---|---|---|---|
| MG-01 | critical | A single KO was a mandatory `signature` of its complex. KEGG assigns neither FwdB (K00201) nor HdrA (K03388) to *Methanocaldococcus jannaschii* -- both UGA/selenoprotein-affected -- so a textbook CO2-reducing hydrogenotroph was called `none_detected`. A single-subunit signature reintroduces the single point of failure the relaxation exists to remove. | `signature` entries are now SETS of interchangeable core subunits (at least one must be detected). | `test_component_signature_tolerates_a_missing_selenoprotein_subunit`, `test_reference_methanogens_known_answers` |
| MG-02 | critical | Two DOIs in `RULES.md` were written from memory and never resolved: `10.1126/science.aaf0812` does not exist (Mayumi 2016 is `10.1126/science.aaf8821`) and `10.1111/j.1432-1033.1997.00670.x` resolves to a malaria glutathione paper (Sauer 1997 is `...t01-1-00670.x`). Aramaki is 2019, not 2020; the selenium paper is Quitzke et al., not Rother & Quitzke. | All 15 DOIs re-resolved against CrossRef and corrected; `citations.tsv` records each resolved first author, year, title and journal, and the leaf gate fails if a cited DOI is missing from it, unverified, or listed but uncited (negative control: a planted fake DOI fails the check). | `tests/methanogenesis/check_owned.py` |
| MG-03 | major | `Complex.__repr__` did not parenthesise a nested alternative, so `K00200+K00201+(K00205,K11260)` re-printed as `...+K00205,K11260`, which re-parses as a top-level alternative with different semantics. | `_wrap` applied to complex subunits; the round-trip is now checked semantically on 400 random KO subsets per module rather than by string equality (KEGG's redundant parentheses are dropped by design). | `test_kegg_definitions_reparse_identically_and_have_kegg_step_counts` |
| MG-04 | major | KOfam thresholds are optimised across all taxa and can sit just above a divergent archaeal orthologue: K00625 (Pta) scores 379.0 in *Methanosarcina acetivorans* against a threshold of 387.57, so both *M. acetivorans* and *M. barkeri* -- canonical acetoclastic methanogens -- were called `acetoclastic=absent` on real genomes. | `near_threshold` evidence (>= 0.9 x threshold): never presence, but it withdraws a confident absence (the component becomes `equivocal`). `ackpta` additionally relaxed on AckA, documented in RULES.md section 5. | `test_near_threshold_hits_block_a_false_absence_without_creating_presence`, `test_real_genome_ko_calls_reproduce_published_physiology` |
| MG-05 | major | `results_registry` took `methanogenesis_stats` as an unconditional input, so with the module off (the default) the run would fail: the master table has no `mgen_*` columns. | One `MGEN_ENABLED` flag in `common.smk`; the registry input, the `all` targets and the Stage-3 rules all gate on it. DAG back to the 103-job baseline. | `tests/check_dag.sh`, `tests/check_methanogenesis_contract.py` |
| MG-06 | minor | The master-tree subtype ring recycled three colours for unrecognised labels, so two different subtypes could share a colour. | Deterministic palette; composites beyond the reserved colours are pooled into one honestly named `other combination` level, with the pooling logged. | `tests/run_leaf.sh synthesis` |

**Judgement calls.** Absence is gated on completeness (`absence_alpha`, default 0.05), so a
60 %-complete MAG missing one locus is `unresolved`, not a confident zero. An Acs-only
methylotroph is `acetoclastic=unresolved` rather than present or absent, because Acs is also
anabolic. Methylamine calls require the pyrrolysine machinery, so non-Pyl MttB-superfamily hits
set `sub_quatamine` instead. Methoxydotrophy has no KEGG module and is one-sided (presence is
evidence, absence is not).

**Validation.** On 12 cultured methanogens taken from assembly to call (RefSeq/GenBank ->
Prodigal 2.6.3 -> the caller with the curated KOfam subset), 34 of 36 asserted pathway calls are
resolved and correct against published physiology with **no mismatch**; the 2 unresolved are the
Acs-only methylotrophs above. Read-through recovered 11 genes Prodigal had truncated (3 UAG /
pyrrolysine, 8 UGA / selenocysteine), including FwdB and HdrA in *M. jannaschii*. The two known
capacity-vs-phenotype discordances are reported as gene-level positives, which is correct for a
capacity call: *M. smithii* encodes `mtaABC` but does not grow on methanol
(doi:10.1073/pnas.0704189104), and *M. acetivorans* encodes the H2/CO2 route but does not use it
(doi:10.1128/JB.00563-08).

**Open handoff.** The independent reference-physiology curation (an externally cited truth table
and held-out validation set, leaf-1.1.2) was stopped by a platform content-safety filter and is
NOT complete; the benchmark leaf that depends on it (1.1.4) is therefore not run. The validation
above uses KEGG's own annotation and published physiology recorded by the driver, which is
weaker than an independent curated truth set. See `.unlazy/methanogenesis-type/status.log`.

## Documented limitations and judgement calls (not defects)

- **DefenseFinder sensitivity.** Detection is DefenseFinder only. With `--cut_ga`, the ProA
  profile's per-domain GA (600 bits) drops a genuine ProA whose alignment splits into two
  sub-threshold domain hits. UniProt A0A0A6X990, for example, scores 1024.8 bits over the
  full sequence but splits into domains of 441.8 + 574.1. Such loci are invisible
  downstream, so prevalence means "detected by DefenseFinder". The integration test
  asserts this behaviour.
- **ANME coding differs by stage.** The association and ancestral tests code ANME as
  non-methanogens and report the count; the disorder comparison excludes them. Both are
  defensible. **Recommendation:** add an ANME-excluded sensitivity fit of the headline
  association.
- **Genome-size scale.** The association model uses z(log genome size); the ecology models
  use z(genome size). Both are documented in SCHEMA.md.
- **Unregistered tests.** The per-domain AU-test p-values (`rf_distances.tsv`) and the opt-in
  bacterial / cross-domain tests are not entered in `results_registry.tsv`. Treat them as
  exploratory, or register them if they are reported as tests.
- **IQ-TREE threads.** IQ-TREE runs are seeded, but multithreaded IQ-TREE is not guaranteed bit-reproducible
  (iqtree2 issue #228; Shen et al. 2020, *Nat. Commun.*, doi:10.1038/s41467-020-20005-6). In the review's 2-thread replicate test, the
  topology and every support value were identical and only branch lengths differed in the last digits. For exact
  reproducibility, run with 1 thread. Placement supports are reported per locus (`placement_support`), so
  near-threshold calls are visible.
- **Methanogenesis subtype** is not assigned (`unresolved` for every methanogen), so the
  master-tree subtype ring is uninformative.
- **Design items not implemented** (orthogonal detection sweep, context-type clustering,
  compositional HGT signals, Count/BadiRate, methanogen niche resolution, reference-anchored
  tree, and others) are listed in the design note's §19.

## Verification

- 251 synthetic tests in 10 suites (preprocess_qc 27, detection_diversity 26, scaffolds 28,
  distribution 20, evolution 34, ecology 21, synthesis 21, crossdomain 36, disorder 28,
  integration 10).
  - **All 10 suites pass on the final code:** integration 10/10, and every stage suite
    re-run through its gate ledger after the last change (preprocess_qc 27, detection_diversity
    26, scaffolds 28, distribution 20, evolution 34, ecology 21, synthesis 21, crossdomain 36,
    disorder 28), together with each stage's rule-hygiene and env-solve gates.
  - **Evolution flake:** in 10 evolution runs, one run, under concurrent load, reported
    33/34. The failing test could not be identified (the gate checker elides suite output),
    and 9 further runs, several under load, all passed. Treat it as an unreproduced
    intermittent failure.
  - **How integration reached 10/10:** after the niche-union fix (DRV-04) the integration
    test first passed 8/10. The two failures were (a) validation-tree placements that
    flipped between reruns and (b) a new assertion that ignored the tier rule.
    - (a) traced to MAFFT's non-reproducible multithreaded refinement (DRV-19) combined
      with a weak reference set in the test world. With `--threadit 0` and canonical
      references, every candidate places at 97–100 UFBoot across seeds.
    - (b) was corrected in the test. Real tools are run wherever they are installed. The R-loop
  re-implementation reproduced 477/477 reference QmRLFS-finder v2.0 calls on an independent
  160-kb set.
- The integration test runs the real Snakemake DAG (Phase 1, then Phase 2) on a 17-genome
  synthetic world built from real UniProt ProA homologues and paralogs, and checks every
  stage contract against the known truth.
- `tests/check_dag.sh` dry-runs `preprocess`, `all` (before and after Phase 1),
  `crossdomain` and `disorder`. `tests/check_rules.py` passes for all 12 rule modules.
  `tests/check_env_solve.sh` solves all 13 envs for linux-64.

**Not run here** (tests use faithful stand-ins or mocks):
- GTDB-Tk, CheckM2 and BBMap, because their reference data run from tens to hundreds of GB.
- geNomad and CheckV.
- GeneRax 2.1.3, which crashes on the macOS test host.
- ALE.
- IUPred3, which is licence-gated.

The tests ran on osx-arm64 with locally installed tools (e.g. phytools 2.1.1, not the
pinned 2.3-0) and with the `conda:` directives stripped. HyPhy was run on small
alignments only. The pinned linux-64 environments were solved but not executed. A first
production run should use a subset of genomes.

## Process notes

The review ran as nine parallel stage reviews, each owning its files under a written
gate ledger (`.unlazy/prometheus-review/`). Each stage was then independently re-verified
before integration, and its leases were released. The cross-domain review once re-applied
a held patch before parent approval; the change was reviewed and kept at the user's
request.

## All stage findings, by severity

One line per finding: the stage file's own ID, location and the first sentence of its issue text, copied verbatim. Follow the stage link for the primary-source evidence, the fix and the covering test.

### Critical — silently wrong or empty results (37)

| ID | Stage | Location | Issue |
|---|---|---|---|
| PQ-01 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:70-77 (gtdbtk_batch) | a crashed GTDB-Tk was swallowed |
| PQ-02 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:72-76 (gtdbtk_batch) | summary merge read the wrong files |
| PQ-03 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:68-70 + batch_genomes.py | GTDB-Tk aborts on genomes named like GTDB references |
| PQ-04 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:102-108 (checkm2_batch) | a crashed CheckM2 was swallowed |
| PQ-05 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | qc.smk:33-36 (prodigal_genes) | a crashed Prodigal produced an empty proteome |
| PQ-06 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | envs/qc.yaml:7,16 | environment unsolvable |
| DD-01 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | A crashed DefenseFinder run became a "no Prometheus" genome |
| DD-12 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | An unpressed Pfam silently removed every domain |
| DD-17 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Tier A silently required the optional DUF3320 |
| SC-01 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:35,37,40 | `gtdbtk identify ... \|\| true` and `gtdbtk align ... \|\| true`, then an empty gzip placeholder when no MSA was found. |
| SC-06 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:80,83,86,88 | `\|\| true` on both IQ-TREE runs, plus `cp species.treefile ... \|\| echo "();"`. |
| RT-01 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:25-27,104-110,163-166 | every failure (dendropy missing, parse error, rooting exception) wrote the UNROOTED input as the canonical rooted tree. |
| MS-03 | [scaffolds](tests/scaffolds/REVIEW.md) | mcr_search.py:136-138,179-191,216 | An hmmsearch non-zero exit was only a WARNING, so its hits were lost. |
| CM-01 | [scaffolds](tests/scaffolds/REVIEW.md) | classify_methanogen.py:24-27,110-173,295 | every ANME-check failure (missing or empty references, mafft or IQ-TREE failure, dendropy missing, McrA missing from the proteome) became `anme_check_skipped` with is_anme False, so ANME were silently counted as methanogens. |
| BM-01 | [scaffolds](tests/scaffolds/REVIEW.md) | build_master_table.py:66-81,91-96,120-125,145-150,170-175,204-225,240-250,353-365 | silent empty joins: a missing or unreadable side table became `{}` (e.g. a broken catalog meant no Prometheus anywhere); |
| D1 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:121,166,224 | headline phyloglm never ran: the master was read with `readr::read_tsv` (a tibble). |
| D2 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:89-91,221,258,303,327,356-365 | every failure was silent: the whole script sat in `tryCatch(... header-only ...)` followed by `quit(status = 0)`; each test had its own `tryCatch` writing a `failed:` row; |
| C1 | [distribution](tests/distribution/REVIEW.md) | cooccurrence.R:162,197-201,70,237 | silent method substitution: when `fitPagel` failed, or phytools was missing, the Fisher p was used with method `fisher(fallback)` and entered the registry's co-occurrence family next to Pagel p-values. |
| EVO-01 | [evolution](tests/evolution/REVIEW.md) | codon_align.py:121-125 | PAL2NAL ran with the default codon table 1. Table 1 does not accept GTG as an initiator Met, and Prodigal translates every start codon as M, so any GTG-start gene was "inconsistent". |
| EVO-05 | [evolution](tests/evolution/REVIEW.md) | run_reconciliation.py:129,146-160,196 | GeneRax output parsing was wrong. `<family>_eventCounts.txt` holds `S: SL: D: T: TL: L: Leaf:`, but `key_map` ignored SL and TL, so speciations and transfers were undercounted. |
| EVO-06 | [evolution](tests/evolution/REVIEW.md) | run_reconciliation.py:204-260 | GeneRax discards the whole family, with only a log warning and exit 0, when any gene maps to a species absent from the species tree, or when alignment and tree labels differ. |
| EVO-11 | [evolution](tests/evolution/REVIEW.md) | ancestral_states.R:188-191,209 | Presence came from `prometheus_present` (any non-rejected tier). Under the primary set {A}, B/C-only genomes were confident absences, against contract v2 item 6. |
| EVO-17 | [evolution](tests/evolution/REVIEW.md) | run_selection.py:113-136,505-576 | HyPhy failures were swallowed ("skipped/failed" rows, exit 0), including a missing hyphy. |
| EVO-34 | [evolution](tests/evolution/REVIEW.md) | evolution.smk:60-61 | `output.tree`/`output.ufboot` ARE IQ-TREE's own `--prefix {params.workdir}/family` outputs, so `cp {params.workdir}/family.treefile {output.tree}` copies a file onto itself. cp exits 1, and the `\|\| echo "();" > {output.tree}` fallback then REPLACED the ML gene tree and the UFBoot trees with `();` on EVERY successful run. |
| ECO-01 | [ecology](tests/ecology/REVIEW.md) | ecology.smk pristine:27,30 (now `rule genomad` L50-84) | `genomad end-to-end ... \|\| true` swallowed every geNomad failure (missing DB, crash, OOM), then wrote a 5-column placeholder (`seq_name,length,topology,n_genes,virus_score`) that is not geNomad's header, so a crashed genome became "0 proviruses" |
| ECO-02 | [ecology](tests/ecology/REVIEW.md) | ecology.smk pristine:53-56 (now `rule checkv` L86-126) | `checkv end_to_end ... \|\| true` + 4-column placeholder masked crashes and a missing DB; the placeholder header did not match CheckV's; |
| ECO-03 | [ecology](tests/ecology/REVIEW.md) | aggregate_viral.py pristine:16-39,42-62,65-94,106-114 | every read problem returned 0 / None (`except Exception: return 0, 0.0`), an unreadable samples table gave a header-only summary, and samples were read without `dtype=str` (an id like `000123` became `123` -> its files "missing" -> silently 0) |
| ECO-10 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:77,140-144,248-249,270-276 | a top-level `tryCatch` wrote a header-only table and exited 0 on ANY error; a missing `phylolm` gave a header-only table; the whole ecology analysis could disappear silently |
| ECO-11 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:191-198,233-240 (now `phyloglm_once`/`fit_phyloglm` L203-263) | non-converged phyloglm fits were reported as valid estimates (no `convergence` check). |
| S1 | [synthesis](tests/synthesis/REVIEW.md) | plot_tree_master.R:70-268 | The whole script ran inside one `tryCatch` that drew a placeholder and `quit(status = 0)` on ANY error (unreadable input, missing column, package missing, render failure): the master figure could silently be an "error" placeholder while the rule succeeded (Contract v2 item 4). |
| CD-01 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | cross-domain tree replaced by `();` (cp onto itself) |
| CD-02 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | tool failures swallowed (`\|\| true`, `\|\| cp`, Pfam guard) |
| CD-03 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | bacterial validation tree has no supports → tier A unreachable in bacteria |
| D-01 | [disorder](tests/disorder/REVIEW.md) | workflow/envs/disorder.yaml:15 | `metapredict==2.6.5` does not exist on PyPI, so the env cannot be built and the module can never run. |
| D-02 | [disorder](tests/disorder/REVIEW.md) | workflow/envs/disorder.yaml:20 | `alphaPredict==1.2` is unsatisfiable next to every metapredict 2.6x (2.63 needs `==1.0`, 2.65.x needs `>=1.3.0`), so even a corrected 2.6x pin would not resolve. |
| D-05 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:118-131, 280-286, 312-313 | A missing metapredict produced header-only outputs (the rule "succeeded"), and any per-sequence exception dropped that locus silently. metapredict raises on X/B/Z/U/`*` (verified), so every protein carrying Prodigal's trailing `*` or an ambiguous `X` would have vanished without a trace (see U4). |
| D-22 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R:47, 86, 136-139 | Top-level `tryCatch` wrote a header-only table and `quit(status = 0)` on any error; a missing `phylolm` was "handled" the same way. |

### Major — wrong biology or statistics, data loss, or failure on realistic inputs (156)

| ID | Stage | Location | Issue |
|---|---|---|---|
| PQ-07 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | run_bbmap_stats.py:69-77 | BBMap failures were warnings |
| PQ-08 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | merge_annotations.py:106-113 | missing BBMap column only warned |
| PQ-09 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | merge_annotations.py:36-46,77-102 | inputs located by glob, read errors skipped |
| PQ-10 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | merge_annotations.py:83, 97, 115 | duplicates and unknown ids silently resolved |
| PQ-11 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:147 | merge rule default `length_col: scaf_bp` |
| PQ-12 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | batch_genomes.py:33-38,54,59-67 | genome set and ids not the ones the DAG counts |
| PQ-13 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py:184-191 | staging failures ignored, real files deleted |
| PQ-14 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py:132 | n_genomes_in_species depended on the aggregation mode |
| PQ-15 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py:86 | invalid niche values silently converted |
| PQ-16 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py:94-107 | keep_domain fragile, filter could be skipped |
| PQ-17 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | merge_annotations.py:128 / filter_collapse.py:119-125 | threshold decided by float noise |
| PQ-18 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | rename_prodigal.py:30-48 | gene ids not Contract v2, no coordinate table, silent CDS drops |
| PQ-19 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py / merge_annotations.py | genomes GTDB-Tk/CheckM2 did not report were invisible |
| DD-02 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Multi-contig genomes were searched as one replicon |
| DD-03 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Every DefenseFinder run depended on the GitHub API |
| DD-04 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Floating detector pins; documented update path breaks the run |
| DD-05 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | The preflight could pass without any usable model set |
| DD-06 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Ingest mode silently replaced a missing provided table with an empty one |
| DD-10 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Genome universe from the filesystem; silent id and read fallbacks |
| DD-11 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Missing per-genome proteomes were skipped silently |
| DD-13 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Sequence count broke on zero matches |
| DD-14 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | ClipKIT and IQ-TREE failures produced "no tree", and every candidate became Unplaced |
| DD-15 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Placement tree had no branch supports and no seed |
| DD-18 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Substring domain matching, an obsolete category, and no order check |
| DD-19 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | "Catalytically intact" was inferred from domain presence alone |
| DD-20 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Placement by nearest reference on an unrooted, unsupported tree |
| DD-22 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Defense context was genome-level, and the window parameter was ignored |
| DD-23 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | "Fragmented" meant "short", and contig-edge fragments were Rejected |
| DD-27 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Subfamily column names broke the rarefaction figure |
| DD-28 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | MMseqs2 failures produced blank or header-only subfamily tables |
| DD-29 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Clustering relied on implicit MMseqs2 defaults |
| DD-30 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Typology used a hard-coded, outdated domain list |
| DD-31 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Rarefaction over carriers only, and subfamily ids labelled as architectures |
| SC-02 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:39 | `find -name '*user_msa.fasta.gz' \| head -n1` can pick `gtdbtk.bac120.user_msa.fasta.gz`. |
| SC-03 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:19,33 | `tree.scope: domain` was a silent no-op: dropping `--skip_gtdb_refs` does not put references in the user MSA. |
| SC-04 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:34-40 | genomes that GTDB-Tk leaves out of the MSA disappeared from every phylogenetic analysis without a trace. |
| SC-05 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:34 | `--genome_dir` processes every `*.fna` in the staging directory, including stale genomes from an earlier Phase-1 run. |
| SC-07 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:75 | `n=$(grep -c '>' $msa \|\| echo 0)` yields `0\n0` for an empty MSA (grep -c prints 0 and exits 1), which breaks `[ "$n" -ge 4 ]`; it also counted `>` anywhere on a line. |
| RT-02 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:141,147,151,154 | SH-aLRT/UFBoot supports ended up on the wrong bipartitions after rerooting. |
| RT-03 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:124-145 | the outgroup MRCA was computed on IQ-TREE's arbitrary root. A monophyletic outgroup straddling that root was declared unrootable (midpoint fallback). |
| RT-04 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:72 | `outgroup_taxon` matched by substring (`tx in classification.lower()`), so `p__Altiarchaeota` also matched `p__Altiarchaeota_A`, which GTDB treats as a different lineage. |
| RT-05 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:158 | the dendropy writer emitted a `[&R] ` prefix and single-quoted labels containing `/` (`'100/100'`). ape keeps the quote characters in node labels (verified), so the supports were non-numeric in R. |
| RT-06 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:133,147,151,154 | `reroot_at_midpoint` mishandles a midpoint that falls exactly on a node: it turns a tip into the root node: `(A:1,B:1,C:1)` became `(A:1.0,C:1.0)B;`, losing tip B; |
| RT-07 | [scaffolds](tests/scaffolds/REVIEW.md) | (not handled) | IQ-TREE rewrites characters outside `[A-Za-z0-9_-.\|/]` to `_` (25/25 characters predicted on a real run). |
| MS-01 | [scaffolds](tests/scaffolds/REVIEW.md) | mcr_search.py:15-24,214; config `methanogen.kofam_score_fraction: 0.5` | hits were accepted at half the KOfam adaptive threshold. |
| MS-02 | [scaffolds](tests/scaffolds/REVIEW.md) | mcr_search.py:81-107,214 | the ko_list `score_type` was ignored (the full score was always used). |
| MS-04 | [scaffolds](tests/scaffolds/REVIEW.md) | mcr_search.py:45-56 | the subunit came from the file stem only. A TIGRFAM profile named `TIGR03256.hmm` became gene `TIGR03256`, which classify ignores, so no methanogens were called (resources/README offers TIGRFAM as an option). |
| CM-02 | [scaffolds](tests/scaffolds/REVIEW.md) | classify_methanogen.py:220,284-292 | only the exact label `ANME` counted. An `other` nearest reference (e.g. an alkyl-CoM reductase) and subgroup labels such as `ANME-1` / `ANME-2a` both led to methanogen = True. |
| CM-03 | [scaffolds](tests/scaffolds/REVIEW.md) | classify_methanogen.py:176-221 | assignment was by the patristically nearest reference alone, which is sensitive to long branches. |
| BM-02 | [scaffolds](tests/scaffolds/REVIEW.md) | build_master_table.py:212-213 | the focal system was excluded by the substring `prometheus` in `type`. |
| BM-03 | [scaffolds](tests/scaffolds/REVIEW.md) | (not handled) | DefenseFinder 2.0 with `-a` reports AntiDefenseFinder systems in the same table (`activity = Antidefense`), and they were counted as defense systems. |
| BM-04 | [scaffolds](tests/scaffolds/REVIEW.md) | build_master_table.py:178-192 | catalog tiers were not validated. A blank or unknown tier counted as present with tier `''`, which tier-based consumers read as absent. |
| BM-05 | [scaffolds](tests/scaffolds/REVIEW.md) | Contract v2 item 6 (rev 2) | tier consumers coded B/C-only genomes as confident absences in the tier-A analysis. |
| BM-06 | [scaffolds](tests/scaffolds/REVIEW.md) | build_master_table.py:256,305-321 | blank or NA ecology cells, and genomes absent from the species table, were imputed as 0 (absence). |
| D3 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:330 | D statistic always failed on real species trees |
| D4 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:131-135 | B/C-only genomes coded as absences |
| D5 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:224-226,272-275 | a single phyloglm fit is not reliable |
| D6 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:226,238-241 | bootstrap CI built on unchecked refits |
| D7 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:197-215 | naive baseline wrong under separation; no Fisher test |
| D8 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:112-113`, `cooccurrence.R:87-88 | silent re-rooting: an unrooted tree was passed to `phangorn::midpoint` inside `tryCatch`. phangorn was never loaded (`requireNamespace` absent), so any failure silently kept the unrooted tree. |
| C2 | [distribution](tests/distribution/REVIEW.md) | cooccurrence.R:50-56,144 | focal Prometheus vector was the raw DefenseFinder call |
| C3 | [distribution](tests/distribution/REVIEW.md) | cooccurrence.R:175-182 | perfect co-occurrence labelled "ambiguous" |
| C4 | [distribution](tests/distribution/REVIEW.md) | cooccurrence.R:190-195 | runtime: phytools' default `lik.func = "lik"`, run serially, took ~20 s per system for 120 tips and scales linearly (500 tips: 52 s). |
| M1 | [distribution](tests/distribution/REVIEW.md) | build_defense_matrix.py:22,43 | numeric-looking ids lost their systems |
| M2 | [distribution](tests/distribution/REVIEW.md) | build_defense_matrix.py:21-25,41-50 | swallowed failures: an unreadable master gave a header-only output; an unreadable or 0-byte systems table hit `except: pass` and gave zeros; a missing table gave zeros. |
| P1 | [distribution](tests/distribution/REVIEW.md) | compute_prevalence.py:81-90 | B/C-only genomes were confident absences |
| P2 | [distribution](tests/distribution/REVIEW.md) | `compute_prevalence.py` (design 4.1) | only one tier set reported: the design requires prevalence "at tier A and tier A+B". |
| P3 | [distribution](tests/distribution/REVIEW.md) | compute_prevalence.py:51-78 | silent failures and silent denominators |
| PL1 | [distribution](tests/distribution/REVIEW.md) | plot_distribution.py:169 | exclusion drawn as co-occurrence |
| PE1 | [distribution](tests/distribution/REVIEW.md) | plot_distribution_extra.py:47-52 | enrichment flags compared mismatched tier sets |
| R1 | [distribution](tests/distribution/REVIEW.md) | distribution.smk:5,27,105 | rules without logs: `defense_matrix`, `compute_prevalence` and `plot_distribution` had no `log:`. |
| EVO-02 | [evolution](tests/evolution/REVIEW.md) | codon_align.py:124 | `-nogap` removes every column with a gap in ANY sequence. For a multidomain family (DUF3320 is absent in 54%) that discards most informative columns, and it breaks the codon-column ↔ protein-column correspondence the site→domain map relies on: the residue counter drifts. |
| EVO-03 | [evolution](tests/evolution/REVIEW.md) | codon_align.py:61-150 | Silent failures: a missing pal2nal, a non-zero exit or an exception all wrote an empty file with exit 0. |
| EVO-04 | [evolution](tests/evolution/REVIEW.md) | codon_align.py:73 | `dict(read_fasta(all_cds.ffn))` loaded every CDS of every genome (millions of records) into an 8 GB job |
| EVO-07 | [evolution](tests/evolution/REVIEW.md) | run_reconciliation.py:224-228 | The species tree went to GeneRax unchanged. Its IQ-TREE support values (`100/100`, non-unique) become species-node NAMES, so per-species counts merge. |
| EVO-08 | [evolution](tests/evolution/REVIEW.md) | run_reconciliation.py:208-243,267-330 | Every failure was swallowed (missing generax/ALE, non-zero exit, exceptions → count-0 rows, exit 0) |
| EVO-10 | [evolution](tests/evolution/REVIEW.md) | run_reconciliation.py:263-371 | ALE: the parser never matched the `Total` line of `.uml_rec` (its regexes need a number after "duplications"). |
| EVO-12 | [evolution](tests/evolution/REVIEW.md) | ancestral_states.R:266 | `rate_gain`/`rate_loss` were always NA: Q was read from `attr(s, "Q")`, but phytools stores it as the element `s$Q` |
| EVO-13 | [evolution](tests/evolution/REVIEW.md) | ancestral_states.R:246,340 | Model fixed at ARD with the implicit root prior "equal". `n_independent_acquisitions` was simply the mean number of gains, which gain-loss-regain histories and transient gains inflate. |
| EVO-14 | [evolution](tests/evolution/REVIEW.md) | ancestral_states.R:296-338 | The structured null had four defects: (a) the `sim.Mk` root state was drawn uniformly (`anc=NULL`); |
| EVO-15 | [evolution](tests/evolution/REVIEW.md) | ancestral_states.R:178-424 | Every error was swallowed by a top-level `tryCatch` (header-only output, exit 0). |
| EVO-18 | [evolution](tests/evolution/REVIEW.md) | run_selection.py:466-492 | The gene tree was pruned to the codon records only in intact-only mode; otherwise extra tips break HyPhy. |
| EVO-19 | [evolution](tests/evolution/REVIEW.md) | run_selection.py:150-246 | MEME β+ was never parsed: headers are HTML (`&beta;<sup>+</sup>`). Substring header matching was fragile (2.5.103 inserts `LRT MEME vs FEL`). |
| EVO-20 | [evolution](tests/evolution/REVIEW.md) | run_selection.py:343-405 | Site→domain used one "reference" sequence and SUBSTRING domain matching |
| EVO-21 | [evolution](tests/evolution/REVIEW.md) | run_selection.py:298-335 | RELAX labelled only TERMINAL branches. ANME hosts were tagged Reference. |
| EVO-22 | [evolution](tests/evolution/REVIEW.md) | run_selection.py (whole) | Recombination was ignored: every HyPhy test assumes one tree, and domain shuffling is the stage's own hypothesis |
| EVO-24 | [evolution](tests/evolution/REVIEW.md) | domain_trees.py:35-48,77-91,104,125 | SUBSTRING domain matching with its own domain logic; ALI coordinates; |
| EVO-25 | [evolution](tests/evolution/REVIEW.md) | evolution.smk:22-28 | `mafft --auto` (FFT-NS-2 for ~2.2-kaa sequences) contradicted the documented L-INS-i, and `evolution.align_tool` was never read. |
| EVO-26 | [evolution](tests/evolution/REVIEW.md) | evolution.smk:56-63 | `iqtree2 … \|\| true` plus `cp … \|\| echo "();"` turned an IQ-TREE crash into an empty tree, which every downstream rule then accepted |
| EVO-33 | [evolution](tests/evolution/REVIEW.md) | run_selection.py:113-136 | HyPhy was never seeded. It seeds its RNG from time+pid, and BUSTED draws a random initial grid, so reruns gave different results. |
| ECO-04 | [ecology](tests/ecology/REVIEW.md) | aggregate_viral.py pristine:16-37,65-92 (now `summarize_genome` L107-176) | the "proviral load" counted every geNomad viral sequence, i.e. whole viral contigs (free/extrachromosomal viruses or assembly-split prophages) as proviruses; the R-loop count used all sequences |
| ECO-05 | [ecology](tests/ecology/REVIEW.md) | aggregate_viral.py pristine:91 | `frac_rloop_proviruses = 0.0` for a genome with no provirus: "no provirus" was coded as "proviruses without R-loops" in the association test |
| ECO-07 | [ecology](tests/ecology/REVIEW.md) | rloop_predict.py pristine:49-122 (now L69-232) | the "QmRLFS-style" detector was not QmRLFS: RIZ = >=3 G3-runs inside a 100-bp window (no m2 model, no 1-10 nt run linkers, no 50 %-G RIZ rule); |
| ECO-12 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:124 | response = `prometheus_present` (any validated tier): B/C-only genomes were confident presences, and Contract v2 item 6 (rev 2) requires B/C-only genomes to be EXCLUDED from the tier-A analysis |
| ECO-13 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:68-69,146-154; ecology.smk pristine:116 | config `ecology.predictors` (biome_envo, oxygen, ...) are not master columns, and the script silently ignored `--predictors` (autodetect, undocumented in the CLI) |
| ECO-14 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:212-242 | binary niche predictors were unadjusted (only the viral terms had covariates), raw completeness was unscaled, and blank niche cells were coded 0 |
| ECO-15 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:165-211 | continuous viral predictors entered on their raw scale (log-odds per provirus; heavy-tailed counts) |
| ECO-17 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:91-94 | an unrooted species tree was silently midpoint-rooted |
| ECO-19 | [ecology](tests/ecology/REVIEW.md) | plot_ecology.py pristine:76-82,163-171 | intervals drawn as estimate +/- 1 SE (68 %) on a forest labelled as effects, and the viral panel printed a naive Mann-Whitney p (not phylogeny-controlled) as its inference; presence = any tier; non-pass_primary genomes included |
| ECO-20 | [ecology](tests/ecology/REVIEW.md) | plot_ecology_extra.py pristine:80-88,94-104,123-160,33-40 | the niche forest mixed per-SD viral effects with binary niche log-odds; |
| S2 | [synthesis](tests/synthesis/REVIEW.md) | plot_tree_master.R:87-92, 179, 190 | ggtree/ggplot2 missing produced a placeholder (exit 0); ggnewscale was optional: without it every ring's `scale_fill_*` replaces the previous one, so all rings silently share the last palette/legend. |
| S3 | [synthesis](tests/synthesis/REVIEW.md) | plot_tree_master.R:129-135, 146-158 | Tips without a master row were painted tier "none", methanogen "no", subtype "none", niche richness 0, group "no": missing data shown as biological absence. |
| S4 | [synthesis](tests/synthesis/REVIEW.md) | plot_tree_master.R:81, 129 | `ape::read.tree` keeps the single quotes of quoted Newick labels, and tips were matched on `genome_id` only (`tree_tip` ignored): quoted labels (dendropy writes them for labels with special characters) silently never match. |
| S5 | [synthesis](tests/synthesis/REVIEW.md) | plot_tree_master.R:132-135 | ANME genomes (methanogen=False, is_anme=True; classify_methanogen.py:284-290) were painted plain "Methanogen: no", although README §10 says the master tree paints "methanogen and ANME status". |
| S8 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:242-253 | `load_table` returned `[]` for a missing file, a 0-byte file or ANY parse exception (log warning only): a crashed producer silently removed its whole family from the registry and the FDR. |
| S9 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:24-38 | Confirmatory tagging by substring: `"phyloglm" in test and "methanogen" in term` made BOTH phyloglm models confirmatory (two tests of one pre-registered hypothesis); |
| S10 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:288-293 | The BH family included the independence baselines (naive_glm; the revised producer adds naive_fisher). |
| S11 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:265, 288-293 | `--fdr` was accepted and ignored (always BH): a configured BY/holm silently became BH. |
| S12 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:143-157 | Selection estimate = `value`; for BUSTED run_selection.py writes `value = p` (statistic `episodic_diversifying_p`), so the registry's estimate held a p-value. |
| S13 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:181-204 | run_reconciliation.py writes count-0 rows when a tool did not run ("ALE not on PATH (skipped)", "generax failed (exit N)", ...). |
| S18 | [synthesis](tests/synthesis/REVIEW.md) | plot_summary_forest.py:38-101 | The "forest of confirmatory results" put every confirmatory row on ONE x-axis with a null line at 0: log-odds (null 0), RELAX K (a ratio with null 1; |
| CD-04 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | default config (no bacterial CDS) aborts / archaea-only "cross-domain" dN/dS |
| CD-05 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | `python` run inside the R-only conda env (niche) |
| CD-06 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | validation universe taken from the QC table |
| CD-07 | [crossdomain](tests/crossdomain/REVIEW.md) | extract_bacterial_prometheus.py | substring system matching |
| CD-08 | [crossdomain](tests/crossdomain/REVIEW.md) | extract_bacterial_prometheus.py | genes table never split; no tables for zero-system genomes; stale tables |
| CD-09 | [crossdomain](tests/crossdomain/REVIEW.md) | extract_bacterial_prometheus.py | AntiDefenseFinder rows counted as defense |
| CD-10 | [crossdomain](tests/crossdomain/REVIEW.md) | extract_bacterial_prometheus.py | missing proteome/protein silently skipped; gzip read as text |
| CD-11 | [crossdomain](tests/crossdomain/REVIEW.md) | extract_bacterial_prometheus.py | NCBI `cds_from_genomic` CDS never matched |
| CD-12 | [crossdomain](tests/crossdomain/REVIEW.md) | extract / build_bacterial_table | genome ids not validated (blank, duplicate, `__`, zero overlap) |
| CD-13 | [crossdomain](tests/crossdomain/REVIEW.md) | prune_tree_to_taxa.py | GTDB tip matching / no relabel / silent trivial tree |
| CD-14 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain_analysis.py | rooting on a non-monophyletic outgroup; supports misassigned on reroot |
| CD-15 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain_analysis.py | direction from the ML Fitch root state alone |
| CD-16 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain_analysis.py | archaeal acquisitions counted on unsupported splits |
| CD-17 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain_analysis.py | dendropy missing / unclassified tips handled silently |
| CD-18 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | reference ids not sanitised → label lookups fail after IQ-TREE |
| CD-19 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain_analysis.py | no synonymous-saturation control for cross-domain dN/dS |
| CD-20 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | HyPhy gene tree ≠ codon-alignment taxa |
| CD-21 | [crossdomain](tests/crossdomain/REVIEW.md) | build_bacterial_table.py | `pass_primary` ignores contamination |
| CD-22 | [crossdomain](tests/crossdomain/REVIEW.md) | build_bacterial_table.py | taxonomy formats, duplicates, universe; coupling to build_master_table internals |
| CD-23 | [crossdomain](tests/crossdomain/REVIEW.md) | build_bacterial_niche_master.py | unknown niche coded 0; missing niche file passed through |
| CD-24 | [crossdomain](tests/crossdomain/REVIEW.md) | plot_crossdomain_extra.py, plot_viral_crossdomain.py | Contract-v2 presence rule ignored (tier-agnostic prevalence) |
| CD-25 | [crossdomain](tests/crossdomain/REVIEW.md) | plot_viral_crossdomain.py | 0/0 R-loop fraction plotted as 0; unlabelled p |
| CD-26 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | geNomad/CheckV placeholders with wrong headers mask crashes |
| CD-27 | [crossdomain](tests/crossdomain/REVIEW.md) | extract / crossdomain.smk | bacterial fragmentation by length proxy, archaeal by Prodigal flags |
| CD-28 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | bacterial validation lacks the archaeal evidence/thresholds (motif alignments, focal system, placement cut-offs) |
| CD-29 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | tier sets/thresholds not passed to tier-based responses; no A+B sensitivity fits |
| CD-41 | [crossdomain](tests/crossdomain/REVIEW.md) | build_bacterial_table.py | genome_size blank for compressed assemblies (decode error swallowed) |
| CD-43 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | geNomad resumed stale intermediates (no `--restart`); no `--splits`; DB check = any directory |
| CD-44 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | CheckV resumed a stale tmp/ (no `--restart`); DB check = any directory; no post-run table/header check |
| D-03 | [disorder](tests/disorder/REVIEW.md) | workflow/envs/disorder.yaml:12 | `blast=2.15.0` (libzlib <1.3) makes the env unsolvable together with any current conda-forge PyTorch. |
| D-04 | [disorder](tests/disorder/REVIEW.md) | workflow/envs/disorder.yaml (whole file) | metapredict has no linux wheel and its sdist must compile a required Cython extension (no compiler in the env → install fails); pip would also pull the multi-GB CUDA torch wheel; |
| D-06 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:145-186 | IUPred3: `res[0] if isinstance(res,(list,tuple))` took the first *score* of a bare list (TypeError → silently switched to the CLI path); a failure for one locus returned None, so that locus became metapredict-only while still labelled `metapredict+iupred3`; library and CLI paths could use different smoothing defaults. |
| D-07 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:189-213 | segmasker missing or failing → None → no LCR masking while `frac_disordered_nonlcr` was still reported as the LCR-excluded "preferred" metric. |
| D-08 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:217-234 | Missing/unreadable catalog or missing columns returned None → the intact filter was silently disabled; a missing `fragmented` value was kept. |
| D-09 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:61-71 | Consensus silently truncated to the shorter track (`min(n, len)`) and treated a missing second score as "ordered". |
| D-11 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:101 | `frac_disordered_nonlcr` divided by all residues, so masked LCR residues still sat in the denominator as "ordered" — not "disorder excluding LCR". |
| D-12 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:237-251; disorder.smk:34; disorder_stats.R (whole) | Background column auto-picked (first name containing "background" or "frac"), any error → {}; the file was not a rule input; partial coverage mixed normalised and raw loci; and the statistics never used the background at all. |
| D-15 | [disorder](tests/disorder/REVIEW.md) | map_disorder_to_alignment.py:24-52, 100-101 | Own substring name matching (Contract v2 item 3) and a domain track taken from a single "reference" (longest) locus: when that locus lacks DUF3320 (≈54 % of Prometheus loci) the DUF3320 columns are labelled linker for everyone. |
| D-16 | [disorder](tests/disorder/REVIEW.md) | map_disorder_to_alignment.py:98-122 | No check that aligned rows are the scored sequences: a stale/re-ordered alignment, a residue edited by the aligner, or a scored locus missing from the alignment shifted or dropped data silently; ragged rows (family_align copies the raw FASTA for < 3 sequences) were treated as an MSA. |
| D-19 | [disorder](tests/disorder/REVIEW.md) | plot_disorder.py:26-27 | Hard-coded `DOMAIN_COLORS` including the non-existent `nuclease` category. |
| D-20 | [disorder](tests/disorder/REVIEW.md) | plot_disorder.py:96-118 | Left join + `_bool(NaN)` coded genomes absent from the master table as non-methanogens; |
| D-23 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R:118-123 | Temperature covariates entered only the niche models, but methanogen status is itself confounded with growth temperature (many methanogens are (hyper)thermophiles) and disordered content tracks temperature. |
| D-24 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R:70-75 (Contract v2 item 6) | Per-genome disorder averaged every scored locus regardless of tier; no exclusion or count. |
| D-25 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R:94-96 | ANME (methanogen = False by construction) sat in the non-methanogen reference group, so mcr-carrying methanotrophs diluted the methanogen contrast ("ANME separated" in the design). |
| D-30 | [disorder](tests/disorder/REVIEW.md) | disorder.smk:52-68, 135-157, 13, 131 | `disorder_map` and `plot_disorder` had no `log:`; the host background was not an input; tier sets, threads, network version and plot formats were not passed; helper names `_d` / `_bact_disorder_locus` live in Snakemake's shared namespace. |

### Minor — robustness, reproducibility, documentation (108)

| ID | Stage | Location | Issue |
|---|---|---|---|
| PQ-20 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | qc.smk:30-32 | sub-20-kb fallback test did not match Prodigal's count |
| PQ-21 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | qc.smk:33 | translation table implicit; meta mode ignores it |
| PQ-22 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py:40 | empty representatives.txt began with "genome_id" |
| PQ-23 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py:181 | quality_tier hard-coded to "A" |
| PQ-24 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | filter_collapse.py:140-147 | representative tie-break |
| PQ-25 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | merge_annotations.py:130-149 | "KB" heuristic instead of a real check |
| PQ-26 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | batch_genomes.py / merge_annotations.py:72 | a mid-line '#' truncates a metadata row |
| PQ-27 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | plot_preprocess.py:40,55,59-61 | "pass" did not mean pass |
| PQ-28 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | plot_qc.py:70-73,97-114 | hard-coded thresholds, forbidden label, unit guesses |
| PQ-29 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | plot_qc_extra.py:47-63 | degenerate groups and a mislabelled fallback |
| PQ-30 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:220, qc.smk:44,62 | rules without log: |
| PQ-31 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:220-241, qc.smk:62-78 | plot outputs vs PLOT_FORMATS |
| PQ-32 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | aggregate_qc.py | dead code removed |
| PQ-33 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk gtdbtk_batch | stale GTDB-Tk state reused |
| PQ-34 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:96-108 | CheckM2 intermediates |
| PQ-35 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | preprocess.smk:125 | bbmap_batch reserved 8 threads |
| PQ-36 | [preprocess_qc](tests/preprocess_qc/REVIEW.md) | batch_genomes.py | batch paths with spaces |
| DD-07 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Placeholder headers did not match DefenseFinder's tables |
| DD-08 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | DefenseFinder index written into another stage's directory |
| DD-09 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Substring system matching, and every gene of a system taken as ProA |
| DD-16 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Reference labels did not survive IQ-TREE tip renaming |
| DD-21 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Missing dendropy or an unparsable tree silently made everything Unplaced |
| DD-24 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Unused `--reps`, stop codons counted, invalid spec accepted |
| DD-25 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Missing tiers were plotted as Tier C |
| DD-26 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Validation-tree figure always exited 0 |
| DD-32 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Extra diversity figure had a hard-coded domain list |
| DD-33 | [detection_diversity](tests/detection_diversity/REVIEW.md) | — | Rules without `log:` |
| SC-08 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:79 | the PMSF guide tree used `-m LG+G -fast`; the IQ-TREE docs obtain the guide tree under the simpler LG+F+G model. |
| SC-09 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:49-50,98 | `rooting_method.txt` (cited in the README) was written but not declared as an output. |
| SC-10 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:161-192 | `build_master_table` had no `log:`. |
| SC-11 | [scaffolds](tests/scaffolds/REVIEW.md) | scaffolds.smk:110-112,135-141 | `methanogen.mcr_kos` was unused, and `mcr_required` was passed but ignored by the script. |
| RT-08 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:155 | an unknown `--method` was silently midpoint-rooted. |
| RT-09 | [scaffolds](tests/scaffolds/REVIEW.md) | root_tree.py:35-42,94 | the rooting record was a single free-text line, and log write errors were swallowed. |
| MS-05 | [scaffolds](tests/scaffolds/REVIEW.md) | mcr_search.py:59-78 | GA mode compared only GA1 with the full score. HMMER 3.4 `--cut_ga` lists in `--tblout` a sequence that passes GA1 but has no domain >= GA2 (N = 0, verified). |
| MS-06 | [scaffolds](tests/scaffolds/REVIEW.md) | mcr_search.py:196-222 | one protein could count as two subunits, inflating `n_mcr`. |
| MS-07 | [scaffolds](tests/scaffolds/REVIEW.md) | mcr_search.py:218 | `--pmap` was ignored and the genome was taken from `genome_of`. |
| MS-08 | [scaffolds](tests/scaffolds/REVIEW.md) | (missing) | copy number, fragment status and the thresholds used were not reported. |
| CM-04 | [scaffolds](tests/scaffolds/REVIEW.md) | classify_methanogen.py:149 | IQ-TREE ran without `--seed`. |
| CM-05 | [scaffolds](tests/scaffolds/REVIEW.md) | classify_methanogen.py:118 | the whole combined proteome was loaded into a dict in an 8 GB rule. |
| CM-06 | [scaffolds](tests/scaffolds/REVIEW.md) | classify_methanogen.py:230-231 | `--required` was parsed but McrA was hard-coded. |
| CM-07 | [scaffolds](tests/scaffolds/REVIEW.md) | classify_methanogen.py:245-251 | only genomes with mcr hits were written. |
| BM-07 | [scaffolds](tests/scaffolds/REVIEW.md) | build_master_table.py:19-20,376 | `--envo-map` was silently ignored, and config `ecology.envo_map` / `predictors` refer to free-text columns that the species table does not have. |
| BM-08 | [scaffolds](tests/scaffolds/REVIEW.md) | (missing) | the README lists contig count and N50 as quality covariates, but they were absent from the master table. |
| D9 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:181-191,271-275 | covariate coding and a biased fallback |
| D10 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:333-339 | D p-values: only in a note, and can be 0 |
| D11 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:309-312`, `cooccurrence.R:186-195 | fitPagel optimiser failures: `lik.ratio = 2(logL_dep - logL_indep)` and `P = pchisq(LR, df = 4)` are taken as given. |
| D12 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:87`, `cooccurrence.R:68 | RNG and locale: one `set.seed` at start, so bootstrap and permutation draws depended on which earlier tests ran. |
| D13 | [distribution](tests/distribution/REVIEW.md) | association_tests.R:121`, `cooccurrence.R:95,106 | numeric-looking ids: readr type guessing turns zero-padded ids such as `000123` into numbers that no longer match tree tips. |
| D14 | [distribution](tests/distribution/REVIEW.md) | `association_tests.R` (new) | zero-length branches: `phylo.d` stops on any zero-length branch ("cannot currently simulate"), and the Brownian VCV is singular. |
| C5 | [distribution](tests/distribution/REVIEW.md) | `cooccurrence.R` (design 4.6) | no genome-size correction: the design requires co-occurrence "correcting for genome size". |
| C6 | [distribution](tests/distribution/REVIEW.md) | cooccurrence.R:169 | prevalence window: the window [3, n-3] admits systems seen in only 3 genomes. |
| C7 | [distribution](tests/distribution/REVIEW.md) | cooccurrence.R:139-143,220-222 | strictness: `bin()` turned NA or garbage into 0, and an unknown `--fdr` reached `p.adjust` inside the blanket `tryCatch` (header-only). |
| M3 | [distribution](tests/distribution/REVIEW.md) | build_defense_matrix.py:44-48 | anti-defense systems counted as defense |
| P4 | [distribution](tests/distribution/REVIEW.md) | compute_prevalence.py:25-35 | Wilson edge cases: n = 0 returned (0, 0), a CI for an undefined proportion; z was rounded to 1.96. |
| P5 | [distribution](tests/distribution/REVIEW.md) | compute_prevalence.py | no overall row; unassigned genomes dropped unannounced |
| PL2 | [distribution](tests/distribution/REVIEW.md) | plot_distribution.py:100-133 | the forest mixed scales: every row with an estimate was drawn on one "log-odds" axis, including D (a different statistic) and the separated GLM (25 ± 7475). |
| PL3 | [distribution](tests/distribution/REVIEW.md) | plot_distribution.py:45-48,95-98,153-156`, `plot_distribution_extra.py | placeholders masked interface errors; formats hard-coded |
| PE2 | [distribution](tests/distribution/REVIEW.md) | plot_distribution_extra.py:101-118 | strongest results dropped from the OR figure |
| EVO-09 | [evolution](tests/evolution/REVIEW.md) | run_reconciliation.py:224-229 | The pipeline seed was not passed (GeneRax defaults to 123), `--strategy SPR` and `LG+G` were hard-coded, and the model was not validated (`-r Auto` is invalid in 2.1.3) |
| EVO-16 | [evolution](tests/evolution/REVIEW.md) | ancestral_states.R:96-130,246 | Methanogen edges were computed on the input edge order, while `make.simmap` reorders the tree cladewise. |
| EVO-23 | [evolution](tests/evolution/REVIEW.md) | run_selection.py | Stop codons in the codon alignment only failed inside HyPhy; aBSREL used HyPhy's fixed 0.05; |
| EVO-27 | [evolution](tests/evolution/REVIEW.md) | evolution.smk:97,272 | `combined_cds` and `plot_selection` had no `log:` |
| EVO-28 | [evolution](tests/evolution/REVIEW.md) | evolution.smk:136-162 | Reconciliation asked for T_MED threads, but GeneRax with one family is a single process (it parallelises over families via MPI) |
| EVO-29 | [evolution](tests/evolution/REVIEW.md) | plot_tanglegram.R:60-150 | Errors were swallowed (placeholder, exit 0). The species tree was not pruned to carriers. |
| EVO-30 | [evolution](tests/evolution/REVIEW.md) | plot_selection.py:47-58 | Raw positive-site counts per domain (not normalised by domain size, no uncertainty). |
| EVO-31 | [evolution](tests/evolution/REVIEW.md) | plot_evolution_extra.py:31,77,126 | Undefined dN/dS was drawn at 1.0 (`fillna(1.0)`), i.e. as "neutral". The domain list was hard-coded (incl. |
| EVO-32 | [evolution](tests/evolution/REVIEW.md) | common.smk (driver) | `evolution/domains/rf_distances.tsv` is not in `all_targets()`, so `domain_trees` never runs in `all` (PER_DOMAIN_TREES is unused) |
| EVO-35 | [evolution](tests/evolution/REVIEW.md) | evolution.smk:22-28; domain_trees.py | MAFFT needs an absolute TMPDIR. The 7.526 wrapper sets `MAFFT_TMPDIR=$TMPDIR` (L1081-1082), creates `mktemp -dt` (L1085) and `pushd`es into it (L2295), so a relative TMPDIR (the driver reports `./.tmp` in its sandbox) resolves against a changed cwd. |
| ECO-06 | [ecology](tests/ecology/REVIEW.md) | aggregate_viral.py pristine:42-62 | mean CheckV completeness averaged over all CheckV rows (incl. whole viral contigs), matched columns by guesswork |
| ECO-08 | [ecology](tests/ecology/REVIEW.md) | ecology.smk pristine:60-81; config `ecology.rloop` | `min_rlfs_score` was never passed; `window` / `riz_min_g` semantics did not correspond to QmRLFS |
| ECO-09 | [ecology](tests/ecology/REVIEW.md) | rloop_predict.py pristine:148-152 | CRLF files kept `\r` inside sequences (`utils.read_fasta` strips only `\n`), breaking the RIZ regex; duplicate ids and non-FASTA input passed silently |
| ECO-16 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:217-230 | (quasi-)separation not flagged; phyloglm warnings not reported; min-data rules inconsistent (n >= 5; no response-class minimum for binary terms) |
| ECO-18 | [ecology](tests/ecology/REVIEW.md) | ecology_stats.R pristine:120-121,153 | without `--viral` (bacterial niche call) two meaningless "skipped" viral rows were emitted; readr/dplyr were required for no reason |
| ECO-21 | [ecology](tests/ecology/REVIEW.md) | ecology.smk pristine:60-156 | `rloop_predict`, `viral_summary`, `plot_ecology` had no `log:`; `ecology_stats` reserved T_MED threads for single-threaded R; plot formats not passed |
| ECO-22 | [ecology](tests/ecology/REVIEW.md) | ecology.smk | CheckV/geNomad DB paths unquoted; ecology tier sets and covariates not configurable |
| S6 | [synthesis](tests/synthesis/REVIEW.md) | plot_tree_master.R:236-239 | R-loop load not shown (README §10 claims "proviral / R-loop load"). aggregate_viral.py writes `frac_rloop_proviruses = 0` when a genome has no provirus, which would read as "no R-loop-prone proviruses". |
| S7 | [synthesis](tests/synthesis/REVIEW.md) | plot_tree_master.R:160-172, 199-239 | Layout: rings started at offset 0 under tip labels drawn at offset 0.001, so long IDs were overprinted; the height depended on the tip count only, so legends were clipped (the proviral legend was cut off in the 9-tip baseline); |
| S14 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:106-121 | Co-occurrence estimate = `statistic` = Pagel LR (or OR on the Fisher fallback): mixed and not an effect size. |
| S15 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:89-178 | No intervals for ecology rows (SE available), no record of the estimate scale/null, and the D and naive rows were not labelled except by `hypothesis`. |
| S16 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py (no handling) | The same test at several tier sets (the A vs A+B sensitivity of design 7.3) would be double-counted in one BH family and could be double-tagged confirmatory. |
| S17 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py:41-51 | p-values outside [0, 1] were accepted. |
| S19 | [synthesis](tests/synthesis/REVIEW.md) | plot_summary_forest.py:48 | Confirmatory rows with a missing estimate were silently dropped, so a failed headline fit vanished from the figure. |
| S20 | [synthesis](tests/synthesis/REVIEW.md) | plot_summary_forest.py:22-28, 101 | A missing registry produced a placeholder (silent); a registry without the columns gave a placeholder; only png/svg were written even if `plotting.formats` asks for more (rule outputs `expand(PLOT_FORMATS)`). |
| S21 | [synthesis](tests/synthesis/REVIEW.md) | synthesis.smk:58-109 | `results_registry` and `plot_summary_forest` had no `log:`. Output-directory params were hard-coded (snakemake --lint). |
| S22 | [synthesis](tests/synthesis/REVIEW.md) | results_registry.py (new `--confirmatory`) | Found while testing the fix: the disorder family's test label `phylolm(disorder~pred, lambda)` contains a comma, so a comma-separated key list cannot express it. |
| CD-30 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | `grep -c … \|\| echo 0` |
| CD-31 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | 11 rules without `log:` |
| CD-32 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | GTDB tree passed as `params` (no dependency/missing-file check) |
| CD-33 | [crossdomain](tests/crossdomain/REVIEW.md) | plot_crossdomain_extra.py | hard-coded domain list incl. `nuclease` |
| CD-34 | [crossdomain](tests/crossdomain/REVIEW.md) | plot scripts | `plotting.formats` ignored (png+svg hard-coded) |
| CD-35 | [crossdomain](tests/crossdomain/REVIEW.md) | plot_crossdomain_tree.py | any error → placeholder; unrooted tree drawn |
| CD-36 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | viral summary `--samples` = QC table; no log |
| CD-37 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | hmmscan alignments dumped into the log |
| CD-38 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | no empty-output checks after MAFFT/ClipKIT/IQ-TREE (ClipKIT exits 0 on errors) |
| CD-39 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | missing per-genome systems table counted as zero systems |
| CD-40 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | codon-alignment interface: bacterial-only `--proteins`, whole archaeal CDS copied |
| CD-42 | [crossdomain](tests/crossdomain/REVIEW.md) | crossdomain.smk | cross-domain selection without `--seed` (HyPhy RNG seed) |
| CD-45 | [crossdomain](tests/crossdomain/REVIEW.md) | tests/crossdomain/run_tests.py | parity canary did not cover the geNomad/CheckV invocations |
| D-10 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:68 | IUPred threshold applied as `>=`; IUPred's convention is score **> 0.5**. |
| D-13 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:126; disorder.smk:36-48 | Network version implicit (metapredict's default moved V2→V3 in 3.0), device auto-selected (CUDA/MPS scores differ slightly), torch threads unbounded on cluster nodes, one call per sequence. |
| D-14 | [disorder](tests/disorder/REVIEW.md) | predict_disorder.py:319-337 | The per-residue table (~10⁷ rows for thousands of ~2,200-aa loci) was built in memory; no residue letters for downstream identity checks. |
| D-17 | [disorder](tests/disorder/REVIEW.md) | map_disorder_to_alignment.py:92-96 | Bad JSON → `{}` and an empty/unmatched domtbl silently produced an all-NA domain track (and a "need both domain and linker" placeholder downstream). |
| D-18 | [disorder](tests/disorder/REVIEW.md) | map_disorder_to_alignment.py:135-136 | `mean_disorder` and `frac_loci_disordered` are identical by construction (binary consensus). |
| D-21 | [disorder](tests/disorder/REVIEW.md) | plot_disorder.py:77, 124, 148; disorder.smk:146-149 | Always wrote png+svg while the rule expects `plotting.formats`; one SVG rectangle per alignment column; gappy low-occupancy columns plotted at full weight. |
| D-26 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R:90, 99, 107 | No minimum group size for methanogen; misspelt temperature columns silently ignored (`intersect`); non-binary niche values coerced; the phylolm error text and λ-at-bound warnings lost. |
| D-27 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R (missing test) | The README promises comparison across lineages, but no lineage / phylogenetic-signal test existed. |
| D-28 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R:107-112 | No interval estimates. |
| D-29 | [disorder](tests/disorder/REVIEW.md) | disorder_stats.R:53-58 | Unrooted tree left unrooted when phangorn was absent (phylolm accepts it silently); branch lengths never validated. |
