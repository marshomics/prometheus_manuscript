# Table schemas (the I/O contract every script adheres to)

## Contract v2 (pipeline review) — cross-cutting rules every stage follows

1. **Gene ids are contig-aware and gembase-compatible.** Per-genome protein and CDS
   ids are `{genome}_{contig#}_{gene#}` (`utils.make_gene_id`; contig# = 1-based order of
   the contig in the assembly, gene# = 1-based order of the gene on that contig). The
   replicon is everything before the last `_` (`{genome}_{contig#}`), so DefenseFinder is
   run with `--db-type gembase` and never assembles a system across a contig junction
   (with a protein FASTA it otherwise defaults to `ordered_replicon` = one replicon per
   genome). Gene adjacency = same contig#, consecutive gene#. `utils.split_gene_id`
   parses ids; combined ids stay `genome__{gene id}`.
2. **Per-gene coordinates travel in `genes/{genome}/{genome}.genes.tsv`** (prodigal_genes):
   `protein_id, contig, contig_index, gene_index, start, end, strand (+/-), partial,
   contig_length, length_aa`, then (additive) `transl_table` — `partial` is Prodigal's
   two-character code, kept as a string (`00` complete; `10`/`01`/`11` = runs off the
   left/right/both contig ends); `length_aa` excludes the terminal `*`; `transl_table` is the
   table Prodigal used (qc.prodigal_translation_table, 11, in single mode; the chosen model's
   table, possibly 4, in the < 20,000-bp meta fallback). This is what makes the
   fragmentation flag and the ±N-ORF defense-island test direct rather than proxies.
3. **One domain specification, exact matching.** config `prometheus.domains` (key order =
   canonical N->C order DUF4011, DNA2 [SF1 helicase AAA_11/AAA_12], REase [MTES_1575],
   DUF3320; `required: false` for DUF3320, absent in ~54% of Prometheus) is consumed only
   through `utils.load_domain_spec / hit_category / architecture_tokens /
   is_canonical_order / domain_intervals / other_domains`. No script hardcodes a domain
   list; plots derive categories from the spec or from the `has_*` columns present.
   There is no `nuclease` category.
4. **No silent failure.** A tool that fails makes its rule fail. `|| true`, header-only
   placeholder outputs and `();` trees are allowed only for a genuinely empty *input*
   (e.g. zero candidates), never to mask a crashed tool or a missing database, and every
   such line carries a `# nonfatal: <reason>` comment. A missing required R/Python package
   is an error, not a silently skipped test.
5. **Interfaces are additive.** Table columns listed below may be added to, not renamed or
   removed; CLI flags of a script invoked from another stage's rule file stay
   backward-compatible (new flags must be optional).
6. **(rev 2) One presence rule for every tier-based response.** For a tier set T (the
   primary `stats.confidence_tiers_primary`, or the A+B sensitivity set) and a genome whose
   best validated tier is `prometheus_tier` (A > B > C; `""` when it has no non-rejected
   locus): **present** (1) if `prometheus_tier` is in T; **absent** (0) if `prometheus_tier`
   is `""`; **excluded** (NA, dropped from that analysis, counted and reported) if it has a
   validated locus whose best tier is outside T. A B/C-only genome is therefore never coded
   as a confident absence in the tier-A analysis. Prevalence denominators, phylogenetic
   regressions, Pagel/D, co-occurrence, ancestral reconstruction and ecology all apply this
   rule.

IDs: combined protein/CDS ids are `genome__origProteinId` (see utils.ID_SEP).
Species-tree tips = `genome_id`. Gene-tree tips = combined protein id. Recover the
genome from a protein id with `utils.genome_of()`. All tables are TSV with a header.
Every script must tolerate empty inputs and still write a header-only output.

## preprocess/species_metadata.tsv  (filter_collapse.py) — the SINGLE filtered table (samples + metadata)
genome_id, filename, source, classification, gtdb_domain..gtdb_species, assembly_path,
completeness, contamination, gc, genome_size, n_contigs, n50, quality_score,
n_genomes_in_species, n_passing_in_species, group_* and niche_* ({1, 0, NA}, unioned per
species over niche_aggregation 'all' or 'passing' genomes). This one table is both
`config.samples` and `config.metadata`. Only genomes GTDB-Tk places in d__Archaea are
retained (preprocess.keep_domain); novel/unassigned GTDB species are dropped (or kept as
singletons, preprocess.unclassified_species). `assembly_path` is the staged genome
(`{config genomes_dir}/{genome_id}.{config genome_extension}`) each analysis step reads.
BBMap columns: gc = gc_avg (a FRACTION 0-1), genome_size = contig_bp (gap Ns excluded),
n_contigs = n_contigs, n50 = ctg_L50 (BBMap's N50 length in bp). quality_score =
completeness - 5 x contamination, rounded to 1e-6 so thresholds apply exactly as stated.
n_genomes_in_species counts every genome of the species (both aggregation modes). Integer
columns are written without '.0'. preprocess/representatives.txt: one id per line (empty
file when there is none); taxonomy.tsv: user_genome, classification.

## preprocess/annotated_metadata.tsv  (merge_annotations.py)
Every input genome with taxonomy, CheckM2 and BBMap values, plus (additive) batch,
source_path, gtdbtk_status {classified, unclassified, missing}, gtdbtk_marker_set,
gtdbtk_classification_method, gtdbtk_closest_genome_reference, gtdbtk_closest_genome_ani,
gtdbtk_msa_percent, gtdbtk_red_value, gtdbtk_translation_table, gtdbtk_note,
gtdbtk_warnings, checkm2_status {ok, missing}, checkm2_model, checkm2_translation_table,
checkm2_contig_n50, checkm2_genome_size, checkm2_notes, bbmap_status {ok, missing}.
Duplicate or unknown tool ids are errors; a genome a tool ran on but did not report stays
in the table with its status (never silently dropped).

## preprocess/filter_report.tsv  (filter_collapse.py)
genome_id, gtdbtk_status, gtdb_domain, gtdb_species, species_key, completeness,
contamination, quality_score, n_contigs, n50, passing, status (e.g. representative,
failed_quality, excluded_domain, no_species, ...), reasons, species_representative: one
row per INPUT genome saying why it was kept or dropped.

## preprocess/batches/manifest.tsv and batch_{b}/gtdbtk_batchfile.tsv  (batch_genomes.py)
manifest: genome_id, batch, source_path, staged_path, gtdbtk_id. The GTDB-Tk batchfile is
`absolute link path<TAB>user_{genome_id}`: GTDB-Tk 2.4.0 aborts when a user id equals a GTDB
reference id (incl. bare GCA_/GCF_ accessions), so it only ever sees prefixed ids; the
prefix is stripped from every output (batch_{b}/gtdbtk_summary.tsv = the ar53 + bac120
classify summaries plus a marker_set column). Input ids must match
[A-Za-z0-9][A-Za-z0-9_.-]* without '__' (checked at DAG build and by batch_genomes.py);
only whole-line '#' comments are allowed in the input table.

## genes/{genome}/{genome}.{faa,fna}  (prodigal_genes + rename_prodigal.py)
Prodigal proteins (.faa, terminal `*` kept), CDS (.fna) and GFF, with every id rewritten
to `{genome}_{contig#}_{gene#}` (Contract item 1; contig# = 1-based FASTA record index,
gene# = Prodigal's per-record gene index), IDENTICAL for the matching protein and CDS and
in the GFF `ID=` attribute; plus `{genome}.genes.tsv` (Contract item 2). A crashed Prodigal,
an empty input or a genome with 0 genes fails the rule (never an empty proteome). CDS
records without a matching protein are errors. The .faa feeds DefenseFinder (gembase),
the McrA search and the validation.

## qc/genome_metadata.tsv  (cp of preprocess/genome_metadata.tsv)
genome_id, source, completeness, contamination, genome_size, n_contigs, n50, gc,
quality_tier {A,B,C}, pass_primary {True,False}. quality_tier A = MIMAG high quality
(qc.tier_high_completeness / tier_high_contamination); pass_primary = completeness >=
qc.completeness_floor AND contamination <= qc.contamination_ceiling. With the Phase-1
filter (>= 90 / <= 5) every representative is tier A and pass_primary. No MAG/isolate
genome_type: quality is completeness / contamination / n_contigs / n50.

## prometheus/protein_map.tsv  (combine_proteomes.py)
protein_id (`genome__{gene id}`), genome_id, orig_id. A missing per-genome proteome, an
unreadable sample table and duplicate genome or protein ids are errors. Reference mode
(validation trees: `--candidates`, `--ref GROUP=FASTA` repeatable, `--out-labels`) writes
ids sanitised with `utils.sanitize` ([^A-Za-z0-9_.-] -> `_`), which equal IQ-TREE's tip
names; duplicates after sanitising are errors.

## prometheus/candidate_sources.tsv  (parse_defensefinder.py)
protein_id, genome_id, detect_defensefinder {0,1}, df_score, df_subtype, then (additive)
df_sys_id, df_gene_name, df_i_eval, df_profile_cov, df_seq_cov.
Candidates are the `--gene` hits (default `Prometheus__ProA`) of DefenseFinder systems
whose type or subtype EQUALS the configured system (case-insensitive; no substring
match). Detection is DefenseFinder-only, so detect_defensefinder is 1 for every
candidate. Missing per-genome tables, unknown protein ids and unreadable files are errors.

## defense/defensefinder/{genome}/{genome}_defense_finder_{systems,genes}.tsv  (defensefinder rule)
DefenseFinder 3.0.0 (macsyfinder 2.1.4) output on the gembase-named proteome, run with
`--db-type gembase --skip-model-version-check` against pinned models (defense-finder-models
3.1.0, CasFinder 3.1.0: config prometheus.models_pin / casfinder_pin / models_dir).
Outputs are checked after the run (check_defensefinder.py --check-outputs). The sentinel
results/defense/.defensefinder_ok holds key<TAB>value provenance: versions, model paths
and md5s, the Prometheus profile header, and the positive-control result.

## prometheus/prometheus_catalog.tsv  (validate_prometheus.py)
First 13 columns (original contract, order unchanged): locus_id (=protein_id), genome_id,
protein_id, length_aa (a terminal `*` is not counted), detect_defensefinder, architecture
(spec categories N->C via `utils.architecture_tokens`, repeats kept, e.g.
"DUF4011-DNA2-REase-DUF3320"), n_domains (distinct spec categories), catalytic_status
{intact,partial,degenerate,unknown}, defense_context {True,False}, tree_assignment
{Prometheus,Paralog,Unplaced}, fragmented {True,False}, confidence_tier {A,B,C,Rejected},
reason. Additive: missing_required, canonical_order, other_domains, has_DUF4011, has_DNA2,
has_REase, has_DUF3320, helicase_module, nuclease_module, catalytic_evidence,
placement_support, placement_clade, placement_note, context_method, context_n_systems,
context_systems, partial, fragment_evidence, length_in_bounds.
- catalytic_status is residue-level, read from hmmscan's own alignments (`--hmmscan-out`):
  helicase module = SF1 Walker A/motif I [AG]x4GK[ST] and Walker B/motif II DE; nuclease
  module (REase_MTES_1575, PF18741) = every consensus D/E match state carries D or E.
  intact = both modules intact; partial = one compromised; degenerate = both; unknown =
  not assessable. Without alignments a labelled domain proxy is used
  (catalytic_evidence=domain_proxy).
- tree_assignment: the validation tree is rooted on the DNA2/NurA/HerA paralog outgroup;
  a locus is Prometheus (Paralog) only inside a clade with a Prometheus (paralog)
  reference whose SH-aLRT >= prometheus.placement_min_alrt (80) and UFBoot >=
  placement_min_ufboot (95); otherwise Unplaced, with placement_note.
- defense_context: another (non-focal, non-Antidefense) DefenseFinder system within
  +/-N genes on the same contig (genes.tsv adjacency). It can raise a tier, never reject.
- fragmented: Prodigal partial flag != 00 (fragment_evidence=prodigal_partial); a length
  proxy only where no genes table exists (labelled). Fragments are never length-rejected.
- Tiers. Rejected = supported placement with a cellular paralog, or a complete ORF
  < 0.5 x min_len_aa or > 2 x max_len_aa. A = DefenseFinder call + every REQUIRED category
  (DUF4011, DNA2, REase; DUF3320 is optional) in canonical N->C order + supported
  Prometheus placement + intact catalytic residues + not fragmented + length within
  [min_len_aa, max_len_aa]. B = supported Prometheus placement, or required core + defense
  context, with any other shortfall. C = every other retained locus. `reason` lists the
  basis and every shortfall.

## prometheus/detection_sensitivity.tsv  (validate_prometheus.py)
genome_id, detect_defensefinder, prometheus_present {True,False}

## prometheus/prometheus_proteins.faa
Validated (tier A/B/C, not Rejected) proteins; header = locus_id.

## diversity/subfamilies.tsv  (cluster_diversity.py)
protein_id, genome_id, then one column per diversity identity, always named
`cluster_<identity %.2f>` (cluster_0.30, cluster_0.50, cluster_0.70, cluster_0.90 by
default). MMseqs2 failures or a missing binary are errors; a single sequence is its own
cluster. diversity/subfamily_summary.tsv: identity, n_proteins, n_clusters, n_singletons,
largest_cluster. diversity/subfamily_rarefaction.tsv (plot_diversity.py): n_genomes,
expected_subfamilies, perm_lo95, perm_hi95 -- exact incidence-based rarefaction over all
master-table genomes (Chao2 annotated on the figure).

## diversity/architecture_typology.tsv  (architecture_typology.py)
architecture, n_loci, example_locus, has_DUF4011, has_DUF3320, has_DNA2, has_nuclease,
has_REase, then (additive) has_<other spec categories>, core_complete, canonical_order,
rearranged, variant_class {canonical_full, canonical_core, rearranged, core_loss, none},
missing_required, n_with_fusion, fusion_domains, n_tier_A, n_tier_B, n_tier_C.
has_nuclease is kept for backward compatibility and is always False (there is no nuclease
category, Contract item 3). Architecture 'none' = loci without any spec domain. A
catalog/domtbl id mismatch is an error.

## scaffolds/metabolism/mcr_hits.tsv  (mcr_search.py)
genome_id, gene {mcrA,mcrB,mcrG}, protein_id, score, evalue, then (additive)
best_domain_score, criterion_score, threshold, score_type, threshold_source {ko_list,GA},
model, call {full_threshold, partial_rescue}, partial, n_copies. A protein is called when
its KOfam score reaches methanogen.kofam_score_fraction (1.0) x the curated KOfam
threshold (ko_list); a Prodigal-partial gene (genes.tsv) can be rescued at
methanogen.partial_score_fraction (0.5) x threshold (call=partial_rescue). At most one
subunit per protein. A subunit without a profile or curated threshold, or a failed
hmmsearch, is an error. scaffolds/metabolism/mcr_hits_all.tsv lists every hit (audit).

## scaffolds/metabolism/methanogen_calls.tsv  (classify_methanogen.py)
One row per genome of protein_map.tsv: genome_id, mcrA {True,False}, n_mcr, methanogen
{True,False}, is_anme {True,False}, methanogenesis_subtype {unresolved, NA}, evidence, then
(additive) mcr_group {methanogen, ANME, other, ambiguous, unplaced, not_placed, NA},
mcrA_protein, mcrA_call, n_mcrA_copies, nearest_ref, nearest_ref_group.
With methanogen.anme_via_mcrA_tree, each McrA is placed among labelled reference McrA
(scaffolds/metabolism/mcrA_placement.tsv). methanogen = mcrA called AND placed with the
methanogen references; ANME, 'other' (e.g. alkane-oxidising Mcr), ambiguous and unplaced
McrA give methanogen False. The subtype is not inferred here (unresolved for methanogens); with
methanogenesis.enabled it is resolved from pathway genes by methanogenesis_pathways.py.
How each analysis uses ANME: distribution and ancestral code them as non-methanogens (the
count is reported in the note); the disorder contrast excludes them.

## scaffolds/metabolism/methanogenesis_ko_hits.tsv  (methanogenesis_search.py)
Only when methanogenesis.enabled. One row per accepted protein of a methanogen genome:
genome_id, protein_id, ko, score, best_domain_score, criterion_score, threshold, score_type,
threshold_source {ko_list,GA}, call {full_threshold, partial_rescue, readthrough_TAG,
readthrough_TGA, near_threshold}, partial, contig_index, gene_index, start, end, strand, then
(additive) readthrough_length_aa, merged_proteins. Acceptance is the KofamScan rule (kofam.py,
identical to mcr_search). A gene Prodigal truncated at an in-frame UAG (pyrrolysine) or UGA
(selenocysteine) is re-translated through the codon and must pass the FULL threshold
(rules/RULES.md 3). call=near_threshold rows are EVIDENCE ONLY: never presence, but they
withdraw a confident absence. At most one KO per protein. A rules KO without a profile or a
curated threshold is an error. methanogenesis_ko_hits_all.tsv lists every hit (audit);
methanogenesis_readthrough.faa holds every attempted read-through translation.

## scaffolds/metabolism/methanogenesis_pathways.tsv  (methanogenesis_pathways.py)
Only when methanogenesis.enabled. One row per genome of protein_map.tsv: genome_id, methanogen,
completeness, then for each pathway P in co2_reducing, acetoclastic, methylotrophic,
methyl_dismutating, methyl_reducing: P {present, absent, unresolved, not_applicable} with
P_p_lo, P_p_hi (bounds on the probability that the pathway is truly encoded), i.e.
co2_reducing, co2_reducing_p_lo, co2_reducing_p_hi, acetoclastic, acetoclastic_p_lo,
acetoclastic_p_hi, methylotrophic, methylotrophic_p_lo, methylotrophic_p_hi,
methyl_dismutating, methyl_dismutating_p_lo, methyl_dismutating_p_hi, methyl_reducing,
methyl_reducing_p_lo, methyl_reducing_p_hi; the two-sided
flags sub_h2, sub_formate, sub_methanol, sub_methylamine (same four values); the one-sided
flags sub_methylthiol, sub_quatamine, sub_methoxyl {present, not_detected, not_applicable};
methanogenesis_subtype (composite label, unresolved, none_detected or NA), subtype_status
{resolved, partial, unresolved, not_applicable}, components, evidence, taxonomy_expected,
taxonomy_concordance {concordant, discordant, no_expectation, NA}, taxonomy_discordant_terms.
not_applicable = not a methanogen (ANME included). An absence needs completeness: a pathway is
absent only when the upper bound on a false absence is <= methanogenesis.absence_alpha.
methanogenesis_reconciliation.tsv lists every gene-call / taxon-expectation disagreement for
manual review (design 3.3).

## scaffolds/master_table.tsv  (build_master_table.py)  — the analysis backbone
genome_id, representative {True,False}, tree_tip (genome_id or ""),
gtdb_domain, gtdb_phylum, gtdb_class, gtdb_order, gtdb_family, gtdb_genus, gtdb_species,
completeness, contamination, quality_tier, pass_primary, genome_size,
mcrA, n_mcr, methanogen, is_anme, methanogenesis_subtype,
prometheus_present, prometheus_tier, prometheus_n_loci,
detect_defensefinder, n_defense_systems, n_defense_systems_excl_prometheus,
group_* and niche_* (species-level {1, 0, NA}; group_* sorted then niche_*), then
(additive, after the ecology block) prometheus_status_primary, prometheus_status_sensitivity
{present, absent, excluded} (Contract item 6 for stats.confidence_tiers_primary and the
A+B sensitivity set), n_contigs, n50, gc, tree_status, mcr_group, and -- only with
methanogenesis.enabled -- mgen_co2_reducing, mgen_acetoclastic, mgen_methylotrophic,
mgen_methyl_dismutating, mgen_methyl_reducing, mgen_sub_h2, mgen_sub_formate, mgen_sub_methanol,
mgen_sub_methylamine {1, 0, NA}, mgen_sub_methylthiol, mgen_sub_quatamine, mgen_sub_methoxyl
{1, NA}, subtype_status, subtype_taxonomy_concordance. NA covers unresolved and not-a-methanogen;
a one-sided flag is never 0. With the pathways table, methanogenesis_subtype carries the
composite pathway label instead of 'unresolved'.
- prometheus_tier = best non-rejected tier (A > B > C; "" if none); prometheus_present =
  any non-rejected locus (ANY tier -- analyses use the tier rule, i.e. the status columns
  or prometheus_tier, never prometheus_present alone); prometheus_n_loci = non-rejected
  loci; an invalid tier is an error.
- n_defense_systems excludes AntiDefenseFinder rows (activity Antidefense);
  n_defense_systems_excl_prometheus also excludes type == prometheus.defensefinder_system
  (exact, case-insensitive). A missing DefenseFinder table is an error, so defense counts
  are always integers; completeness/genome_size are empty only when unknown.
- Ecology: a blank group_*/niche_* cell stays NA (never imputed 0); a universe genome
  absent from the species table is an error. scaffolds/ecology_missingness.tsv: column,
  n_genomes, n_1, n_0, n_NA.

## scaffolds/species_tree/species_tree.treefile  (species_tree + root_tree.py)
Rooted with a bifurcating root (ape::is.rooted TRUE); tips are genome_ids (IQ-TREE name
rewrites reverted; ids needing Newick quoting are refused); SH-aLRT/UFBoot labels kept on
their bipartitions and written unquoted; no `[&R]`, no root branch length. A rooting
failure is an error unless tree.outgroup_on_failure: midpoint (recorded). tree.scope:
domain infers the tree with the GTDB references, roots it, then prunes it to the analysed
genomes (species_tree.full.treefile keeps the rooted tree before pruning).
rooting_method.txt (line 1 `method=...; tips=N`, then key=value lines),
species_tree_membership.tsv (genome_id, in_msa, in_tree, status) and
scaffolds/gtdbtk_align/gtdbtk_dropped.tsv (genome_id<TAB>GTDB-Tk reason, no header)
account for every genome that is not a tip.

## scaffolds/ecology_harmonized.tsv  (build_master_table.py)
genome_id, group_* and niche_* (species-level {1, 0, NA}; blank = NA, never 0): the
harmonized ecology subset of the master table.

## distribution/defense_matrix.tsv  (build_defense_matrix.py)
genome_id, one 0/1 column per DefenseFinder system type (includes Prometheus).
AntiDefenseFinder rows (activity Antidefense) are excluded. The archaeal rule passes
--require-all (a missing per-genome table is an error); unreadable or 0-byte tables are
errors.

## distribution/prevalence_by_taxon.tsv  (compute_prevalence.py)
rank, taxon, n_genomes, n_present, prevalence, ci_low, ci_high, median_completeness,
tier_set, then (additive) n_excluded, n_genomes_total, tier_set_sensitivity,
n_genomes_sensitivity, n_present_sensitivity, n_excluded_sensitivity,
prevalence_sensitivity, ci_low_sensitivity, ci_high_sensitivity. The first row is
rank=overall, taxon=all. n_genomes = tier-rule denominator (present + absent among
pass_primary genomes, Contract item 6); n_genomes_total = all pass_primary genomes of the
taxon (the old meaning). Wilson 95% CI (z = 1.95996); prevalence/CI are NA when
n_genomes = 0; median_completeness is over the denominator genomes.

## distribution/association_tests.tsv  (association_tests.R)
test, model, term, estimate, std_error, statistic, p_value, ci_low, ci_high, n, note.
Six rows are always written (NA plus the reason in `note` when not estimable):
naive_glm (logistic GLM; likelihood-ratio p; estimate NA under separation), naive_fisher
(log conditional-MLE odds ratio, exact CI and p), phyloglm_methanogen,
phyloglm_confounders (the headline model: covariates z(completeness), z(log genome_size),
z(log1p(n_defense_systems_excl_prometheus)); n = complete-covariate subset), fitPagel
(correlated evolution of presence and methanogen status; LR test) and D_statistic
(Fritz & Purvis D; p = (#{D_random <= D_obs} + 1)/(N + 1), a test of phylogenetic signal
against a random pattern; stats.phylo_d_permutations).
Response = presence under stats.confidence_tiers_primary (Contract item 6; excluded
genomes counted in `note`); methanogen from methanogen_calls (ANME coded non-methanogen,
n_anme in `note`). phyloglm rows are multi-start logistic_MPLE fits (best penalised
likelihood over four starts incl. the logistic_IG10 solution, warm-restarted); ci_* are
parametric-bootstrap percentile CIs (stats.phyloglm_boot; Wald when 0). Errors exit
non-zero; an unrooted tree is an error (no midpoint fallback).

## distribution/cooccurrence.tsv  (cooccurrence.R)
system_a (Prometheus), system_b, n_a, n_b, n_both, method, statistic, p_value, p_adj,
direction, then (additive) n, odds_ratio, or_ci_low, or_ci_high, fisher_p,
logL_independent, logL_dependent, phyloglm_logodds, phyloglm_se, phyloglm_p, note.
Prometheus presence = prometheus_tier under --primary-tiers (Contract item 6) over
pass_primary genomes on the species tree. A system is tested only with >=
stats.cooccurrence_min_count genomes with and without it. method is always fitPagel;
p_value is the Pagel LR test (the registry test); odds_ratio/fisher_p are the naive
tip-level Fisher description; phyloglm_* is a supplementary multi-start MPLE fit
(Prometheus ~ system + z(log genome size)). Perfect co-occurrence is labelled
'co-occurrence (OR=Inf)'. Rows are sorted by system_b.

## evolution/align/family.{aln,trim.aln,align_method.txt}  (family_align)
family.aln is the UNTRIMMED alignment of the validated proteins (terminal `*` stripped,
internal `*` -> X): MAFFT L-INS-i up to evolution.linsi_max_seqs (200), FAMSA above that or
when align_tool=famsa; family.trim.aln = ClipKIT (evolution.clipkit_mode) for the tree;
family.align_method.txt records aligner, n, thresholds and trimmer. Unknown
align_tool/trim_tool values fail. 0 sequences -> empty outputs (nonfatal); 1 is copied.

## evolution/tree/family.{treefile,ufboot}  (family_tree)
The IQ-TREE ML gene tree (evolution.gene_tree_model) with SH-aLRT/UFBoot labels, and the
UFBoot replicate trees. These ARE IQ-TREE's own --prefix outputs (no copy; the pristine
self-copy replaced both with '();' on every run). The rule fails if IQ-TREE fails or the
tip count differs from the number of aligned sequences; '();' only for < 4 sequences.

## evolution/domains/rf_distances.tsv  (domain_trees.py)
domain, n_tips, rf_to_family, normalized_rf, note, then (additive) n_seqs,
n_splits_domain, n_supported_conflicts (domain splits with UFBoot >= 95 incompatible with
family splits with UFBoot >= 95), au_pvalue_family_topology (AU test of the family topology
on the domain alignment), model, n_rearranged. Domains are cut by exact category matching
on envelope coordinates (utils.domain_intervals); a domain tree needs >=
evolution.domain_min_taxa sequences. Requested by `all` when evolution.per_domain_trees.

## evolution/align/family.codon.fasta and family.codon_report.tsv  (codon_align.py)
PAL2NAL codon alignment (NCBI table evolution.codon_table, default 11; gaps KEPT, so codon
column k = protein alignment column k). Every pep/CDS pair is validated first; the report
(id, status {kept, no_cds, cds_not_codon_multiple, length_mismatch, translation_mismatch,
internal_stop, protein_mismatch}, in_proteins, n_aa, n_codons, note) accounts for every
record. More than evolution.max_cds_drop_frac of the expected loci unusable, or any PAL2NAL
error or output mismatch, is an error.

## evolution/reconciliation/reconciliation_events.tsv  (run_reconciliation.py)
tool, event {Duplication,Transfer,Loss,Speciation,Origination}, count, fraction, note.
GeneRax (evolution.reconciliation.generax_model UndatedDTL, generax_strategy, subst_model,
seed) on the rooted, binary species tree (internal nodes relabelled sp_n<k>; support labels
never become node names). Event accounting from `<family>_eventCounts.txt`: Speciation =
S + SL, Transfer = T + TL, Loss = L + SL + TL, Duplication = D; `fraction` on the
Speciation and Transfer rows = vertical fraction (S+SL)/(S+SL+T+TL). Genes whose genome is
not a species-tree tip are pruned and counted in `note`. GeneRax/ALE failures, discarded
families and empty scenarios are errors. ALE (run_ale) is optional and exploratory.
evolution/reconciliation/reconciliation_transfers.tsv: tool, donor, recipient, count,
recipient_methanogen_clade, recipient_n_genomes.

## evolution/ancestral/gains_losses.tsv  (ancestral_states.R)
metric {n_gains, n_losses, rate_gain, rate_loss, n_independent_acquisitions,
gains_on_methanogen_branches, expected_gains_methanogen, enrichment_p, n_tips_analysed,
n_tips_excluded_tier, n_tips_present}, value, note, then (additive) tier_set. One block per
tier set (stats.confidence_tiers_primary and _sensitivity; Contract item 6 presence).
phytools make.simmap (evolution.ancestral.nsim maps; Mk model ER or ARD chosen by AIC;
Q = q_method; root prior = root_prior, FitzJohn by default) on the ROOTED species tree
(an unrooted tree is an error). n_independent_acquisitions = distinct seeding gain events;
methanogen branches = mcrA-defined methanogens (ANME and 'other' mcr are not methanogens).
enrichment_p comes from a prevalence-matched structured null (evolution.ancestral.nnull
accepted simulations within prevalence_tol of the observed prevalence);
gains_losses_null.tsv holds the null draws (with tier_set).

## evolution/selection/selection_summary.tsv  (run_selection.py)
method {FEL,MEME,BUSTED,aBSREL,RELAX,GARD}, statistic, value, p_value, n_sites_or_branches,
note; rows include FEL/sites_negative, GARD/breakpoints (evolution.selection.gard) and
'not run: <reason>' rows. HyPhy runs seeded (ENV=SetParameter(RANDOM_SEED, seed)) with a
CPU limit on the codon alignment (intact loci when evolution.selection.intact_only; no
silent fallback to all loci); the gene tree is always pruned to the codon records (the
note gives the number of pruned tips). RELAX (relax_partition methanogen,
relax_labelling clade): Test = branches whose descendants are all methanogens, Reference
= all-non-methanogen branches, ANME/ambiguous/unplaced hosts Unclassified; statistic K is
the confirmatory row. A missing or failing hyphy, or bad JSON, is an error.

## evolution/selection/site_table.tsv  (run_selection.py)
site, alignment_col, domain, beta_alpha_or_dnds, p_value, classification
{positive,negative,neutral}, then (additive) fel_alpha, fel_beta, meme_p_value,
meme_classification, n_seqs_with_codon, domain_support. beta_alpha_or_dnds may be 'inf'
(alpha = 0, beta > 0). domain = the per-column category from every sequence's exact domain
envelopes (not one reference sequence), with its support. MEME beta+ is parsed from
HyPhy's HTML column headers.

## ecology/viral/viral_summary.tsv  (aggregate_viral.py)
genome_id, n_proviruses, total_provirus_kb, mean_checkv_completeness, n_rloop_proviruses,
frac_rloop_proviruses, then (additive) n_viral_contigs, total_viral_contig_kb,
n_viral_sequences, n_rloop_viral_sequences, frac_rloop_viral_sequences,
n_proviruses_checkv_mq, mean_provirus_kb, mean_provirus_gc.
n_proviruses / total_provirus_kb count geNomad topology == 'Provirus' only (integrated
proviruses; whole viral contigs are reported separately); mean_checkv_completeness is over
proviruses; frac_rloop_proviruses is EMPTY (undefined) when n_proviruses == 0. Missing,
empty, malformed or mutually inconsistent per-genome geNomad/CheckV/R-loop files are
errors (never zeros).

## ecology/viral/{genome}/rloop.tsv and rloop_rlfs_calls.tsv  (rloop_predict.py)
rloop_predict.py is a call-for-call re-implementation of QmRLFS-finder v2.0 (models m1/m2,
both strands; ecology.rloop). rloop.tsv: genome_id, provirus_id, length, n_rlfs (QmRLFS
call count, all models, both strands), rloop_prone (merged RLFS loci >=
ecology.rloop.min_rlfs_score), n_rlfs_loci, n_rlfs_plus, n_rlfs_minus, n_rlfs_m1,
n_rlfs_m2, gc, rlfs_loci_per_kb. rloop_rlfs_calls.tsv: genome_id, provirus_id, model,
location, start_RIZ, end_RIZ, length_RIZ, G_RIZ, 3Gs_RIZ, 4Gs_RIZ, perc_G_RIZ, Linker,
start_REZ, end_REZ, length_REZ, G_REZ, 3Gs_REZ, 4Gs_REZ, perc_G_REZ, strand (the reference
table columns without the sequences; '-' calls in plus-strand coordinates). R-loop
proneness is a sequence-composition proxy (transcription direction is not modelled).

## ecology/ecology_associations.tsv  (ecology_stats.R)
predictor, level, model, estimate, std_error, p_value, p_adj, n, note, then (additive)
ci_low, ci_high (Wald 95%), n_present, n_absent, n_level1, n_excluded_tier, alpha,
converged, sens_model, sens_estimate, sens_std_error, sens_p_value, sens_n.
One phyloglm per predictor: binary species-level group_/niche_ columns (level "1";
ecology.predictors, default auto) and, with --viral, the continuous viral predictors
(level "continuous"; n_proviruses as z(log1p(count)), frac_rloop_proviruses as z(fraction):
log-odds per SD); level "" = a requested predictor absent from the master. Covariates
ecology.covariates (z-scored) are in every model, e.g. 'phyloglm(present ~ niche_x +
z(completeness) + z(genome_size), logistic_MPLE)'; '... logistic_IG10 fallback)' flags a
fit taken from the IG10 fallback (weak phylogenetic signal). Response = Contract item 6
under the primary tiers; NA written as 'NA'; p_adj = BH across reported terms. sens_* =
methanogen-adjusted partial effect (not in the BH family). Non-converged fits are not
reported as valid. ecology/ecology_associations_sensitivity.tsv = the same models with
stats.confidence_tiers_sensitivity (A+B). An unrooted tree is an error.

## ecology/methanogenesis_associations.tsv  (methanogenesis_stats.R)
Only when methanogenesis.enabled (design 6.2: within methanogens, does methanogenesis TYPE
predict Prometheus presence, with type and habitat fitted jointly?).
predictor, analysis, role, estimate, std_error, ci_low, ci_high, p_value, p_adj, n, n_present,
n_absent, n_level1, n_changes_mean, n_changes_lo, n_changes_hi, replicated, delta_aic,
best_model, covariates, vif, model, converged, note, then (additive) tier_set,
n_excluded_tier, n_unresolved, n_missing_covariate, n_events, n_terms, epv, habitat_terms,
n_habitat_eligible, alpha, mk_model, n_changes_parsimony, n_gains_mean, n_losses_mean and the
corHMM aic_*/aicw_* columns.
Analysis set: master rows on the rooted species tree with methanogen True, not ANME and
pass_primary; response = Contract item 6 under the primary tiers; rows with an unresolved
predictor ('NA') are dropped and counted (n_unresolved). Eight rows per predictor:
`marginal`/secondary and `habitat_adjusted`/**primary** phyloglm (habitat terms chosen by
|phi(x, habitat)| under the events-per-variable budget floor(n_events / epv_min), Peduzzi et
al. 1996; VIF(x) > vif_max -> the row is reported unfitted, p_value NA), `simmap`/descriptive
(the replication guard: Mk ER vs ARD by AIC then make.simmap; replicated = n_changes_mean >=
min_changes -- when FALSE every p-value of that predictor is withheld and the row leaves the
FDR family, per Maddison & FitzJohn 2015), `corhmm`/descriptive (fitCorrelationTest,
delta_aic = min independent AIC - min correlated AIC; no p-value) and four
`sens_*`/sensitivity re-fits (unresolved-as-absent, completeness >= complete_threshold,
leave-out the gtdb_order with most x = 1 tips, A+B tiers). p_adj = BH over the PRIMARY rows
with a p-value only. The registry family `methanogenesis` is always exploratory.

## --- Bacterial / cross-domain branch (opt-in; bacteria.enabled: true) ---------

## bacteria/prometheus/bacterial_prometheus_catalog.tsv  (validate_prometheus.py, reused)
Same columns and tier rules as prometheus_catalog.tsv, for bacterial loci from the provided
DefenseFinder tables. extract_bacterial_prometheus.py splits BOTH combined tables
(systems and genes) into bacteria/defense/{g}/{g}_defense_finder_{systems,genes}.tsv for
every universe genome (header-only when empty; AntiDefenseFinder rows dropped and
counted), matches bacteria.system exactly (type or subtype, case-insensitive), reads
gzip-aware, and fails on blank/duplicate/unsafe ids, missing proteomes (unless
bacteria.allow_missing_proteins) or zero overlap; bacteria/prometheus/extraction_report.tsv
reports every genome. The bacterial validation tree mirrors the archaeal one (MAFFT,
ClipKIT, IQ-TREE -m prometheus.validation_tree_model with SH-aLRT/UFBoot, --seed), so
bacterial loci can reach tier A. Fragmentation comes from Prodigal partial flags where the
proteome has Prodigal headers (bacteria/genes/{g}/{g}.genes.tsv), else the labelled
length proxy (genes_table=none in the extraction report).

## bacteria/bacterial_master_table.tsv  (build_bacterial_table.py)
genome_id, tree_tip, gtdb_domain..gtdb_species, completeness, contamination, pass_primary,
prometheus_present, prometheus_tier, prometheus_n_loci, n_defense_systems,
n_defense_systems_excl_prometheus, genome_size, then (additive) detect_defensefinder,
prometheus_status_primary, prometheus_status_sensitivity (Contract item 6).
Rows = bacteria/bacterial_samples.tsv (the analysis universe). pass_primary =
completeness >= bacteria.min_completeness AND contamination <= bacteria.max_contamination
(5.0), matching the archaeal filter; unknown QC does not pass. genome_size comes from the
QC table or is streamed (gzip-aware) from the assemblies. Duplicate ids and missing
per-genome defense tables are errors. Slim master so compute_prevalence.py and
cooccurrence.R are reused as-is.

## bacteria/distribution/bacterial_prevalence_by_taxon.tsv  (compute_prevalence.py, reused)
Same columns as the archaeal prevalence table (tier rule and A+B sensitivity columns
included), for bacteria. Report side by side; matched sampling is the user's
responsibility (bacterial GTDB dwarfs archaeal).

## crossdomain/crossdomain_summary.tsv  (crossdomain_analysis.py)
metric, value, note. Topology-based (Fitch parsimony on the outgroup-rooted ML gene tree;
no DTL model; not a significance test). The tree is rooted on the edge separating every
Paralog_ref tip from the ingroup (rooting, outgroup_monophyletic, outgroup_split_support;
a non-monophyletic outgroup means no direction is called). Metrics: status,
n_archaeal_tips, n_bacterial_tips, n_reference_tips, n_crossdomain_transitions (Fitch
minimum A<->B changes) and n_crossdomain_transitions_supported (after collapsing branches
below bacteria.crossdomain.min_alrt / min_ufboot), n_archaeal_clades,
n_archaeal_clades_supported, n_archaeal_clades_nested_supported, archaeal_monophyletic,
ingroup_root_state, collapsed_root_state, n_bootstrap_trees, bootstrap_root_B/A_fraction,
n_subsample, subsample_root_B/A_fraction (the larger domain thinned to the smaller,
seeded), midpoint_root_state, root_state_excluding_references, inferred_direction_ml
(ML-tree root state alone) and inferred_direction: bacteria->archaea or archaea->bacteria
only when ALL hold -- clean outgroup split, the ML and support-collapsed root states agree,
and the state holds in >= bacteria.crossdomain.robust_fraction of the UFBoot replicate
trees and of the balanced subsamples -- otherwise 'ambiguous' with the failed criteria;
direction_method; cross_domain_dS_saturated_fraction and dnds_saturation_warning. Also
written: the tip table (adds archaeal_clade), crossdomain_archaeal_clades.tsv and
crossdomain/tree/crossdomain_rooted.treefile.

## crossdomain/selection/crossdomain_selection_summary.tsv  (run_selection.py, reused)
Same columns as the archaeal selection summary; FEL/MEME/BUSTED/aBSREL (no methanogen
RELAX partition across domains) on the cross-domain codon alignment of the validated
archaeal + bacterial loci, with the HyPhy seed, on the gene tree pruned to exactly the
codon-aligned taxa (crossdomain_selection.nwk). Not run when bacteria.cds_dir is unset
(note 'not run: bacteria.cds_dir not set ...'). EXPLORATORY: read with
crossdomain/selection/crossdomain_saturation.tsv (pair_class, n_pairs, n_saturated,
frac_saturated, median_dS, median_dN, median_codons, note; NG86; a pair is saturated if dS
is undefined or > bacteria.crossdomain.saturation_ds) and dnds_saturation_warning (> 50% of
cross-domain pairs saturated). Bacteria with a non-11 genetic code are dropped by
codon_align as translation mismatches.

## bacteria/distribution/bacterial_cooccurrence.tsv  (cooccurrence.R, reused)
Same columns and rules as the archaeal co-occurrence table (same --primary-tiers,
--min-count, --btol), on the bacterial defense matrix (--require-all; Antidefense excluded)
and the GTDB bac120 tree pruned to the bacterial genomes by prune_tree_to_taxa.py: tiered
tip matching (exact, prefix-, version-stripped, paired GCA/GCF accession, species
representative via bacteria.gtdb_bac_taxonomy), one genome per tip, tips relabelled to
genome_id, an error below bacteria.min_tree_match_fraction;
bacteria/scaffolds/bacterial_species_tree_tips.tsv reports every genome.

## bacteria/ecology/viral/bacterial_viral_summary.tsv  (aggregate_viral.py, reused)
Same columns and semantics as the archaeal viral_summary (integrated proviruses only;
frac_rloop_proviruses empty when n_proviruses == 0), from geNomad + CheckV +
rloop_predict.py on the provided bacterial assemblies (bacteria.genomes_dir), run with the
same safeguards as the archaeal rules (geNomad --restart and a version.txt database check;
CheckV --restart --remove_tmp, database and output-header checks).

## bacteria/ecology/bacterial_viral_associations.tsv  (ecology_stats.R, reused)
Same columns as ecology_associations.tsv: bacterial Prometheus presence (Contract item 6)
vs the continuous viral predictors, phylogeny-controlled on the pruned bacterial species
tree and adjusted for bacteria.covariates (default ecology.covariates).
bacterial_niche_associations.tsv does the same for the niche labels of
bacteria.niche_metadata (bacteria/ecology/bacterial_niche_master.tsv holds annotated
genomes only; a missing label is NA, never 0). Each has an A+B *_sensitivity.tsv. The
viral_rloop_crossdomain figure visualises the mechanism across both domains.

## --- Intrinsic-disorder module (opt-in; disorder.enabled: true) ---------------

## disorder/disorder_per_locus.tsv  (predict_disorder.py)
locus_id, genome_id, length_aa (a trailing `*` excluded), intact ('' when no catalog is
given), predictors {metapredict+iupred3 | metapredict+alphapredict | metapredict_only},
frac_disordered (consensus), frac_disordered_metapredict, frac_disordered_second, n_idrs,
longest_idr, frac_in_idrs, frac_lcr, frac_disordered_nonlcr, disorder_minus_background,
then (additive) confidence_tier, n_nonstandard, host_background, predictor_detail.
metapredict 3.1.1 (disorder network V3); the second predictor is IUPred3 or metapredict's
pLDDT network (V1 = the alphaPredict network; no separate alphaPredict package).
frac_disordered_nonlcr = disordered non-LCR residues / non-LCR residues ('' if every
residue is LCR) -- the preferred comparative metric. Intact, non-fragmented loci only by
default; any tool failure is an error.

## disorder/disorder_per_residue.tsv  (predict_disorder.py)
locus_id, position (1-based), meta_score, second_score, consensus {0,1}, lcr {0,1}, then
(additive) aa (sanitised residue). second_score is the IUPred3 score (0-1) or the predicted
pLDDT (0-100), per the locus' predictors value.

## disorder/disorder_by_column.tsv  (map_disorder_to_alignment.py)
alignment_col, n_loci, mean_disorder, frac_loci_disordered, domain, then (additive)
occupancy, domain_support, mean_meta_score. domain = the per-column plurality spec
category across the annotated aligned loci (support >= disorder.min_domain_support, 0.5;
ties NA; NA = linker/tail). mean_disorder equals frac_loci_disordered by construction.
Alignment/disorder disagreements are errors.

## disorder/disorder_associations.tsv  (disorder_stats.R)
predictor, level, model, estimate, std_error, p_value, p_adj, n, note, then (additive)
response, tier_set ('A' / 'A,B'), n_level1, n_excluded_tier, lambda, ci_low, ci_high.
phylolm (Pagel lambda) of per-genome consensus disorder on methanogen status (ANME
excluded from the contrast) and each niche, with the configured temperature columns (and
host_background when present) as covariates; predictor 'phylogenetic_signal' (level
'lambda') reports the signal in disorder itself. One block of rows per tier set (primary +
A+B sensitivity); model strings carry the covariates and '[tiers X]'. EXPLORATORY
(sequence-predicted): enters results_registry under family "disorder", never confirmatory.

## synthesis/results_registry.tsv  (results_registry.py)
family, test, term, estimate, ci_low, ci_high, p_value, p_adj, hypothesis
{confirmatory,exploratory}, then (additive) role {primary, secondary, baseline,
descriptive, sensitivity}, estimate_scale {log_odds, odds_ratio, ratio, count,
coefficient, D_statistic, statistic, none}, null_value, ci_method {producer, wald_z,
none}, statistic, std_error, n, fdr_method, fdr_m, tier_set, model, note.
Families: association, ancestral, reconciliation, selection, cooccurrence, ecology, and
disorder when enabled (always exploratory). Multiple testing (stats.fdr_method: BH, BY,
holm or bonferroni) is applied within each family over rows with a p-value; the naive
baselines (naive_glm, naive_fisher) and non-primary tier_set rows (sensitivity) are
reported with an empty p_adj, outside the family. Confirmatory = exactly one pre-registered
primary test per hypothesis (stats.confirmatory overrides; default):
association:phyloglm_confounders:methanogen;
ancestral:methanogen_gain_enrichment:methanogen_independent_acquisition_enrichment;
selection:RELAX:K; reconciliation:generax:transfer. Everything else is exploratory
(including phyloglm_methanogen and ALE transfer). The rule runs with
--strict-confirmatory: a non-empty producer table whose confirmatory label is missing
fails the rule. Missing, empty or unparseable inputs are errors. The bacterial and
cross-domain tables are not registered (opt-in, exploratory; each carries its own p_adj).
