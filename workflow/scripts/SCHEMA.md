# Table schemas (the I/O contract every script adheres to)

IDs: combined protein/CDS ids are `genome__origProteinId` (see utils.ID_SEP).
Species-tree tips = `genome_id`. Gene-tree tips = combined protein id. Recover the
genome from a protein id with `utils.genome_of()`. All tables are TSV with a header.
Every script must tolerate empty inputs and still write a header-only output.

## preprocess/species_metadata.tsv  (filter_collapse.py) — the SINGLE filtered table (samples + metadata)
genome_id, filename, source, classification, gtdb_domain..gtdb_species, assembly_path,
completeness, contamination, gc, genome_size, n_contigs, n50, quality_score,
n_genomes_in_species, n_passing_in_species, group_* and niche_* (0/1, unioned per
species). This one table is both `config.samples` and `config.metadata`. Only
genomes GTDB-Tk places in d__Archaea are retained (config preprocess.keep_domain).
`assembly_path` is the staged genome each analysis step (Prodigal first) reads.
BBMap columns used: gc=gc_avg, genome_size=contig_bp, n_contigs=n_contigs, n50=ctg_L50.

## genes/{genome}/{genome}.{faa,fna}  (prodigal_genes + rename_prodigal.py)
Prodigal proteins (.faa) and CDS nucleotides (.fna). Headers are rewritten to
`>{genome}_{N}` (N from 1), IDENTICAL for the matching protein and CDS so the two
files stay consistent gene-by-gene. The .faa feeds DefenseFinder + downstream.

## qc/genome_metadata.tsv  (cp of preprocess/genome_metadata.tsv)
genome_id, source, completeness, contamination, genome_size, n_contigs, n50, gc,
quality_tier {A}, pass_primary {True}.  No MAG/isolate genome_type — quality is
represented by completeness / contamination / n_contigs / n50.

## prometheus/protein_map.tsv  (combine_proteomes.py)
protein_id, genome_id, orig_id

## prometheus/candidate_sources.tsv  (parse_defensefinder.py)
protein_id, genome_id, detect_defensefinder {0,1}, df_score, df_subtype
(Detection is DefenseFinder-only; detect_defensefinder is 1 for every candidate.)

## prometheus/prometheus_catalog.tsv  (validate_prometheus.py)
locus_id (=protein_id), genome_id, protein_id, length_aa, detect_defensefinder,
architecture (e.g. "DUF4011-DUF3320-DNA2-REase"), n_domains,
catalytic_status {intact,partial,degenerate,unknown},
defense_context {True,False}, tree_assignment {Prometheus,Paralog,Unplaced},
fragmented {True,False}, confidence_tier {A,B,C,Rejected}, reason
(Tier A = DefenseFinder call + full canonical architecture (>=4 domains) +
Prometheus placement + intact catalytic modules + full length.)

## prometheus/detection_sensitivity.tsv  (validate_prometheus.py)
genome_id, detect_defensefinder, prometheus_present {True,False}

## prometheus/prometheus_proteins.faa
Validated (tier A/B/C, not Rejected) proteins; header = locus_id.

## diversity/subfamilies.tsv  (cluster_diversity.py)
protein_id, genome_id, cluster_0.30, cluster_0.50, cluster_0.70, cluster_0.90

## diversity/architecture_typology.tsv  (architecture_typology.py)
architecture, n_loci, example_locus, has_DUF4011, has_DUF3320, has_DNA2, has_nuclease, has_REase

## scaffolds/metabolism/mcr_hits.tsv  (mcr_search.py)
genome_id, gene {mcrA,mcrB,mcrG}, protein_id, score, evalue

## scaffolds/metabolism/methanogen_calls.tsv  (classify_methanogen.py)
genome_id, mcrA {True,False}, n_mcr, methanogen {True,False}, is_anme {True,False},
methanogenesis_subtype {hydrogenotrophic,methylotrophic,acetoclastic,methyl_reducing,unresolved,NA},
evidence

## scaffolds/master_table.tsv  (build_master_table.py)  — the analysis backbone
genome_id, representative {True,False}, tree_tip (genome_id or ""),
gtdb_domain, gtdb_phylum, gtdb_class, gtdb_order, gtdb_family, gtdb_genus, gtdb_species,
completeness, contamination, quality_tier, pass_primary, genome_size,
mcrA, n_mcr, methanogen, is_anme, methanogenesis_subtype,
prometheus_present, prometheus_tier (best of A>B>C or ""), prometheus_n_loci,
detect_defensefinder, n_defense_systems,
group_* and niche_* (0/1, one column per group/niche; species-level
presence-absence from the metadata table, group_* sorted then niche_*; absent
genome -> 0)

## scaffolds/ecology_harmonized.tsv  (build_master_table.py)
genome_id, group_* and niche_* (0/1, one column per group/niche; the harmonized
ecology subset of the master table)

## distribution/defense_matrix.tsv  (build_defense_matrix.py)
genome_id, <one column per defense system type, 0/1> ... (includes Prometheus)

## distribution/prevalence_by_taxon.tsv  (compute_prevalence.py)
rank, taxon, n_genomes, n_present, prevalence, ci_low, ci_high, median_completeness, tier_set

## distribution/association_tests.tsv  (association_tests.R)
test, model, term, estimate, std_error, statistic, p_value, ci_low, ci_high, n, note
(rows: naive_glm, phyloglm_methanogen, phyloglm_methanogen+confounders, fitPagel, D_statistic)

## distribution/cooccurrence.tsv  (cooccurrence.R)
system_a (Prometheus), system_b, n_a, n_b, n_both, method, statistic, p_value, p_adj, direction

## evolution/domains/rf_distances.tsv  (domain_trees.py)
domain, n_tips, rf_to_family, normalized_rf, note

## evolution/reconciliation/reconciliation_events.tsv  (run_reconciliation.py)
tool, event {Duplication,Transfer,Loss,Speciation,Origination}, count, fraction, note

## evolution/ancestral/gains_losses.tsv  (ancestral_states.R)
metric {n_gains,n_losses,rate_gain,rate_loss,n_independent_acquisitions,gains_on_methanogen_branches,expected_gains_methanogen,enrichment_p}, value, note

## evolution/selection/selection_summary.tsv  (run_selection.py)
method {FEL,MEME,BUSTED,aBSREL,RELAX}, statistic, value, p_value, n_sites_or_branches, note

## evolution/selection/site_table.tsv  (run_selection.py)
site, alignment_col, domain, beta_alpha_or_dnds, p_value, classification {positive,negative,neutral}

## ecology/viral/viral_summary.tsv  (aggregate_viral.py)
genome_id, n_proviruses, total_provirus_kb, mean_checkv_completeness,
n_rloop_proviruses, frac_rloop_proviruses

## ecology/ecology_associations.tsv  (ecology_stats.R)
predictor, level, model, estimate, std_error, p_value, p_adj, n, note
(predictor is a binary species-level group_/niche_ column [level = "1"] or a
continuous viral predictor {n_proviruses, frac_rloop_proviruses}; one phyloglm
term per predictor; p_adj = BH across all reported terms.)

## synthesis/results_registry.tsv  (results_registry.py)
family, test, term, estimate, ci_low, ci_high, p_value, p_adj, hypothesis {confirmatory,exploratory}
