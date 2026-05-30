# Prometheus in archaea — Snakemake pipeline

A reproducible workflow for the diversity, distribution, evolution, taxonomy, and
ecology of the Prometheus anti-phage defense system across archaea, and its
association with methanogens. It implements
[`Prometheus_archaea_pipeline_design.md`](Prometheus_archaea_pipeline_design.md);
read that first for the scientific rationale. This README covers installation,
data setup, running, the decisions baked in, a per-step program reference, and
the outputs.

Prometheus is a single ~2,221-aa ORF of the bacterial DNA2-like (Bad) family with
four domains (DUF4011, DUF3320, DNA2 helicase/nuclease, REase MTES) that acts on
phage transcription R-loops (van den Berg et al., *Cell Host & Microbe* 2024). In
archaea it has only been detected, never tested; this pipeline characterizes it.

## Layout

```
config/            config.yaml, samples.tsv, metadata.tsv, envo_map.tsv
resources/         genomes/, db/, refs/ (you provide; see resources/README.md)
workflow/
  Snakefile        includes the 8 stage modules + provenance
  rules/*.smk      common + qc, prometheus, diversity, scaffolds,
                   distribution, evolution, ecology, synthesis
  envs/*.yaml      10 per-stage conda environments
  scripts/py/*.py  27 analysis + plotting scripts (matplotlib, Illustrator-ready)
  scripts/R/*.R    8 phylogenetic-stats + tree-figure scripts (ggtree/phytools)
  schemas/         config + samples JSON schemas (validated at start)
  profiles/sge/    Grid Engine submission profile + wrapper
results/           created on run (see "Outputs")
```

## 1. Install

You need one small controller environment with Snakemake; every rule then builds
its own pinned environment automatically via `--use-conda`. Use mamba (fast solver).

```bash
# controller env
mamba create -n prometheus-smk -c conda-forge -c bioconda \
    "snakemake-minimal>=7.32" "snakemake-executor-plugin-cluster-generic" mamba
conda activate prometheus-smk
```

Per-rule environments in `workflow/envs/` are materialized on first run by adding
`--use-conda`. You do not install them by hand. To build them ahead of time (e.g.
on a login node with internet, before going to compute nodes):

```bash
snakemake --use-conda --conda-create-envs-only --cores 4
```

A few tools need one-time post-install data fetches (the pipeline does not do these
for you, by design — pin the versions):

```bash
# in the relevant materialized env, or any env that has the tool:
defense-finder update                  # DefenseFinder models (the sole Prometheus detector)
checkm2 database --download --path resources/db/checkm2
genomad download-database resources/db
checkv download_database resources/db
# GTDB-Tk + KOfam + Pfam: download the release; see resources/README.md
```

## 2. Provide data

- Put assemblies at `resources/genomes/{genome_id}.fna` and list them in
  `config/samples.tsv` (genome_id, source, assembly_level, genome_type).
- Put ecological metadata in `config/metadata.tsv` (keyed by genome_id; blanks are
  fine — missingness is tracked, never imputed).
- Stage databases and reference sequences and point `config/config.yaml` at them.
  **`resources/README.md` lists every file, what it is, and how to get it.** The
  reference protein sets (`prometheus_reference.faa`, `paralog_reference.faa`, the
  mcr HMMs, and the labelled McrA reference alignment) are the
  defensibility-critical inputs — supply them carefully.

If you already have DefenseFinder outputs and prefer to ingest rather than re-run,
set `prometheus.run_defensefinder: false` and drop them under
`resources/defensefinder/{genome_id}/`. Re-running under a pinned version (the
default) is more reproducible.

## 3. Run

Always dry-run first to see the DAG:

```bash
snakemake -n -p                      # nothing executes; checks wiring
```

Local workstation:

```bash
snakemake --use-conda --cores 24 --resources mem_mb=120000
```

Grid Engine ("Sun"/SGE) cluster:

```bash
snakemake --workflow-profile workflow/profiles/sge --use-conda --jobs 200
```

The profile (`workflow/profiles/sge/`) translates each rule's `threads`,
`resources.mem_mb`, and `resources.runtime` into a `qsub` submission via
`sge-submit.sh`. Two site-specific things to check before launching:

- **Parallel environment name.** The wrapper requests `-pe smp {threads}`. Run
  `qconf -spl` and, if your PE isn't called `smp`, export `SGE_PE_NAME=<name>`.
- **Memory semantics.** The wrapper divides total memory by threads because most
  SGE installs enforce `h_vmem` *per slot*. If yours treats `h_vmem` per job, edit
  the one line in `sge-submit.sh` (it's commented).

Target a single stage or file by name, e.g.
`snakemake --use-conda --cores 8 results/scaffolds/master_table.tsv`. The design
doc's phased rollout (foundation → distribution → evolution/ecology) maps to
requesting those stage outputs in order.

## 4. Decisions baked in (resolved for defensibility)

Where the design left a choice, the pipeline takes the more defensible option and
records it here:

- **"Methanogen" is defined from genes, not taxonomy.** Presence of McrA
  (`mcrABG` HMMs), with anaerobic methanotrophs (ANME) separated by placing each
  McrA against a labelled reference. This keeps the trait independent of the GTDB
  taxonomy it's tested against (methanogens are polyphyletic).
- **Detection is DefenseFinder-only** (no PADLOC, no custom HMM sweep), then
  curated. DefenseFinder already applies its mandatory-gene architecture rule
  before calling Prometheus; nothing is counted before the validation below.
- **Validation separates Prometheus from housekeeping DNA2/NurA-HerA** by domain
  architecture, a catalytic-domain proxy, defense-island context, and phylogenetic
  placement against references. Candidates that place with cellular paralogs are
  rejected. Loci are tiered A/B/C; primary analyses use tier A, with A+B as a
  sensitivity check (`stats.confidence_tiers_primary`).
- **Absence is modeled, not assumed.** Genome completeness, a detection-sensitivity
  flag, genome size, and total defense-system count enter the association models as
  covariates. The headline methanogen estimate is the phylogenetic-logistic term
  *after* this conditioning.
- **Every association is phylogenetically controlled** (`phyloglm`, Pagel's
  correlated-evolution test, the D statistic; co-occurrence via `fitPagel`). The
  naive contingency-table result is reported only as a labelled baseline.
- **The species tree is rooted before any comparative or ancestral use**, at a
  position you set in `config` (`tree.outgroup` genome_ids, or `tree.outgroup_taxon`
  resolved via the taxonomy), with midpoint as a labelled fallback. This matters
  because ancestral gain/loss (and the independent-acquisition count), `phyloglm`,
  `fitPagel`, `phylo.d`, and GeneRax's UndatedDTL all assume a correct root. The
  rooting actually applied is written to `results/scaffolds/species_tree/rooting_method.txt`.
  Set an outgroup appropriate to your genome set; midpoint is a heuristic, not a
  defensible default for rate-variable archaeal trees.
- **Vertical inheritance is never assumed** — it's an outcome of GeneRax
  reconciliation (ALE optional). Gains/losses come from stochastic mapping with a
  methanogen-enrichment permutation test.
- **dN/dS is restricted to catalytically-intact loci** so pseudogene relaxation
  can't masquerade as signal; RELAX tests whether selection intensity differs on
  methanogen branches.
- **One FDR regime** (Benjamini-Hochberg, per test family) is applied centrally in
  the results registry, which also tags confirmatory vs exploratory tests.

Change any of these in `config/config.yaml`; the schema is validated at start.

## 5. Pipeline steps and the programs each runs

Each box in the DAG is one Snakemake rule. The tables list every rule, what it
does, the program(s) it invokes, and its main output. **Bold** names are external
programs (from the rule's conda env in `workflow/envs/`); plain `name.py`/`name.R`
are in-house steps under `workflow/scripts/`. Versions are pinned in the env files;
per-step parameters live in `config/config.yaml`. Per-genome rules fan out over
`config/samples.tsv`; aggregation rules gather.

### Stage 0 — genome set, QC, gene calls  (`rules/qc.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `prodigal_genes` | Predict ORFs/proteins in one coordinate space | **Prodigal** | `genes/{g}/{g}.{faa,ffn,gff}` |
| `checkm2` | Completeness & contamination (ML, archaea-aware) | **CheckM2** | `qc/checkm2/quality_report.tsv` |
| `seqkit_stats` | Assembly stats: size, N50, contigs, GC | **seqkit** `stats` | `qc/seqkit_stats.tsv` |
| `barrnap_rrna` | rRNA genes (MIMAG tiering) | **barrnap** (`--kingdom arc`) | `qc/rrna/{g}.rrna.gff` |
| `trnascan_trna` | tRNA genes (MIMAG tiering) | **tRNAscan-SE** (`-A`) | `qc/trna/{g}.trna.txt` |
| `genome_metadata` | Aggregate QC, assign A/B/C/fail tiers, dRep genomeInfo | `aggregate_qc.py` | `qc/genome_metadata.tsv`, `qc/genomeInfo.csv` |
| `dereplicate` | Species-level dereplication (optional) | **dRep** (fastANI) | `qc/representatives.txt` |
| `plot_qc` | Completeness–contamination + assembly-stat diagnostics | `plot_qc.py` (**matplotlib**) | `qc/plots/qc_*.{png,svg}` |

### Stage 1 — Prometheus identification & validation  (`rules/prometheus.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `defensefinder` | Detect anti-phage systems; Prometheus is taken from these calls (run pinned, or ingest) | **DefenseFinder** (MacSyFinder + **HMMER**) | `defense/defensefinder/{g}/{g}_defense_finder_{systems,genes}.tsv` |
| `combined_proteome` | Concatenate proteomes with newick-safe `genome__orig` IDs | `combine_proteomes.py` | `prometheus/all_proteins.faa`, `protein_map.tsv` |
| `assemble_candidates` | Pull the proteins of DefenseFinder's Prometheus systems | `parse_defensefinder.py` | `prometheus/candidates.faa`, `candidate_sources.tsv` |
| `prometheus_domain_scan` | Annotate domain architecture | **HMMER** `hmmscan` vs Pfam-A (`--cut_ga`) | `prometheus/domains.domtbl` |
| `prometheus_validation_tree` | Place candidates against Prometheus + DNA2/NurA-HerA references | **MAFFT** + **ClipKIT** + **IQ-TREE2** | `prometheus/validation_tree/validation.treefile` |
| `validate_prometheus` | Integrate architecture + placement + context → A/B/C/Rejected catalog | `validate_prometheus.py` (**dendropy**) | `prometheus/prometheus_catalog.tsv`, `prometheus_proteins.faa`, `detection_sensitivity.tsv` |
| `plot_prometheus` | Detection→validation funnel + architecture | `plot_prometheus_validation.py` | `prometheus/plots/*.{png,svg}` |

### Stage 2 — diversity  (`rules/diversity.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `mmseqs_cluster` | Cluster Prometheus proteins into subfamilies at graded identities | **MMseqs2** `easy-cluster` (via `cluster_diversity.py`) | `diversity/subfamilies.tsv` |
| `architecture_typology` | Enumerate distinct domain architectures | `architecture_typology.py` | `diversity/architecture_typology.tsv` |
| `plot_diversity` | Subfamily rarefaction + architecture-by-phylum | `plot_diversity.py` | `diversity/plots/*.{png,svg}` |

### Stage 3 — reference scaffolds  (`rules/scaffolds.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `gtdbtk_classify` | Standardized GTDB taxonomy | **GTDB-Tk** `classify_wf` (Prodigal, HMMER, pplacer, skani) | `scaffolds/taxonomy_gtdbtk.tsv` |
| `gtdbtk_align` | Concatenated ar53 marker MSA (focal or domain-wide) | **GTDB-Tk** `identify` + `align` | `scaffolds/gtdbtk_align/user_msa.fasta.gz` |
| `species_tree` | ML archaeal species tree, two-step PMSF for deep nodes | **IQ-TREE2** | `scaffolds/species_tree.treefile` |
| `mcr_search` | Detect `mcrA/B/G` (defines methanogens from genes) | **HMMER** `hmmsearch` vs mcr HMMs (via `mcr_search.py`) | `scaffolds/metabolism/mcr_hits.tsv` |
| `methanogen_classify` | Methanogen vs ANME by McrA phylogenetic placement | `classify_methanogen.py` (**MAFFT** `--add`, **IQ-TREE2**, **dendropy**) | `scaffolds/metabolism/methanogen_calls.tsv` |
| `build_master_table` | Join QC, taxonomy, tree tip, metabolism, ecology, Prometheus, defense load | `build_master_table.py` | `scaffolds/master_table.tsv`, `ecology_harmonized.tsv` |

### Stage 4 — distribution & taxonomic association  (`rules/distribution.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `defense_matrix` | Genome × defense-system presence/absence | `build_defense_matrix.py` | `distribution/defense_matrix.tsv` |
| `compute_prevalence` | Prevalence by GTDB rank with Wilson CIs | `compute_prevalence.py` | `distribution/prevalence_by_taxon.tsv` |
| `association_tests` | Methanogen↔Prometheus: naive + phylogenetic + correlated-evolution + signal | **R**: `phylolm` (phyloglm), `phytools` (fitPagel), `caper` (phylo.d), `ape` | `distribution/association_tests.tsv` |
| `cooccurrence` | Prometheus vs each other system (correlated evolution) | **R**: `phytools::fitPagel` (+ Fisher direction) | `distribution/cooccurrence.tsv` |
| `plot_distribution` | Prevalence bars, association forest, co-occurrence heat | `plot_distribution.py` | `distribution/plots/*.{png,svg}` |

### Stage 5 — molecular evolution  (`rules/evolution.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `family_align` | Align + trim Prometheus proteins | **MAFFT** (L-INS-i) + **ClipKIT** | `evolution/align/family.trim.aln` |
| `family_tree` | ML gene tree (ModelFinder, UFBoot2 + SH-aLRT) | **IQ-TREE2** | `evolution/tree/family.treefile`, `.ufboot` |
| `domain_trees` | Per-domain trees; Robinson–Foulds vs family (domain shuffling) | `domain_trees.py` (**MAFFT**, **ClipKIT**, **IQ-TREE2**, **dendropy**) | `evolution/domains/rf_distances.tsv` |
| `combined_cds` | Concatenate CDS with matching IDs | `combine_proteomes.py` (`--suffix ffn`) | `evolution/all_cds.ffn` |
| `codon_align` | In-frame codon alignment | **PAL2NAL** (via `codon_align.py`) | `evolution/align/family.codon.fasta` |
| `reconciliation` | Gene-tree/species-tree DTL reconciliation | **GeneRax** (ALE optional) via `run_reconciliation.py` | `evolution/reconciliation/reconciliation_events.tsv` |
| `ancestral_states` | Ancestral gain/loss + methanogen-enrichment permutation | **R**: `phytools::make.simmap` | `evolution/ancestral/gains_losses.tsv` (+ figure) |
| `selection` | Site/branch dN/dS on catalytically-intact loci | **HyPhy** FEL, MEME, BUSTED, aBSREL, RELAX (via `run_selection.py`) | `evolution/selection/selection_summary.tsv`, `site_table.tsv` |
| `plot_tanglegram` | Gene-tree vs species-tree tanglegram | **R**: `phytools::cophylo` | `evolution/plots/genetree_speciestree_tanglegram.{png,svg}` |
| `plot_selection` | Positively-selected sites per domain | `plot_selection.py` | `evolution/plots/selection_by_domain.{png,svg}` |

### Stage 6 — ecology  (`rules/ecology.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `genomad` | Predict proviruses / viral contigs (viral-pressure proxy) | **geNomad** `end-to-end` | `ecology/viral/{g}/{g}_summary/{g}_virus_summary.tsv` |
| `checkv` | Provirus completeness / quality | **CheckV** `end_to_end` | `ecology/viral/{g}/checkv/quality_summary.tsv` |
| `rloop_predict` | R-loop-forming-sequence (RLFS) prediction in proviruses | `rloop_predict.py` (RIZ/REZ model) | `ecology/viral/{g}/rloop.tsv` |
| `viral_summary` | Aggregate proviral load + CheckV + R-loop per genome | `aggregate_viral.py` | `ecology/viral/viral_summary.tsv` |
| `ecology_stats` | Environment/host/viral association, phylogeny-controlled | **R**: `phylolm` phyloglm | `ecology/ecology_associations.tsv` |
| `plot_ecology` | Environment enrichment + viral/R-loop figures | `plot_ecology.py` | `ecology/plots/*.{png,svg}` |

### Stage 7 — synthesis  (`rules/synthesis.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `plot_validation_tree` | Validation tree (candidates vs references, by tier) | **R**: `ggtree` | `prometheus/plots/prometheus_validation_tree.{png,svg}` |
| `master_tree_figure` | Annotated species tree (presence/tier, methanogen, ecology, proviral rings) | **R**: `ggtree` + `gheatmap` + `ggnewscale` | `synthesis/figures/master_tree.{png,svg}` |
| `results_registry` | Collect every test; one Benjamini–Hochberg FDR regime; tag confirmatory/exploratory | `results_registry.py` | `synthesis/results_registry.tsv` |
| `plot_summary_forest` | Forest of confirmatory results | `plot_summary_forest.py` | `synthesis/figures/summary_forest.{png,svg}` |

`provenance` (in `Snakefile`) records resolved tool/DB versions to `results/provenance/versions.txt`.

### External programs (versions pinned in `workflow/envs/`)

| Program | Role | Env(s) |
|---|---|---|
| Prodigal | ORF / gene prediction | qc, defensefinder, gtdbtk |
| CheckM2 | genome completeness / contamination | qc |
| seqkit | assembly statistics | qc |
| barrnap | rRNA detection | qc |
| tRNAscan-SE | tRNA detection | qc |
| dRep | species dereplication (fastANI) | qc |
| DefenseFinder (MacSyFinder) | anti-phage system detection (Prometheus source) | defensefinder |
| HMMER (`hmmscan`/`hmmsearch`) | profile-HMM search (Pfam domains, mcr genes) | annotation, phylo |
| MMseqs2 | sequence clustering (subfamilies) | annotation |
| MAFFT | multiple sequence alignment | phylo |
| FAMSA | large-set alignment fallback | phylo |
| ClipKIT | alignment trimming | phylo |
| IQ-TREE2 | ML phylogenetics (ModelFinder, UFBoot2, PMSF) | phylo |
| GTDB-Tk | taxonomy + concatenated marker MSA | gtdbtk |
| PAL2NAL | protein→codon alignment | selection |
| GeneRax | species-tree-aware DTL reconciliation | reconciliation |
| HyPhy | selection: FEL / MEME / BUSTED / aBSREL / RELAX | selection |
| geNomad | provirus / virus prediction | viral |
| CheckV | virus genome quality | viral |
| R `phylolm` | phylogenetic logistic regression | rstats |
| R `phytools` | fitPagel, make.simmap, cophylo | rstats |
| R `caper` | phylogenetic D statistic | rstats |
| R `ape` / `geiger` | tree handling, tip pruning + ordering | rstats |
| R `ggtree` / `treeio` | tree figures | rstats |
| R `svglite` / `ragg` | editable SVG + 600-dpi PNG | rstats |
| matplotlib | Python figures (PNG + SVG) | plotting |

In-house scripts (`workflow/scripts/py/*.py`, `workflow/scripts/R/*.R`) are named per rule above; the data contracts between them are documented in `workflow/scripts/SCHEMA.md`.

## 6. Outputs

Key tables (all TSV; schema in `workflow/scripts/SCHEMA.md`):

```
results/qc/genome_metadata.tsv               quality, tiers, pass flags
results/prometheus/prometheus_catalog.tsv    curated, tiered loci (the catalog)
results/prometheus/detection_sensitivity.tsv per-method detection
results/scaffolds/master_table.tsv           the analysis backbone (one row/genome)
results/distribution/association_tests.tsv   naive + phylogenetic + Pagel + D
results/distribution/cooccurrence.tsv        Prometheus vs other systems
results/evolution/reconciliation/...         DTL events, transfer fraction
results/evolution/ancestral/gains_losses.tsv independent acquisitions, enrichment
results/evolution/selection/selection_summary.tsv  FEL/MEME/BUSTED/aBSREL/RELAX
results/ecology/ecology_associations.tsv     environment/viral, phylo-controlled
results/synthesis/results_registry.tsv       every test, one FDR regime
results/provenance/versions.txt              tool/db versions for the methods
```

Figures — **every figure is written as both `.png` (600 dpi) and `.svg`**:

```
qc/plots/qc_completeness_contamination, qc_assembly_stats     (diagnostic)
prometheus/plots/prometheus_detection_summary, _architecture, _validation_tree
diversity/plots/subfamily_rarefaction, architecture_by_taxon
distribution/plots/prevalence_by_taxon, association_forest, cooccurrence_heatmap
evolution/plots/genetree_speciestree_tanglegram, ancestral_gain_loss, selection_by_domain
ecology/plots/ecology_enrichment, viral_rloop_association
synthesis/figures/master_tree, summary_forest
```

### Editing figures in Adobe Illustrator

The SVGs keep **text as editable `<text>`** (matplotlib `svg.fonttype='none'`;
R via `svglite`), so labels stay selectable and editable in Illustrator rather than
being outlined to paths. Figures use a standard sans family (Arial, falling back to
DejaVu Sans, which ships with matplotlib) so glyphs are present everywhere; if a
machine lacks Arial, Illustrator substitutes the family without rasterizing. The
600-dpi PNG is the always-identical raster companion. (If you also want a
font-embedded vector that renders identically on any machine, open the SVG and
`File ▸ Save As ▸ PDF`, or ask and I'll add a PDF export rule — matplotlib embeds
TrueType as Type 42 in PDF.)

## 7. Verification status

- `snakemake -n` builds the full 63-job DAG with correct per-genome fan-out, no
  missing-input/cyclic/ambiguity errors; config + samples schemas validate.
- All 27 Python scripts pass `python -m py_compile`; each tolerates empty inputs
  (writes header-only outputs) so the DAG never breaks on sparse data.
- The 8 R scripts were structurally reviewed (flag parity, bracket balance, output
  headers vs schema). They have **not** been executed against a live R 4.3 stack in
  this environment — run the `rstats` rules on a small test set first.
- The first real run should be a subset of genomes to confirm tool/database paths
  and resource sizing before scaling out.
- A full robustness/accuracy audit (tool-format verification against docs,
  synthetic runtime tests of the Python steps, statistical review, and the bugs
  fixed) is documented in [`AUDIT.md`](AUDIT.md).
