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
config/            config.yaml (input + preprocess + analysis settings)
resources/         db/, refs/ (you provide; see resources/README.md)
input_genomes_metadata_with_niches.tsv   your raw input table (Phase-1 input)
workflow/
  Snakefile        includes preprocess + 8 stage modules + provenance
  rules/*.smk      common + preprocess, qc, prometheus, diversity, scaffolds,
                   distribution, evolution, ecology, synthesis
  envs/*.yaml      11 per-stage conda environments
  scripts/py/*.py  33 preprocessing + analysis + plotting scripts (Illustrator-ready)
  scripts/R/*.R    8 phylogenetic-stats + tree-figure scripts (ggtree/phytools)
  schemas/         config + samples JSON schemas (validated at start)
  profiles/sge/    Grid Engine submission profile + wrapper
results/           preprocess/ (Phase 1) + analysis stages (Phase 2)
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

The pipeline starts from a metadata table over **all** your input genomes plus the
genome FASTAs; Phase-1 preprocessing turns those into the species-level,
quality-filtered set the analysis uses. You do not prepare a samples sheet or a
genomes folder by hand, and dRep is not used.

- **Input metadata** (`config: input.metadata`, e.g. `input_genomes_metadata_with_niches.tsv`):
  one row per genome with `genome_id`, `filename`, `source`, and one binary (0/1)
  column per environment named `group_*` / `niche_*` (1 = the genome was observed
  in that group/niche). These niche columns become the ecology predictors.
- **Input genome FASTAs**: point `config: input.genomes_dir` at the directory that
  holds the files named in the metadata's `filename` column (your store, not the
  working folder).
- **Databases + reference sequences**: stage them and point `config/config.yaml` at
  them. **`resources/README.md` lists every file.** The reference protein sets
  (`prometheus_reference.faa`, `paralog_reference.faa`, the mcr HMMs, and the
  labelled McrA reference alignment) are the defensibility-critical inputs.

Phase 1 then writes **one** filtered table,
`results/preprocess/species_metadata.tsv` (the genome list + taxonomy + QC +
unioned niches: both the samples sheet and the metadata), and stages the
representative genomes. That single table is what Phase 2 consumes. Input genomes
are all `*.fasta` (`config: input.genome_extension`); only genomes GTDB-Tk places in
`d__Archaea` are kept; novel / unassigned GTDB species are dropped; MAG/isolate type
is not recorded or used anywhere (quality is completeness, contamination, contig
count, and N50). The table carries an `assembly_path` column naming the assembly each
step reads. Gene calling (Prodigal) then writes per-genome `.faa` proteins and `.fna`
CDS whose headers are rewritten to `{genome}_{N}` (identical between the two files);
the `.faa` is the input to DefenseFinder and the downstream analyses.

If you already have DefenseFinder outputs and prefer to ingest rather than re-run,
set `prometheus.run_defensefinder: false` and drop them under
`resources/defensefinder/{genome_id}/`. Re-running under a pinned version (the
default) is more reproducible.

## 3. Run

Two phases, because the analysis genome set is decided by preprocessing. Dry-run
each first (`-n`).

**Phase 1 — preprocess** (classify, QC, score, quality-filter, collapse to species):

```bash
snakemake preprocess -n                                              # check the preprocess DAG
snakemake preprocess --workflow-profile workflow/profiles/sge --use-conda --jobs 200
```

Batches GTDB-Tk / CheckM2 / BBMap (~20k genomes/batch), scores and filters, and
collapses to one representative per GTDB species (unioning niche/group presence).
Inspect `results/preprocess/plots/preprocess_quality.{png,svg}` and
`species_metadata.tsv` before launching the expensive analysis.

**Phase 2 — analysis** (runs on the species representatives Phase 1 produced):

```bash
snakemake -n                                                        # now sees the species set
snakemake --workflow-profile workflow/profiles/sge --use-conda --jobs 200
```

Local workstation (either phase): `snakemake [preprocess] --use-conda --cores 24
--resources mem_mb=120000`.

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

- **Genome set = quality-filtered GTDB species representatives.** Preprocessing
  classifies all input genomes (GTDB-Tk, batched), QC's them (CheckM2 + BBMap),
  scores them (completeness − 5×contamination), keeps completeness ≥ 90,
  contamination ≤ 5, < 1000 contigs, contig-N50 > 5 kb, score > 50, and collapses to
  the highest-scoring genome per GTDB species. Niche/group presence is unioned across
  each species' genomes (by default over all of them, so quality filtering doesn't
  erase where a species occurs). This replaces dRep; no separate QC runs before
  analysis. Novel/unassigned GTDB species are dropped, and MAG/isolate type is not
  used anywhere — where models need quality it enters as completeness, contamination,
  contig count, and N50, not a genome-type label.
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

### Phase 1 — preprocessing  (`rules/preprocess.smk`; `snakemake preprocess`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `batch_genomes` | Split inputs into ~20k batches; stage `{genome_id}.fna` symlinks | `batch_genomes.py` | `preprocess/batches/batch_*/.ready` |
| `gtdbtk_batch` | GTDB taxonomy per batch | **GTDB-Tk** `classify_wf` | `.../gtdbtk_summary.tsv` |
| `checkm2_batch` | Completeness / contamination per batch | **CheckM2** | `.../checkm2_quality.tsv` |
| `bbmap_batch` | Assembly stats (GC, length, contigs, N50) per batch | **BBMap** `statswrapper.sh` (format=6) | `.../bbmap_stats.tsv` |
| `merge_annotations` | Join taxonomy + CheckM2 + BBMap; add quality score | `merge_annotations.py` | `preprocess/annotated_metadata.tsv` |
| `filter_collapse` | Keep `d__Archaea`; quality filter → 1 rep per GTDB species (novel dropped); union niches; add `assembly_path` | `filter_collapse.py` | `preprocess/species_metadata.tsv` (the single table) + staged genomes |
| `plot_preprocess` | Quality scatter + input→pass→species funnel | `plot_preprocess.py` | `preprocess/plots/preprocess_quality.{png,svg}` |

### Stage 0 — gene calls + per-genome QC table  (`rules/qc.smk`)

| Rule | What it does | Program(s) run | Main output |
|---|---|---|---|
| `prodigal_genes` | Predict genes; rewrite headers to `{g}_{N}` (consistent .faa/.fna) | **Prodigal** + `rename_prodigal.py` | `genes/{g}/{g}.{faa,fna,gff}` |
| `genome_metadata` | Per-genome QC table (from preprocessing; no QC re-run) | `cp` | `qc/genome_metadata.tsv` |
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
| `gtdbtk_align` | Concatenated ar53 marker MSA of the species reps (taxonomy itself comes from Phase 1) | **GTDB-Tk** `identify` + `align` | `scaffolds/gtdbtk_align/user_msa.fasta.gz` |
| `species_tree` | ML archaeal species tree, two-step PMSF, then rooted (outgroup/midpoint) | **IQ-TREE2** + `root_tree.py` | `scaffolds/species_tree.treefile` |
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
| `combined_cds` | Concatenate CDS with matching IDs | `combine_proteomes.py` (`--suffix fna`) | `evolution/all_cds.ffn` |
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
| CheckM2 | completeness / contamination (Phase-1, batched) | qc |
| BBMap (`statswrapper.sh`) | assembly stats: GC, length, contigs, N50 (Phase-1, batched) | bbmap |
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

- `snakemake preprocess -n` builds the Phase-1 DAG (one batch set per ~20k genomes;
  20 jobs for the 92k-genome input), and `snakemake -n` builds the full preprocess +
  analysis DAG, both with no missing-input/cyclic/ambiguity errors and validated schemas.
- All 33 Python scripts pass `python -m py_compile`; each tolerates empty inputs
  (writes header-only outputs) so the DAG never breaks on sparse data. The
  preprocessing merge/filter/species-collapse/niche-union logic is unit-tested.
- The 8 R scripts were structurally reviewed (flag parity, bracket balance, output
  headers vs schema). They have **not** been executed against a live R 4.3 stack in
  this environment — run the `rstats` rules on a small test set first.
- The first real run should be a subset of genomes to confirm tool/database paths
  and resource sizing before scaling out.
- A full robustness/accuracy audit (tool-format verification against docs,
  synthetic runtime tests of the Python steps, statistical review, and the bugs
  fixed) is documented in [`AUDIT.md`](AUDIT.md).
