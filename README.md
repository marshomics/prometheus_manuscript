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

## Scientific and computational plan

The aim is to characterize Prometheus across archaea (how diverse it is, where it sits
taxonomically and phylogenetically, how it evolves, and which traits and habitats it
tracks) and to test whether it is associated with methanogens. Optional branches bring
bacteria in as the baseline the archaeal findings are read against (origin and direction
of transfer, a matched comparison, and a validation of the R-loop mechanism) and add a
comparison of ProA intrinsic disorder.

The workflow is a sequence of questions, each answered by a stage whose outputs feed the
next. Phase 1 turns the raw genome set into quality-filtered, species-level
representatives; Phase 2 runs the analysis on that set. The per-rule reference is §5; the
table below maps each question to what the pipeline does and the results that answer it.

| Scientific question | What the pipeline does | Results that answer it |
|---|---|---|
| Are the archaeal hits genuine Prometheus, not cellular DNA2 / NurA-HerA paralogs? | DefenseFinder detection, then validation by domain architecture, phylogenetic placement against Prometheus vs paralog references, defense-island context, and length/fragmentation; loci tiered A/B/C/Rejected | `prometheus_catalog.tsv`, `detection_sensitivity.tsv`, `prometheus_proteins.faa`, validation-tree + detection-funnel figures |
| Where does it sit on the archaeal tree, and is it clade-restricted or scattered? | prevalence by GTDB rank with Wilson 95% CIs; the Fritz–Purvis D statistic for phylogenetic clumping | `prevalence_by_taxon.tsv`, prevalence-caterpillar figure, D term in `association_tests.tsv` |
| How diverse is the system, in sequence and in architecture? | MMseqs2 subfamily clustering at graded identities; enumeration of domain architectures (losses, fusions, reorderings) | `subfamilies.tsv`, `architecture_typology.tsv`, subfamily + architecture figures |
| Is it associated with methanogens, beyond phylogeny and confounders? | methanogens defined from genes (McrA, with ANME separated); phylogenetic logistic regression adjusted for completeness, genome size, and defense load excluding Prometheus; Pagel correlated-evolution test; a naive test kept only as a labelled baseline | `association_tests.tsv`, association-forest figure, `master_table.tsv` |
| Does it co-occur with, or exclude, other defense systems? | correlated-evolution test (fitPagel) per system with a Fisher odds-ratio direction, BH-corrected | `cooccurrence.tsv`, co-occurrence odds-ratio + heatmap figures |
| Is it vertically inherited or transferred, how many times was it acquired, and are gains enriched in methanogen lineages? | GeneRax DTL reconciliation (transfer fraction); stochastic mapping for independent acquisitions and gains/losses, with a phylogenetically structured methanogen-enrichment null; gene-tree vs species-tree tanglegram; per-domain trees + Robinson–Foulds for domain shuffling | `reconciliation_events.tsv`, `gains_losses.tsv`, `rf_distances.tsv`, tanglegram / ancestral / gain-enrichment figures |
| Is it under diversifying selection, and on which domains? | HyPhy FEL/MEME/BUSTED/aBSREL/RELAX on catalytically-intact loci, with RELAX contrasting methanogen branches | `selection_summary.tsv`, `site_table.tsv`, dN/dS-along-protein + selection-by-domain figures |
| Does it track environment and viral pressure / R-loops? | environment and niche phylogenetic regression (temperature-aware); geNomad proviruses and the QmRLFS R-loop-prone fraction tested against presence, completeness- and size-adjusted | `ecology_associations.tsv`, `viral_summary.tsv`, ecology + viral/R-loop figures |
| Which results are confirmatory, and do they survive one multiple-testing regime? | every test collected, BH-corrected per family, tagged confirmatory vs exploratory; an annotated species tree | `results_registry.tsv`, `master_tree`, `summary_forest` |
| (opt-in) Where did the archaeal copies come from, and how do they compare with bacteria? | cross-domain phylogeny with a topology-based transfer direction and independent-acquisition count; matched prevalence and architecture; cross-domain dN/dS and co-occurrence | `crossdomain_summary.tsv`, bacterial catalog/prevalence/typology, cross-domain figures |
| (opt-in) Does the R-loop mechanism hold in the domain where it was demonstrated? | geNomad/CheckV/R-loop on bacterial assemblies; bacterial Prometheus vs proviral R-loop load, phylogeny-controlled | `bacterial_viral_associations.tsv`, `viral_rloop_crossdomain` figure |
| (opt-in) Do ProA disordered regions vary across niches, lineages, and domains? | consensus disorder (two independent predictors) on intact loci, low-complexity-aware, alignment-anchored, compared with phylolm | disorder tables + along-protein / by-group figures |

Two choices run through every row. Absence is modeled, not assumed: genome completeness,
a detection-sensitivity flag, genome size, and defense load all enter as covariates, so a
missing call is treated as possible non-detection rather than confirmed absence. And
every association is phylogenetically controlled, because methanogens are polyphyletic
and habitat tracks clade in archaea; the naive contingency result is reported only as a
labelled baseline. The conventions that make each step defensible (genome set,
methanogen-from-genes, single-detector validation, rooting, one FDR regime) are listed
in §4.

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

## 2. Required inputs and configuration

Phase 1 starts from a metadata table over **all** your input genomes plus the genome
FASTAs, and turns them into the species-level, quality-filtered set Phase 2 analyzes.
You do not prepare a samples sheet, a genomes folder, or run dRep by hand: Phase 1
writes `results/preprocess/species_metadata.tsv`, the single table both downstream
phases consume. Everything you supply is set in `config/config.yaml`; paths are
relative to the working directory unless absolute. `resources/README.md` gives the
download command for each database and reference file.

### 2a. What you must provide (core archaeal run)

| Input | config key(s) | What it is / how to get it |
|---|---|---|
| Input genome metadata | `input.metadata` | one row per genome: `genome_id`, `filename`, `source`, and one 0/1 column per environment named `group_*` / `niche_*` (1 = observed there). The niche columns become the ecology predictors. |
| Input genome assemblies | `input.genomes_dir`, `input.genome_extension` (default `fasta`) | **must set** `input.genomes_dir` — the directory holding the FASTAs named in the metadata `filename` column. |
| CheckM2 database | `db.checkm2` | completeness/contamination, Phase 1. `checkm2 database --download`. |
| GTDB-Tk reference data | `db.gtdbtk_data` | taxonomy (Phase 1) + ar53 marker MSA (Phase 2). Download the release matching your GTDB-Tk version. |
| Pfam-A HMMs (pressed) | `db.pfam_hmm` | domain architecture. `hmmpress Pfam-A.hmm`. |
| geNomad database | `db.genomad` | provirus / virus prediction (ecology). `genomad download-database`. |
| CheckV database | `db.checkv` | provirus quality (ecology). `checkv download_database`. |
| Prometheus reference proteins | `refs.prometheus_faa` | validated Prometheus incl. the characterized *P. aeruginosa* protein; roots the validation tree. |
| Paralog outgroup proteins | `refs.paralog_faa` | DNA2 / NurA / HerA and other helicase-nucleases; the outgroup that rejects cellular paralogs. |
| mcr profile HMMs | `refs.mcr_hmm_dir` | mcrA / mcrB / mcrG HMMs; defines methanogens from genes. |
| McrA reference alignment + labels | `refs.mcrA_ref_aln`, `refs.mcrA_ref_labels` | curated McrA alignment with ANME / methanogen labels; separates methanogens from ANME. |
| DefenseFinder models | (post-install, not a path) | run `defense-finder update` to fetch models **≥ v2.0.0** (the version that first contains Prometheus). The env pins the software at ≥ 2.0.0 and the `defensefinder_check` rule fails fast if either is too old. |

The four reference sets (`prometheus_faa`, `paralog_faa`, the mcr HMMs, and the labelled
McrA alignment) are the defensibility-critical inputs — detection and validation are only
as good as these, so curate them deliberately. `db.kofam_profiles` / `db.kofam_ko_list`
and `db.interproscan` are optional and left `""` by default.

### 2b. Config keys to review before launching

| Key | Set it to | Why |
|---|---|---|
| `input.genomes_dir` | your assemblies directory | required; empty by default |
| `db.*` and `refs.*` | the staged paths in 2a | required for the stages that use them |
| `tree.outgroup` **or** `tree.outgroup_taxon` | a lineage appropriate to your genome set | gives a defensible root; otherwise the species tree falls back to midpoint, a heuristic for rate-variable archaeal trees |
| `prometheus.run_defensefinder` | `true` to run, `false` to ingest | if `false`, drop existing DefenseFinder output under `prometheus.defensefinder_provided_dir/{genome_id}/` |
| `preprocess.quality.*`, `qc.*`, `stats.confidence_tiers_primary` | defaults are defensible | change only deliberately (see §4) |

Everything else has a defensible default (§4).

### 2c. Optional modules (extra inputs, off by default)

Each adds nothing to the core run until enabled, and is built with its own target
(`crossdomain`, `disorder`); see §8 and §9 for the science.

**Bacterial / cross-domain branch** — `bacteria.enabled: true`, target `crossdomain`:

| Input | config key(s) | Required when enabled? |
|---|---|---|
| Combined DefenseFinder systems + genes tables | `bacteria.df_systems`, `bacteria.df_genes` | yes |
| Per-genome proteomes + CDS | `bacteria.faa_dir`, `bacteria.cds_dir` | yes (CDS needed for cross-domain dN/dS) |
| GTDB taxonomy + QC (completeness) | `bacteria.taxonomy`, `bacteria.qc` | yes |
| GTDB bac120 reference tree | `bacteria.gtdb_bac_tree` | for the phylogeny-controlled tests (co-occurrence, niche, viral) |
| Genome-column names in the DF tables | `bacteria.systems_genome_col`, `bacteria.genes_genome_col` | adjust to match your tables |
| Oxygen / biome / host labels | `bacteria.niche_metadata` | optional — activates the niche specificity control |
| Bacterial assemblies | `bacteria.genomes_dir`, `bacteria.genome_extension` | optional — activates the provirus / R-loop mechanism test |

*Format of the bacterial DefenseFinder inputs.* `df_systems` and `df_genes` are your
per-genome DefenseFinder `*_defense_finder_systems.tsv` and `*_defense_finder_genes.tsv`
concatenated into one of each, with a **genome column added** so every row records which
genome it came from. Only a few columns are read; the rest are carried through unchanged.

- `df_systems` must contain: the genome column (`bacteria.systems_genome_col`, default
  `genome_id`; `genome` / `user_genome` / `sample` / `assembly` / `replicon` are also
  recognized), `sys_id`, and `type` and/or `subtype`. The system is matched
  case-insensitively against `bacteria.system` (default `Prometheus`) in `type`/`subtype`;
  `type` is also read downstream to build the bacterial defense matrix, so keep it.
- `df_genes` must contain: the genome column (`bacteria.genes_genome_col`), `sys_id`
  (matching the systems table for the same genome), `hit_id` (the protein identifier),
  and optionally `hit_score` (or `score`).
- `faa_dir` holds `{genome_id}.faa` — the **same proteomes DefenseFinder was run on**, so
  each `hit_id` equals a protein's FASTA header (first whitespace token). `cds_dir`
  (optional) holds the matching `{genome_id}.fna` CDS used for cross-domain dN/dS.

`extract_bacterial_prometheus.py` consumes these: it splits `df_systems` into the
per-genome layout the archaeal scripts already understand
(`results/bacteria/defense/{g}/{g}_defense_finder_systems.tsv`), so `validate_prometheus`,
`build_defense_matrix`, and the defense-load count reuse it unchanged, and it pulls the
Prometheus systems' member proteins (via the `df_genes` `hit_id`s) from `faa_dir` to form
the bacterial candidates that feed validation and the cross-domain tree. The `genome_id`
values must be identical across `df_systems`, `df_genes`, the `.faa`/`.fna` filenames,
`bacteria.taxonomy`, and `bacteria.qc`, and must not contain `__`.

**Intrinsic-disorder module** — `disorder.enabled: true`, target `disorder`:

| Input | config key | Required when enabled? |
|---|---|---|
| IUPred3 install | `disorder.iupred3_path` | recommended — the independent second predictor (free academic, license-gated; download once). Without it the module uses the metapredict + alphaPredict fallback and labels it as the weaker consensus. |
| Temperature-proxy niche columns | `disorder.temperature_cols` | optional — covariates for the niche tests |
| Host-proteome background table | `disorder.host_background` | optional — normalizes disorder against lineage baseline |

metapredict, alphaPredict, and segmasker install with the module's conda env; only
IUPred3 is staged by hand.

### 2d. What Phase 1 produces, and the ingest shortcut

Phase 1 writes `results/preprocess/species_metadata.tsv` (genome list + taxonomy + QC +
unioned niches), which serves as both the samples sheet and the metadata for Phase 2,
and stages one representative genome per GTDB species. Only genomes GTDB-Tk places in
`d__Archaea` are kept; novel / unassigned GTDB species are dropped; quality is
completeness, contamination, contig count, and N50 (MAG/isolate type is never used).
Gene calling (Prodigal) then writes per-genome `.faa` proteins and `.fna` CDS with
matching `{genome}_{N}` headers; the `.faa` feeds DefenseFinder and everything
downstream. Inspect that table and the Phase-1 quality plot before launching Phase 2.

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
  Prometheus only exists in **defense-finder-models v2.0.0+**, which needs the
  **software at ≥ 2.0.0** (an older binary detects zero Prometheus loci silently);
  the env pins this and a `defensefinder_check` preflight rule fails fast unless the
  detector can actually call Prometheus. `prometheus.strict_detection` additionally
  fails the run if DefenseFinder emits zero systems of any type over a non-empty input.
- **Validation separates Prometheus from housekeeping DNA2/NurA-HerA** by domain
  architecture, a catalytic-domain proxy, defense-island context, and phylogenetic
  placement against references. Candidates that place with cellular paralogs are
  rejected. Loci are tiered A/B/C; primary analyses use tier A, with A+B as a
  sensitivity check (`stats.confidence_tiers_primary`).
- **Absence is modeled, not assumed.** Genome completeness, a detection-sensitivity
  flag, genome size, and defense-system count (**excluding the focal Prometheus
  system**, so the model never conditions on its own response) enter the association
  models as covariates. The headline methanogen estimate is the phylogenetic-logistic
  term *after* this conditioning.
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
  reconciliation (ALE optional). Gains/losses come from stochastic mapping; the
  methanogen-enrichment test uses a *phylogenetically structured* null (the methanogen
  trait is simulated on the tree under a fitted Mk model, preserving its clustering),
  not free tip-label shuffling.
- **dN/dS is restricted to catalytically-intact loci** so pseudogene relaxation
  can't masquerade as signal; RELAX tests whether selection intensity differs on
  methanogen branches.
- **One FDR regime** (Benjamini-Hochberg, per test family) is applied centrally in
  the results registry over **every** test family — association, co-occurrence,
  ecology, selection, reconciliation, and ancestral (so the methanogen-gain-enrichment
  p and the transfer signal are corrected too). Confirmatory vs exploratory is tagged
  by (family, test, term): the phylogenetic methanogen tests (not the naive baseline),
  RELAX, transfer, and the gain-enrichment test are confirmatory.

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

Extra figures (`rules/plots_extra.smk` + the cross-domain rule; added alongside the
above, same PNG+SVG editable-text convention). Each shows effect size AND uncertainty
AND sample size, and degrades to a labelled placeholder on empty input:

```
qc/plots/prometheus_covariates                 size/defense-load/completeness by presence (the confound the model adjusts for)
distribution/plots/prevalence_caterpillar      prevalence ± Wilson CI by clade; CI vs overall flags enriched/depleted
distribution/plots/cooccurrence_oddsratio      directional log2 OR per system + BH significance
diversity/plots/architecture_combinations      UpSet of which domains co-occur, and how many loci
diversity/plots/architecture_by_phylum         architecture × phylum count heatmap
diversity/plots/subfamily_identity_curve       subfamily count vs clustering identity + rank-abundance
ecology/plots/viral_pressure_prometheus        proviral / R-loop load by presence + prevalence vs R-loop-pressure (Wilson CI)
ecology/plots/niche_effect_forest              phyloglm environment effects ± CI, BH-coloured
ecology/plots/prevalence_niche_phylum          prevalence bubble grid (niche × phylum): habitat vs clade
evolution/plots/dnds_along_protein             per-site dN/dS along the alignment + domain track (conflict signature)
evolution/plots/reconciliation_events          D/T/L/S counts per tool + vertical fraction
evolution/plots/gain_enrichment_null           methanogen-gain enrichment: structured null + observed + p
crossdomain/figures/crossdomain_comparison     bacteria vs archaea: prevalence (Wilson CI) + domain-architecture fractions
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

## 8. Bacterial / cross-domain branch (opt-in)

Prometheus was characterised in bacteria (van den Berg 2024); the archaeal pipeline
is the extension. This branch brings bacteria in as the baseline the archaeal numbers
are read against. It's **opt-in and additive**: it never changes the archaeal
`all`/`preprocess` runs, it only consumes the archaeal Prometheus proteins and CDS as
read-only inputs. Build it with the `crossdomain` target after setting `bacteria:` in
`config/config.yaml`.

You supply pre-computed bacterial DefenseFinder output (combined systems + genes
tables), per-genome `.faa` (and `.fna` CDS), GTDB taxonomy, and completeness. The
branch then:

- **Validates bacterial Prometheus with the SAME scripts as archaea.**
  `extract_bacterial_prometheus.py` splits the combined systems table into the
  per-genome layout the pipeline already understands, so `validate_prometheus.py`,
  `build_defense_matrix.py`, `compute_prevalence.py`, `architecture_typology.py` and
  `cooccurrence.R` are reused unchanged. Bacterial loci are tiered A/B/C exactly as
  archaeal ones.
- **Builds the cross-domain phylogeny** (archaeal + bacterial validated Prometheus +
  the Pseudomonas reference + DNA2/NurA-HerA outgroup) and answers the question the
  archaea-only design can't: where the archaeal copies came from.
  `crossdomain_analysis.py` reports, from topology alone (no DTL model — a single
  rooted species tree spanning both domains is contested), the minimum number of
  cross-domain transitions (Fitch on the paralog-rooted ingroup), whether the archaeal
  sequences are monophyletic, the count of independent archaeal acquisitions, and a
  polarity call (bacteria→archaea when the ingroup root is bacterial and the archaeal
  clades are derived/nested).
- **Compares prevalence and architecture** across domains, runs cross-domain dN/dS
  (FEL/MEME/BUSTED/aBSREL on the codon alignment; no methanogen RELAX partition across
  domains; per-domain homology caveats apply), and tests co-occurrence against the
  GTDB-bac120 reference tree pruned to your genomes.

Matched sampling is your call: bacterial GTDB dwarfs archaeal, so report prevalence
side-by-side and use species representatives with comparable QC rather than pooling.

Two further modules activate only when you supply the extra data. The
**anaerobe/niche specificity control** (`bacteria.niche_metadata`) is a clean control,
since bacteria can't be methanogens, for whether the archaeal methanogen signal is
partly a niche effect. The bacterial **provirus / R-loop mechanism test** turns on when
you point `bacteria.genomes_dir` at the bacterial assemblies (and set
`bacteria.genome_extension`): it runs geNomad → CheckV → the QmRLFS R-loop detector →
aggregation on each assembly (reusing the archaeal scripts unchanged), then tests
bacterial Prometheus presence against proviral load and R-loop-prone fraction in the
same phylogeny-controlled, completeness/size-adjusted `phyloglm` framework. This matters
because R-loop targeting was demonstrated in a bacterium (Pseudomonas), so reproducing
the signal in bacteria validates the proxy the archaeal side leans on. The
`viral_rloop_crossdomain` figure shows the R-loop-prone fraction by presence in both
domains and the prevalence-vs-R-loop-pressure dose-response (Wilson CIs) side by side.
The phylogeny-controlled test needs the bacterial species tree, so set
`bacteria.gtdb_bac_tree` too; without it you still get the descriptive figure and the
per-genome viral summary.

```bash
# after editing config: bacteria.enabled: true + paths
snakemake crossdomain -n                                            # check the cross-domain DAG
snakemake crossdomain --workflow-profile workflow/profiles/sge --use-conda --jobs 200
```

## 9. Intrinsic-disorder module (opt-in)

Large multidomain proteins like ProA carry disordered linkers and tails, and that
disorder can vary across niches and lineages. This module predicts it and compares it
in the same phylogeny-controlled framework as everything else. It's **opt-in**
(`disorder.enabled: true`, `disorder` target) and reuses the archaeal validated
proteins, catalog, family alignment, and species tree read-only; the default run is
untouched.

Disorder is called by **consensus of two independent methods** — a residue counts as
disordered only where both agree. metapredict (deep-learning) is always one; the second
is IUPred3 (biophysical, genuinely independent) when you stage it
(`disorder.iupred3_path`), otherwise alphaPredict (predicted pLDDT), which installs
freely but is more correlated with metapredict and is labelled as the weaker fallback.
IUPred3 is free for academic use but license-gated, so you download it once and point
the config at it; the recommended consensus is metapredict + IUPred3.

The controls that keep it defensible are baked in: it scores **intact, non-fragmented
loci only** (truncation fakes disorder at broken ends), reports disorder **excluding
low-complexity** regions (segmasker; LCR isn't disorder), and anchors disorder to the
**family alignment** so lineages are compared on homologous positions rather than by raw
whole-protein percentage. Comparisons use `phylolm` (Pagel's λ) on the species tree;
because disordered content tracks growth temperature, any columns you list in
`disorder.temperature_cols` are added as covariates to the niche tests. An optional
host-proteome background (`disorder.host_background`, a `genome_id, background_frac`
table) normalises each locus so the metric is "unusually disordered for its host," not
lineage baseline. All disorder results enter the registry tagged **exploratory** —
sequence prediction, not experiment.

Outputs: `disorder/disorder_per_locus.tsv` (IDR count/length, disorder fractions),
`disorder/disorder_by_column.tsv` (per-alignment-column disorder + domain), and figures
`disorder_along_protein` (with the domain track — disorder should pile up in linkers),
`disorder_by_group` (methanogen, and bacteria-vs-archaea when the bacterial branch is
on), and `disorder_linker_vs_domain`.

```bash
# stage IUPred3 (recommended) and set disorder.enabled: true + disorder.iupred3_path
snakemake disorder -n
snakemake disorder --workflow-profile workflow/profiles/sge --use-conda --jobs 200
```

## 10. Synthesis: what the results show, together

Three outputs pull the per-stage results into one picture. `master_table.tsv` carries
one row per genome (QC, taxonomy, species-tree tip, gene-based methanogen status,
Prometheus presence and tier, total and Prometheus-excluded defense load, and the
species-level niche columns), so every downstream test draws from the same backbone.
`results_registry.tsv` collects every test across families, applies one
Benjamini–Hochberg regime per family, and tags each result confirmatory or exploratory.
The `master_tree` figure paints presence and tier, methanogen and ANME status, ecology,
and proviral / R-loop load onto the rooted species tree, so the distribution and its
correlates are read in phylogenetic context rather than as a list.

Read together, they support a layered account. First, a confidence-tiered catalog of
genuine archaeal Prometheus, separated from the cellular DNA2 / NurA-HerA enzymes it
resembles. Second, its distribution: which lineages carry it, how clumped it is, and how
much its domain architecture varies. Third, the methanogen question, answered not by one
test but by agreement across the confounder-adjusted phylogenetic regression, the Pagel
correlated-evolution test, and the gains-on-methanogen-branches enrichment, three
independent angles held to one central FDR. Fourth, its evolutionary mode: the balance of
vertical inheritance against transfer, the number of independent acquisitions, and
whether the nuclease and REase domains carry the diversifying-selection signature
expected of a protein locked in conflict with phages. Fifth, its ecology: whether
presence tracks viral pressure and R-loop-prone proviruses, the prediction that follows
from the characterized mechanism.

With the bacterial branch on, two further layers join. The cross-domain phylogeny places
the archaeal sequences among bacterial Prometheus and polarizes the direction of
transfer (the expected reading being a bacterial origin with one or more acquisitions
into archaea), while the matched comparison asks whether archaeal Prometheus is rarer or
architecturally reduced relative to the experimentally characterized bacterial form. The
bacterial provirus/R-loop test re-runs the ecological prediction in the domain where the
mechanism was shown, so a positive result there validates the R-loop proxy the archaeal
side relies on. The disorder module, if run, adds where ProA carries intrinsic disorder
and whether that varies across niches and lineages.

The throughline is that no claim rests on a single number. Detection sensitivity and
absence are modeled, every association is phylogenetically controlled, confounders enter
explicitly rather than by assumption, and the confirmatory hypotheses are each
triangulated across independent methods and corrected together. What the pipeline is
built to deliver is a defensible, evolutionarily framed account of an anti-phage system
only recently recognized in archaea: what it is, where it lives, how it evolves, what it
travels with, and (with the optional branches) where it came from and whether its
mechanism generalizes beyond the bacterium in which it was found.
