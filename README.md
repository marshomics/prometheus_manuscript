# Reference data the pipeline expects

None of this is committed (see `.gitignore`). Stage it once on shared or scratch storage
and point `config/config.yaml` at it. **Pin every release.** The `provenance` rule writes
resolved tool and database versions to `results/provenance/versions.txt`; quote that file
in the methods. The pipeline never downloads anything itself.

## Genomes

Phase 1 reads the assemblies from `input.genomes_dir`: one FASTA per row of
`input.metadata`, named in its `filename` column (extension `input.genome_extension`).
`filter_collapse` stages the chosen species representatives under
`results/preprocess/genomes/`, and records each path in the `assembly_path` column of
`results/preprocess/species_metadata.tsv`. Phase 2 reads only those staged copies. There is
no samples sheet to write by hand.

**DefenseFinder ingest (optional).** With `prometheus.run_defensefinder: false`, the
`defensefinder` rule copies, for each genome, exactly one `*_defense_finder_systems.tsv` and
one `*_defense_finder_genes.tsv` from `prometheus.defensefinder_provided_dir/{genome_id}/`. A
missing or ambiguous pair is an error. Every `hit_id` is then checked against this pipeline's
own Prodigal proteome, `results/genes/{genome_id}/{genome_id}.faa`, whose ids are
`{genome_id}_{contig#}_{gene#}`. The tables must therefore come from DefenseFinder run on
those proteins: tables from another gene caller fail the check. There is no merged-table
mode.

## Databases (`config: db`)

| Key | What | How to get it |
|---|---|---|
| `checkm2` | the CheckM2 DIAMOND database **file** (`.dmnd`) | from the `checkm2` env: `checkm2 database --download --path resources/db/checkm2` writes `resources/db/checkm2/CheckM2_database/uniref100.KO.1.dmnd` (the config default). CheckM2 1.0.2 accepts database versions ≤ 2 (Zenodo 10.5281/zenodo.5571251); `checkm2 predict` sha256-checks the file and fails on version 3 (Zenodo 14897628, for CheckM2 ≥ 1.1.0). |
| `gtdbtk_data` | GTDB-Tk reference package | **R220** for GTDB-Tk 2.4.0: `gtdbtk_r220_data.tar.gz` (MD5 `5aafa1b9c27ceda003d75adf238ed9e0`) from the GTDB data server; unpack and point `db.gtdbtk_data` at the release directory. |
| `pfam_hmm` | `Pfam-A.hmm`, pressed | from the Pfam/InterPro FTP; record the release; run `hmmpress Pfam-A.hmm`. An unpressed file fails `prometheus_domain_scan`. |
| `kofam_ko_list` | KOfam `ko_list` (curated per-KO score thresholds) | from the KEGG KOfam FTP (`ko_list.gz`, gunzip). **Required:** KOfam profiles carry no GA lines, so this table is the only acceptance threshold for the mcr search. |
| `genomad` | geNomad database directory (contains `version.txt`) | `genomad download-database resources/db` → `resources/db/genomad_db` |
| `checkv` | CheckV database directory (contains `genome_db/` and `hmm_db/`) | `checkv download_database resources/db` → `resources/db/checkv-db-v1.5` (the config default) |
| `prometheus.defensefinder_models_dir` | DefenseFinder models | the pinned releases, installed into a project directory: `macsydata install --target resources/db/defense-finder-models --org mdmparis 'defense-finder-models==3.1.0'` and `macsydata install --target resources/db/defense-finder-models 'CasFinder==3.1.0'`. Do **not** use `defense-finder update`: it installs the latest CasFinder, whose grammar macsyfinder 2.1.4 rejects. The `defensefinder_check` preflight fails on any other model version, and runs a positive control first. |

## Reference sequences (`config: refs`) — the defensibility-critical inputs

| File | What it is | Why it matters |
|---|---|---|
| `prometheus_reference.faa` | validated Prometheus proteins, **including the experimentally characterised *Pseudomonas aeruginosa* Prometheus** (van den Berg et al. 2024, *Cell Host & Microbe* 32:1427, PMID 39094584) | the ingroup of the validation tree: candidates that place here with support (every reported value meets its threshold: SH-aLRT ≥ `prometheus.placement_min_alrt`, UFBoot ≥ `prometheus.placement_min_ufboot`) can reach tier A |
| `paralog_reference.faa` | housekeeping outgroup: eukaryotic and bacterial DNA2, archaeal NurA and HerA, other helicase-nucleases | roots the validation tree; a candidate that places with these is **Rejected** |
| `mcr_hmm/` | one `*.hmm` per `methanogen.mcr_kos`: copy `K00399.hmm` (McrA), `K00401.hmm` (McrB), `K00402.hmm` (McrG) from KOfam `profiles/` (TIGRFAM / NCBIfam models with GA lines also work) | defines "methanogen" from genes, not taxonomy; a listed KO without a model is an error |
| `mcrA_reference.afa` + `mcrA_reference_labels.tsv` | curated McrA alignment; every sequence labelled `methanogen`, `other` (e.g. alkyl-CoM reductases) or `ANME*` | separates true methanogens from anaerobic methanotrophs and divergent McrA-like enzymes; both ANME and methanogen labels must be present |

Prometheus detection comes only from DefenseFinder. These reference proteins are used for
validation (phylogenetic placement against cellular paralogs), never for an HMM detection
sweep.

**Curating `prometheus_reference.faa`.** Use canonical, full-length, validated ProA proteins,
ideally several spanning the known diversity. Divergent homologues that fall below
DefenseFinder's ProA threshold (e.g. kinase-fused variants) make poor references. In the
review's integration test, two such proteins grouped with the paralogs and left the paralog
root edge at ~77 UFBoot. Genuine candidates then stalled at 91–95 UFBoot ("Unplaced", so
not tier A). With the consensus, a four-domain ProA and a haloarchaeal ProA as references, the
same candidates placed at 98–100. Check `placement_support` in the catalog and the root-edge
support in the `validate_prometheus` log.

## A note on Pfam accessions

Domains are matched **exactly**, first on the version-stripped Pfam accession and then on
the case-insensitive model name; there is no substring matching (`config:
prometheus.domains`). DUF4011 = PF13195; the DNA2 helicase core = AAA_11 PF13086 / AAA_12
PF13087 / AAA_19 PF13245 / AAA_30 PF13604; REase = MTES_1575 PF18741; DUF3320 = PF11784
(optional: absent from ~54% of Prometheus, van den Berg et al. 2024). Verify any accession
you rely on against your installed Pfam release, since accessions and names can drift
between releases.

## Bacterial / cross-domain branch (`config: bacteria`, opt-in)

Stage your pre-computed bacterial inputs and point `bacteria:` at them. The archaeal run
ignores all of this unless `bacteria.enabled: true`.

| Key | What |
|---|---|
| `df_systems`, `df_genes` | combined bacterial DefenseFinder systems and genes tables, each with a genome column (`systems_genome_col` / `genes_genome_col`) |
| `faa_dir`, `cds_dir` | per-genome `{genome_id}.faa` proteomes (the ones DefenseFinder was run on), and optionally `{genome_id}.fna` CDS |
| `taxonomy`, `qc` | GTDB classification and completeness / contamination per bacterial genome |
| `gtdb_bac_tree` (+ optional `gtdb_bac_taxonomy`) | GTDB **bac120** reference tree (newick) from the same release as the archaeal data, pruned to your genomes for the co-occurrence, niche and viral tests |
| `niche_metadata` (optional) | oxygen / biome / host labels; activates the anaerobe / niche specificity control |
| `genomes_dir` + `genome_extension` (optional) | bacterial **assemblies** (`{genome_id}.<ext>`); activates the provirus / R-loop mechanism test (geNomad, CheckV and R-loop on each assembly). Needs the `genomad` / `checkv` databases and, for the phylogeny-controlled test, `gtdb_bac_tree`. |

The archaeal `prometheus_reference.faa` and `paralog_reference.faa` also anchor the
cross-domain tree, so the same references serve both domains.

## Intrinsic-disorder module (`config: disorder`, opt-in)

metapredict 3.1.1 installs with `workflow/envs/disorder.yaml` (its predicted-pLDDT network,
the "alphaPredict" network, is built in), and segmasker comes with BLAST+ 2.17.0. The
recommended second predictor, **IUPred3**, is free for academic use but licence-gated, so
it is not installable here. Download it once (register at the IUPred site), unpack it, and
set `disorder.iupred3_path` to the folder containing `iupred3.py` / `iupred3_lib.py`. With
that set, the consensus is metapredict + IUPred3, two independent methods. Without it, the
module uses metapredict + its pLDDT network and labels the result as the weaker,
correlated consensus.
