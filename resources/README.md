# Reference data the pipeline expects

None of this is committed (see `.gitignore`); stage it once on shared/scratch and
point `config/config.yaml` at it. **Pin every release** and record the version in
`results/provenance/` (the pipeline writes tool/db versions there).

## Genomes
`resources/genomes/{genome_id}.fna` — one nucleotide assembly per row of
`config/samples.tsv`. Your existing DefenseFinder outputs go in
`resources/defensefinder/{genome_id}/` (the standard DefenseFinder per-genome
output folder). If you only have a merged DefenseFinder table, drop it at
`resources/defensefinder/all_systems.tsv` and set `prometheus.defensefinder_merged: true`.

## Databases (`config: db`)
| Key | What | How to get it |
|---|---|---|
| `checkm2` | CheckM2 DIAMOND db | `checkm2 database --download --path resources/db/checkm2` |
| `gtdbtk_data` | GTDB-Tk reference | download the release matching your gtdbtk version; `GTDBTK_DATA_PATH` |
| `pfam_hmm` | `Pfam-A.hmm` (run `hmmpress`) | InterPro/Pfam FTP; `hmmpress Pfam-A.hmm` |
| `interproscan` | optional `interproscan.sh` | EBI InterProScan (richer domain calls; Pfam-only works without it) |
| `kofam_profiles`, `kofam_ko_list` | KOfam HMMs + thresholds | KEGG KOfam ftp (`ko_list`, `profiles/`) |
| `genomad` | geNomad db | `genomad download-database resources/db` |
| `checkv` | CheckV db | `checkv download_database resources/db` |

## Reference sequences (`config: refs`) — the defensibility-critical inputs
| File | What it is | Why it matters |
|---|---|---|
| `prometheus_reference.faa` | Validated Prometheus proteins, **including the experimentally characterized *P. aeruginosa* Prometheus** (van den Berg et al. 2024, PMID 39094584) | Roots the validation tree: candidates that place here (not with paralogs) are kept |
| `paralog_reference.faa` | Housekeeping outgroup: eukaryotic/bacterial DNA2, archaeal NurA + HerA, other helicase-nucleases | Stage 1 places candidates against these to reject cellular paralogs |
| `mcr_hmm/` | mcrA/mcrB/mcrG profile HMMs (KOfam K00399/K00402/K00401 or TIGRFAM) | Defines "methanogen" from genes, not taxonomy |
| `mcrA_reference.afa` + `mcrA_reference_labels.tsv` | Curated McrA alignment with ANME and methanogen references, labelled | Separates true methanogens from anaerobic methanotrophs |

Detection of Prometheus comes only from DefenseFinder; these reference proteins are
used for validation (phylogenetic placement of candidates vs cellular paralogs),
not for an HMM detection sweep.

## A note on accessions
Pfam accessions drift between releases (e.g. DUF4011 = PF13195). The pipeline
detects domains by **name** from InterProScan/Pfam output by default and treats
accessions as optional hints in `config: prometheus.domains`. Verify any accession
you rely on against your installed Pfam release.

## Bacterial / cross-domain branch (`config: bacteria`, opt-in)
Stage your pre-computed bacterial inputs and point `bacteria:` at them (the archaeal
run ignores all of this unless `bacteria.enabled: true`):
| Key | What |
|---|---|
| `df_systems`, `df_genes` | combined bacterial DefenseFinder systems + genes tables (a genome_id column each) |
| `faa_dir`, `cds_dir` | per-genome `{genome_id}.faa` proteomes and `{genome_id}.fna` CDS |
| `taxonomy`, `qc` | GTDB classification + completeness per bacterial genome |
| `gtdb_bac_tree` | GTDB **bac120** reference tree (newick), pruned to your genomes for the co-occurrence/niche tests |
| `niche_metadata` (optional) | oxygen/biome/host labels — activates the anaerobe/niche specificity control |
| `genomes_dir` + `genome_extension` (optional) | directory of bacterial **assemblies** (`{genome_id}.<ext>`) — activates the provirus/R-loop mechanism test (geNomad + CheckV + R-loop on each assembly). Needs the `genomad`/`checkv` DBs (config `db`) and, for the phylogeny-controlled test, `gtdb_bac_tree`. |
The reused `prometheus_reference.faa` + `paralog_reference.faa` also anchor the
cross-domain tree, so the same references serve both domains.

## Intrinsic-disorder module (`config: disorder`, opt-in)
metapredict and alphaPredict install automatically (`workflow/envs/disorder.yaml`), and
segmasker comes with BLAST+. The gold-standard second predictor, **IUPred3**, is free
for academic use but license-gated, so it is NOT installable here — download it once
(register at the IUPred site), unpack it, and set `disorder.iupred3_path` to the folder
containing `iupred3.py`/`iupred3_lib.py`. With that set, the consensus is
metapredict + IUPred3 (two independent methods); without it the module falls back to
metapredict + alphaPredict and labels the result as the weaker (correlated) consensus.
