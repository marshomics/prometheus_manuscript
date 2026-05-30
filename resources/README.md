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
