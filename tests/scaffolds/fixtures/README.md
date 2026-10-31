# tests/scaffolds fixtures (provenance)

| file | source | retrieved |
|---|---|---|
| `uniprot_mcrA.faa` | UniProt REST `gene_exact:mcrA AND reviewed:true AND taxonomy_id:2157` (9 Swiss-Prot archaeal McrA) | 2026-10-30 |
| `uniprot_mcrB.faa` | same query for `mcrB` (7 sequences) | 2026-10-30 |
| `uniprot_mcrG.faa` | same query for `mcrG` (7 sequences) | 2026-10-30 |
| `kofam_ko_list_mcr_subset.tsv` | KOfam `ko_list` (https://www.genome.jp/ftp/db/kofam/ko_list.gz): header + K00399, K00401, K00402, K03421, K03422 rows | 2026-10-30 |

The suite builds small profile HMMs from these sequences at test time (MAFFT + hmmbuild,
deterministic) and runs real `hmmsearch`; thresholds used in the tests are derived from
independently parsed (Biopython SearchIO) hmmsearch scores, so the assertions test the
acceptance logic rather than particular bit scores. The ko_list subset documents the real
KOfam thresholds (K00399 McrA 741.47 full; K00401 McrB 506.80 full; K00402 McrG 359.37
full) and is parsed by the ko_list reader test.
