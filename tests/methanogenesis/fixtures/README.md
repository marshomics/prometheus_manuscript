# tests/methanogenesis fixtures (provenance)

| file | source | retrieved |
|---|---|---|
| `uniprot_mttB.faa` | UniProtKB REST `gene:mttB AND taxonomy_id:2157 AND reviewed:true` (4 Swiss-Prot trimethylamine methyltransferases) | 2026-10-31 |
| `uniprot_mtbA.faa` | UniProtKB REST `gene:mtbA AND taxonomy_id:2157 AND reviewed:true` (3 Swiss-Prot methylamine-corrinoid:CoM methyltransferases) | 2026-10-31 |
| `uniprot_mtaB.faa` | UniProtKB REST `gene:mtaB AND taxonomy_id:2157` (first 8; one gene-name collision, a tRNA methylthiotransferase, removed) | 2026-10-31 |
| `uniprot_mtrA.faa` | UniProtKB REST `protein_name:"tetrahydromethanopterin S-methyltransferase" AND taxonomy_id:2157 AND reviewed:true`, subunit-A records, first 8 | 2026-10-31 |
| `uniprot_mtrH.faa` | the same query, subunit-H records, first 8 | 2026-10-31 |

Every Swiss-Prot MttB carries pyrrolysine (`O`) at position 334 (335 in *Methanococcoides
burtonii* Q12TR2), i.e. a single in-frame UAG in the gene. The suite back-translates one of these
proteins with a real TAG at the O position, runs Prodigal's truncation by construction (the gene
record stops at the TAG) and checks that `methanogenesis_search.py` recovers the full-length
protein by read-through — the fixture therefore tests the pyrrolysine logic against a real
pyrrolysine protein rather than a synthetic one.

Profile HMMs are built at test time (MAFFT + hmmbuild, deterministic: `--threadit 0`), and the
thresholds used in the assertions are derived from independently parsed hmmsearch scores, so the
tests exercise the acceptance and pathway logic rather than particular bit scores.

| `kegg_reference_ko_sets.json` | KEGG REST `/link/ko/<org>` for 12 cultured methanogens (msi, mmp, mja, mth, mhu, mcj, mba, mac, mbu, mst, max, mhz); `{org: {KO: n_genes}}` | 2026-10-31 |

The reference-genome test uses these KEGG assignments as an INDEPENDENT annotation: it checks
the pathway rules, not the HMM search (the search is checked separately on real HMMER runs).
