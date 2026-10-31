# Test fixtures (real data, small)

- `proA_homologs.faa` — 10 UniProtKB ProA-like proteins (DUF4011 + REase_MTES_1575;
  6 also carry DUF3320; 2 Pkinase fusions; 2 haloarchaeal ResIII/Helicase_C fusions)
  plus the `hmmemit -c` consensus of the DefenseFinder `Prometheus__ProA` profile
  (defense-finder-models, master).
- `pfam_proA_subset.hmm` — 12 Pfam-A HMMs fetched from InterPro (PF13195 DUF4011,
  PF13086 AAA_11, PF13087 AAA_12, PF18741 MTES_1575, PF11784 DUF3320, PF08696 Dna2,
  PF12705 PDDEXK_1, PF00069 Pkinase, PF13245 AAA_19, PF13604 AAA_30, PF00271
  Helicase_C, PF04851 ResIII).
- `proA_homologs.pfam_subset.domtbl` — `hmmscan --cut_ga --domtblout` of the proteins
  against that subset (HMMER 3.4). Real hmmscan output for parser/architecture tests.
- `paralogs.faa` — 12 reviewed UniProtKB/Swiss-Prot housekeeping paralogs used as the
  validation-tree outgroup in tests: DNA2 (human P51530, rat D3ZG52, mouse Q6ZQJ5, yeast
  P38859), archaeal NurA (Q8U1N8, Q97WH1, Q96YR4, Q58658) and archaeal HerA (Q97WG8,
  Q96YR7, F2Z5Z6, Q8U1P0). Headers are UniProt-style (`sp|ACC|NAME`), which also exercises
  reference-id sanitisation (IQ-TREE rewrites `|`).
