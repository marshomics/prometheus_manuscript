# Stage 6 (ecology) test fixtures

## R-loop (QmRLFS) reference data

| file | what | provenance |
|---|---|---|
| `rloop_real_viruses.fna` | 3 real virus genomes: Enterobacteria phage lambda (NC_001416.1, 48,502 bp), Sulfolobus spindle-shaped virus 1 (NC_001338.1, 15,465 bp, integrates as a provirus) and Haloarcula phage SH1 (NC_007217.1, 30,889 bp, 68 % GC) | NCBI nuccore efetch (`rettype=fasta`), 2026-10-30 |
| `rloop_reference_input_v1.fna` | 12 synthetic sequences (seeded constructs + random) | `python tests/ecology/make_rloop_fixtures.py` (`build(1)`) |
| `rloop_synthetic.fna` | the hand-derived constructs used for known-answer assertions (`build(2)`: a 'T' spacer keeps the REZ's first G out of the RIZ) | same script |
| `qmrlfs_ref_real_viruses.table.tsv` | reference QmRLFS-finder output for `rloop_real_viruses.fna` (26 RLFS: lambda 1, SSV1 0, SH1 25) | see below |
| `qmrlfs_ref_synthetic_v1.table.tsv` | reference QmRLFS-finder output for `rloop_reference_input_v1.fna` (57 RLFS) | see below |

Reference tool: `QmRLFS-finder.py` v2.0 (github.com/piroonj/QmRLFS-finder, master at commit
40b551e, 2023-12-10; file sha256 `253b91bbe58c4d0438b74c8bc7008980f366e08f23be0d7b6f563abed7902b08`),
default models `m1,m2`, run with Biopython 1.83 / Python 3.11 as

    python QmRLFS-finder.py -o <prefix> < <input.fna>      # table written to <prefix>.out.table.txt

(stdin redirection: the script's `open(..., "rU")` is invalid on Python >= 3.11). The tables
are the unmodified tool output (`#Name` = full FASTA description). sha256:
`qmrlfs_ref_real_viruses.table.tsv` 9af0aeee34712e9f1fb5a85eb053ab0cd740ef785a875312f49dd30f907cbbf0,
`qmrlfs_ref_synthetic_v1.table.tsv` ad59e200352a275b8f66e7968758d14d8815715a144f6b08038a2a2e55e75c7c.
Every RIZ/REZ sequence in `qmrlfs_ref_synthetic_v1.table.tsv` was checked to occur at its reported
coordinates in the regenerated `rloop_reference_input_v1.fna`, which proves the regenerated file is
the input the tool saw. The script itself is not vendored (its repository carries no licence text).

`workflow/scripts/py/rloop_predict.py` must reproduce both tables call-for-call (model, strand,
RIZ/REZ coordinates, lengths, G counts, %G, linker, location) and in the same order:
`test_rloop_matches_reference_qmrlfs`. The lambda result (one RLFS) agrees with van den Berg et al.
(Cell Host Microbe 2024), who report that phage lambda carries one predicted R-loop.
