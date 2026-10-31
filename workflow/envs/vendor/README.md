# Vendored source packages

| File | Source | SHA-256 | Licence | Installed by |
|---|---|---|---|---|
| `corHMM_2.8.tar.gz` | CRAN, <https://cran.r-project.org/src/contrib/corHMM_2.8.tar.gz> (Date/Publication 2022-06-13; NeedsCompilation: no) | `d6e46641f5b2b15720c8bb6c51d715220519dd8a6ba44171d13ef75bfa2eb6df` | GPL (>= 2) | `workflow/envs/rstats.post-deploy.sh` |

corHMM (Beaulieu, O'Meara, Oliver & Boyko) is not packaged on conda-forge. Version 2.8
provides `fitCorrelationTest` (Boyko & Beaulieu 2023, Syst Biol 72:476,
doi:10.1093/sysbio/syac066), the hidden-rate correlation test used by
`workflow/scripts/R/methanogenesis_stats.R`. Its dependencies (Depends: ape, nloptr,
GenSA; Imports: expm, numDeriv, corpcor, MASS, nnet, phangorn, parallel, viridis, Rmpfr,
igraph, phytools) are pinned in `workflow/envs/rstats.yaml`; Snakemake runs the
post-deploy script once per env creation, which checks the SHA-256 above and runs
`R CMD INSTALL --no-docs --no-html` into the env's R library (any mismatch or error fails
the deployment).

To refresh: download the new tarball from CRAN, replace the file, and update the file
name, `EXPECTED_SHA256`/`PKG_VERSION` in the post-deploy script and this table.
