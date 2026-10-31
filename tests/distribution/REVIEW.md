# Stage 4 (distribution) review — leaf-1.1.4

Scope: `workflow/rules/distribution.smk`, `build_defense_matrix.py`, `compute_prevalence.py`,
`plot_distribution.py`, `plot_distribution_extra.py`, `association_tests.R`, `cooccurrence.R`.
Line numbers refer to the pristine copy (commit 56f5437); every owned file was rewritten.
Tests: `tests/distribution/run_tests.py` (20 tests; `bash tests/run_leaf.sh distribution`), with
`sim_data.R` (simulated trees and traits with a known generative model) and
`verify_direct.R` (independent reference computations that call the tools directly).

Primary sources used (installed versions in `prom-r` unless stated):
- phylolm 2.6.5: `phyloglm` source printed from the installed namespace and its Rd (`tools::Rd_db`).
- phytools 2.1.1 (installed) and 2.3-0 (pinned in `rstats.yaml`; `github.com/cran/phytools@2.3-0`
  `R/fitPagel.R`, `R/fitMk.R`): `fitPagel` return object and test are identical in both versions.
- caper 1.0.4 (installed) and 1.0.3 (pinned; `github.com/cran/caper@1.0.3` `R/phylo.d.R`,
  `R/comparative.data.R`): `phylo.d` identical between the two versions.
- DefenseFinder v2.0.0 (`github.com/mdmparis/defense-finder@v2.0.0`: `defense_finder_posttreat/
  df_systems.py`, `best_solution.py`, `defense_finder_cli/main.py`); defense-finder-models 2.0.0
  repository tree (`definitions/DefenseFinder/Prometheus/Prometheus.xml`, `definitions/ADF/...`).
- scipy `binomtest(...).proportion_ci(method="wilson")` as the Wilson reference.

Severity: critical = the headline or a primary result is wrong or silently missing; major = a
reported number or label is wrong, or a failure is masked; minor = robustness, clarity, runtime.

## Findings

### D1 — critical — `association_tests.R:121,166,224` — headline phyloglm never ran
- Issue: the master was read with `readr::read_tsv` (a tibble). `rownames(dat) <- dat$tip` does not
  survive on a tibble, so `phyloglm` matched no data row to any tip, dropped every tip and stopped with
  "only 0 or 1 leaf with data on all variables". Both phyloglm rows (including the confirmatory
  HEADLINE) were therefore `NA` with a `failed:` note on every input, and the script exited 0.
- Evidence: phylolm 2.6.5 source — `taxa_without_data = setdiff(phy$tip.label, rownames(mf))`,
  `phy = drop.tip(phy, taxa_without_data)`, then `stop("only 0 or 1 leaf ...")`; Rd Note: "tip labels in
  the tree are matched to the data names (row names in the data frame)". Reproduced: the pristine script
  on a simulated 120-tip dataset wrote `failed: only 0 or 1 leaf with data on all variables` for both rows.
- Fix: every table is read with base `read.delim(colClasses = "character")` into a data.frame whose
  row names are the tree tips, in tree order.
- Test: `test_association_signal_against_direct_calls` (estimates present, positive and significant
  for a simulated beta_methanogen = 1.5; n matches an independent Python recomputation).

### D2 — critical — `association_tests.R:89-91,221,258,303,327,356-365` — every failure was silent
- Issue: the whole script sat in `tryCatch(... header-only ...)` followed by `quit(status = 0)`; each
  test had its own `tryCatch` writing a `failed:` row; `requireNamespace()` turned a missing package
  into a skipped test. Violates Contract v2 item 4.
- Evidence: SCHEMA.md Contract v2 item 4; negative control: the pristine script exited 0 on the simulated data while writing `failed:` rows for both phyloglm fits and for D.
- Fix: packages are loaded with `library()` (a missing one is an error); errors propagate and Rscript
  exits non-zero. Header-only output only for a genuinely empty input (`();` tree, zero-row master),
  marked `# nonfatal:`. Data that exist but cannot support a test give that test's row with `NA` and the
  reason in `note` (e.g. "no variation in methanogen status").
- Test: `test_association_empty_not_estimable_and_failures` (unrooted tree, missing column,
  inconsistent master, no branch lengths and a bad tier set all exit non-zero; `();` gives a header).

### D3 — major — `association_tests.R:330` — D statistic always failed on real species trees
- Issue: `caper::comparative.data` stops when node labels are duplicated. The species tree comes from
  IQ-TREE with `-B 1000 -alrt 1000`, so many nodes carry the same `'100/100'` support label.
- Evidence: caper 1.0.3/1.0.4 `comparative.data`: `if (any(duplicated(c(phy$tip.label, phy$node.label))))
  stop("Labels duplicated between tips and nodes in phylogeny")`. Reproduced with a 30-tip tree.
  `root_tree.py` writes `'100/100'` labels and a `[&R]` prefix (checked with dendropy 4.6).
- Fix: node labels are dropped before all comparative fits (they carry no information for them).
- Test: `test_association_clumped_D_and_node_labels` (simulated tree with duplicated quoted labels;
  D < 0.3 and p ≤ 0.01 for a single-clade trait).

### D4 — major — `association_tests.R:131-135` — B/C-only genomes coded as absences
- Issue: presence was `prometheus_present & tier in T`, so a genome whose best locus is tier B or C
  was a confident absence in the tier-A analysis (Contract v2 item 6, rev 2).
- Evidence: SCHEMA.md Contract v2 item 6 (rev 2); pristine lines 131-135; negative control: pristine n = 95 on the tier-A simulated set, versus the independently recomputed tier-rule n.
- Fix: `tier_response()`: present if `prometheus_tier` in T, absent if `""`, excluded (dropped and
  counted as `n_excluded_tier` in every note) otherwise. `prometheus_present` must agree with
  `prometheus_tier` or the script stops.
- Test: `test_association_B_only_excluded_under_A_present_under_AB` (B-only genomes are absent from
  the tier-A set and present under A+B; counts match an independent recomputation).

### D5 — major — `association_tests.R:224-226,272-275` — a single phyloglm fit is not reliable
- Issue: `phyloglm(logistic_MPLE)` optimises the penalised likelihood with L-BFGS-B at
  `factr = 1e12` behind a hard 1e10 barrier at the alpha bounds. Its convergence code is not an
  optimality check. In a 120-dataset simulation (n = 120 tips): a default fit reported code 0 at a point
  18 log-units below the optimum (beta = -0.47, p = 0.59 vs beta = 4.98, p = 0.0015 at the optimum), and
  code-52 fits stopped with beta far from optimal (null dataset: beta = 1.69, p = 2.5e-4, a false positive;
  at the best optimum found, beta = 0.72, p = 0.093). Warnings (alpha at a bound, btol reached) went only
  to the log.
- Evidence: phylolm 2.6.5 source (`npllh` returns 1e10 outside the bounds; `optim(..., method =
  "L-BFGS-B", control = list(factr = 1e+12))`; `penlogLik = logLik + log(det(info))/2`, the MPLE
  objective). Simulation study: plain MPLE had 0/58 null rejections and power 0.78; multi-start MPLE
  power 0.82 and always reached a penalised likelihood at least as high as IG10's (minimum gain +0.005);
  logistic_IG10 power 0.65 and failed its own convergence criterion in 48% of null datasets.
- Fix: `fit_phyloglm_ms()` fits logistic_MPLE from four starts (alpha*T = exp(-3), 1, exp(3), and the
  logistic_IG10 solution), polishes each with warm restarts until the gain is < 1e-4 log-units, and
  reports the fit with the highest penalised likelihood. The note records the winning start, the gap to
  the runner-up, the optim code, alpha*T, and alpha-bound and btol flags; a fit still improving at the
  restart limit is flagged CAUTION. The method stays logistic_MPLE (phylolm's default, Firth-penalised,
  finite under separation).
- Test: `test_association_signal_against_direct_calls` (multi-start penlogLik ≥ the default single fit
  computed independently; estimate equals it when that fit converged at the same optimum).

### D6 — major — `association_tests.R:226,238-241` — bootstrap CI built on unchecked refits
- Issue: phyloglm's own `boot` refits every replicate once, starting from the original start values, with
  no convergence check (see D5). The percentile CI therefore mixed unconverged fits.
- Evidence: phylolm 2.6.5 source (`boot_model` calls `optim(...)` once per replicate inside `try`).
- Fix: `boot_phyloglm()` simulates `--boot` replicates under the fitted model with
  `phylolm::rbinTrait` (as phyloglm does) and refits each from the estimate with the same warm-restart
  polishing. Replicates are run with `parallel::mclapply` (`--threads`; deterministic because the data are
  drawn up front from a fixed seed). The number fitted is reported; with fewer than max(10, B/2) fitted the
  CI is withheld. With `--boot 0` the CI is Wald and is labelled so.
- Test: `test_association_signal_against_direct_calls` (note "replicates fitted"; CI brackets the
  estimate); `test_association_complete_separation` (`--boot 0` gives a Wald CI label).

### D7 — major — `association_tests.R:197-215` — naive baseline wrong under separation; no Fisher test
- Issue: with a zero cell in the 2x2 table the Wald GLM reported 25.36 ± 7475 (p = 0.997) and a profile CI
  of (-695.9, NA) for an association whose exact p is about 2e-9. The design's Fisher exact baseline (4.3)
  was missing.
- Evidence: pristine run on simulated data (2x2 = 113/0 vs 1/6); design §8, 4.3.
- Fix: `naive_glm` reports the likelihood-ratio p (valid under separation) and the profile CI; when a cell is
  zero, the estimate is `NA` and the note explains why. A new `naive_fisher` row reports the log conditional-MLE
  odds ratio, exact CI and exact p. Both rows contain "naive", so the registry keeps them exploratory.
- Test: `test_association_complete_separation`; `test_association_signal_against_direct_calls` (both rows
  equal direct `glm`/`fisher.test` computations).

### D8 — major — `association_tests.R:112-113`, `cooccurrence.R:87-88` — silent re-rooting
- Issue: an unrooted tree was passed to `phangorn::midpoint` inside `tryCatch`. phangorn was never loaded
  (`requireNamespace` absent), so any failure silently kept the unrooted tree. All four methods are
  root-dependent, and rooting belongs to `root_tree.py` (README §4).
- Evidence: pristine lines 112-113 (`tryCatch(phangorn::midpoint(tree), error = function(e) tree)`); README §4 ("species tree is rooted before any comparative or ancestral use").
- Fix: an unrooted tree (basal polytomy) is an error that names `root_tree.py`. Other checks that stop the
  script: missing or negative branch lengths, duplicated tips, a multi-tree file.
- Test: `test_association_empty_not_estimable_and_failures`;
  `test_cooccurrence_bacterial_cli_empty_and_failures` (unrooted → non-zero exit, "not rooted").

### D9 — minor — `association_tests.R:181-191,271-275` — covariate coding and a biased fallback
- Issue: completeness entered raw (0-100) and genome size scaled but not log-transformed; the defense count
  entered raw. When `n_defense_systems_excl_prometheus` was absent, the script fell back to the total count,
  which includes the response.
- Evidence: pristine lines 181-191 (raw covariates; total-count fallback) and 271-275 (`scale(genome_size)` only); SCHEMA association_tests note (headline adjusts for defense load EXCLUDING Prometheus).
- Fix: z(completeness), z(log genome size), z(log1p(defense load excluding Prometheus)), on complete cases
  (the number dropped is reported). Zero-variance covariates are dropped and named in the note. The
  excl-Prometheus column is required. `detect_defensefinder` is deliberately not a covariate: it is 1
  exactly when DefenseFinder called Prometheus, i.e. close to the response.
- Test: `test_association_empty_not_estimable_and_failures` (missing column → exit 1); signal test (HEADLINE
  note names the covariates).

### D10 — minor — `association_tests.R:333-339` — D p-values: only in a note, and can be 0
- Issue: `phylo.d`'s `Pval1 = sum(random < obs)/permut` (strict inequality, no +1) can be exactly 0, and it
  appeared only in `note`, so the registry never saw the test.
- Evidence: caper 1.0.3 `phylo.d` source (`soccpval1 <- sum(ransocc < obssocc) / permut`).
- Fix: `p_value = (#{D_random ≤ D_obs} + 1)/(N + 1)` (test of D < 1, i.e. phylogenetic signal). The
  Brownian-model p (same correction) and caper's raw values are in `note`. `--d-permutations` defaults to 1000.
- Test: signal test (equality with a direct seeded `phylo.d`; p lies on the (k+1)/(N+1) lattice).

### D11 — minor — `association_tests.R:309-312`, `cooccurrence.R:186-195` — fitPagel optimiser failures
- Issue: `lik.ratio = 2(logL_dep - logL_indep)` and `P = pchisq(LR, df = 4)` are taken as given. A dependent
  fit stuck below the nested independent model gives LR < 0 and P = 1 silently. An under-optimised
  independent fit inflates LR.
- Evidence: phytools 2.1.1 and 2.3-0 `fitPagel` source (identical); df = k.dQ - k.iQ = 8 - 4 for ARD.
- Fix: `pagel_test()`. The independent model is the Kronecker sum of two 2-state ARD chains, so its MLE equals
  the sum of two separate 2-state `fitMk` fits (verified numerically). That sum guards the independent logL.
  A dependent fit below the independent one is refitted from the independent optimum
  (`pagel_refit_dependent`), and LR is truncated at 0 with a flag if needed.
- Test: `test_pagel_guard_helpers` (LR equals `fitPagel`'s; refit ≥ independent logL; tier helper);
  signal test (factorisation equality to 1e-3).

### D12 — minor — `association_tests.R:87`, `cooccurrence.R:68` — RNG and locale
- Issue: one `set.seed` at start, so bootstrap and permutation draws depended on which earlier tests ran.
  The session locale leaked German warnings into notes.
- Evidence: phylolm 2.6.5 `plogregfunct` calls `rnorm()` for IG10 restarts; caper `phylo.d` uses `sample()`/`rmvnorm()`; pristine notes contained German warnings (e.g. 'Algorithmus konvergierte nicht').
- Fix: seeds are set immediately before each stochastic step (bootstrap draws, `phylo.d`, IG10 starts; per
  system in co-occurrence); `Sys.setenv(LANGUAGE = "en")`.
- Test: `test_association_reproducible` (byte-identical output across runs and thread counts, bootstrap
  included); `test_cooccurrence_thread_count_does_not_change_results`; D equality with an independently
  seeded `phylo.d`.

### D13 — minor — `association_tests.R:121`, `cooccurrence.R:95,106` — numeric-looking ids
- Issue: readr type guessing turns zero-padded ids such as `000123` into numbers that no longer match tree
  tips.
- Evidence: readr guesses column types (zero-padded digits parse as double); the simulated ids `000001`... exercise this, and the counts match an independent recomputation in every test.
- Fix: all-character reads. Simulated genome ids are `000001`...
- Test: every association and co-occurrence test (ids are zero-padded).

### D14 — minor — `association_tests.R` (new) — zero-length branches
- Issue: `phylo.d` stops on any zero-length branch ("cannot currently simulate"), and the Brownian VCV is
  singular.
- Evidence: caper `phylo.d` source.
- Fix: exact zeros are set to 1e-6 x tree height (logged); this is negligible for every model used here.
- Test: `test_association_zero_length_branch`.

### C1 — critical — `cooccurrence.R:162,197-201,70,237` — silent method substitution
- Issue: when `fitPagel` failed, or phytools was missing, the Fisher p was used with method
  `fisher(fallback)` and entered the registry's co-occurrence family next to Pagel p-values. Any other error
  wrote a header-only file and exited 0.
- Evidence: pristine lines 162 (`requireNamespace("phytools")`), 197-201 (`method <- "fisher(fallback)"`), 70 and 237 (blanket `tryCatch` and `quit(status = 0)`); results_registry.py `normalize_cooccurrence` takes `p_value` of every row regardless of method.
- Fix: phytools is loaded with `library()`; fitPagel is the only `p_value`, never replaced; errors propagate.
  `mclapply` worker failures stop the script.
- Test: `test_cooccurrence_known_answer_against_direct` (every row is `fitPagel`; LR and p equal direct
  `phytools::fitPagel` calls to 1e-3 / 1e-4).

### C2 — major — `cooccurrence.R:50-56,144` — focal Prometheus vector was the raw DefenseFinder call
- Issue: Prometheus presence came from the defense-matrix column (any DefenseFinder call, including loci
  later Rejected or tiered B/C), over all genomes rather than pass_primary ones. This was inconsistent with
  the association tests and Contract v2 item 6.
- Evidence: SCHEMA.md Contract v2 item 6; negative control on the same simulated data: pristine n_a = 9 (raw DefenseFinder vector, B/C/Rejected included) versus the tier-A count used by the new script.
- Fix: presence from the master's `prometheus_tier` under the new optional `--primary-tiers` (default A);
  pass_primary genomes on the tree only. Excluded genomes are counted in `note`. Without a master (or with
  no `prometheus_tier`), the raw matrix column is used and this is logged.
- Test: `test_cooccurrence_known_answer_against_direct` (n, n_a and n_excluded match an independent
  recomputation).

### C3 — major — `cooccurrence.R:175-182` — perfect co-occurrence labelled "ambiguous"
- Issue: `direction` was computed only for finite OR, so OR = Inf (every Prometheus genome also carries the
  system: the strongest co-occurrence) became `ambiguous (OR=Inf)`.
- Evidence: pristine lines 175-182 (`if (is.finite(orval))`); negative control on the simulated data: CoSys (n_both = n_a) labelled `ambiguous (OR=Inf)`.
- Fix: `direction_of()`: OR > 1 (including Inf) → co-occurrence, OR < 1 (including 0) → exclusion. The label
  format is unchanged, and new `odds_ratio`/`or_ci_*`/`fisher_p` columns carry the numbers.
- Test: `test_cooccurrence_known_answer_against_direct` (CoSys → co-occurrence, ExSys → exclusion).

### C4 — major — `cooccurrence.R:190-195` — runtime
- Issue: phytools' default `lik.func = "lik"`, run serially, took ~20 s per system for 120 tips and scales
  linearly (500 tips: 52 s). With ~100 testable systems on thousands of tips this is hours.
- Evidence: benchmark on 500 tips: `lik` 52.0 s, `pruning` 34.9 s, identical logL to 6 decimals.
  `method = "ace"` took 2.0 s with the same LR, but emits "imaginary parts discarded" and "model fit
  suspicious" warnings (eigendecomposition of a non-symmetric Q), so the `expm`-based `fitMk` is kept.
- Fix: `lik.func = "pruning"`; systems run in parallel (`--threads`, `mclapply`) with a fixed seed per system;
  BLAS pinned to 1 thread in the rule.
- Test: `test_cooccurrence_thread_count_does_not_change_results`.

### C5 — minor — `cooccurrence.R` (design 4.6) — no genome-size correction
- Issue: the design requires co-occurrence "correcting for genome size". Fisher direction and Pagel ignore it.
- Evidence: design §8 step 4.6 ("correcting for genome size and FDR across the system pairs").
- Fix: supplementary columns `phyloglm_logodds`, `phyloglm_se`, `phyloglm_p` from multi-start
  `phyloglm(Prometheus ~ system + z(log genome size))`: a phylogenetically controlled, size-adjusted signed
  effect. When the master has no genome_size or the model cannot be fitted for a system, they are NA and
  `note` gives the reason. `p_value` stays the Pagel test.
- Test: `test_cooccurrence_known_answer_against_direct` (columns filled).

### C6 — minor — `cooccurrence.R:169` — prevalence window
- Issue: the window [3, n-3] admits systems seen in only 3 genomes. With so few occurrences the χ²(4)
  approximation of the LR is doubtful and a single gain in one clade can drive the test. Such systems add
  little power and inflate the BH family.
- Evidence: design §8 step 4.6; pristine line 169 (`if (nb < 3 || nb > (n - 3)) next`); phytools `fitPagel` tests LR against `pchisq(df = 4)`.
- Fix (judgement call): `--min-count` (default 5) for present and absent genomes, applied to Prometheus as
  well. The number of skipped systems is logged.
- Test: `test_cooccurrence_known_answer_against_direct` (the 2-genome RareSys is not tested).

### C7 — minor — `cooccurrence.R:139-143,220-222` — strictness
- Issue: `bin()` turned NA or garbage into 0, and an unknown `--fdr` reached `p.adjust` inside the blanket
  `tryCatch` (header-only).
- Evidence: pristine lines 139-143 (`bin()` maps NA to 0) and 220-222 (`p.adjust(method = opt$fdr)` inside the blanket `tryCatch`).
- Fix: the matrix must be strictly 0/1; `--fdr` is validated (bh/fdr/by/bonferroni/holm/hochberg/hommel/none);
  matrix genomes missing from the master are an error.
- Test: `test_cooccurrence_bacterial_cli_empty_and_failures`.

### M1 — major — `build_defense_matrix.py:22,43` — numeric-looking ids lost their systems
- Issue: `pd.read_csv` without `dtype=str` turns `000123` into `123`. The per-genome path is then wrong, the
  file "does not exist", and the genome silently gets all zeros.
- Evidence: pandas `read_csv` infers numeric dtype for digit-only ids; pristine run of the suite's fixture failed (`test_defense_matrix_known_answer`, negative control).
- Fix: `dtype=str, keep_default_na=False` throughout.
- Test: `test_defense_matrix_known_answer` (genome `000123`).

### M2 — major — `build_defense_matrix.py:21-25,41-50` — swallowed failures
- Issue: an unreadable master gave a header-only output; an unreadable or 0-byte systems table hit
  `except: pass` and gave zeros; a missing table gave zeros.
- Evidence: pristine lines 21-25 (`except Exception: write_empty`) and 41-50 (`except Exception: pass`); negative control: the missing-table case exited 0.
- Fix: read errors and a missing `type` column exit 1. New optional `--require-all` (passed by the archaeal
  rule) makes a missing table an error. Without it (the bacterial branch, whose splitter writes tables only for
  genomes with systems), a missing table means zero systems and is counted in the log.
- Test: `test_defense_matrix_missing_empty_and_failures`.

### M3 — minor — `build_defense_matrix.py:44-48` — anti-defense systems counted as defense
- Issue: with `defense-finder run -a`, or user-supplied bacterial tables, AntiDefenseFinder systems appear
  as rows with `activity = Antidefense` and became matrix columns.
- Evidence: DefenseFinder v2.0.0 `best_solution.py` (`activity = 'Antidefense'` when `model_fqn` contains
  "ADF"); defense-finder-models 2.0.0 tree (`definitions/ADF/...`).
- Fix: rows with activity "Antidefense" are excluded (counted in the log). Tables without an `activity`
  column (the prometheus.smk placeholder header) are read unchanged.
- Test: `test_defense_matrix_known_answer` (`Anti_CBASS` absent).

### P1 — major — `compute_prevalence.py:81-90` — B/C-only genomes were confident absences
- Issue: the denominator counted every pass_primary genome, so B/C-only genomes were absences in the tier-A
  prevalence (Contract v2 item 6).
- Evidence: SCHEMA.md Contract v2 item 6 (rev 2); pristine `is_present` (lines 81-90) codes B/C-only genomes as 0 in the denominator.
- Fix: they are excluded from that tier set's denominator and counted in the new `n_excluded` column.
  `n_genomes` is now the analysed denominator (present + absent); the new `n_genomes_total` keeps every
  pass_primary genome of the taxon. `median_completeness` is computed over the denominator genomes.
- Test: `test_prevalence_tier_rule_known_answer` (P1: n = 8, present = 3, excluded = 2, total = 10; Wilson CI
  equal to scipy's).

### P2 — major — `compute_prevalence.py` (design 4.1) — only one tier set reported
- Issue: the design requires prevalence "at tier A and tier A+B".
- Evidence: design §8 step 4.1 ("report it at tier A and tier A+B").
- Fix: additive `*_sensitivity` columns for `--sensitivity-tiers` (optional, default A,B; rule key
  `stats.confidence_tiers_sensitivity`).
- Test: `test_prevalence_tier_rule_known_answer` (A+B: n = 9, present = 4, excluded = 1, so the B-only
  genome is present).

### P3 — major — `compute_prevalence.py:51-78` — silent failures and silent denominators
- Issue: an unreadable master and a missing `prometheus_present` gave a header-only output; a missing
  `pass_primary` silently used every genome as the denominator; unknown tiers counted as absent; a
  disagreement between `prometheus_present` and `prometheus_tier` was ignored.
- Evidence: pristine lines 51-78; negative control: `test_prevalence_empty_and_failures` got exit 0 where a failure was required.
- Fix: required columns, strict booleans, known tiers and present/tier agreement are enforced (exit 1).
  Only a zero-row master gives a header-only output.
- Test: `test_prevalence_empty_and_failures`.

### P4 — minor — `compute_prevalence.py:25-35` — Wilson edge cases
- Issue: n = 0 returned (0, 0), a CI for an undefined proportion; z was rounded to 1.96.
- Evidence: scipy `binomtest(...).proportion_ci(method="wilson")`; negative control: the pristine `wilson_ci` differed by 8e-6 at k = 0, n = 7 (z = 1.96 rounding), and returned (0, 0) for n = 0.
- Fix: NaN for n = 0; exact z = Φ⁻¹(0.975); exact endpoints at k = 0 and k = n. The pristine formula was
  otherwise correct: it matches scipy's Wilson interval once z is exact.
- Test: `test_wilson_matches_scipy` (matches `scipy.stats.binomtest(...).proportion_ci("wilson")` to 1e-12).

### P5 — minor — `compute_prevalence.py` — no overall row; unassigned genomes dropped unannounced
- Evidence: design §8 step 4.1; the pristine output had no overall row, which forced plot_distribution_extra to recompute it (PE1).
- Fix: a `rank = overall, taxon = all` row (all pass_primary genomes). The number of genomes unassigned at
  each rank is logged.
- Test: `test_prevalence_tier_rule_known_answer` (overall row; genus rows omit genus-less genomes).

### PL1 — major — `plot_distribution.py:169` — exclusion drawn as co-occurrence
- Issue: the heatmap sign required `direction` to be exactly `exclusion`, but cooccurrence.R writes
  `exclusion (OR=0.2)`, so every system was drawn with a positive sign.
- Evidence: pristine line 169 compares the whole label to 'exclusion'; cooccurrence.R writes `exclusion (OR=...)` (pristine cooccurrence.R line 205).
- Fix: `direction_sign()` matches the label prefix; neutral and ambiguous are 0.
- Test: `test_plots_real_inputs_placeholders_and_helpers`.

### PL2 — minor — `plot_distribution.py:100-133` — the forest mixed scales
- Issue: every row with an estimate was drawn on one "log-odds" axis, including D (a different statistic) and
  the separated GLM (25 ± 7475).
- Evidence: pristine lines 100-133 plot every finite estimate; the pristine naive GLM row on separated data was 25.36 with SE 7475.
- Fix: only the four methanogen log-odds rows are drawn (naive grey, phylogenetic blue, HEADLINE larger);
  Pagel and D are printed as text. A missing CI falls back to ±1.96 SE.
- Test: `test_plots_real_inputs_placeholders_and_helpers`.

### PL3 — minor — `plot_distribution.py:45-48,95-98,153-156`, `plot_distribution_extra.py` — placeholders masked interface errors; formats hard-coded
- Issue: a table missing its required columns produced a "no data" placeholder. Figures were always PNG+SVG
  while the rule declares `PLOT_FORMATS`.
- Evidence: pristine plot_distribution.py lines 45-48, 95-98 and 153-156 (placeholders on missing columns); `plotting_style.save_fig(formats=("png","svg"))` default while the rule expands `PLOT_FORMATS`.
- Fix: missing columns exit 1; placeholders only for header-only inputs; optional `--formats` (the
  distribution rule passes `PLOT_FORMATS`).
- Test: `test_plots_real_inputs_placeholders_and_helpers`.

### PE1 — major — `plot_distribution_extra.py:47-52` — enrichment flags compared mismatched tier sets
- Issue: the "overall" prevalence counted any validated tier (A/B/C) via `prometheus_present`, while the
  per-taxon prevalence was tier A. Taxa were therefore biased toward "depleted".
- Evidence: pristine plot_distribution_extra.py line 52 (overall = mean of `prometheus_present` over pass_primary genomes, any tier) versus compute_prevalence tier-set prevalence.
- Fix: the overall row of the prevalence table (same tier set and denominator). For older tables, it is
  recomputed from the master with the tier rule. Flags are labelled descriptive (no multiplicity or phylogeny
  control).
- Test: `test_plots_real_inputs_placeholders_and_helpers` (`overall_prevalence` = 1/3 with a B-only genome
  excluded).

### PE2 — minor — `plot_distribution_extra.py:101-118` — strongest results dropped from the OR figure
- Issue: the OR was parsed from the label, so OR = Inf and OR = 0 (perfect co-occurrence or exclusion) were
  filtered out.
- Evidence: pristine lines 117-118 drop non-finite parsed ORs; pristine cooccurrence wrote `OR=Inf` and `OR=0` for perfect co-occurrence or exclusion.
- Fix: log2 OR from the 2x2 counts, with the Haldane–Anscombe +0.5 when a cell is 0; falls back to label
  parsing for old tables.
- Test: `test_plots_real_inputs_placeholders_and_helpers` (known values).

### R1 — major — `distribution.smk:5,27,105` — rules without logs
- Issue: `defense_matrix`, `compute_prevalence` and `plot_distribution` had no `log:`. Output went to the
  terminal and was lost on a cluster.
- Evidence: `tests/check_rules.py` on the pristine rule file: "rule defense_matrix / compute_prevalence / plot_distribution has no log: directive" (RULES HYGIENE FAILED, 3 findings).
- Fix: all five rules log stdout and stderr; new flags are read through `cfg()` with defaults; `threads` are
  used (`--threads`); BLAS is pinned to 1 thread per worker.
- Test: gate G2 (`check_rules.py` → RULES HYGIENE OK); `test_rules_hygiene_and_logs`.

## Verified correct, no change
- DefenseFinder v2.0.0 systems table columns: `sys_id, type, subtype, activity, sys_beg, sys_end,
  protein_in_syst, genes_count, name_of_profiles_in_sys`. With no hits it is header-only, and it is named
  `{input basename}_defense_finder_systems.tsv` (`df_systems.py`, `main.py`), matching the rule's paths.
- `type` is the model family (`model_fqn.split('/')[-2]`); Prometheus is
  `definitions/DefenseFinder/Prometheus/Prometheus.xml`, so DefenseFinder reports type "Prometheus".
  Family level (`type`) is the right unit for co-occurrence: subtypes of one family are near-exclusive variants.
- `phyloglm` summary columns (`Estimate, StdErr, z.value, [lowerbootCI, upperbootCI], p.value`): p is Wald,
  conditional on alpha; `bootconfint95` is 2 x (p+1) with column names = coefficients + "alpha" (phylolm 2.6.5).
- `fitPagel` (2.1.1 = 2.3-0): returns `independent.logL, dependent.logL, lik.ratio, P`; ARD df = 4; flat root
  prior (`pi = "equal"`); deterministic start values. `lik.func = "pruning"` gives identical likelihoods.
- caper `phylo.d` (1.0.3 = 1.0.4): `binvar`/`names.col` use non-standard evaluation (unquoted names work
  inside functions); D = (obs - mean_Brownian)/(mean_random - mean_Brownian).
- ape 5.8 `read.tree` parses dendropy's `[&R]` prefix; the result is rooted and binary.
- `rbinTrait(n > 1)` returns a tips x n matrix with tip-label row names in tree order.
- Wilson formula in the pristine `compute_prevalence.py` is correct: with exact z it matches `scipy.stats.binomtest(...).proportion_ci(method="wilson")` to 1e-12 (`test_wilson_matches_scipy`).
- Snakemake 8.30 dry runs (private copy, `.test` config): `all` = 102 jobs and `crossdomain` = 86 jobs both
  resolve. The rendered commands carry the new flags; the bacterial rules still call the scripts with their
  unchanged CLIs.

## Judgement calls (documented defaults)
- ANME genomes keep the master's coding (methanogen = False, as `classify_methanogen.py` defines it); their
  number is reported in every note.
- Co-occurrence minimum count 5 (configurable, `stats.cooccurrence_min_count`).
- Zero-length branches are set to 1e-6 x tree height rather than collapsed.
- D's registry p is the test against random (phylogenetic signal); the Brownian test is in `note`.
