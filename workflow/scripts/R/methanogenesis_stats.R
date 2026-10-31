#!/usr/bin/env Rscript
# =====================================================================
# Design 6.2 -- within methanogens, does methanogenesis TYPE predict Prometheus
# presence, with type and habitat fitted jointly (partial effects)?
# SCHEMA ecology/methanogenesis_associations.tsv; registry family `methanogenesis`
# (always exploratory).
#
# Predictors (--predictors): the gene-based pathway-capacity indicators of the master
# table (mgen_co2_reducing, mgen_acetoclastic, mgen_methylotrophic,
# mgen_methyl_dismutating, mgen_methyl_reducing; '1' / '0' / 'NA' = unresolved).
#
# Analysis set: master rows whose tree_tip is a tip of the ROOTED species tree, with
# methanogen = True, not ANME (is_anme, when present) and pass_primary (when present;
# the same genome universe as association_tests.R / ecology_stats.R).
# Response (SCHEMA Contract v2 item 6, as ecology_stats.R): for the tier set T
# (--primary-tiers) present (1) = best validated tier in T; absent (0) = no validated
# locus (tier ""); EXCLUDED (dropped, counted in n_excluded_tier) = a validated locus
# whose best tier lies outside T. Rows with an unresolved predictor ('NA') are dropped
# and counted (n_unresolved); rows with a missing covariate likewise (n_missing_covariate).
# Covariates (--covariates, default completeness,genome_size) are z-scored within every
# model (estimates are invariant to their units); a constant covariate is not adjusted.
#
# Rows per predictor (columns `analysis` / `role`), in this order:
#  marginal / secondary
#     phylolm::phyloglm(present ~ x + z(covariates), logistic_MPLE): Firth-penalised
#     phylogenetic logistic regression (Ives & Garland 2010 Syst Biol 59:9,
#     doi:10.1093/sysbio/syp074) fitted with the multi-start MPLE of association_tests.R
#     (fit_phyloglm_ms via run_phyloglm). Wald z p-value and 95 % CI conditional on
#     alpha; estimates are log-odds. The VIF of x on the covariates is reported (no guard).
#  habitat_adjusted / primary  (the partial effect of design 6.2)
#     present ~ x + z(covariates) + habitat terms. Candidate habitat terms (--habitat;
#     auto = every ^(group|niche)_ column) are binary columns with >= --min-class genomes
#     in each level among the complete cases; they are added in order of |phi(x, h)|
#     (the most confounded habitat first) while the events-per-variable rule holds:
#     n_events = min(n_present, n_absent) and total model terms (x + covariates + habitat)
#     <= floor(n_events / --epv-min) (Peduzzi et al. 1996 J Clin Epidemiol 49:1373,
#     doi:10.1016/S0895-4356(96)00236-3). A habitat column with NA cells is admitted only
#     if its complete-case subset still meets the EPV rule and the class minima (else it
#     is skipped and named in the note). VIF(x) = 1 / (1 - R^2) of the linear regression
#     of x on the other model terms; VIF > --vif-max -> the row is reported unfitted with
#     p_value NA ("partial effect not identifiable"). 'covariates' lists every adjustment
#     term, 'habitat_terms' the habitat terms alone, 'vif' the VIF. The selection never
#     looks at the x-response association (|phi| is between x and the habitat; the EPV
#     budget uses only the response counts), so it cannot steer the p-value. A habitat
#     column redundant with the terms already included (R^2 = 1) is skipped.
#  simmap / descriptive  (replication guard)
#     Maddison & FitzJohn 2015 (Syst Biol 64:127, doi:10.1093/sysbio/syu070) show that
#     phylogenetic correlation tests are falsely significant when a trait arose once: the
#     effective sample size is the number of independent origins, not the number of tips.
#     The history of x on the species tree pruned to the analysed tips is mapped: Mk
#     models ER and ARD (phytools::fitMk, flat root prior), the one with the lower AIC
#     kept, then phytools::make.simmap(nsim = --nsim) with Q fixed at that model's MLE
#     (phytools' "empirical" Q) under a fixed seed. Changes per stochastic map are
#     counted (n_changes_mean, 2.5 % / 97.5 % quantiles n_changes_lo / n_changes_hi);
#     n_changes_parsimony is the Fitch-Hartigan minimum, a model-free lower bound that
#     every map attains. replicated = n_changes_mean >= --min-changes. With k independent changes a fully
#     concordant pattern has sign-test p = 2^-k; k = 5 is the smallest k with
#     2^-k < 0.05, hence the default 5. When replicated is FALSE (or the history could not
#     be mapped) the regression rows of that predictor (marginal, habitat_adjusted and the
#     sensitivity rows) keep estimate and CI but get p_value NA, a 'replication guard: ...'
#     note, and stay out of the FDR family.
#  corhmm / descriptive
#     corHMM::fitCorrelationTest (Boyko & Beaulieu 2023 Syst Biol 72:476,
#     doi:10.1093/sysbio/syac066): independent vs correlated evolution of (x, present),
#     each with one or two rate classes (the hidden-rate independent model is the null
#     that removes the single-origin false positives of Pagel's test). delta_aic =
#     min(AIC independent, AIC hidden-rate independent) - min(AIC correlated, AIC
#     hidden-rate correlated): > 0 favours correlated evolution. No p-value. The four
#     AICs and Akaike weights are in aic_* / aicw_*; best_model = lowest AIC.
#  sens_unresolved_as_absent / sens_complete_only / sens_leave_clade_out / sens_tiers_AB
#     role sensitivity (never in the FDR family): the habitat_adjusted procedure
#     (eligibility, |phi| order, EPV, VIF guard) re-applied with unresolved x coded 0;
#     on genomes with completeness >= --complete-threshold; without the gtdb_order that
#     carries the most x = 1 tips (named in the note; ties broken alphabetically); and
#     with the response defined by --sensitivity-tiers. A sensitivity data set identical
#     to the primary one (e.g. no unresolved genome) repeats the primary values (noted).
# p_adj: --fdr-method over the primary (habitat_adjusted) rows that carry a p-value; NA
# elsewhere. results_registry.py applies the same rule to its methanogenesis family.
# A row is not fitted (written with the reason in `note`) when n < --min-n or a response
# or predictor class has fewer than --min-class genomes (simmap: predictor classes only).
# Output: the SCHEMA columns above in order, then (additive) tier_set, n_excluded_tier,
# n_unresolved, n_missing_covariate, n_events, n_terms, epv, habitat_terms,
# n_habitat_eligible, alpha, mk_model, n_changes_parsimony, n_gains_mean, n_losses_mean,
# aic_* and aicw_* (independent, correlated, hmm_independent, hmm_correlated).
# --threads analyses predictors in parallel; every random step (IG10 start, stochastic
# maps, corHMM) is seeded per predictor, so the output does not depend on it.
#
# Failure policy (SCHEMA Contract v2 item 4): an empty species tree ("();"), a master
# table without rows, or no methanogen tree tip is a genuinely empty input -> header-only
# output, exit 0. A missing input, a missing required column, a malformed value, an
# unrooted tree, a missing R package (corHMM when --corhmm TRUE), a crashed Mk / stochastic
# map / corHMM fit, or every attempted phyloglm fit failing -> error (non-zero exit).
# =====================================================================

Sys.setenv(LANGUAGE = "en")
suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phylolm)
  library(phytools)
})

DEFAULT_PREDICTORS <- c("mgen_co2_reducing", "mgen_acetoclastic", "mgen_methylotrophic",
                        "mgen_methyl_dismutating", "mgen_methyl_reducing")
OUT_COLS <- c("predictor", "analysis", "role", "estimate", "std_error", "ci_low", "ci_high",
              "p_value", "p_adj", "n", "n_present", "n_absent", "n_level1", "n_changes_mean",
              "n_changes_lo", "n_changes_hi", "replicated", "delta_aic", "best_model",
              "covariates", "vif", "model", "converged", "note")
EXTRA_COLS <- c("tier_set", "n_excluded_tier", "n_unresolved", "n_missing_covariate",
                "n_events", "n_terms", "epv", "habitat_terms", "n_habitat_eligible", "alpha",
                "mk_model", "n_changes_parsimony", "n_gains_mean", "n_losses_mean", "aic_independent",
                "aic_correlated", "aic_hmm_independent", "aic_hmm_correlated",
                "aicw_independent", "aicw_correlated", "aicw_hmm_independent",
                "aicw_hmm_correlated")
ALL_COLS <- c(OUT_COLS, EXTRA_COLS)
ANALYSES <- c(marginal = "secondary", habitat_adjusted = "primary", simmap = "descriptive",
              corhmm = "descriptive", sens_unresolved_as_absent = "sensitivity",
              sens_complete_only = "sensitivity", sens_leave_clade_out = "sensitivity",
              sens_tiers_AB = "sensitivity")
CORHMM_MODELS <- c(independent = "independent_model_fit",
                   correlated = "correlated_model_fit",
                   hmm_independent = "hidden_Markov_independent_model_fit",
                   hmm_correlated = "hidden_Markov_correlated_model_fit")
TRUE_S <- c("1", "true", "t", "yes", "y", "1.0")
FALSE_S <- c("0", "false", "f", "no", "n", "0.0")
NA_S <- c("", "na", "nan", "none", "null")

say <- function(...) message("[methanogenesis_stats] ", ...)
clean <- function(x) gsub("[\t\r\n]+", " ", x)
join_notes <- function(...) {
  v <- unlist(list(...))
  v <- v[!is.na(v) & nzchar(trimws(v))]
  clean(paste(v, collapse = "; "))
}
fmt <- function(x, d = 4) ifelse(is.na(x), "NA", as.character(signif(x, d)))

# ---- shared implementation: association_tests.R (sourced, main() not run) ---------- #
# fit_phyloglm_ms / run_phyloglm (multi-start Firth-penalised phyloglm), the strict
# parsers, the Contract v2 item 6 tier rule, read_species_tree / prepare_tree.
script_dir <- function() {
  # source()d (unit tests, possibly from another Rscript): the source() frame records
  # this file in `ofile`; under Rscript the path is the --file= argument.
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (is.character(of) && length(of) == 1L && file.exists(of))
      return(dirname(normalizePath(of, mustWork = TRUE)))
  }
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f) == 1L && nzchar(f) && file.exists(f)) return(dirname(normalizePath(f, mustWork = TRUE)))
  "workflow/scripts/R"                                  # repository root as working directory
}
ASSOC_R <- file.path(script_dir(), "association_tests.R")
if (!file.exists(ASSOC_R)) stop("association_tests.R not found next to methanogenesis_stats.R: ", ASSOC_R)
A <- new.env()
sys.source(ASSOC_R, envir = A)
for (fn in c("fit_phyloglm_ms", "run_phyloglm", "read_species_tree", "prepare_tree",
             "tier_response", "parse_bool", "parse_num", "parse_tiers", "read_tsv_chr", "zscore"))
  if (!is.function(A[[fn]])) stop("association_tests.R does not define ", fn, "()")

# ---- small helpers ---------------------------------------------------------------- #
parse_list <- function(x) {
  if (is.null(x)) return(character(0))
  v <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  v[nzchar(v)]
}

parse_flag <- function(x, what) {
  s <- tolower(trimws(as.character(x)))
  if (s %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (s %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(what, " must be TRUE or FALSE (got '", x, "')")
}

# 1 / 0 / NA; anything else is an error naming the column (fail loudly).
parse_indicator <- function(x, col) {
  s <- tolower(trimws(as.character(x)))
  v <- rep(NA_integer_, length(s))
  v[s %in% TRUE_S] <- 1L
  v[s %in% FALSE_S] <- 0L
  bad <- !(s %in% c(TRUE_S, FALSE_S, NA_S))
  if (any(bad))
    stop(sprintf("column '%s': %d value(s) not in {1, 0, NA}, e.g. '%s'", col, sum(bad),
                 as.character(x)[bad][1]))
  v
}

header_only <- function(path, why) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(ALL_COLS, collapse = "\t"), path)
  say("header-only output (", why, ") -> ", path)
  invisible(NULL)
}

new_row <- function(predictor, analysis, tier_set) {
  r <- as.list(rep(NA, length(ALL_COLS)))
  names(r) <- ALL_COLS
  r$predictor <- predictor
  r$analysis <- analysis
  r$role <- ANALYSES[[analysis]]
  r$tier_set <- tier_set
  r$note <- ""
  r$model <- ""
  r$covariates <- ""
  r$habitat_terms <- ""
  r$best_model <- ""
  r$mk_model <- ""
  r
}

counts_into <- function(r, d) {
  r$n <- nrow(d)
  r$n_present <- sum(d$y == 1L)
  r$n_absent <- sum(d$y == 0L)
  r$n_level1 <- sum(d$x == 1L)
  r
}

gate <- function(d, opt, response = TRUE) {
  why <- character(0)
  n <- nrow(d)
  if (n < opt$min_n) why <- c(why, sprintf("n=%d < min_n=%d", n, opt$min_n))
  if (response) {
    np <- sum(d$y == 1L)
    if (min(np, n - np) < opt$min_class)
      why <- c(why, sprintf("present=%d absent=%d (min %d each)", np, n - np, opt$min_class))
  }
  n1 <- sum(d$x == 1L)
  if (min(n1, n - n1) < opt$min_class)
    why <- c(why, sprintf("predictor 1=%d 0=%d (min %d each)", n1, n - n1, opt$min_class))
  why
}

# R^2 of the linear regression of column `resp` on `terms` (aliased terms do not change
# the fitted span, so R^2 stays valid when two adjustment terms coincide).
r2_on <- function(d, resp, terms) {
  if (!length(terms)) return(0)
  f <- stats::lm(stats::reformulate(terms, response = resp), data = d)
  r2 <- suppressWarnings(summary(f)$r.squared)        # "essentially perfect fit" is the point
  if (is.finite(r2)) r2 else 1
}

vif_of_x <- function(d, terms) {
  r2 <- r2_on(d, "x", terms)
  if (r2 >= 1 - 1e-10) Inf else 1 / (1 - r2)
}

prune_to <- function(tree, tips) {
  tr <- if (length(tips) < length(tree$tip.label)) ape::keep.tip(tree, tips) else tree
  if (!ape::is.rooted(tr)) stop("pruned species tree lost its root")
  tr
}

# ---- data assembly for one analysis ----------------------------------------------- #
# base: the analysis-set frame; x / y: predictor and response vectors aligned to base;
# keep: optional row filter (sensitivity subsets). Returns complete cases on
# (y, x, covariates) with z-scored covariates zc1.. and the drop counts.
assemble <- function(base, x, y, keep, ctx) {
  if (is.null(keep)) keep <- rep(TRUE, nrow(base))
  b <- base[keep, , drop = FALSE]; x <- x[keep]; y <- y[keep]
  n_excl <- sum(is.na(y))
  ok <- !is.na(y)
  n_unres <- sum(ok & is.na(x))
  ok <- ok & !is.na(x)
  cov_ok <- rep(TRUE, nrow(b))
  for (cv in ctx$covariates) cov_ok <- cov_ok & !is.na(b[[cv]])
  n_cov <- sum(ok & !cov_ok)
  ok <- ok & cov_ok
  d <- data.frame(tip = b$tip[ok], y = as.integer(y[ok]), x = as.integer(x[ok]),
                  order = b$gtdb_order[ok], completeness = b$completeness[ok],
                  stringsAsFactors = FALSE)
  for (h in ctx$habitat_cols) d[[ctx$habitat_id[[h]]]] <- b[[h]][ok]
  terms <- character(0); labels <- character(0); notes <- character(0)
  for (k in seq_along(ctx$covariates)) {
    cv <- ctx$covariates[k]
    if (nrow(d) < 2L) next                      # nothing to fit; the row reports why
    z <- A$zscore(b[[cv]][ok])
    if (is.null(z)) { notes <- c(notes, paste0(cv, " constant: not adjusted")); next }
    id <- paste0("zc", k)
    d[[id]] <- z
    terms <- c(terms, id); labels <- c(labels, paste0("z(", cv, ")"))
  }
  rownames(d) <- d$tip
  list(d = d, terms = terms, labels = labels, notes = notes, n_excl = n_excl,
       n_unres = n_unres, n_cov = n_cov)
}

# ---- habitat selection: eligibility, |phi| order, EPV, complete cases ------------- #
select_habitat <- function(dat, ctx, opt) {
  d <- dat$d
  base_k <- 1L + length(dat$terms)
  elig <- character(0); phi <- numeric(0); inelig <- character(0)
  for (h in ctx$habitat_cols) {
    id <- ctx$habitat_id[[h]]
    hv <- d[[id]]; okh <- !is.na(hv)
    n1 <- sum(hv[okh] == 1L); n0 <- sum(hv[okh] == 0L)
    if (min(n1, n0) < opt$min_class) { inelig <- c(inelig, h); next }
    xs <- d$x[okh]
    ph <- if (stats::sd(xs) > 0) stats::cor(xs, hv[okh]) else NA_real_
    if (!is.finite(ph)) { inelig <- c(inelig, h); next }
    elig <- c(elig, h); phi <- c(phi, ph)
  }
  ord <- if (length(elig)) order(-abs(phi), elig, method = "radix") else integer(0)  # C-locale ties
  cur <- d; inc <- character(0); inc_phi <- numeric(0); skipped <- character(0)
  for (j in ord) {
    h <- elig[j]; id <- ctx$habitat_id[[h]]
    nxt <- cur[!is.na(cur[[id]]), , drop = FALSE]
    np <- sum(nxt$y == 1L); ne <- min(np, nrow(nxt) - np)
    budget <- floor(ne / opt$epv_min)
    k <- base_k + length(inc) + 1L
    if (k > budget) {
      skipped <- c(skipped, sprintf("%s (|phi|=%.2f): EPV %d terms > floor(%d/%g)=%d%s", h,
                                    abs(phi[j]), k, ne, opt$epv_min, budget,
                                    if (nrow(nxt) < nrow(cur)) sprintf(" on its %d complete cases", nrow(nxt)) else ""))
      next
    }
    why <- gate(nxt, opt)
    hv <- nxt[[id]]
    if (min(sum(hv == 1L), sum(hv == 0L)) < opt$min_class)
      why <- c(why, sprintf("%s level < min_class on the complete cases", h))
    # a habitat column that duplicates the adjustment terms already in the model adds no
    # adjustment and would make the design singular (x itself is NOT among these terms:
    # a habitat that duplicates x is the non-identifiable case the VIF guard reports)
    if (!length(why) && r2_on(nxt, id, c(dat$terms, unname(unlist(ctx$habitat_id[inc])))) >= 1 - 1e-8)
      why <- c(why, "redundant with the adjustment terms already included")
    if (length(why)) {
      skipped <- c(skipped, sprintf("%s: %s", h, paste(why, collapse = ", ")))
      next
    }
    cur <- nxt; inc <- c(inc, h); inc_phi <- c(inc_phi, phi[j])
  }
  np <- sum(cur$y == 1L); ne <- min(np, nrow(cur) - np)
  list(d = cur, habitat = inc, phi = inc_phi, ids = unname(unlist(ctx$habitat_id[inc])),
       n_eligible = length(elig), n_ineligible = length(inelig), skipped = skipped,
       n_events = ne, budget = floor(ne / opt$epv_min), base_k = base_k,
       n_dropped_na = nrow(d) - nrow(cur))
}

# ---- one phyloglm fit (shared multi-start MPLE) ----------------------------------- #
fit_logit <- function(d, terms, tree, opt) {
  tr <- prune_to(tree, d$tip)
  d <- d[match(tr$tip.label, d$tip), , drop = FALSE]
  rownames(d) <- d$tip
  fml <- stats::as.formula(paste("y ~", paste(c("x", terms), collapse = " + ")))
  res <- tryCatch(A$run_phyloglm(fml, d, tr, "x", opt$btol, 0L, opt$seed, 1L),
                  error = function(e) e)
  if (inherits(res, "error")) return(list(ok = FALSE, error = conditionMessage(res)))
  res$ok <- is.finite(res$est) && is.finite(res$se)
  if (!res$ok) res$error <- "non-finite estimate or standard error"
  res
}

model_label <- function(pred, labels) {
  sprintf("phyloglm(present ~ %s, logistic_MPLE, multi-start)",
          paste(c(pred, labels), collapse = " + "))
}

# Fit one regression row: marginal (no habitat selection) or the habitat_adjusted
# procedure (habitat = TRUE; also used by the sensitivity rows).
regression_row <- function(r, dat, pred, tree, opt, ctx, habitat, extra_note = character(0)) {
  r$n_excluded_tier <- dat$n_excl; r$n_unresolved <- dat$n_unres
  r$n_missing_covariate <- dat$n_cov
  notes <- c(extra_note, dat$notes)
  d <- dat$d; terms <- dat$terms; labels <- dat$labels
  sel <- NULL
  if (habitat) {
    why0 <- gate(d, opt)
    if (length(why0)) {
      r <- counts_into(r, d)
      r$model <- model_label(pred, labels)
      r$covariates <- paste(labels, collapse = ",")
      r$note <- join_notes(notes, paste("not fitted:", paste(why0, collapse = "; ")))
      r$.status <- "skipped"
      return(r)
    }
    sel <- select_habitat(dat, ctx, opt)
    d <- sel$d
    terms <- c(terms, sel$ids); labels <- c(labels, sel$habitat)
    r$habitat_terms <- paste(sel$habitat, collapse = ",")
    r$n_habitat_eligible <- sel$n_eligible
    hnote <- if (!length(ctx$habitat_cols)) "no habitat column in the master table" else
      sprintf("habitat: %d of %d column(s) eligible (>= %d per level), %d included by |phi(x,h)|%s",
              sel$n_eligible, length(ctx$habitat_cols), opt$min_class, length(sel$habitat),
              if (length(sel$habitat)) paste0(" (", paste(sprintf("%s %.2f", sel$habitat, sel$phi),
                                                         collapse = ", "), ")") else "")
    notes <- c(notes, hnote)
    if (sel$base_k >= sel$budget && sel$n_eligible > 0L && !length(sel$habitat)) {
      notes <- c(notes, sprintf(paste0("EPV: x + %d covariate(s) = %d term(s) already %s floor(%d events/%g)=%d; ",
                                       "none of the %d eligible habitat column(s) can be added"),
                                length(dat$terms), sel$base_k, if (sel$base_k > sel$budget) "exceed" else "reach",
                                sel$n_events, opt$epv_min, sel$budget, sel$n_eligible))
    } else if (length(sel$skipped)) {
      notes <- c(notes, paste0("not added: ", paste(sel$skipped, collapse = " | ")))
    }
    if (sel$n_dropped_na) notes <- c(notes, sprintf("%d genome(s) with NA in an included habitat column dropped",
                                                    sel$n_dropped_na))
  }
  r <- counts_into(r, d)
  np <- sum(d$y == 1L); ne <- min(np, nrow(d) - np)
  r$n_events <- ne; r$n_terms <- 1L + length(terms)
  r$epv <- ne / r$n_terms
  r$covariates <- paste(labels, collapse = ",")
  r$model <- model_label(pred, labels)
  why <- gate(d, opt)
  if (length(why)) {
    r$note <- join_notes(notes, paste("not fitted:", paste(why, collapse = "; ")))
    r$.status <- "skipped"
    return(r)
  }
  if (!habitat && r$epv < opt$epv_min)
    notes <- c(notes, sprintf("EPV=%.1f < epv_min=%g (interpret with caution)", r$epv, opt$epv_min))
  r$vif <- vif_of_x(d, terms)
  if (habitat && r$vif > opt$vif_max) {
    r$note <- join_notes(notes, sprintf(paste0("partial effect not identifiable: VIF(x)=%s > vif_max=%g ",
                                               "(x is nearly a linear function of the adjustment terms); not fitted"),
                                        fmt(r$vif), opt$vif_max))
    r$.status <- "vif"
    return(r)
  }
  tab <- table(factor(d$x, levels = 0:1), factor(d$y, levels = 0:1))
  if (any(tab == 0))
    notes <- c(notes, sprintf(paste0("quasi-separation (x0: absent=%d present=%d; x1: absent=%d present=%d): ",
                                     "Firth-penalised estimate, interpret with caution"),
                              tab[1, 1], tab[1, 2], tab[2, 1], tab[2, 2]))
  f <- fit_logit(d, terms, tree, opt)
  if (!isTRUE(f$ok)) {
    r$note <- join_notes(notes, paste("fit failed:", f$error))
    r$converged <- FALSE
    r$.status <- "failed"
    say("WARNING: ", pred, " ", r$analysis, ": phyloglm failed: ", f$error)
    return(r)
  }
  r$estimate <- f$est; r$std_error <- f$se; r$p_value <- f$p
  r$ci_low <- f$lo; r$ci_high <- f$hi; r$alpha <- f$alpha
  r$converged <- isTRUE(f$converged)
  r$note <- join_notes(notes, sprintf("log-odds; Wald z/p and CI conditional on alpha=%s (alpha*T=%s)",
                                      fmt(f$alpha), fmt(f$alpha * f$tmax)), f$flags)
  r$.status <- "fitted"
  r
}

# ---- replication guard: Mk (ER vs ARD by AIC) + stochastic character maps --------- #
# Minimum number of state changes of a binary tip trait (Fitch-Hartigan parsimony,
# polytomies allowed): a model-free lower bound that every stochastic map attains.
fitch_changes <- function(tr, x) {
  tr <- ape::reorder.phylo(tr, "postorder")
  n <- ape::Ntip(tr); N <- n + tr$Nnode
  xs <- as.integer(x[tr$tip.label])
  s0 <- s1 <- logical(N); s0[seq_len(n)] <- xs == 0L; s1[seq_len(n)] <- xs == 1L
  k0 <- k1 <- nc <- integer(N); cost <- 0L
  settle <- function(v) {
    K <- max(k0[v], k1[v])
    s0[v] <<- k0[v] == K; s1[v] <<- k1[v] == K
    cost <<- cost + nc[v] - K
  }
  for (e in seq_len(nrow(tr$edge))) {
    a <- tr$edge[e, 1]; b <- tr$edge[e, 2]
    if (b > n) settle(b)                     # all of b's child edges precede (postorder)
    k0[a] <- k0[a] + s0[b]; k1[a] <- k1[a] + s1[b]; nc[a] <- nc[a] + 1L
  }
  settle(n + 1L)
  cost
}

mk_q <- function(fit) {
  st <- fit$states
  Q <- matrix(0, length(st), length(st), dimnames = list(st, st))
  idx <- fit$index.matrix
  off <- !is.na(idx) & idx > 0
  Q[off] <- fit$rates[idx[off]]
  diag(Q) <- -rowSums(Q)
  Q
}

simmap_row <- function(r, dat, pred, tree, opt) {
  d <- dat$d
  r <- counts_into(r, d)
  r$n_excluded_tier <- dat$n_excl; r$n_unresolved <- dat$n_unres
  r$n_missing_covariate <- dat$n_cov
  r$model <- sprintf("phytools fitMk ER vs ARD (AIC, flat root prior) + make.simmap(nsim=%d, Q fixed at the MLE)",
                     opt$nsim)
  why <- gate(d, opt, response = FALSE)
  if (length(why)) {
    r$note <- join_notes(dat$notes, paste("not mapped:", paste(why, collapse = "; ")))
    return(r)
  }
  tr <- prune_to(tree, d$tip)
  xs <- setNames(factor(as.character(d$x[match(tr$tip.label, d$tip)]), levels = c("0", "1")),
                 tr$tip.label)
  fits <- list(ER = phytools::fitMk(tr, xs, model = "ER", pi = "equal"),
               ARD = phytools::fitMk(tr, xs, model = "ARD", pi = "equal"))
  ll <- vapply(fits, function(f) as.numeric(f$logLik), numeric(1))
  k <- vapply(fits, function(f) max(f$index.matrix, na.rm = TRUE), numeric(1))
  aic <- 2 * k - 2 * ll
  if (any(!is.finite(aic))) stop(pred, ": Mk fit returned a non-finite likelihood")
  best <- if (aic[["ARD"]] < aic[["ER"]]) "ARD" else "ER"
  Q <- mk_q(fits[[best]])
  set.seed(opt$seed)
  maps <- phytools::make.simmap(tr, xs, model = best, nsim = opt$nsim, Q = Q, pi = "equal",
                                message = FALSE)
  if (opt$nsim == 1L) maps <- list(maps)
  per <- t(vapply(maps, function(m) {
    ch <- 0; g <- 0; l <- 0
    for (seg in m$maps) {
      s <- names(seg)
      if (length(s) > 1L) {
        a <- s[-length(s)]; b <- s[-1L]
        ch <- ch + sum(a != b); g <- g + sum(a == "0" & b == "1"); l <- l + sum(a == "1" & b == "0")
      }
    }
    c(ch, g, l)
  }, numeric(3)))
  nch <- per[, 1]
  r$n_changes_parsimony <- fitch_changes(tr, setNames(d$x[match(tr$tip.label, d$tip)], tr$tip.label))
  if (any(nch < r$n_changes_parsimony))
    stop(pred, ": a stochastic map has fewer changes than the parsimony minimum (internal error)")
  r$n_changes_mean <- mean(nch)
  q <- stats::quantile(nch, c(0.025, 0.975), names = FALSE, type = 7)
  r$n_changes_lo <- q[1]; r$n_changes_hi <- q[2]
  r$n_gains_mean <- mean(per[, 2]); r$n_losses_mean <- mean(per[, 3])
  r$replicated <- r$n_changes_mean >= opt$min_changes
  r$mk_model <- best
  r$note <- join_notes(dat$notes,
    sprintf("Mk AIC: ER=%.3f ARD=%.3f -> %s (Q01=%s, Q10=%s)", aic[["ER"]], aic[["ARD"]], best,
            fmt(Q["0", "1"]), fmt(Q["1", "0"])),
    sprintf("%d stochastic maps, seed %d: changes mean %.2f (95%% %s-%s; gains %.2f, losses %.2f)",
            opt$nsim, opt$seed, r$n_changes_mean, fmt(q[1]), fmt(q[2]), r$n_gains_mean, r$n_losses_mean),
    sprintf("parsimony (Fitch) minimum %d change(s)", r$n_changes_parsimony),
    if (r$n_changes_mean > nrow(d))
      "more changes than tips: x carries little phylogenetic signal (Mk rate on a flat likelihood ridge)",
    sprintf("replicated = mean changes >= min_changes=%g (sign test: k concordant changes give p=2^-k): %s",
            opt$min_changes, if (r$replicated) "TRUE" else "FALSE"))
  r$.status <- "mapped"
  r
}

# ---- corHMM hidden-rate correlation test ------------------------------------------ #
corhmm_row <- function(r, dat, pred, tree, opt) {
  d <- dat$d
  r <- counts_into(r, d)
  r$n_excluded_tier <- dat$n_excl; r$n_unresolved <- dat$n_unres
  r$n_missing_covariate <- dat$n_cov
  r$model <- "corHMM::fitCorrelationTest(x, present): independent vs correlated, 1 or 2 rate classes"
  if (!opt$corhmm) {
    r$note <- "corHMM test disabled (--corhmm FALSE)"
    return(r)
  }
  why <- gate(d, opt)
  if (length(why)) {
    r$note <- join_notes(dat$notes, paste("not fitted:", paste(why, collapse = "; ")))
    return(r)
  }
  tr <- prune_to(tree, d$tip)
  dd <- d[match(tr$tip.label, d$tip), , drop = FALSE]
  cd <- data.frame(species = dd$tip, x = as.character(dd$x), present = as.character(dd$y),
                   stringsAsFactors = FALSE)
  set.seed(opt$seed)
  out <- NULL
  log_txt <- utils::capture.output(
    out <- withCallingHandlers(corHMM::fitCorrelationTest(phy = tr, data = cd),
                               warning = function(w) {
                                 say(pred, " corHMM warning: ", conditionMessage(w))
                                 invokeRestart("muffleWarning")
                               }))
  if (is.null(out) || !all(CORHMM_MODELS %in% names(out)))
    stop(pred, ": corHMM::fitCorrelationTest returned no four-model list (got ",
         paste(names(out), collapse = ","), ")")
  aic <- vapply(CORHMM_MODELS, function(m) as.numeric(out[[m]]$AIC), numeric(1))
  ll <- vapply(CORHMM_MODELS, function(m) as.numeric(out[[m]]$loglik), numeric(1))
  np <- vapply(CORHMM_MODELS, function(m) max(out[[m]]$index.mat, na.rm = TRUE), numeric(1))
  names(aic) <- names(ll) <- names(np) <- names(CORHMM_MODELS)
  if (any(!is.finite(aic))) stop(pred, ": corHMM returned a non-finite AIC")
  w <- exp(-0.5 * (aic - min(aic))); w <- w / sum(w)
  r$delta_aic <- min(aic[c("independent", "hmm_independent")]) - min(aic[c("correlated", "hmm_correlated")])
  r$best_model <- names(aic)[which.min(aic)]
  for (m in names(aic)) {
    r[[paste0("aic_", m)]] <- aic[[m]]
    r[[paste0("aicw_", m)]] <- w[[m]]
  }
  r$note <- join_notes(dat$notes,
    paste0("delta_aic = min AIC(independent models) - min AIC(correlated models); > 0 favours ",
           "correlated evolution; no p-value"),
    paste(sprintf("%s: logL=%.3f k=%d AIC=%.3f w=%.3f", names(aic), ll, as.integer(np), aic, w),
          collapse = " | "))
  r$.status <- "fitted"
  r
}

# ---- one predictor ----------------------------------------------------------------- #
leave_out_order <- function(d) {
  o <- trimws(as.character(d$order[d$x == 1L]))
  o <- o[nzchar(o) & !tolower(o) %in% NA_S]
  if (!length(o)) return(NULL)
  tab <- table(o)
  cand <- names(tab)[tab == max(tab)]
  list(order = sort(cand, method = "radix")[1], n = as.integer(max(tab)), tied = length(cand) > 1L)
}

copy_fit <- function(r, prim) {     # values of an identical-data fit
  for (k in c("estimate", "std_error", "ci_low", "ci_high", "p_value", "n", "n_present",
              "n_absent", "n_level1", "covariates", "vif", "model", "converged", "n_events",
              "n_terms", "epv", "habitat_terms", "n_habitat_eligible", "alpha",
              "n_excluded_tier", "n_unresolved", "n_missing_covariate"))
    r[[k]] <- prim[[k]]
  r$.status <- prim$.status
  r
}

analyse_predictor <- function(p, ctx, opt) {
  base <- ctx$base
  rows <- list()
  xp <- base[[p]]
  T1 <- ctx$tier_lab; T2 <- ctx$tier_lab_sens
  secs <- numeric(0)
  timed <- function(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    v <- expr
    secs[[label]] <<- proc.time()[["elapsed"]] - t0
    v
  }
  prim <- assemble(base, xp, base$y_prim, NULL, ctx)
  rows$marginal <- timed("marginal", regression_row(new_row(p, "marginal", T1), prim, p, ctx$tree, opt,
                                                   ctx, habitat = FALSE))
  rows$habitat_adjusted <- timed("habitat_adjusted", regression_row(new_row(p, "habitat_adjusted", T1), prim,
                                                                   p, ctx$tree, opt, ctx, habitat = TRUE))
  rows$simmap <- timed("simmap", simmap_row(new_row(p, "simmap", T1), prim, p, ctx$tree, opt))
  rows$corhmm <- timed("corhmm", corhmm_row(new_row(p, "corhmm", T1), prim, p, ctx$tree, opt))

  # ---- sensitivity rows (habitat_adjusted procedure re-applied) ---- #
  hp <- rows$habitat_adjusted
  sens <- function(analysis, dat, note, identical_to_primary = FALSE) {
    r <- new_row(p, analysis, if (analysis == "sens_tiers_AB") T2 else T1)
    if (identical_to_primary) {
      r <- copy_fit(r, hp)
      r$note <- join_notes(note, "identical data to habitat_adjusted: its values are repeated", hp$note)
      return(r)
    }
    r <- regression_row(r, dat, p, ctx$tree, opt, ctx, habitat = TRUE, extra_note = note)
    if (identical(r$.status, "fitted") && !identical(r$habitat_terms, hp$habitat_terms))
      r$note <- join_notes(r$note, sprintf("habitat terms differ from habitat_adjusted (%s)",
                                           if (nzchar(hp$habitat_terms)) hp$habitat_terms else "none"))
    r
  }
  n_unres <- prim$n_unres
  x0 <- xp; x0[is.na(x0)] <- 0L
  t_sens <- proc.time()[["elapsed"]]
  rows$sens_unresolved_as_absent <- sens("sens_unresolved_as_absent",
    if (n_unres) assemble(base, x0, base$y_prim, NULL, ctx) else NULL,
    sprintf("unresolved (NA) %s coded 0 (%d genome(s))", p, n_unres), identical_to_primary = n_unres == 0L)
  comp <- base$completeness
  keep_c <- !is.na(comp) & comp >= opt$complete_threshold
  n_drop_c <- sum(!keep_c)
  rows$sens_complete_only <- sens("sens_complete_only",
    if (n_drop_c) assemble(base, xp, base$y_prim, keep_c, ctx) else NULL,
    sprintf("completeness >= %g only (%d genome(s) of the analysis set removed)", opt$complete_threshold, n_drop_c),
    identical_to_primary = n_drop_c == 0L)
  lo <- leave_out_order(prim$d)
  if (is.null(lo)) {
    r <- new_row(p, "sens_leave_clade_out", T1)
    r$note <- "not fitted: no x = 1 genome with a gtdb_order in the analysis set"
    r$.status <- "skipped"
    rows$sens_leave_clade_out <- r
  } else {
    keep_o <- trimws(as.character(base$gtdb_order)) != lo$order
    rows$sens_leave_clade_out <- sens("sens_leave_clade_out", assemble(base, xp, base$y_prim, keep_o, ctx),
      sprintf("left out gtdb_order %s (%d of %d x=1 genome(s)%s; %d genome(s) of the analysis set removed)",
              lo$order, lo$n, sum(prim$d$x == 1L), if (lo$tied) ", tie broken alphabetically" else "",
              sum(!keep_o)))
  }
  same_resp <- identical(is.na(base$y_prim), is.na(base$y_sens)) &&
    isTRUE(all(base$y_prim == base$y_sens, na.rm = TRUE))
  rows$sens_tiers_AB <- sens("sens_tiers_AB",
    if (!same_resp) assemble(base, xp, base$y_sens, NULL, ctx) else NULL,
    sprintf("response by tier set %s", T2), identical_to_primary = same_resp)

  secs[["sensitivity"]] <- proc.time()[["elapsed"]] - t_sens
  say(sprintf("%s done in %.1f s (%s)", p, sum(secs),
              paste(sprintf("%s %.1f", names(secs), secs), collapse = ", ")))

  # ---- replication guard ---- #
  rep_ok <- rows$simmap$replicated
  for (a in setdiff(names(rows), c("simmap", "corhmm"))) {
    r <- rows[[a]]
    r$replicated <- rep_ok
    if (!is.na(r$p_value) && !isTRUE(rep_ok)) {
      why <- if (is.na(rep_ok)) "history of x not mapped (see the simmap row)" else
        sprintf("mean %.2f state changes of %s < min_changes=%g (Maddison & FitzJohn 2015)",
                rows$simmap$n_changes_mean, p, opt$min_changes)
      r$p_value <- NA_real_
      r$note <- join_notes(paste0("replication guard: ", why, "; p-value withheld, estimate and CI kept"), r$note)
    }
    rows[[a]] <- r
  }
  rows$corhmm$replicated <- rep_ok
  rows[names(ANALYSES)]
}

# --------------------------------------------------------------------------------------- #
main <- function() {
  opt <- parse_args(OptionParser(option_list = list(
    make_option("--master", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--out", type = "character"),
    make_option("--predictors", type = "character", default = paste(DEFAULT_PREDICTORS, collapse = ",")),
    make_option("--habitat", type = "character", default = "auto",
                help = "auto (every group_/niche_ column), none, or a comma list of columns"),
    make_option("--covariates", type = "character", default = "completeness,genome_size"),
    make_option("--primary-tiers", type = "character", default = "A", dest = "primary_tiers"),
    make_option("--sensitivity-tiers", type = "character", default = "A,B", dest = "sensitivity_tiers"),
    make_option("--min-n", type = "integer", default = 10L, dest = "min_n"),
    make_option("--min-class", type = "integer", default = 3L, dest = "min_class"),
    make_option("--min-changes", type = "double", default = 5, dest = "min_changes"),
    make_option("--nsim", type = "integer", default = 100L),
    make_option("--corhmm", type = "character", default = "TRUE"),
    make_option("--epv-min", type = "double", default = 10, dest = "epv_min"),
    make_option("--vif-max", type = "double", default = 5, dest = "vif_max"),
    make_option("--complete-threshold", type = "double", default = 95, dest = "complete_threshold"),
    make_option("--fdr-method", type = "character", default = "BH", dest = "fdr_method"),
    make_option("--btol", type = "double", default = 10),
    make_option("--seed", type = "integer", default = 12345L),
    make_option("--threads", type = "integer", default = 1L,
                help = "predictors analysed in parallel (results do not depend on it)"))))
  for (k in c("master", "tree", "out"))
    if (is.null(opt[[k]]) || !nzchar(opt[[k]])) stop("--", k, " is required")
  opt$corhmm <- parse_flag(opt$corhmm, "--corhmm")
  if (is.na(opt$min_n) || opt$min_n < 3L) stop("--min-n must be >= 3")
  if (is.na(opt$min_class) || opt$min_class < 1L) stop("--min-class must be >= 1")
  if (!is.finite(opt$min_changes) || opt$min_changes <= 0) stop("--min-changes must be > 0")
  if (is.na(opt$nsim) || opt$nsim < 1L) stop("--nsim must be >= 1")
  if (!is.finite(opt$epv_min) || opt$epv_min <= 0) stop("--epv-min must be > 0")
  if (!is.finite(opt$vif_max) || opt$vif_max < 1) stop("--vif-max must be >= 1")
  if (!is.finite(opt$complete_threshold)) stop("--complete-threshold must be a number")
  if (!is.finite(opt$btol) || opt$btol <= 0) stop("--btol must be > 0")
  if (is.na(opt$threads) || opt$threads < 1L) stop("--threads must be >= 1")
  fdr <- if (toupper(opt$fdr_method) %in% c("BH", "FDR")) "BH" else opt$fdr_method
  if (!fdr %in% stats::p.adjust.methods) stop("--fdr-method '", opt$fdr_method, "' is not a p.adjust method")
  tiers <- A$parse_tiers(opt$primary_tiers, "--primary-tiers")
  tiers_s <- A$parse_tiers(opt$sensitivity_tiers, "--sensitivity-tiers")
  preds <- parse_list(opt$predictors)
  if (!length(preds)) stop("--predictors is empty")
  if (anyDuplicated(preds)) stop("--predictors lists a column twice")
  covars <- parse_list(opt$covariates)
  if (opt$corhmm && !requireNamespace("corHMM", quietly = TRUE))
    stop("R package 'corHMM' is not installed (required with --corhmm TRUE; the rstats env ",
         "installs the vendored corHMM 2.8 in rstats.post-deploy.sh)")
  if (opt$corhmm && !is.function(tryCatch(getExportedValue("corHMM", "fitCorrelationTest"),
                                          error = function(e) NULL)))
    stop("installed corHMM ", as.character(utils::packageVersion("corHMM")),
         " lacks fitCorrelationTest (corHMM >= 2.8 required)")

  # ---- inputs ---------------------------------------------------------------- #
  if (!file.exists(opt$tree)) stop("species tree not found: ", opt$tree)
  if (!file.exists(opt$master)) stop("master table not found: ", opt$master)
  if (file.info(opt$master)$size == 0) stop("master table is a 0-byte file: ", opt$master)
  tree0 <- A$read_species_tree(opt$tree)
  if (is.null(tree0)) return(header_only(opt$out, "species tree is empty ('();')"))  # nonfatal: empty input
  tree0 <- A$prepare_tree(tree0)
  master <- A$read_tsv_chr(opt$master)
  need <- unique(c("genome_id", "tree_tip", "methanogen", "prometheus_tier", "gtdb_order",
                   "completeness", covars, preds))
  miss <- setdiff(need, names(master))
  if (length(miss)) stop("master table lacks required column(s): ", paste(miss, collapse = ", "))
  if (nrow(master) == 0) return(header_only(opt$out, "master table has zero genomes"))  # nonfatal: empty input
  gid <- trimws(master$genome_id)
  if (any(!nzchar(gid))) stop("master table has blank genome_id values")
  if (anyDuplicated(gid)) stop("master table has duplicated genome_id: ", gid[duplicated(gid)][1])
  tip <- trimws(master$tree_tip)
  if (anyDuplicated(tip[nzchar(tip)])) stop("master table has duplicated tree_tip values")
  meth <- A$parse_bool(master$methanogen, "methanogen")
  anme <- if ("is_anme" %in% names(master)) A$parse_bool(master$is_anme, "is_anme", blank = FALSE) else
    rep(FALSE, nrow(master))
  pass <- if ("pass_primary" %in% names(master)) A$parse_bool(master$pass_primary, "pass_primary", blank = FALSE) else
    rep(TRUE, nrow(master))
  y1 <- A$tier_response(master$prometheus_tier, tiers)
  y2 <- A$tier_response(master$prometheus_tier, tiers_s)
  on_tree <- nzchar(tip) & tip %in% tree0$tip.label
  sel <- on_tree & meth & !anme & pass
  say(sprintf(paste0("master=%d genomes; methanogen=%d; methanogen tree tips=%d; ANME dropped=%d; ",
                     "not pass_primary dropped=%d; analysis set=%d"),
              nrow(master), sum(meth), sum(meth & on_tree), sum(meth & on_tree & anme),
              sum(meth & on_tree & !anme & !pass), sum(sel)))
  if (!any(sel)) return(header_only(opt$out, "no methanogen tip of the species tree"))  # nonfatal: empty input

  # habitat columns (binary 1/0/NA; any other value is an error)
  hab_spec <- parse_list(opt$habitat)
  if (!length(hab_spec) || identical(tolower(hab_spec), "auto")) {
    hab <- grep("^(group|niche)_", names(master), value = TRUE)
  } else if (identical(tolower(hab_spec), "none")) {
    hab <- character(0)
  } else {
    hm <- setdiff(hab_spec, names(master))
    if (length(hm)) stop("--habitat column(s) not in the master table: ", paste(hm, collapse = ", "))
    hab <- unique(hab_spec)
  }
  hab <- setdiff(hab, c(preds, covars))
  base <- data.frame(tip = tip[sel], y_prim = y1[sel], y_sens = y2[sel],
                     gtdb_order = trimws(master$gtdb_order[sel]),
                     completeness = A$parse_num(master$completeness, "completeness")[sel],
                     stringsAsFactors = FALSE)
  for (cv in covars) base[[cv]] <- A$parse_num(master[[cv]], cv)[sel]
  for (h in hab) base[[h]] <- parse_indicator(master[[h]], h)[sel]
  for (p in preds) base[[p]] <- parse_indicator(master[[p]], p)[sel]
  ctx <- list(base = base, tree = prune_to(tree0, base$tip), covariates = covars,
              habitat_cols = hab, habitat_id = as.list(setNames(sprintf("h%d", seq_along(hab)), hab)),
              tier_lab = paste(tiers, collapse = ","), tier_lab_sens = paste(tiers_s, collapse = ","))
  say(sprintf("analysis set: %d methanogen genomes; tiers %s: present=%d absent=%d excluded=%d; %d habitat column(s)",
              nrow(base), ctx$tier_lab, sum(base$y_prim %in% 1L), sum(base$y_prim %in% 0L),
              sum(is.na(base$y_prim)), length(hab)))

  # ---- per predictor (parallel over predictors; every random step is seeded) ---- #
  one <- function(p) {
    say("predictor ", p)
    analyse_predictor(p, ctx, opt)
  }
  res <- if (opt$threads > 1L && length(preds) > 1L)
    parallel::mclapply(preds, one, mc.cores = min(opt$threads, length(preds)), mc.preschedule = FALSE)
  else lapply(preds, one)
  bad <- vapply(res, function(r) inherits(r, "try-error") || !is.list(r), logical(1))
  if (any(bad)) stop("analysis of ", preds[which(bad)[1]], " failed: ",
                     paste(as.character(res[[which(bad)[1]]]), collapse = " "))
  rows <- unlist(res, recursive = FALSE)
  st <- vapply(rows, function(r) if (is.null(r$.status)) "" else r$.status, character(1))
  reg <- vapply(rows, function(r) !r$analysis %in% c("simmap", "corhmm"), logical(1))
  attempted <- sum(reg & st %in% c("fitted", "failed")); failed <- sum(reg & st == "failed")
  if (attempted > 0 && failed == attempted)
    stop("every attempted phyloglm fit failed (", failed, "); see the notes above")

  out <- do.call(rbind, lapply(rows, function(r) {
    r$.status <- NULL
    as.data.frame(r[ALL_COLS], stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out$p_adj <- NA_real_
  fam <- which(out$role == "primary" & !is.na(out$p_value))
  if (length(fam)) out$p_adj[fam] <- stats::p.adjust(out$p_value[fam], method = fdr)
  for (cc in c("converged", "replicated"))
    out[[cc]] <- ifelse(is.na(out[[cc]]), NA, ifelse(as.logical(out[[cc]]), "TRUE", "FALSE"))
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  say(sprintf("%d row(s) for %d predictor(s): %d regression fit(s), %d failed, %d with a p-value in the %s family -> %s",
              nrow(out), length(preds), attempted - failed, failed, length(fam), fdr, opt$out))
  invisible(NULL)
}

if (sys.nframe() == 0L) main()
