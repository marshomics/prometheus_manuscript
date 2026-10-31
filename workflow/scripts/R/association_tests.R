#!/usr/bin/env Rscript
# =====================================================================
# Stage 4.3-4.5 -- Prometheus x methanogen association in a phylogenetic
# comparative framework (design Stage 4; SCHEMA distribution/association_tests.tsv).
#
# Methanogens are polyphyletic, so genomes are not independent observations.
# The naive tests (logistic GLM with a likelihood-ratio p, Fisher's exact test)
# are reported only as labelled baselines. The defensible estimates are:
#   * phylolm::phyloglm (Ives & Garland 2010; method logistic_MPLE = Firth-
#     penalised likelihood, which stays finite under separation), unadjusted
#     and adjusted for z(completeness), z(log genome size) and
#     z(log1p(defense load EXCLUDING Prometheus)) -> the HEADLINE estimate.
#     Fitted from four starts with warm restarts, keeping the highest penalised
#     likelihood (see fit_phyloglm_ms), because a single phyloglm call can stop
#     far from the optimum while reporting success. Wald z/p (conditional on
#     alpha) plus a parametric-bootstrap percentile CI (--boot replicates, each
#     refitted with the same polishing; failures counted in `note`);
#   * phytools::fitPagel (Pagel 1994) correlated evolution, ARD, LR test (df 4);
#   * caper::phylo.d (Fritz & Purvis 2010) phylogenetic signal of presence.
#
# Presence (SCHEMA Contract v2 item 6): for tier set T (--primary-tiers) a genome
# is PRESENT if prometheus_tier is in T, ABSENT if prometheus_tier is "", and
# EXCLUDED (dropped, counted, reported in `note`) when its best validated tier
# lies outside T. Only pass_primary genomes that are tips of the species tree
# are analysed.
#
# Failure policy (Contract v2 item 4): an empty species tree ("();", fewer than
# four genomes upstream) or a master table with zero rows is a genuinely empty
# input and yields a header-only TSV with exit status 0. Data that exist but
# cannot support a test (too few genomes, no variation) yield that test's row
# with NA values and the reason in `note`. Everything else -- unreadable or
# unrooted tree, missing columns, malformed values, a missing R package, a
# crashed fit -- stops with a non-zero exit status.
#
# Output columns: test, model, term, estimate, std_error, statistic, p_value,
# ci_low, ci_high, n, note. Rows (fixed order): naive_glm, naive_fisher,
# phyloglm_methanogen, phyloglm_confounders, fitPagel, D_statistic.
# Estimates of the first four rows are log-odds for methanogen. The single
# cross-family FDR regime is applied in results_registry.py, so this table
# carries raw p-values; --fdr is accepted for interface parity only.
# =====================================================================

Sys.setenv(LANGUAGE = "en")  # English warnings/notes regardless of the user's locale
suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phylolm)
  library(phytools)
  library(caper)
})

OUT_COLS <- c("test", "model", "term", "estimate", "std_error", "statistic",
              "p_value", "ci_low", "ci_high", "n", "note")
REQUIRED_COLS <- c("genome_id", "tree_tip", "pass_primary", "prometheus_present",
                   "prometheus_tier", "methanogen", "completeness", "genome_size",
                   "n_defense_systems_excl_prometheus")
TEST_ORDER <- c("naive_glm", "naive_fisher", "phyloglm_methanogen",
                "phyloglm_confounders", "fitPagel", "D_statistic")
MODEL_LABEL <- c(
  naive_glm = "glm(present~methanogen, binomial); LRT p, profile CI",
  naive_fisher = "fisher.test(methanogen x present); conditional MLE odds ratio",
  phyloglm_methanogen = "phyloglm(present~methanogen, logistic_MPLE)",
  phyloglm_confounders = paste0("phyloglm(present~methanogen+z(completeness)+z(log(genome_size))",
                                "+z(log1p(n_defense_systems_excl_prometheus)), logistic_MPLE)"),
  fitPagel = "fitPagel(methanogen, present; ARD, dep.var=xy, fitMk)",
  D_statistic = "caper::phylo.d(present)")
TERM_LABEL <- c(naive_glm = "methanogen", naive_fisher = "methanogen",
                phyloglm_methanogen = "methanogen", phyloglm_confounders = "methanogen",
                fitPagel = "correlated_evolution", D_statistic = "D")

log_msg <- function(...) message("[association_tests] ", ...)
clean_note <- function(x) gsub("[\t\r\n]+", " ", paste(x, collapse = "; "))
fmt <- function(x, d = 4) ifelse(is.na(x), "NA", as.character(signif(x, d)))

write_header_only <- function(path, reason) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(OUT_COLS, collapse = "\t"), path)
  log_msg("header-only output: ", reason)
}

na_row <- function(test, n, note) {
  data.frame(test = test, model = MODEL_LABEL[[test]], term = TERM_LABEL[[test]],
             estimate = NA_real_, std_error = NA_real_, statistic = NA_real_,
             p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
             n = as.integer(n), note = clean_note(note), stringsAsFactors = FALSE)
}

make_row <- function(test, n, note, estimate = NA_real_, std_error = NA_real_,
                     statistic = NA_real_, p_value = NA_real_,
                     ci_low = NA_real_, ci_high = NA_real_) {
  r <- na_row(test, n, note)
  r$estimate <- as.numeric(estimate); r$std_error <- as.numeric(std_error)
  r$statistic <- as.numeric(statistic); r$p_value <- as.numeric(p_value)
  r$ci_low <- as.numeric(ci_low); r$ci_high <- as.numeric(ci_high)
  r
}

# ---- strict value parsers (no silent coercion) ---------------------------- #
# `blank`: value for empty cells (NULL = an empty cell is an error).
parse_bool <- function(x, col, blank = NULL) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(s))
  out[s %in% c("true", "t", "1", "yes", "y", "1.0")] <- TRUE
  out[s %in% c("false", "f", "0", "no", "n", "0.0")] <- FALSE
  if (!is.null(blank)) out[s == ""] <- blank
  bad <- is.na(out)
  if (any(bad))
    stop(sprintf("column '%s': %d non-boolean value(s), e.g. '%s'", col, sum(bad), s[bad][1]))
  out
}

parse_num <- function(x, col) {
  s <- trimws(as.character(x))
  miss <- s == "" | tolower(s) %in% c("na", "nan", "none", "null")
  v <- suppressWarnings(as.numeric(s))
  bad <- !miss & is.na(v)
  if (any(bad))
    stop(sprintf("column '%s': %d non-numeric value(s), e.g. '%s'", col, sum(bad), s[bad][1]))
  v[miss] <- NA_real_
  v
}

parse_tiers <- function(s, flag) {
  t <- unique(trimws(strsplit(s, ",", fixed = TRUE)[[1]]))
  t <- t[nzchar(t)]
  if (!length(t) || !all(t %in% c("A", "B", "C")))
    stop(sprintf("%s must be a comma-separated subset of A,B,C (got '%s')", flag, s))
  t
}

# Contract v2 item 6: 1 present, 0 absent, NA excluded (validated tier outside T).
tier_response <- function(tier, tier_set) {
  tier <- trimws(as.character(tier))
  ok <- tier %in% c("", "A", "B", "C")
  if (!all(ok))
    stop(sprintf("prometheus_tier: unexpected value(s) %s (allowed: '', A, B, C)",
                 paste(unique(tier[!ok]), collapse = ",")))
  out <- rep(NA_integer_, length(tier))
  out[tier == ""] <- 0L
  out[tier %in% tier_set] <- 1L
  out
}

read_tsv_chr <- function(path) {
  utils::read.delim(path, colClasses = "character", na.strings = character(0),
                    check.names = FALSE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE)
}

# ---- species tree --------------------------------------------------------- #
# Returns NULL for a genuinely empty tree; stops on anything unusable.
read_species_tree <- function(path) {
  txt <- trimws(paste(readLines(path, warn = FALSE), collapse = ""))
  if (!nzchar(txt) || grepl("^\\(\\s*\\)\\s*;?$", txt)) return(NULL)
  tr <- ape::read.tree(text = txt)
  if (is.null(tr) || !inherits(tr, "phylo"))
    stop("could not parse the species tree as a single Newick tree: ", path)
  # ape keeps Newick single quotes in labels; genome_ids never carry them.
  tr$tip.label <- sub("^'(.*)'$", "\\1", tr$tip.label)
  if (anyDuplicated(tr$tip.label)) stop("species tree has duplicated tip labels")
  if (is.null(tr$edge.length) || anyNA(tr$edge.length))
    stop("species tree lacks branch lengths (required by phyloglm, fitPagel and phylo.d)")
  if (any(tr$edge.length < 0)) stop("species tree has negative branch lengths")
  if (!ape::is.rooted(tr))
    stop("species tree is not rooted (basal polytomy); root it upstream (root_tree.py) -- ",
         "phyloglm, fitPagel and phylo.d are root-dependent")
  tr
}

# Support values (IQ-TREE '100/100' labels) are irrelevant to the comparative
# models and duplicated labels make caper::comparative.data stop. Exactly-zero
# branches make caper::phylo.d stop and give a singular VCV, so they are set to
# 1e-6 x tree height (numerically negligible for every model used here).
prepare_tree <- function(tr) {
  tr$node.label <- NULL
  h <- max(ape::node.depth.edgelength(tr))
  if (!is.finite(h) || h <= 0) stop("species tree has no positive branch lengths")
  z <- tr$edge.length == 0
  if (any(z)) {
    tr$edge.length[z] <- 1e-6 * h
    log_msg(sum(z), " zero-length branch(es) set to 1e-6 x tree height (", signif(1e-6 * h, 3), ")")
  }
  tr
}

# ---- Pagel (1994) correlated-evolution test, guarded against optimiser failure
# The independent model is the Kronecker sum of two 2-state ARD chains, so its
# MLE equals the sum of two separate 2-state fits (same flat root prior); that
# sum guards against an under-optimised independent fit (which would inflate
# LR). A dependent fit below the independent one (LR < 0, impossible at the
# optimum of nested models) is refitted from the independent solution.
# Dependent (8-rate) Pagel model refitted from the independent optimum: the
# independent Q is a point of the dependent parameter space, so the refit can only
# reach a likelihood >= the independent one. Index matrix as in phytools::fitPagel.
pagel_refit_dependent <- function(tree, fx, fy, indep_Q) {
  lev <- c("0", "1")
  dQ <- matrix(c(0, 1, 2, 0, 3, 0, 0, 4, 5, 0, 0, 6, 0, 7, 8, 0), 4, 4, byrow = TRUE)
  lv <- as.vector(sapply(lev, paste, lev, sep = "|"))
  rownames(dQ) <- colnames(dQ) <- lv
  q0 <- vapply(1:8, function(k) indep_Q[which(dQ == k)], numeric(1))
  xy <- setNames(factor(paste(fx, fy, sep = "|"), levels = lv), names(fx))
  fd <- phytools::fitMk(tree, phytools::to.matrix(xy, lv), model = dQ,
                        q.init = q0, lik.func = "pruning")
  as.numeric(fd$logLik)
}

pagel_test <- function(tree, x, y) {
  lev <- c("0", "1")
  fx <- setNames(factor(as.character(x), levels = lev), names(x))
  fy <- setNames(factor(as.character(y), levels = lev), names(y))
  fp <- phytools::fitPagel(tree, x = fx, y = fy, model = "ARD", dep.var = "xy",
                           method = "fitMk", lik.func = "pruning")
  li <- as.numeric(fp$independent.logL)
  ld <- as.numeric(fp$dependent.logL)
  flags <- character(0)
  sx <- phytools::fitMk(tree, fx, model = "ARD", lik.func = "pruning")
  sy <- phytools::fitMk(tree, fy, model = "ARD", lik.func = "pruning")
  li_sep <- as.numeric(sx$logLik) + as.numeric(sy$logLik)
  if (li_sep > li + 1e-6) {
    flags <- c(flags, sprintf("independent logL raised %.6g->%.6g from separate 2-state fits", li, li_sep))
    li <- li_sep
  }
  if (ld < li - 1e-8) {
    ld2 <- pagel_refit_dependent(tree, fx, fy, fp$independent.Q)
    flags <- c(flags, sprintf("dependent model refitted from the independent optimum (logL %.6g->%.6g)", ld, ld2))
    ld <- max(ld, ld2)
    if (ld < li - 1e-8) flags <- c(flags, "dependent logL still below independent; LR truncated at 0")
  }
  lr <- max(0, 2 * (ld - li))
  list(lr = lr, p = stats::pchisq(lr, df = 4, lower.tail = FALSE),
       logL_independent = li, logL_dependent = ld,
       p_fitpagel = as.numeric(fp$P), flags = flags)
}

# ---- phylogenetic logistic regression: multi-start MPLE ------------------- #
# phyloglm(logistic_MPLE) optimises the Firth-penalised likelihood with L-BFGS-B at a
# loose tolerance (factr = 1e12) behind a hard barrier at the alpha bounds. Its
# convergence code is therefore not a reliable optimality check: on simulated data a
# default fit can report code 0 at a point many log-units below the optimum, or stop
# with code 52 (line search failed at the alpha barrier) with beta far from optimal.
# The estimate reported here is the highest penalised log-likelihood (penlogLik, the
# MPLE objective, comparable across fits on the same data and tree) over four starts
# -- alpha*T = exp(-3), 1 (phyloglm's default), exp(3), and the logistic_IG10
# solution -- each polished by warm restarts until the gain is < 1e-4 log-units.
pglm_once <- function(fml, d, tr, btol, method = "logistic_MPLE", sb = NULL, sa = NULL) {
  w <- character(0)
  f <- withCallingHandlers(
    phylolm::phyloglm(fml, data = d, phy = tr, method = method, btol = btol,
                      start.beta = sb, start.alpha = sa, boot = 0),
    warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  f$warn_msgs <- unique(w)
  f
}

restartable <- function(f, btol) {
  lp <- as.vector(f$X %*% f$coefficients)
  all(is.finite(c(lp, f$alpha, f$penlogLik))) && max(abs(lp)) < 0.98 * btol
}

clamp_alpha <- function(a, tmax, b = 3.9) min(max(a, exp(-b) / tmax), exp(b) / tmax)

polish <- function(f, fml, d, tr, btol, max_restarts) {
  cur <- f; k <- 0L; gain <- NA_real_
  while (k < max_restarts && restartable(cur, btol)) {
    k <- k + 1L
    nf <- tryCatch(pglm_once(fml, d, tr, btol, sb = as.numeric(cur$coefficients),
                             sa = clamp_alpha(cur$alpha, cur$mean.tip.height)),
                   error = function(e) NULL)
    if (is.null(nf) || !is.finite(nf$penlogLik)) break
    gain <- nf$penlogLik - cur$penlogLik
    if (gain > 0) cur <- nf
    if (gain < 1e-4) break
  }
  list(fit = cur, restarts = k, gain = gain)
}

fit_phyloglm_ms <- function(fml, d, tr, btol, max_restarts = 10L) {
  notes <- character(0)
  f0 <- pglm_once(fml, d, tr, btol)                  # an error here is a genuine failure
  tmax <- f0$mean.tip.height
  cands <- list(`aT=1` = polish(f0, fml, d, tr, btol, max_restarts))
  for (lab in c("aT=exp(-3)", "aT=exp(3)")) {
    a0 <- if (lab == "aT=exp(3)") exp(3) / tmax else exp(-3) / tmax
    f <- tryCatch(pglm_once(fml, d, tr, btol, sa = a0), error = function(e) {
      notes <<- c(notes, sprintf("start %s failed: %s", lab, conditionMessage(e))); NULL })
    if (!is.null(f)) cands[[lab]] <- polish(f, fml, d, tr, btol, max_restarts)
  }
  fi <- tryCatch(pglm_once(fml, d, tr, btol, method = "logistic_IG10"), error = function(e) {
    notes <<- c(notes, paste0("IG10 start unavailable: ", conditionMessage(e))); NULL })
  if (!is.null(fi) && restartable(fi, btol)) {
    f <- tryCatch(pglm_once(fml, d, tr, btol, sb = as.numeric(fi$coefficients),
                            sa = clamp_alpha(fi$alpha, tmax)), error = function(e) NULL)
    if (!is.null(f)) cands[["IG10"]] <- polish(f, fml, d, tr, btol, max_restarts)
  }
  pens <- vapply(cands, function(cn) cn$fit$penlogLik, numeric(1))
  ib <- which.max(pens)
  best <- cands[[ib]]
  srt <- sort(pens, decreasing = TRUE)
  list(fit = best$fit, start = names(cands)[ib], pens = pens, restarts = best$restarts,
       gain = best$gain, runner_up_gap = if (length(srt) > 1) srt[1] - srt[2] else NA_real_,
       notes = notes)
}

# Parametric bootstrap under the fitted model (as phyloglm's own `boot`, but every
# replicate is refitted from the estimate with the same warm-restart polishing; the
# number of replicates that could not be fitted is reported).
boot_phyloglm <- function(ms, fml, d, tr, btol, term, B, seed, threads) {
  best <- ms$fit
  set.seed(seed)
  Y <- phylolm::rbinTrait(n = B, phy = tr, beta = as.numeric(best$coefficients),
                          alpha = best$alpha, X = best$X, model = "LogReg")
  Y <- as.matrix(Y)[match(rownames(d), rownames(as.matrix(Y))), , drop = FALSE]
  resp <- all.vars(fml)[1]
  one <- function(j) {
    dj <- d; dj[[resp]] <- as.integer(Y[, j])
    if (length(unique(dj[[resp]])) < 2) return(NA_real_)
    f <- tryCatch(pglm_once(fml, dj, tr, btol, sb = as.numeric(best$coefficients),
                            sa = clamp_alpha(best$alpha, best$mean.tip.height)),
                  error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    unname(polish(f, fml, dj, tr, btol, 3L)$fit$coefficients[term])
  }
  res <- if (threads > 1) parallel::mclapply(seq_len(B), one, mc.cores = threads)
         else lapply(seq_len(B), one)
  crashed <- vapply(res, function(r) inherits(r, "try-error") || is.null(r) || length(r) != 1, logical(1))
  if (any(crashed)) stop("bootstrap worker crashed: ", paste(unlist(res[crashed])[1], collapse = ""))
  v <- unlist(res)
  ok <- is.finite(v)
  list(n_ok = sum(ok), ci = if (sum(ok) >= max(10, ceiling(0.5 * B))) stats::quantile(v[ok], c(0.025, 0.975), names = FALSE)
                            else c(NA_real_, NA_real_))
}

run_phyloglm <- function(fml, d, tr, term, btol, boot, seed, threads) {
  set.seed(seed)   # the logistic_IG10 start draws random numbers (plogregfunct)
  ms <- fit_phyloglm_ms(fml, d, tr, btol)
  fit <- ms$fit
  co <- summary(fit)$coefficients
  if (!term %in% rownames(co)) stop("phyloglm: term '", term, "' missing from coefficients")
  est <- co[term, "Estimate"]; se <- co[term, "StdErr"]
  z <- co[term, "z.value"]; p <- co[term, "p.value"]
  if (boot > 0) {
    bt <- boot_phyloglm(ms, fml, d, tr, btol, term, boot, seed, threads)
    lo <- bt$ci[1]; hi <- bt$ci[2]
    ci_type <- if (is.na(lo)) sprintf("bootstrap CI not reported: only %d/%d replicates fitted", bt$n_ok, boot)
               else sprintf("parametric-bootstrap percentile CI (%d/%d replicates fitted)", bt$n_ok, boot)
  } else {
    q <- stats::qnorm(0.975)
    lo <- est - q * se; hi <- est + q * se
    ci_type <- "Wald CI (conditional on alpha; --boot 0)"
  }
  stationary <- isTRUE(fit$convergence == 0) || (is.finite(ms$gain) && ms$gain < 1e-4)
  flags <- c(sprintf("multi-start MPLE: best of %d starts = %s (penlogLik %.6g; gap to runner-up %s); %d warm restart(s); optim code %s",
                     length(ms$pens), ms$start, max(ms$pens), fmt(ms$runner_up_gap), ms$restarts,
                     fit$convergence), ms$notes)
  if (!stationary)
    flags <- c(flags, "CAUTION: best fit still improving at the restart limit (not numerically stationary)")
  if (identical(as.integer(fit$alphaWarn), 1L))
    flags <- c(flags, "alpha at LOWER bound (maximal phylogenetic correlation; flat likelihood)")
  if (identical(as.integer(fit$alphaWarn), 2L))
    flags <- c(flags, "alpha at UPPER bound (negligible phylogenetic correlation)")
  lp <- as.vector(fit$X %*% fit$coefficients)
  if (max(abs(lp)) > 0.99 * btol)
    flags <- c(flags, sprintf("linear predictor at btol=%g (separation-like; estimate is bounded)", btol))
  for (w in fit$warn_msgs) log_msg("phyloglm warning (best fit): ", gsub("\\s+", " ", w))
  list(est = est, se = se, z = z, p = p, lo = lo, hi = hi, ci_type = ci_type,
       alpha = fit$alpha, tmax = fit$mean.tip.height, flags = flags, n = fit$n,
       converged = stationary)
}

zscore <- function(v) {
  s <- stats::sd(v)
  if (!is.finite(s) || s == 0) return(NULL)
  (v - mean(v)) / s
}

# --------------------------------------------------------------------------- #
main <- function() {
  opt <- parse_args(OptionParser(option_list = list(
    make_option("--master", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--primary-tiers", type = "character", default = "A", dest = "primary_tiers"),
    make_option("--fdr", type = "character", default = "BH",
                help = "accepted for interface parity; FDR is applied in results_registry.py"),
    make_option("--btol", type = "double", default = 10),
    make_option("--boot", type = "integer", default = 200),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--d-permutations", type = "integer", default = 1000, dest = "d_permutations"),
    make_option("--threads", type = "integer", default = 1),
    make_option("--out", type = "character"))))
  for (k in c("master", "tree", "out"))
    if (is.null(opt[[k]]) || !nzchar(opt[[k]])) stop("--", k, " is required")
  if (opt$boot < 0) stop("--boot must be >= 0")
  if (opt$d_permutations < 99) stop("--d-permutations must be >= 99")
  if (!is.finite(opt$btol) || opt$btol <= 0) stop("--btol must be > 0")
  tiers <- parse_tiers(opt$primary_tiers, "--primary-tiers")
  tier_lab <- paste(tiers, collapse = "+")
  if (opt$threads < 1) stop("--threads must be >= 1")

  # ---- inputs ------------------------------------------------------------ #
  tree0 <- read_species_tree(opt$tree)
  if (is.null(tree0))
    return(write_header_only(opt$out, "species tree is empty ('();'): fewer than 4 genomes upstream"))  # nonfatal: empty input
  tree0 <- prepare_tree(tree0)
  master <- read_tsv_chr(opt$master)
  miss <- setdiff(REQUIRED_COLS, names(master))
  if (length(miss)) stop("master table lacks required column(s): ", paste(miss, collapse = ", "))
  if (nrow(master) == 0)
    return(write_header_only(opt$out, "master table has zero genomes"))  # nonfatal: empty input
  gid <- trimws(master$genome_id)
  if (any(!nzchar(gid))) stop("master table has blank genome_id values")
  if (anyDuplicated(gid)) stop("master table has duplicated genome_id: ", gid[duplicated(gid)][1])
  tip <- trimws(master$tree_tip)
  if (anyDuplicated(tip[nzchar(tip)])) stop("master table has duplicated tree_tip values")

  y_all <- tier_response(master$prometheus_tier, tiers)
  pres_flag <- parse_bool(master$prometheus_present, "prometheus_present")
  incons <- pres_flag != (trimws(master$prometheus_tier) != "")
  if (any(incons))
    stop(sprintf("master: %d genome(s) where prometheus_present disagrees with prometheus_tier (e.g. %s)",
                 sum(incons), gid[incons][1]))
  pass <- parse_bool(master$pass_primary, "pass_primary", blank = FALSE)
  meth_all <- parse_bool(master$methanogen, "methanogen")
  anme_all <- if ("is_anme" %in% names(master))
    parse_bool(master$is_anme, "is_anme", blank = FALSE) else rep(FALSE, nrow(master))
  on_tree <- nzchar(tip) & tip %in% tree0$tip.label
  cand <- pass & on_tree
  excl <- cand & is.na(y_all)
  use <- cand & !is.na(y_all)
  n_excl <- sum(excl)
  log_msg(sprintf(paste0("master=%d genomes; pass_primary=%d; pass_primary on tree=%d; ",
                         "EXCLUDED (best tier outside %s)=%d; analysed=%d; tree tips not analysed=%d"),
                  nrow(master), sum(pass), sum(cand), tier_lab, n_excl, sum(use),
                  length(tree0$tip.label) - sum(use)))
  base_note <- sprintf("tier set %s; n_excluded_tier=%d (best tier outside set)", tier_lab, n_excl)

  d <- data.frame(tip = tip[use], present = y_all[use], methanogen = as.integer(meth_all[use]),
                  anme = anme_all[use],
                  completeness = parse_num(master$completeness, "completeness")[use],
                  genome_size = parse_num(master$genome_size, "genome_size")[use],
                  defense_excl = parse_num(master$n_defense_systems_excl_prometheus,
                                           "n_defense_systems_excl_prometheus")[use],
                  stringsAsFactors = FALSE)
  n <- nrow(d)
  if (sum(d$anme) > 0) base_note <- paste0(base_note, sprintf("; n_anme=%d (coded non-methanogen)", sum(d$anme)))

  out_rows <- list()
  if (n < 5) {
    why <- sprintf("not estimable: %d analysable genome(s) (<5)", n)
    for (t in TEST_ORDER) out_rows[[t]] <- na_row(t, n, c(why, base_note))
    return(finish(out_rows, opt$out, n))
  }
  tr <- ape::keep.tip(tree0, d$tip)
  if (!ape::is.rooted(tr)) stop("pruned species tree lost its root")
  d <- d[match(tr$tip.label, d$tip), , drop = FALSE]
  rownames(d) <- d$tip
  n_pres <- sum(d$present); n_meth <- sum(d$methanogen)
  counts <- sprintf("n_present=%d; n_methanogen=%d", n_pres, n_meth)
  log_msg(sprintf("analysis set: n=%d, %s", n, counts))

  if (n_pres == 0 || n_pres == n) {
    why <- sprintf("not estimable: no variation in Prometheus presence (%s)", counts)
    for (t in TEST_ORDER) out_rows[[t]] <- na_row(t, n, c(why, base_note))
    return(finish(out_rows, opt$out, n))
  }
  meth_ok <- n_meth > 0 && n_meth < n

  # ---- (a) naive logistic GLM (assumes independence) --------------------- #
  tab <- table(methanogen = factor(d$methanogen, levels = 0:1),
               present = factor(d$present, levels = 0:1))
  tab_txt <- sprintf("2x2 [non-methanogen: absent=%d present=%d | methanogen: absent=%d present=%d]",
                     tab[1, 1], tab[1, 2], tab[2, 1], tab[2, 2])
  if (!meth_ok) {
    why <- sprintf("not estimable: no variation in methanogen status (%s)", counts)
    for (t in c("naive_glm", "naive_fisher", "phyloglm_methanogen", "phyloglm_confounders", "fitPagel"))
      out_rows[[t]] <- na_row(t, n, c(why, base_note))
  } else {
    gw <- character(0)
    fit <- withCallingHandlers(glm(present ~ methanogen, family = binomial(), data = d),
                               warning = function(w) { gw <<- c(gw, conditionMessage(w)); invokeRestart("muffleWarning") })
    lr <- fit$null.deviance - fit$deviance
    p_lr <- stats::pchisq(lr, df = 1, lower.tail = FALSE)
    if (any(tab == 0)) {
      out_rows$naive_glm <- make_row("naive_glm", n,
        c("assumes independence; baseline only",
          "quasi-complete separation (zero cell): MLE log-odds is infinite, Wald invalid; LRT p reported",
          tab_txt, base_note),
        statistic = lr, p_value = p_lr)
    } else {
      sm <- summary(fit)$coefficients
      ci <- suppressMessages(stats::confint(fit, "methanogen"))
      out_rows$naive_glm <- make_row("naive_glm", n,
        c("assumes independence; baseline only; log-odds; LRT p; profile-likelihood CI",
          tab_txt, base_note, if (length(gw)) paste("glm warnings:", paste(unique(gw), collapse = " | "))),
        estimate = sm["methanogen", "Estimate"], std_error = sm["methanogen", "Std. Error"],
        statistic = lr, p_value = p_lr, ci_low = ci[1], ci_high = ci[2])
    }
    # ---- (a2) Fisher's exact test (design 4.3) ---------------------------- #
    ft <- stats::fisher.test(tab)
    or <- as.numeric(ft$estimate); oci <- as.numeric(ft$conf.int)
    lg <- function(v) if (is.finite(v) && v > 0) log(v) else NA_real_
    out_rows$naive_fisher <- make_row("naive_fisher", n,
      c(sprintf("assumes independence; baseline only; log conditional-MLE odds ratio (OR=%s, exact 95%% CI %s-%s)",
                fmt(or), fmt(oci[1]), fmt(oci[2])),
        if (!is.finite(lg(or))) "OR is 0 or Inf (zero cell): log-OR not finite", tab_txt, base_note),
      estimate = lg(or), p_value = ft$p.value, ci_low = lg(oci[1]), ci_high = lg(oci[2]))

    # ---- (b) phyloglm, methanogen only ------------------------------------ #
    r <- run_phyloglm(present ~ methanogen, d, tr, "methanogen", opt$btol, opt$boot, opt$seed, opt$threads)
    out_rows$phyloglm_methanogen <- make_row("phyloglm_methanogen", n,
      c(sprintf("phylogenetic logistic regression (Ives & Garland 2010, Firth-penalised MPLE); log-odds; Wald z/p conditional on alpha=%s (alpha*T=%s)",
                fmt(r$alpha), fmt(r$alpha * r$tmax)), r$ci_type, r$flags, counts, base_note),
      estimate = r$est, std_error = r$se, statistic = r$z, p_value = r$p, ci_low = r$lo, ci_high = r$hi)

    # ---- (c) phyloglm + confounders: HEADLINE ----------------------------- #
    dc <- d
    dc$log_gs <- ifelse(is.finite(dc$genome_size) & dc$genome_size > 0, log(dc$genome_size), NA_real_)
    dc$log_def <- ifelse(is.finite(dc$defense_excl) & dc$defense_excl >= 0, log1p(dc$defense_excl), NA_real_)
    keep <- stats::complete.cases(dc$completeness, dc$log_gs, dc$log_def)
    n_drop <- sum(!keep)
    dc <- dc[keep, , drop = FALSE]
    cov_note <- sprintf("n_dropped_missing_covariate=%d", n_drop)
    nc <- nrow(dc)
    if (nc < 10 || sum(dc$present) %in% c(0, nc) || sum(dc$methanogen) %in% c(0, nc)) {
      out_rows$phyloglm_confounders <- na_row("phyloglm_confounders", nc,
        c(sprintf("not estimable on complete-covariate subset (n=%d, present=%d, methanogen=%d)",
                  nc, sum(dc$present), sum(dc$methanogen)), cov_note, base_note))
    } else {
      covs <- list(z_completeness = zscore(dc$completeness), z_log_genome_size = zscore(dc$log_gs),
                   z_log1p_defense_excl = zscore(dc$log_def))
      dropped <- names(covs)[vapply(covs, is.null, logical(1))]
      covs <- covs[!vapply(covs, is.null, logical(1))]
      for (k in names(covs)) dc[[k]] <- covs[[k]]
      fml <- stats::as.formula(paste("present ~ methanogen", paste(c("", names(covs)), collapse = " + ")))
      trc <- if (nc < n) ape::keep.tip(tr, dc$tip) else tr
      dc <- dc[match(trc$tip.label, dc$tip), , drop = FALSE]
      rownames(dc) <- dc$tip
      r <- run_phyloglm(fml, dc, trc, "methanogen", opt$btol, opt$boot, opt$seed, opt$threads)
      out_rows$phyloglm_confounders <- make_row("phyloglm_confounders", nc,
        c(sprintf("HEADLINE: methanogen log-odds adjusted for %s (defense load EXCLUDES Prometheus); Firth-penalised MPLE; Wald z/p conditional on alpha=%s (alpha*T=%s)",
                  paste(names(covs), collapse = ", "), fmt(r$alpha), fmt(r$alpha * r$tmax)),
          if (length(dropped)) paste0("zero-variance covariate(s) dropped: ", paste(dropped, collapse = ",")),
          r$ci_type, r$flags, cov_note, base_note),
        estimate = r$est, std_error = r$se, statistic = r$z, p_value = r$p, ci_low = r$lo, ci_high = r$hi)
    }

    # ---- (d) Pagel (1994) correlated evolution ---------------------------- #
    pg <- pagel_test(tr, setNames(d$methanogen, d$tip), setNames(d$present, d$tip))
    out_rows$fitPagel <- make_row("fitPagel", n,
      c(sprintf("Pagel 1994 dependent vs independent (ARD, flat root prior); LR statistic, chi-square df=4; logL_independent=%.6g; logL_dependent=%.6g",
                pg$logL_independent, pg$logL_dependent), pg$flags, counts, base_note),
      statistic = pg$lr, p_value = pg$p)
  }

  # ---- (e) Fritz & Purvis D (needs only presence variation) --------------- #
  cd_df <- data.frame(tip = d$tip, present = d$present, stringsAsFactors = FALSE)
  cdat <- caper::comparative.data(phy = tr, data = cd_df, names.col = tip)
  set.seed(opt$seed)
  pdres <- caper::phylo.d(cdat, binvar = present, permut = opt$d_permutations)
  obs <- as.numeric(pdres$Parameters$Observed)
  ran <- as.numeric(pdres$Permutations$random); bro <- as.numeric(pdres$Permutations$brownian)
  p_ran <- (sum(ran <= obs) + 1) / (length(ran) + 1)
  p_bro <- (sum(bro >= obs) + 1) / (length(bro) + 1)
  out_rows$D_statistic <- make_row("D_statistic", n,
    c(sprintf(paste0("Fritz-Purvis D (0 = Brownian-threshold clumping, 1 = random); p_value = P(D_random <= D_obs), ",
                     "(k+1)/(N+1) with N=%d permutations; P vs Brownian (D>0) = %s; caper Pval1=%s Pval0=%s"),
              opt$d_permutations, fmt(p_bro), fmt(pdres$Pval1), fmt(pdres$Pval0)),
      sprintf("n_present=%d", n_pres), base_note),
    estimate = as.numeric(pdres$DEstimate), p_value = p_ran)

  finish(out_rows, opt$out, n)
}

finish <- function(out_rows, path, n) {
  out <- do.call(rbind, out_rows[TEST_ORDER[TEST_ORDER %in% names(out_rows)]])
  out <- out[, OUT_COLS, drop = FALSE]
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  log_msg(sprintf("wrote %d rows (n=%d genomes) -> %s", nrow(out), n, path))
  invisible(NULL)
}

# Run only under Rscript; source() exposes the helpers for unit tests.
if (sys.nframe() == 0L) main()
