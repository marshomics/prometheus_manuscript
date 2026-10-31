#!/usr/bin/env Rscript
# =====================================================================
# Stage 6 -- ecology: does Prometheus presence track niche / host association or viral
# pressure BEYOND what the phylogeny explains (design Stage 6; confounder #11: biome and
# clade are deeply confounded in archaea)? Every predictor is tested in the same
# phylogenetic logistic framework as Stage 4 (phylolm::phyloglm, Ives & Garland 2010,
# method "logistic_MPLE" = Firth-penalised likelihood, which keeps estimates finite under
# separation), one model per predictor:
#
#   present ~ predictor + z(completeness) + z(genome_size)
#
# Response (SCHEMA Contract v2 item 6): for the tier set T (--primary-tiers, default A)
#   1 = best validated tier (prometheus_tier) in T; 0 = no validated locus (tier "");
#   EXCLUDED = a validated locus whose best tier is outside T (counted in n_excluded_tier
#   and the log, never coded as a confident absence). Genomes must pass_primary and have a
#   tree tip present in the species tree.
#
# Predictors
#   binary  (--predictors; default "auto" = every ^(group|niche)_ column of the master =
#           species-level presence/absence, README "the niche columns become the ecology
#           predictors"). A comma list may name exact columns and/or globs ("niche_*");
#           a named column absent from the master gets an explicit "not found" row. Values
#           1/0 (true/false/yes/no); blank/NA = unknown (dropped), anything else = the
#           predictor is skipped with a note. level = "1".
#   continuous (only with --viral): n_proviruses as z(log1p(count)) and
#           frac_rloop_proviruses as z(fraction) -- estimates are log-odds per SD;
#           frac_rloop_proviruses is undefined (NA) for genomes without proviruses, so that
#           model is fitted on genomes with >= 1 provirus. level = "continuous".
#   Covariates (--covariates, default completeness,genome_size; e.g. add temperature-proxy
#   niche columns for temperature-aware niche tests) are z-scored within each model; a
#   covariate missing from the master is reported, never silently dropped, and a predictor
#   is never its own covariate.
#   frac_rloop_proviruses is tested two-sided on purpose: Prometheus could be retained under
#   R-loop-forming phage pressure (positive) or purge R-loop-forming proviruses (negative).
# Sensitivity columns (sens_*; not part of the BH family): binary -> + methanogen status
#   (partial effect, when the master has a varying methanogen column); n_proviruses ->
#   all geNomad viral sequences (proviruses + whole viral contigs; MAG fragmentation);
#   frac_rloop_proviruses -> + z(mean provirus GC) + z(log mean provirus length)
#   (QmRLFS calls are G-content and length dependent).
# A predictor is fitted only with n >= --min-n genomes and >= --min-class genomes in each
#   response class and (binary) in each predictor class; otherwise its row says why.
#   A zero cell of the 2x2 table is flagged as (quasi-)separation. phyloglm warnings
#   (alpha at a bound, btol reached) are copied into the note.
# Convergence: phyloglm's L-BFGS-B (factr 1e12) can stop abnormally (code 52) at a WRONG
#   beta when the likelihood is flat in alpha (negligible phylogenetic signal) -- on an
#   exactly balanced null (glm log-OR 0) it reported log-OR 1.57, p = 0.0016. So every model
#   is fitted from --mple-starts start values of alpha (exp(k)/T, T = mean root-to-tip depth
#   = phyloglm's Tmax) and only CONVERGED fits are used, the one with the highest penalised
#   log-likelihood. If no start converges, the Ives & Garland (2010) estimator
#   (logistic_IG10) is used if it converges (flagged in model + note); else the fit failed.
# p_adj: --fdr (BH) across every reported term with a p-value.
#
# Output (SCHEMA ecology/ecology_associations.tsv): predictor, level, model, estimate,
# std_error, p_value, p_adj, n, note + ci_low, ci_high (Wald 95 %), n_present, n_absent,
# n_level1, n_excluded_tier, alpha, converged, sens_model, sens_estimate, sens_std_error,
# sens_p_value, sens_n.
# The species tree must be rooted (the pipeline roots it upstream; an unrooted tree is an
# error, never silently midpoint-rooted).
# Header-only output only for a genuinely empty input (empty/"();" tree, empty master, or
# no predictor at all). Missing packages, unreadable inputs, and every fit failing are
# errors (non-zero exit).
# =====================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phylolm)
})

OUT_COLS <- c("predictor", "level", "model", "estimate", "std_error", "p_value", "p_adj",
              "n", "note", "ci_low", "ci_high", "n_present", "n_absent", "n_level1",
              "n_excluded_tier", "alpha", "converged", "sens_model", "sens_estimate",
              "sens_std_error", "sens_p_value", "sens_n")
CONT_PREDICTORS <- c("n_proviruses", "frac_rloop_proviruses")
TRUE_S <- c("1", "true", "yes", "t")
FALSE_S <- c("0", "false", "no", "f")
NA_S <- c("", "na", "nan", "none")

say <- function(...) message("ecology_stats: ", ...)

parse_list <- function(x) {
  if (is.null(x)) return(character(0))
  v <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  v[nzchar(v)]
}

read_tsv_chr <- function(path, what) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) stop(what, " not found: ", path)
  if (file.info(path)$size == 0) stop(what, " is an empty file (no header): ", path)
  utils::read.delim(path, colClasses = "character", na.strings = character(0),
                    check.names = FALSE, quote = "", comment.char = "")
}

flag <- function(x) {                       # TRUE / FALSE / NA
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(s))
  out[s %in% TRUE_S] <- TRUE
  out[s %in% FALSE_S] <- FALSE
  out
}

num <- function(x) {
  s <- trimws(as.character(x))
  s[tolower(s) %in% NA_S] <- NA
  suppressWarnings(as.numeric(s))
}

zscore <- function(x) {
  s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(NULL)
  (x - mean(x)) / s
}

clean <- function(x) gsub("[\t\r\n]+", " ", x)

join_notes <- function(...) {
  v <- unlist(list(...))
  v <- v[!is.na(v) & nzchar(trimws(v))]
  clean(paste(v, collapse = "; "))
}

header_only <- function(path, why) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(OUT_COLS, collapse = "\t"), path)
  say("header-only output (", why, ") -> ", path)
}

blank_row <- function(predictor, level, model = "phyloglm(present ~ predictor, logistic_MPLE)", note = "") {
  r <- as.list(rep(NA, length(OUT_COLS))); names(r) <- OUT_COLS
  r$predictor <- predictor; r$level <- level; r$model <- model; r$note <- note
  r$sens_model <- ""
  r
}

# ---- presence rule (Contract v2 item 6) ---------------------------------------------
presence_rule <- function(m, tiers) {
  if ("prometheus_tier" %in% names(m)) {
    tier <- toupper(trimws(m$prometheus_tier))
    y <- ifelse(tier %in% tiers, 1L, ifelse(tier == "", 0L, NA_integer_))
    if ("prometheus_present" %in% names(m)) {
      odd <- sum(tier == "" & flag(m$prometheus_present) %in% TRUE)
      if (odd) say("WARNING: ", odd, " genome(s) have prometheus_present=True but no tier; ",
                   "coded absent by the tier rule")
    }
    return(list(y = y, basis = paste0("tiers=", paste(tiers, collapse = "+"))))
  }
  if ("prometheus_present" %in% names(m)) {
    say("WARNING: master has no prometheus_tier column; presence = prometheus_present ",
        "(any validated locus) and the tier rule cannot be applied")
    f <- flag(m$prometheus_present)
    return(list(y = ifelse(is.na(f), NA_integer_, as.integer(f)),
                basis = "prometheus_present (no tier column)"))
  }
  stop("master has neither prometheus_tier nor prometheus_present")
}

resolve_predictors <- function(spec, cols) {
  toks <- parse_list(spec)
  if (!length(toks)) toks <- "auto"
  found <- character(0); missing <- character(0)
  for (t in toks) {
    if (tolower(t) == "auto") {
      hit <- grep("^(group|niche)_", cols, value = TRUE)
    } else if (grepl("[*?]", t)) {
      hit <- grep(utils::glob2rx(t), cols, value = TRUE)
    } else {
      hit <- if (t %in% cols) t else character(0)
      if (!length(hit)) missing <- c(missing, t)
    }
    found <- c(found, hit)
  }
  found <- setdiff(unique(found), CONT_PREDICTORS)
  ord <- c(sort(grep("^group_", found, value = TRUE)), sort(grep("^niche_", found, value = TRUE)))
  list(found = c(ord, setdiff(found, ord)), missing = unique(missing))
}

as_binary <- function(x) {
  s <- tolower(trimws(as.character(x)))
  v <- rep(NA_integer_, length(s))
  v[s %in% TRUE_S] <- 1L
  v[s %in% FALSE_S] <- 0L
  bad <- !(s %in% c(TRUE_S, FALSE_S, NA_S))
  list(v = v, bad = unique(as.character(x)[bad]))
}

# Short, stable wording for the phyloglm (phylolm 2.6.5) warnings copied into the notes.
short_warning <- function(w) {
  w <- clean(trimws(w))
  if (grepl("reached the upper bound", w, fixed = TRUE))
    return("alpha at its upper bound (phylogenetic signal negligible)")
  if (grepl("reached the lower bound", w, fixed = TRUE))
    return("alpha at its lower bound (very strong phylogenetic signal)")
  if (grepl("beyond 'btol'", w, fixed = TRUE))
    return("ordinary-GLM start beyond btol (separation); started from beta=0")
  if (grepl("boundary of the linear predictor", w, fixed = TRUE) ||
      grepl("linear predictor reaches its bound", w, fixed = TRUE))
    return("linear predictor reached btol")
  if (grepl("failed to converge", w, fixed = TRUE)) return("phyloglm did not converge")
  substr(w, 1, 160)
}

# ---- one phyloglm model: multi-start MPLE, converged fits only, IG10 fallback ----------
phyloglm_once <- function(f, d, tr, method, btol, start_alpha = NULL) {
  warns <- character(0)
  fit <- tryCatch(withCallingHandlers(
    phylolm::phyloglm(f, data = d, phy = tr, method = method, btol = btol,
                      start.alpha = start_alpha),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) return(list(ok = FALSE, error = conditionMessage(fit)))
  co <- summary(fit)$coefficients
  if (!"pred" %in% rownames(co)) return(list(ok = FALSE, error = "predictor term absent from the fit"))
  est <- unname(co["pred", "Estimate"]); se <- unname(co["pred", "StdErr"])
  list(ok = isTRUE(fit$convergence == 0) && is.finite(est) && is.finite(se) && se > 0,
       est = est, se = se, p = unname(co["pred", "p.value"]), alpha = fit$alpha,
       pll = if (is.null(fit$penlogLik)) NA_real_ else fit$penlogLik,
       code = fit$convergence, warns = warns, error = "")
}

fit_phyloglm <- function(d, terms, tree, btol, n_starts = 5L) {
  tr <- ape::keep.tip(tree, d$tip)
  rownames(d) <- d$tip
  d <- d[tr$tip.label, , drop = FALSE]
  f <- stats::as.formula(paste("y ~", paste(terms, collapse = " + ")))
  tmax <- mean(ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))])   # phyloglm's Tmax
  ks <- c(0, -3, 3, -1.5, 1.5, -2.25, 2.25, -0.75, 0.75)[seq_len(max(1L, min(9L, n_starts)))]
  fits <- lapply(ks, function(k)
    phyloglm_once(f, d, tr, "logistic_MPLE", btol, if (k == 0) NULL else exp(k) / tmax))
  ok <- vapply(fits, function(x) isTRUE(x$ok), logical(1))
  notes <- character(0)
  if (any(ok)) {
    cand <- fits[ok]
    pll <- vapply(cand, function(x) x$pll, numeric(1))
    best <- cand[[if (all(is.na(pll))) 1L else which.max(pll)]]
    near <- cand[!is.na(pll) & pll >= max(pll, na.rm = TRUE) - 0.5]
    ests <- vapply(near, function(x) x$est, numeric(1))
    if (length(ests) > 1 && diff(range(ests)) > best$se)
      notes <- c(notes, sprintf("flat likelihood: equally good starts give log-odds %.2f..%.2f",
                                min(ests), max(ests)))
    best$method <- "logistic_MPLE"
    best$detail <- sprintf("MPLE converged from %d/%d starts", sum(ok), length(fits))
  } else {
    codes <- vapply(fits, function(x) if (is.null(x$code)) NA_real_ else as.numeric(x$code), numeric(1))
    ig <- phyloglm_once(f, d, tr, "logistic_IG10", btol)
    if (!isTRUE(ig$ok)) {
      errs <- unique(c(vapply(fits, function(x) x$error, character(1)), ig$error))
      stop(sprintf("no converged fit (MPLE codes %s; IG10 %s)%s",
                   paste(unique(codes), collapse = ","),
                   if (nzchar(ig$error)) ig$error else paste0("code ", ig$code),
                   if (any(nzchar(errs))) paste0(": ", paste(errs[nzchar(errs)], collapse = " | ")) else ""))
    }
    best <- ig
    best$method <- "logistic_IG10"
    best$detail <- sprintf(paste0("MPLE did not converge from %d starts (codes %s): Ives & ",
                                  "Garland 2010 estimator (logistic_IG10) used"),
                           length(fits), paste(unique(codes), collapse = ","))
  }
  best$warns <- unique(vapply(best$warns, short_warning, character(1), USE.NAMES = FALSE))
  best$warns <- setdiff(best$warns, "phyloglm did not converge")
  best$notes <- notes
  best
}

# Prepare the per-model data: complete cases on y, pred and the usable covariates.
prepare <- function(d, covars, min_n, min_class) {
  d <- d[!is.na(d$y) & !is.na(d$pred), , drop = FALSE]
  used <- character(0); notes <- character(0)
  for (cv in covars) {
    v <- d[[cv]]
    if (all(is.na(v))) { notes <- c(notes, paste0(cv, " all NA: not adjusted")); next }
    d <- d[!is.na(v), , drop = FALSE]
    used <- c(used, cv)
  }
  terms <- character(0)
  for (cv in used) {
    z <- zscore(d[[cv]])
    if (is.null(z)) { notes <- c(notes, paste0(cv, " constant: not adjusted")); next }
    d[[paste0("z_", cv)]] <- z
    terms <- c(terms, paste0("z_", cv))
  }
  list(d = d, terms = terms, notes = notes)
}

run_model <- function(base, pname, level, pred, covars, tree, opt, binary,
                      pred_label = pname, extra_note = "") {
  r <- blank_row(pname, level)
  covars <- setdiff(covars, pname)
  d <- base; d$pred <- pred
  pp <- prepare(d, covars, opt$min_n, opt$min_class)
  d <- pp$d
  n <- nrow(d); np <- sum(d$y == 1L); na_ <- n - np
  r$n <- n; r$n_present <- np; r$n_absent <- na_
  r$model <- sprintf("phyloglm(present ~ %s%s, logistic_MPLE)", pred_label,
                     if (length(pp$terms))
                       paste0(" + ", paste0("z(", sub("^z_", "", pp$terms), ")", collapse = " + "))
                     else "")
  notes <- c(extra_note, pp$notes)
  if (binary) r$n_level1 <- sum(d$pred == 1L)
  why <- character(0)
  if (n < opt$min_n) why <- c(why, sprintf("n=%d < min_n=%d", n, opt$min_n))
  if (min(np, na_) < opt$min_class)
    why <- c(why, sprintf("present=%d absent=%d (min %d each)", np, na_, opt$min_class))
  if (binary) {
    n1 <- sum(d$pred == 1L); n0 <- n - n1
    if (min(n1, n0) < opt$min_class)
      why <- c(why, sprintf("predictor 1=%d 0=%d (min %d each)", n1, n0, opt$min_class))
  } else if (n > 0 && length(unique(d$pred)) < 2) {
    why <- c(why, "predictor has no variation")
  }
  if (length(why)) {
    r$note <- join_notes(notes, paste("not fitted:", paste(why, collapse = "; ")))
    r$status <- "skipped"
    return(r)
  }
  if (!binary) {
    d$pred <- zscore(d$pred)
  } else {
    tab <- table(factor(d$pred, levels = 0:1), factor(d$y, levels = 0:1))
    if (any(tab == 0))
      notes <- c(notes, sprintf(paste0("quasi-separation (pred0: absent=%d present=%d; ",
                                       "pred1: absent=%d present=%d): Firth-penalised ",
                                       "estimate, interpret with caution"),
                                tab[1, 1], tab[1, 2], tab[2, 1], tab[2, 2]))
  }
  res <- tryCatch(fit_phyloglm(d, c("pred", pp$terms), tree, opt$btol, opt$mple_starts),
                  error = function(e) e)
  if (inherits(res, "error")) {
    r$note <- join_notes(notes, paste("fit failed:", conditionMessage(res)))
    r$status <- "failed"
    say("WARNING: ", pname, ": phyloglm failed: ", conditionMessage(res))
    return(r)
  }
  r$estimate <- res$est; r$std_error <- res$se; r$p_value <- res$p
  r$ci_low <- res$est - stats::qnorm(0.975) * res$se
  r$ci_high <- res$est + stats::qnorm(0.975) * res$se
  r$alpha <- res$alpha; r$converged <- TRUE
  if (res$method != "logistic_MPLE") r$model <- sub("logistic_MPLE)$", paste0(res$method, " fallback)"), r$model)
  notes <- c(notes, res$detail, res$notes)
  if (length(res$warns)) notes <- c(notes, paste("phyloglm:", paste(res$warns, collapse = " | ")))
  r$note <- join_notes(notes)
  r$status <- "fitted"
  r
}

attach_sensitivity <- function(r, s, label) {
  r$sens_n <- s$n
  if (identical(s$status, "fitted")) {
    r$sens_estimate <- s$estimate; r$sens_std_error <- s$std_error; r$sens_p_value <- s$p_value
    r$sens_model <- if (nzchar(s$note)) paste0(label, " [", s$note, "]") else label
  } else {
    r$sens_model <- paste0(label, " [", s$note, "]")
  }
  r$sens_model <- clean(r$sens_model)
  r
}

main <- function() {
  option_list <- list(
    make_option("--master", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--viral", type = "character", default = NULL),
    make_option("--predictors", type = "character", default = "auto",
                help = "auto (all group_/niche_ columns) or a comma list of columns/globs"),
    make_option("--primary-tiers", type = "character", default = "A", dest = "primary_tiers"),
    make_option("--covariates", type = "character", default = "completeness,genome_size"),
    make_option("--fdr", type = "character", default = "BH"),
    make_option("--btol", type = "double", default = 10),
    make_option("--min-n", type = "integer", default = 10L, dest = "min_n"),
    make_option("--min-class", type = "integer", default = 3L, dest = "min_class"),
    make_option("--mple-starts", type = "integer", default = 5L, dest = "mple_starts",
                help = "start values of alpha for the MPLE fit (1-9)"),
    make_option("--seed", type = "integer", default = 12345L),
    make_option("--out", type = "character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  for (k in c("master", "tree", "out"))
    if (is.null(opt[[k]]) || !nzchar(opt[[k]])) stop("--", k, " is required")
  set.seed(opt$seed)
  tiers <- toupper(parse_list(opt$primary_tiers))
  if (!length(tiers) || !all(tiers %in% c("A", "B", "C")))
    stop("--primary-tiers must be a non-empty subset of A,B,C (got '", opt$primary_tiers, "')")
  fdr <- if (toupper(opt$fdr) %in% c("BH", "FDR")) "BH" else opt$fdr
  if (!fdr %in% stats::p.adjust.methods) stop("--fdr '", opt$fdr, "' is not a p.adjust method")

  # ---- tree ------------------------------------------------------------------------
  if (!file.exists(opt$tree)) stop("tree not found: ", opt$tree)
  txt <- trimws(paste(readLines(opt$tree, warn = FALSE), collapse = ""))
  if (!nzchar(txt) || grepl("^\\(\\s*\\)\\s*;?$", txt))
    return(header_only(opt$out, "empty species tree"))
  tree <- tryCatch(ape::read.tree(text = txt), error = function(e) NULL)
  if (is.null(tree) || !inherits(tree, "phylo")) stop("could not parse a single Newick tree from ", opt$tree)
  if (anyDuplicated(tree$tip.label)) stop("duplicated tip labels in ", opt$tree)
  if (is.null(tree$edge.length)) stop("species tree has no branch lengths")
  if (anyNA(tree$edge.length)) stop("species tree has missing (NA/NaN) branch lengths")
  if (any(tree$edge.length < 0)) stop("species tree has negative branch lengths")
  if (!ape::is.rooted(tree))
    stop("species tree is unrooted; phyloglm needs the rooted species tree (rooted upstream)")
  tree_note <- ""

  # ---- master + presence rule --------------------------------------------------------
  master <- read_tsv_chr(opt$master, "master table")
  if (!"genome_id" %in% names(master)) stop("master table lacks genome_id")
  if (nrow(master) == 0) return(header_only(opt$out, "empty master table"))
  master$genome_id <- trimws(master$genome_id)
  if (anyDuplicated(master$genome_id)) stop("duplicated genome_id in the master table")
  pres <- presence_rule(master, tiers)
  pass <- if ("pass_primary" %in% names(master)) flag(master$pass_primary) %in% TRUE else rep(TRUE, nrow(master))
  tip <- if ("tree_tip" %in% names(master)) trimws(master$tree_tip) else master$genome_id
  keep <- pass & nzchar(tip) & tip %in% tree$tip.label
  dup_tip <- keep & duplicated(ifelse(keep, tip, NA_character_), incomparables = NA)
  if (any(dup_tip)) {
    say("WARNING: ", sum(dup_tip), " duplicated tree_tip(s); first kept")
    keep <- keep & !dup_tip
  }
  n_excl <- sum(keep & is.na(pres$y))
  say(sprintf("%d genomes pass_primary on the tree; presence %s: %d present, %d absent, %d excluded (validated locus outside the tier set)",
              sum(keep), pres$basis, sum(keep & pres$y %in% 1L), sum(keep & pres$y %in% 0L), n_excl))
  base <- data.frame(genome_id = master$genome_id, tip = tip, y = pres$y,
                     stringsAsFactors = FALSE)
  covars <- parse_list(opt$covariates)
  cov_notes <- character(0)
  for (cv in covars) {
    if (cv %in% names(master)) base[[cv]] <- num(master[[cv]])
    else { cov_notes <- c(cov_notes, paste0("covariate ", cv, " not in master: not adjusted")); covars <- setdiff(covars, cv) }
  }
  if (length(cov_notes)) say("WARNING: ", paste(cov_notes, collapse = "; "))
  has_meth <- "methanogen" %in% names(master)
  if (has_meth) base$methanogen_bin <- as.integer(flag(master$methanogen))
  base <- base[keep & !is.na(pres$y), , drop = FALSE]
  common_note <- join_notes(pres$basis, if (n_excl) sprintf("%d excluded by the tier rule", n_excl),
                            tree_note, cov_notes)

  # ---- predictors -------------------------------------------------------------------
  rp <- resolve_predictors(opt$predictors, names(master))
  rows <- list()
  for (p in rp$missing) {
    say("WARNING: requested predictor '", p, "' is not a master column")
    rows[[length(rows) + 1]] <- blank_row(p, "", note = "requested predictor not found in the master table; not tested")
  }
  viral <- NULL
  if (!is.null(opt$viral) && nzchar(opt$viral)) {
    viral <- read_tsv_chr(opt$viral, "viral summary")
    need <- c("genome_id", CONT_PREDICTORS)
    if (!all(need %in% names(viral)))
      stop("viral summary lacks column(s): ", paste(setdiff(need, names(viral)), collapse = ", "))
    viral$genome_id <- trimws(viral$genome_id)
    if (anyDuplicated(viral$genome_id)) stop("duplicated genome_id in the viral summary")
  }
  if (!length(rp$found) && is.null(viral) && !length(rows))
    return(header_only(opt$out, "no group_/niche_ predictor and no --viral table"))

  fitted_or_failed <- 0L; failed <- 0L
  tally <- function(r) {
    if (r$status %in% c("fitted", "failed")) fitted_or_failed <<- fitted_or_failed + 1L
    if (identical(r$status, "failed")) failed <<- failed + 1L
  }
  for (p in rp$found) {
    b <- as_binary(master[[p]][match(base$genome_id, master$genome_id)])
    if (length(b$bad)) {
      rows[[length(rows) + 1]] <- blank_row(p, "1", note = clean(sprintf(
        "non-binary values (%s): not tested", paste(utils::head(b$bad, 3), collapse = ", "))))
      next
    }
    n_unknown <- sum(is.na(b$v))
    r <- run_model(base, p, "1", b$v, covars, tree, opt, TRUE,
                   extra_note = join_notes(common_note, if (n_unknown)
                     sprintf("%d genome(s) with unknown (NA) %s dropped", n_unknown, p)))
    tally(r)
    if (has_meth && identical(r$status, "fitted")) {
      bm <- base; bm$methanogen <- bm$methanogen_bin
      s <- run_model(bm, p, "1", b$v, c(covars, "methanogen"), tree, opt, TRUE)
      r <- attach_sensitivity(r, s, "+ methanogen status (partial effect)")
    }
    rows[[length(rows) + 1]] <- r
  }
  if (!is.null(viral)) {
    idx <- match(base$genome_id, viral$genome_id)
    n_miss <- sum(is.na(idx))
    vnote <- join_notes(common_note, if (n_miss) sprintf("%d genome(s) absent from the viral summary", n_miss))
    vcol <- function(col) if (col %in% names(viral)) num(viral[[col]][idx]) else rep(NA_real_, nrow(base))
    # n_proviruses: integrated proviruses, log1p then z-scored
    r <- run_model(base, "n_proviruses", "continuous", log1p(vcol("n_proviruses")), covars, tree,
                   opt, FALSE, pred_label = "z(log1p(n_proviruses))",
                   extra_note = join_notes(vnote, "log-odds per SD of log1p(n_proviruses)"))
    tally(r)
    if (identical(r$status, "fitted") && "n_viral_sequences" %in% names(viral)) {
      s <- run_model(base, "n_viral_sequences", "continuous", log1p(vcol("n_viral_sequences")),
                     covars, tree, opt, FALSE)
      r <- attach_sensitivity(r, s, "all geNomad viral sequences (proviruses + viral contigs), per SD of log1p")
    }
    rows[[length(rows) + 1]] <- r
    # frac_rloop_proviruses: genomes with >= 1 provirus (NA otherwise), z-scored
    fr <- vcol("frac_rloop_proviruses")
    r <- run_model(base, "frac_rloop_proviruses", "continuous", fr, covars, tree, opt, FALSE,
                   pred_label = "z(frac_rloop_proviruses)",
                   extra_note = join_notes(vnote, "genomes with >= 1 provirus; log-odds per SD of the R-loop-prone fraction"))
    tally(r)
    if (identical(r$status, "fitted") && all(c("mean_provirus_gc", "mean_provirus_kb") %in% names(viral))) {
      bs <- base; bs$provirus_gc <- vcol("mean_provirus_gc")
      bs$log_provirus_kb <- log(vcol("mean_provirus_kb"))
      s <- run_model(bs, "frac_rloop_proviruses", "continuous", fr,
                     c(covars, "provirus_gc", "log_provirus_kb"), tree, opt, FALSE)
      r <- attach_sensitivity(r, s, "+ z(mean provirus GC) + z(log mean provirus length)")
    }
    rows[[length(rows) + 1]] <- r
  }
  if (fitted_or_failed > 0 && failed == fitted_or_failed)
    stop("every attempted phyloglm fit failed (", failed, "); see the notes above")

  out <- do.call(rbind, lapply(rows, function(r) {
    r$status <- NULL
    r$n_excluded_tier <- n_excl
    as.data.frame(r[OUT_COLS], stringsAsFactors = FALSE)
  }))
  pv <- which(!is.na(out$p_value))
  out$p_adj <- NA_real_
  if (length(pv)) out$p_adj[pv] <- stats::p.adjust(out$p_value[pv], method = fdr)
  for (cc in c("level", "model", "note", "sens_model")) out[[cc]][is.na(out[[cc]])] <- ""
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  say(sprintf("%d term(s) (%d fitted, %d failed) -> %s", nrow(out),
              fitted_or_failed - failed, failed, opt$out))
  invisible(NULL)
}

main()
