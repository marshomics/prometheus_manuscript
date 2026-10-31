#!/usr/bin/env Rscript
# =====================================================================
# Disorder as a comparative trait, phylogeny-controlled (EXPLORATORY).
#
# Response: per-genome mean of the locus-level consensus disorder fraction of the
# non-low-complexity residues (frac_disordered_nonlcr; frac_disordered when that column
# is absent) over the genome's scored (intact, non-fragmented) loci of tier set T.
#
# Tier rule (SCHEMA Contract v2 item 6), for every tier set T in --tier-sets
# (';'-separated, e.g. "A;A,B" = primary + A+B sensitivity): a genome enters iff its
# best validated tier (master prometheus_tier) is in T, and only its tier-T loci are
# averaged; a genome whose best tier is outside T is EXCLUDED and counted
# (n_excluded_tier). Population: pass_primary genomes with a species-tree tip.
#
# Models (phylolm, Pagel's lambda, ML; lambda in [1e-7, 1]):
#   phylogenetic_signal   disorder ~ 1 ; estimate = lambda, LRT against lambda = 0
#                         (0.5*chi2_0 + 0.5*chi2_1 boundary null; lambda = 0 keeps each
#                         tip's root-to-tip variance, exactly phylolm's transform)
#   methanogen            disorder ~ methanogen [+ covariates]; ANME genomes (is_anme)
#                         are excluded from this contrast (neither group)
#   group_* / niche_*     disorder ~ predictor [+ covariates], one model per column
# Covariates, in EVERY model: the configured temperature-proxy columns (--temp-cols;
# each must exist in the master table) -- disordered content tracks growth temperature
# and many methanogens are thermophiles -- and host_background when the per-locus table
# carries a host-proteome background. A binary predictor needs >= --min-group genomes
# in each level; otherwise the row is reported with a 'skipped' note.
#
# Output (disorder_associations.tsv): predictor, level, model, estimate, std_error,
#   p_value, p_adj, n, note  + additive: response, tier_set, n_level1, n_excluded_tier,
#   lambda, ci_low, ci_high (95% t interval, df = n - number of coefficients).
#   p_adj = BH (or --fdr) across every row with a p-value in this table.
# Genuinely empty input (no scored locus; a species tree with < 3 tips) -> header-only.
# A missing package, malformed input or crash -> non-zero exit (no silent fallback).
# =====================================================================

for (p in c("optparse", "ape", "phylolm")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("disorder_stats: required R package '%s' is not installed", p))
}
suppressPackageStartupMessages({ library(optparse); library(ape); library(phylolm) })

OUT_COLS <- c("predictor", "level", "model", "estimate", "std_error", "p_value", "p_adj",
              "n", "note")
EXTRA_COLS <- c("response", "tier_set", "n_level1", "n_excluded_tier", "lambda",
                "ci_low", "ci_high")
ALL_COLS <- c(OUT_COLS, EXTRA_COLS)
TRUE_SET <- c("true", "1", "yes", "t")
FALSE_SET <- c("false", "0", "no", "f")

write_table <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(df[, ALL_COLS, drop = FALSE], path, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "NA")
}
header_only <- function(path, why) {
  message("disorder_stats: ", why, "; header-only output.")
  write_table(as.data.frame(setNames(replicate(length(ALL_COLS), character(0),
                                               simplify = FALSE), ALL_COLS)), path)
}
read_tsv_chr <- function(path, what) {
  if (!file.exists(path)) stop(sprintf("disorder_stats: %s '%s' not found", what, path))
  if (file.info(path)$size == 0) stop(sprintf("disorder_stats: %s '%s' is empty (no header)", what, path))
  utils::read.delim(path, sep = "\t", colClasses = "character", na.strings = character(0),
                    check.names = FALSE, quote = "", comment.char = "")
}
flag01 <- function(x, col) {
  v <- tolower(trimws(as.character(x)))
  bad <- !(v %in% c(TRUE_SET, FALSE_SET, "", "na"))
  if (any(bad)) stop(sprintf("disorder_stats: column '%s' is not binary (e.g. '%s')",
                             col, x[which(bad)[1]]))
  out <- ifelse(v %in% TRUE_SET, 1L, ifelse(v %in% FALSE_SET, 0L, NA_integer_))
  as.integer(out)
}
row_df <- function(predictor, level, model, estimate = NA_real_, std_error = NA_real_,
                   p_value = NA_real_, n = NA_integer_, note = "", response = "",
                   tier_set = "", n_level1 = NA_integer_, n_excluded_tier = NA_integer_,
                   lambda = NA_real_, ci_low = NA_real_, ci_high = NA_real_) {
  data.frame(predictor = predictor, level = level, model = model, estimate = estimate,
             std_error = std_error, p_value = p_value, p_adj = NA_real_, n = n, note = note,
             response = response, tier_set = tier_set, n_level1 = n_level1,
             n_excluded_tier = n_excluded_tier, lambda = lambda, ci_low = ci_low,
             ci_high = ci_high, stringsAsFactors = FALSE)
}

# ML log-likelihood of y ~ 1 under Pagel's lambda = 0 as phylolm defines it: the
# covariance keeps each tip's root-to-tip distance on the diagonal and zeroes the rest.
loglik_lambda0 <- function(y, tr) {
  v <- ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))]
  if (any(v <= 0)) stop("disorder_stats: zero root-to-tip distance; lambda LRT undefined")
  w <- 1 / v; mu <- sum(w * y) / sum(w); n <- length(y)
  s2 <- sum(w * (y - mu)^2) / n
  -n / 2 * log(2 * pi * s2) - 0.5 * sum(log(v)) - n / 2
}

# phylolm with every warning captured (returned, not lost on stderr): lambda at its
# bound is reported in the row's note.
fit_phylolm <- function(fml, d, tr) {
  warns <- character(0)
  fit <- withCallingHandlers(
    tryCatch(phylolm::phylolm(fml, data = d, phy = tr, model = "lambda"), error = function(e) e),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  wn <- unique(trimws(gsub("\\s+", " ", warns)))
  list(fit = fit, note = if (length(wn)) paste0("; phylolm warning: ", paste(wn, collapse = " | ")) else "")
}

main <- function() {
  opt <- parse_args(OptionParser(option_list = list(
    make_option("--master", type = "character"),
    make_option("--disorder", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--temp-cols", type = "character", default = "", dest = "temp_cols"),
    make_option("--fdr", type = "character", default = "BH"),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--out", type = "character"),
    make_option("--tier-sets", type = "character", default = "A", dest = "tier_sets"),
    make_option("--min-group", type = "integer", default = 3L, dest = "min_group"),
    make_option("--min-n", type = "integer", default = 5L, dest = "min_n"))))
  for (k in c("master", "disorder", "tree", "out"))
    if (is.null(opt[[k]]) || !nzchar(opt[[k]])) stop(sprintf("disorder_stats: --%s is required", k))
  set.seed(opt$seed)
  fdr <- if (toupper(opt$fdr) %in% c("BH", "FDR")) "BH" else opt$fdr
  if (!fdr %in% stats::p.adjust.methods) stop(sprintf("disorder_stats: unknown --fdr '%s'", opt$fdr))
  tier_sets <- lapply(strsplit(strsplit(opt$tier_sets, ";", fixed = TRUE)[[1]], ","),
                      function(v) { v <- trimws(v); v[nzchar(v)] })
  tier_sets <- tier_sets[lengths(tier_sets) > 0]
  if (!length(tier_sets)) stop("disorder_stats: --tier-sets is empty")

  # ---- per-locus disorder ------------------------------------------------
  dis <- read_tsv_chr(opt$disorder, "disorder table")
  if (nrow(dis) == 0) return(header_only(opt$out, "no scored loci"))
  resp_col <- if ("frac_disordered_nonlcr" %in% names(dis)) "frac_disordered_nonlcr" else "frac_disordered"
  for (k in c("genome_id", resp_col))
    if (!k %in% names(dis)) stop(sprintf("disorder_stats: disorder table lacks '%s'", k))
  dis$resp <- suppressWarnings(as.numeric(dis[[resp_col]]))
  bad <- is.na(dis$resp) & !(trimws(dis[[resp_col]]) %in% c("", "NA"))
  if (any(bad)) stop(sprintf("disorder_stats: non-numeric %s value '%s'", resp_col, dis[[resp_col]][which(bad)[1]]))
  has_ct <- "confidence_tier" %in% names(dis)
  if (!has_ct) message("disorder_stats: WARNING: no confidence_tier column; the tier rule is applied at the genome level only.")
  use_bg <- "host_background" %in% names(dis) && any(nzchar(trimws(dis$host_background)))

  # ---- species tree -------------------------------------------------------
  if (!file.exists(opt$tree)) stop(sprintf("disorder_stats: tree '%s' not found", opt$tree))
  tt <- trimws(paste(readLines(opt$tree, warn = FALSE), collapse = ""))
  if (!nzchar(tt) || grepl("^\\(\\s*\\)\\s*;?$", tt)) return(header_only(opt$out, "empty species tree"))
  tree <- ape::read.tree(text = tt)
  if (inherits(tree, "multiPhylo")) stop("disorder_stats: the tree file holds several trees")
  if (is.null(tree) || ape::Ntip(tree) < 3) return(header_only(opt$out, "species tree has < 3 tips"))
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0))
    stop("disorder_stats: the species tree needs finite, non-negative branch lengths")
  root_note <- ""
  if (!ape::is.rooted(tree)) {
    if (!requireNamespace("phangorn", quietly = TRUE))
      stop("disorder_stats: species tree is unrooted and phangorn (midpoint) is unavailable")
    tree <- phangorn::midpoint(tree); root_note <- "; species tree was unrooted -> midpoint-rooted here"
    message("disorder_stats: WARNING: species tree unrooted; midpoint-rooted.")
  }

  # ---- master table -----------------------------------------------------
  master <- read_tsv_chr(opt$master, "master table")
  for (k in c("genome_id", "tree_tip", "prometheus_tier", "methanogen"))
    if (!k %in% names(master)) stop(sprintf("disorder_stats: master table lacks '%s'", k))
  temp_cols <- trimws(strsplit(opt$temp_cols, ",", fixed = TRUE)[[1]]); temp_cols <- temp_cols[nzchar(temp_cols)]
  miss_t <- setdiff(temp_cols, names(master))
  if (length(miss_t)) stop(sprintf("disorder_stats: temperature column(s) not in the master table: %s",
                                   paste(miss_t, collapse = ", ")))
  master$genome_id <- trimws(master$genome_id); master$tip <- trimws(master$tree_tip)
  pass <- if ("pass_primary" %in% names(master)) tolower(trimws(master$pass_primary)) %in% TRUE_SET else rep(TRUE, nrow(master))
  pop <- master[pass & nzchar(master$tip) & master$tip %in% tree$tip.label, , drop = FALSE]
  if (any(duplicated(pop$tip))) stop("disorder_stats: several genomes share one tree tip")
  pop$best_tier <- trimws(pop$prometheus_tier)
  pop$methanogen01 <- flag01(pop$methanogen, "methanogen")
  pop$anme01 <- if ("is_anme" %in% names(pop)) flag01(pop$is_anme, "is_anme") else 0L
  for (tc in temp_cols) pop[[tc]] <- flag01(pop[[tc]], tc)
  eco_cols <- c(sort(grep("^group_", names(pop), value = TRUE)), sort(grep("^niche_", names(pop), value = TRUE)))
  eco_cols <- setdiff(eco_cols, temp_cols)
  for (ec in eco_cols) pop[[ec]] <- flag01(pop[[ec]], ec)
  if (use_bg) {
    bgv <- suppressWarnings(as.numeric(dis$host_background))
    bg_by_g <- tapply(bgv, dis$genome_id, function(v) { v <- unique(v[is.finite(v)]); if (length(v) == 1) v else NA_real_ })
  }

  rows <- list()
  for (tset in tier_sets) {
    ts_lab <- paste(tset, collapse = ",")
    loci <- if (has_ct) dis[trimws(dis$confidence_tier) %in% tset, , drop = FALSE] else dis
    pg <- tapply(loci$resp, loci$genome_id, function(v) if (any(is.finite(v))) mean(v[is.finite(v)]) else NA_real_)
    n_excl <- sum(nzchar(pop$best_tier) & !(pop$best_tier %in% tset))
    dat <- pop[pop$best_tier %in% tset, , drop = FALSE]
    dat$disorder <- as.numeric(pg[dat$genome_id])
    n_unscored <- sum(!is.finite(dat$disorder))
    dat <- dat[is.finite(dat$disorder), , drop = FALSE]
    covs <- temp_cols
    if (use_bg) {
      dat$host_background <- as.numeric(bg_by_g[dat$genome_id])
      if (any(!is.finite(dat$host_background)))
        stop("disorder_stats: host_background missing or inconsistent for an analysed genome")
      covs <- c(covs, "host_background")
    }
    base_note <- sprintf("EXPLORATORY (sequence-predicted); tiers %s; %d genome(s) outside the tier set excluded; %d tier-%s genome(s) without an intact scored locus%s",
                         ts_lab, n_excl, n_unscored, ts_lab, root_note)
    message(sprintf("disorder_stats: tiers %s -> %d genomes (%d excluded by tier, %d unscored).",
                    ts_lab, nrow(dat), n_excl, n_unscored))
    common <- list(response = resp_col, tier_set = ts_lab, n_excluded_tier = n_excl)

    # -- phylogenetic signal (lambda) of the disorder trait itself --
    sig_model <- sprintf("phylolm(%s ~ 1, lambda) LRT vs lambda=0 [tiers %s]", resp_col, ts_lab)
    if (nrow(dat) < opt$min_n || length(unique(round(dat$disorder, 8))) < 2) {
      rows[[length(rows) + 1]] <- do.call(row_df, c(list("phylogenetic_signal", "lambda", sig_model,
        n = nrow(dat), note = paste0("skipped: n < ", opt$min_n, " or no variation; ", base_note)), common))
    } else {
      tr <- ape::keep.tip(tree, dat$tip); d <- dat[match(tr$tip.label, dat$tip), , drop = FALSE]
      rownames(d) <- d$tip
      fw <- fit_phylolm(disorder ~ 1, d, tr)
      if (inherits(fw$fit, "error")) stop("disorder_stats: lambda fit failed: ", conditionMessage(fw$fit))
      f1 <- fw$fit
      ll0 <- loglik_lambda0(d$disorder, tr)
      stat <- max(0, 2 * (f1$logLik - ll0))
      p <- if (stat > 0) 0.5 * stats::pchisq(stat, df = 1, lower.tail = FALSE) else 1
      rows[[length(rows) + 1]] <- do.call(row_df, c(list("phylogenetic_signal", "lambda", sig_model,
        estimate = unname(f1$optpar), p_value = p, n = nrow(d), lambda = unname(f1$optpar),
        note = sprintf("Pagel lambda (ML); LRT stat %.4g, boundary chi2 mixture%s; %s", stat, fw$note, base_note)), common))
    }

    fit_binary <- function(pname, d0, extra_note = "") {
      model <- sprintf("phylolm(%s ~ %s%s, lambda) [tiers %s]", resp_col, pname,
                       if (length(covs)) paste0(" + ", paste(covs, collapse = " + ")) else "", ts_lab)
      n_na <- sum(!stats::complete.cases(d0[, c("pred", covs), drop = FALSE]))
      d0 <- d0[stats::complete.cases(d0[, c("pred", covs), drop = FALSE]), , drop = FALSE]
      if (n_na) extra_note <- sprintf("%s; %d genome(s) with a missing predictor/covariate dropped", extra_note, n_na)
      n1 <- sum(d0$pred == 1L); n0 <- sum(d0$pred == 0L)
      mk <- function(...) do.call(row_df, c(list(pname, "1", model, ...), common))
      if (nrow(d0) < opt$min_n || n1 < opt$min_group || n0 < opt$min_group)
        return(mk(n = nrow(d0), n_level1 = n1, note = sprintf(
          "skipped: n=%d, level 1 n=%d, level 0 n=%d (need n >= %d and >= %d per level)%s; %s",
          nrow(d0), n1, n0, opt$min_n, opt$min_group, extra_note, base_note)))
      use <- covs[vapply(covs, function(cv) length(unique(d0[[cv]])) > 1 &&
                           !isTRUE(all(d0[[cv]] == d0$pred)) && !isTRUE(all(d0[[cv]] == 1L - d0$pred)),
                         logical(1))]
      dropped <- setdiff(covs, use)
      tr <- ape::keep.tip(tree, d0$tip); d <- d0[match(tr$tip.label, d0$tip), , drop = FALSE]
      rownames(d) <- d$tip
      fml <- stats::as.formula(paste("disorder ~", paste(c("pred", use), collapse = " + ")))
      fw <- fit_phylolm(fml, d, tr); fit <- fw$fit
      cov_note <- paste0(if (length(use)) paste0("; adjusted for ", paste(use, collapse = "+")) else "",
                         if (length(dropped)) paste0("; covariate(s) constant or collinear here, dropped: ",
                                                     paste(dropped, collapse = "+")) else "",
                         fw$note)
      if (inherits(fit, "error"))
        return(mk(n = nrow(d), n_level1 = n1, note = paste0("fit failed: ", conditionMessage(fit),
                                                          extra_note, cov_note, "; ", base_note)))
      co <- summary(fit)$coefficients
      est <- co["pred", "Estimate"]; se <- co["pred", "StdErr"]; p <- co["pred", "p.value"]
      dfree <- fit$n - fit$d; q <- stats::qt(0.975, dfree)
      mk(estimate = est, std_error = se, p_value = p, n = nrow(d), n_level1 = n1,
         lambda = unname(fit$optpar), ci_low = est - q * se, ci_high = est + q * se,
         note = paste0("phylo regression (Pagel lambda, ML)", extra_note, cov_note, "; ", base_note))
    }

    # -- methanogen (ANME excluded from both groups) --
    dm <- dat[dat$anme01 != 1L, , drop = FALSE]; dm$pred <- dm$methanogen01
    rows[[length(rows) + 1]] <- fit_binary("methanogen", dm,
      sprintf("; %d ANME genome(s) excluded", sum(dat$anme01 == 1L)))
    # -- each ecology column --
    for (ec in eco_cols) {
      de <- dat; de$pred <- de[[ec]]
      rows[[length(rows) + 1]] <- fit_binary(ec, de)
    }
  }

  out <- do.call(rbind, rows)
  pv <- which(is.finite(out$p_value))
  if (length(pv)) out$p_adj[pv] <- stats::p.adjust(out$p_value[pv], method = fdr)
  write_table(out, opt$out)
  message(sprintf("disorder_stats: %d rows (%d with p-values) -> %s", nrow(out), length(pv), opt$out))
}

main()
