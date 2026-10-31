#!/usr/bin/env Rscript
# =====================================================================
# Stage 4.6 -- co-occurrence / mutual exclusion of Prometheus with every other
# DefenseFinder system type, with a phylogeny-aware test (design 4.6). Reused
# unchanged by the bacterial branch (bact_cooccurrence).
#
# Per system (prevalence window: >= --min-count genomes with AND without it):
#   * phytools::fitPagel (Pagel 1994; ARD, dep.var = xy): correlated evolution of
#     Prometheus and the system along the tree. LR statistic, chi-square df = 4 ->
#     `statistic`, `p_value` (the test that enters the results registry). Guarded
#     against optimiser failure (see pagel_test). Never replaced by another test.
#   * Fisher's exact test (naive, tip-level): conditional-MLE odds ratio -> the
#     `direction` label (OR > 1 co-occurrence incl. Inf, OR < 1 exclusion incl. 0)
#     and the odds_ratio / or_ci_* / fisher_p columns.
#   * supplementary: phylolm::phyloglm(Prometheus ~ system + z(log genome size)),
#     multi-start MPLE -> phyloglm_logodds / phyloglm_se / phyloglm_p: a
#     phylogenetically controlled, genome-size-adjusted signed effect (design 4.6
#     "correcting for genome size"). NA with the reason in `note` when the master
#     has no genome_size or the model cannot be fitted for that system.
# p_adj: --fdr (BH by default) across the tested systems.
#
# Prometheus presence (SCHEMA Contract v2 item 6): from the master's
# prometheus_tier and --primary-tiers -- present if the best tier is in the set,
# absent if "", EXCLUDED (dropped, counted) otherwise. Only pass_primary genomes on
# the tree are analysed. Without a master table, the defense matrix's raw
# DefenseFinder Prometheus column is used (logged).
#
# Failure policy (Contract v2 item 4): an empty tree ("();"), an empty matrix, or
# no testable variation yields a header-only TSV (exit 0, reason logged).
# Malformed inputs, an unrooted tree, a missing package or a crashed fitPagel stop
# with a non-zero exit status.
#
# Output columns: system_a (Prometheus), system_b, n_a, n_b, n_both, method,
# statistic, p_value, p_adj, direction, then (additive) n, odds_ratio, or_ci_low,
# or_ci_high, fisher_p, logL_independent, logL_dependent, phyloglm_logodds,
# phyloglm_se, phyloglm_p, note.
# =====================================================================

Sys.setenv(LANGUAGE = "en")
suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(phytools)
  library(phylolm)
})

OUT_COLS <- c("system_a", "system_b", "n_a", "n_b", "n_both", "method", "statistic",
              "p_value", "p_adj", "direction", "n", "odds_ratio", "or_ci_low", "or_ci_high",
              "fisher_p", "logL_independent", "logL_dependent", "phyloglm_logodds",
              "phyloglm_se", "phyloglm_p", "note")
FOCAL <- "Prometheus"

log_msg <- function(...) message("[cooccurrence] ", ...)
clean_note <- function(x) gsub("[\t\r\n]+", " ", paste(x[nzchar(x)], collapse = "; "))
fmt <- function(x, d = 4) ifelse(is.na(x), "NA", as.character(signif(x, d)))

write_header_only <- function(path, reason) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(OUT_COLS, collapse = "\t"), path)
  log_msg("header-only output: ", reason)
}

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

read_species_tree <- function(path) {
  txt <- trimws(paste(readLines(path, warn = FALSE), collapse = ""))
  if (!nzchar(txt) || grepl("^\\(\\s*\\)\\s*;?$", txt)) return(NULL)
  tr <- ape::read.tree(text = txt)
  if (is.null(tr) || !inherits(tr, "phylo"))
    stop("could not parse the species tree as a single Newick tree: ", path)
  tr$tip.label <- sub("^'(.*)'$", "\\1", tr$tip.label)
  if (anyDuplicated(tr$tip.label)) stop("species tree has duplicated tip labels")
  if (is.null(tr$edge.length) || anyNA(tr$edge.length))
    stop("species tree lacks branch lengths (required by fitPagel and phyloglm)")
  if (any(tr$edge.length < 0)) stop("species tree has negative branch lengths")
  if (!ape::is.rooted(tr))
    stop("species tree is not rooted (basal polytomy); root it upstream (root_tree.py) -- ",
         "fitPagel and phyloglm are root-dependent")
  tr
}

prepare_tree <- function(tr) {
  tr$node.label <- NULL
  h <- max(ape::node.depth.edgelength(tr))
  if (!is.finite(h) || h <= 0) stop("species tree has no positive branch lengths")
  z <- tr$edge.length == 0
  if (any(z)) {
    tr$edge.length[z] <- 1e-6 * h
    log_msg(sum(z), " zero-length branch(es) set to 1e-6 x tree height")
  }
  tr
}

# Pagel (1994) test guarded against optimiser failure (identical to association_tests.R):
# the independent model is the Kronecker sum of two 2-state ARD chains, so its MLE is
# the sum of two separate 2-state fits; a dependent fit below the independent one is
# refitted from the independent solution.
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
       logL_independent = li, logL_dependent = ld, flags = flags)
}

# Multi-start phyloglm (identical to association_tests.R; no bootstrap here).
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
  f0 <- pglm_once(fml, d, tr, btol)
  tmax <- f0$mean.tip.height
  cands <- list(`aT=1` = polish(f0, fml, d, tr, btol, max_restarts))
  for (lab in c("aT=exp(-3)", "aT=exp(3)")) {
    a0 <- if (lab == "aT=exp(3)") exp(3) / tmax else exp(-3) / tmax
    f <- tryCatch(pglm_once(fml, d, tr, btol, sa = a0), error = function(e) NULL)
    if (!is.null(f)) cands[[lab]] <- polish(f, fml, d, tr, btol, max_restarts)
  }
  fi <- tryCatch(pglm_once(fml, d, tr, btol, method = "logistic_IG10"), error = function(e) NULL)
  if (!is.null(fi) && restartable(fi, btol)) {
    f <- tryCatch(pglm_once(fml, d, tr, btol, sb = as.numeric(fi$coefficients),
                            sa = clamp_alpha(fi$alpha, tmax)), error = function(e) NULL)
    if (!is.null(f)) cands[["IG10"]] <- polish(f, fml, d, tr, btol, max_restarts)
  }
  pens <- vapply(cands, function(cn) cn$fit$penlogLik, numeric(1))
  ib <- which.max(pens)
  list(fit = cands[[ib]]$fit, start = names(cands)[ib], gain = cands[[ib]]$gain)
}

direction_of <- function(or) {
  if (is.na(or)) "ambiguous" else if (or > 1) "co-occurrence" else if (or < 1) "exclusion" else "neutral"
}

# --------------------------------------------------------------------------- #
main <- function() {
  opt <- parse_args(OptionParser(option_list = list(
    make_option("--matrix", type = "character"),
    make_option("--master", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--fdr", type = "character", default = "BH"),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--primary-tiers", type = "character", default = "A", dest = "primary_tiers"),
    make_option("--min-count", type = "integer", default = 5, dest = "min_count",
                help = "minimum number of analysed genomes with AND without a system (default 5)"),
    make_option("--btol", type = "double", default = 10),
    make_option("--threads", type = "integer", default = 1),
    make_option("--out", type = "character"))))
  for (k in c("matrix", "tree", "out"))
    if (is.null(opt[[k]]) || !nzchar(opt[[k]])) stop("--", k, " is required")
  methods <- c(bh = "BH", fdr = "BH", by = "BY", bonferroni = "bonferroni", holm = "holm",
               hochberg = "hochberg", hommel = "hommel", none = "none")
  fdr <- methods[tolower(opt$fdr)]
  if (is.na(fdr)) stop("--fdr must be one of ", paste(names(methods), collapse = ", "))
  if (opt$min_count < 1) stop("--min-count must be >= 1")
  if (opt$threads < 1) stop("--threads must be >= 1")
  tiers <- parse_tiers(opt$primary_tiers, "--primary-tiers")

  tree0 <- read_species_tree(opt$tree)
  if (is.null(tree0))
    return(write_header_only(opt$out, "species tree is empty ('();')"))  # nonfatal: empty input
  tree0 <- prepare_tree(tree0)

  mat <- read_tsv_chr(opt$matrix)
  if (!"genome_id" %in% names(mat)) stop("defense matrix lacks a genome_id column")
  if (nrow(mat) == 0)
    return(write_header_only(opt$out, "defense matrix has zero genomes"))  # nonfatal: empty input
  mat$genome_id <- trimws(mat$genome_id)
  if (anyDuplicated(mat$genome_id)) stop("defense matrix has duplicated genome_id")
  prom_col <- names(mat)[tolower(names(mat)) == tolower(FOCAL)]
  sys_cols <- setdiff(names(mat), c("genome_id", prom_col))
  for (s in c(prom_col, sys_cols)) {
    v <- trimws(mat[[s]])
    if (!all(v %in% c("0", "1")))
      stop(sprintf("defense matrix column '%s' is not strictly 0/1 (e.g. '%s')", s, v[!v %in% c("0", "1")][1]))
  }

  gid <- mat$genome_id
  tip <- gid
  pass <- rep(TRUE, nrow(mat))
  gsize <- rep(NA_real_, nrow(mat))
  n_excl <- 0L
  if (!is.null(opt$master) && nzchar(opt$master)) {
    master <- read_tsv_chr(opt$master)
    if (!"genome_id" %in% names(master)) stop("master table lacks genome_id")
    master$genome_id <- trimws(master$genome_id)
    if (anyDuplicated(master$genome_id)) stop("master table has duplicated genome_id")
    idx <- match(gid, master$genome_id)
    if (anyNA(idx)) stop(sprintf("%d defense-matrix genome(s) absent from the master (e.g. %s)",
                                 sum(is.na(idx)), gid[is.na(idx)][1]))
    m <- master[idx, , drop = FALSE]
    if ("tree_tip" %in% names(m)) tip <- trimws(m$tree_tip)
    if ("pass_primary" %in% names(m)) pass <- parse_bool(m$pass_primary, "pass_primary", blank = FALSE)
    if ("genome_size" %in% names(m)) gsize <- parse_num(m$genome_size, "genome_size")
    if ("prometheus_tier" %in% names(m)) {
      prom_all <- tier_response(m$prometheus_tier, tiers)
      focal_src <- sprintf("master prometheus_tier, tier set %s", paste(tiers, collapse = "+"))
    } else {
      if (!length(prom_col)) stop("no prometheus_tier in master and no Prometheus column in the matrix")
      prom_all <- as.integer(mat[[prom_col]])
      focal_src <- "defense-matrix Prometheus column (raw DefenseFinder; master has no prometheus_tier)"
    }
  } else {
    if (!length(prom_col)) stop("no --master and no Prometheus column in the defense matrix")
    prom_all <- as.integer(mat[[prom_col]])
    focal_src <- "defense-matrix Prometheus column (raw DefenseFinder; no --master given)"
  }
  log_msg("Prometheus presence from ", focal_src)
  on_tree <- nzchar(tip) & tip %in% tree0$tip.label
  if (anyDuplicated(tip[on_tree])) stop("several genomes map to the same tree tip")
  cand <- pass & on_tree
  n_excl <- sum(cand & is.na(prom_all))
  use <- cand & !is.na(prom_all)
  log_msg(sprintf("matrix=%d genomes; pass_primary on tree=%d; EXCLUDED (best tier outside set)=%d; analysed=%d",
                  nrow(mat), sum(cand), n_excl, sum(use)))
  n <- sum(use)
  prom <- prom_all[use]
  n_a <- sum(prom)
  if (n_a < opt$min_count || (n - n_a) < opt$min_count)
    return(write_header_only(opt$out, sprintf(  # nonfatal: no testable variation in Prometheus
      "Prometheus present in %d of %d analysed genomes (< --min-count %d present or absent)", n_a, n, opt$min_count)))

  if (!length(sys_cols))
    return(write_header_only(opt$out, "defense matrix has no system type other than Prometheus"))  # nonfatal: nothing to test
  tips <- tip[use]
  tr <- ape::keep.tip(tree0, tips)
  if (!ape::is.rooted(tr)) stop("pruned species tree lost its root")
  ord <- match(tr$tip.label, tips)
  tips <- tips[ord]; prom <- prom[ord]
  S <- sapply(sys_cols, function(s) as.integer(mat[[s]][use][ord]))
  if (!is.matrix(S)) S <- matrix(S, nrow = length(tips), dimnames = list(NULL, sys_cols))
  rownames(S) <- tips
  nb_all <- colSums(S)
  testable <- sys_cols[nb_all >= opt$min_count & (n - nb_all) >= opt$min_count]
  log_msg(sprintf("%d system type(s) in the matrix; %d testable (>= %d genomes with and without); %d skipped",
                  length(sys_cols), length(testable), opt$min_count, length(sys_cols) - length(testable)))
  if (!length(testable))
    return(write_header_only(opt$out, "no system type inside the prevalence window"))  # nonfatal: nothing to test

  # genome-size-adjusted phyloglm subset (supplementary)
  gs_note <- ""
  gs_ok <- is.finite(gsize[use][ord]) & gsize[use][ord] > 0
  if (sum(gs_ok) >= 10) {
    d_gs <- data.frame(prom = prom[gs_ok], z_log_gs = as.numeric(scale(log(gsize[use][ord][gs_ok]))),
                       row.names = tips[gs_ok])
    tr_gs <- if (all(gs_ok)) tr else ape::keep.tip(tr, tips[gs_ok])
    d_gs <- d_gs[tr_gs$tip.label, , drop = FALSE]
    if (sd(d_gs$z_log_gs) == 0 || !is.finite(sd(d_gs$z_log_gs))) d_gs$z_log_gs <- NULL
  } else {
    d_gs <- NULL
    gs_note <- "phyloglm not run: genome_size unavailable for >= 10 analysed genomes"
  }

  one_system <- function(i) {
    s <- testable[i]
    set.seed(opt$seed + i)   # IG10 starts draw random numbers: fix per system for reproducibility
    sysv <- S[, s]
    nb <- sum(sysv); n_both <- sum(prom == 1 & sysv == 1)
    tab <- table(factor(prom, levels = 0:1), factor(sysv, levels = 0:1))
    ft <- stats::fisher.test(tab)
    or <- as.numeric(ft$estimate)
    pg <- pagel_test(tr, setNames(prom, tips), setNames(sysv, tips))
    pl_est <- pl_se <- pl_p <- NA_real_; pl_note <- gs_note
    if (!is.null(d_gs)) {
      dd <- d_gs; dd$system <- as.integer(sysv[match(rownames(dd), tips)])
      if (length(unique(dd$system)) < 2 || length(unique(dd$prom)) < 2) {
        pl_note <- "phyloglm not estimable: no variation in the genome-size subset"
      } else {
        fml <- if ("z_log_gs" %in% names(dd)) prom ~ system + z_log_gs else prom ~ system
        ms <- tryCatch(fit_phyloglm_ms(fml, dd, tr_gs, opt$btol), error = function(e) e)
        if (inherits(ms, "error")) {
          pl_note <- paste0("phyloglm not estimable: ", conditionMessage(ms))
          log_msg("WARNING: ", s, ": ", pl_note)
        } else {
          co <- summary(ms$fit)$coefficients
          pl_est <- co["system", "Estimate"]; pl_se <- co["system", "StdErr"]; pl_p <- co["system", "p.value"]
          pl_note <- sprintf("phyloglm(Prometheus~system+z(log genome size), multi-start MPLE, best start %s, alpha*T=%s%s)",
                             ms$start, fmt(ms$fit$alpha * ms$fit$mean.tip.height),
                             if (identical(as.integer(ms$fit$alphaWarn), 2L)) ", alpha at upper bound" else "")
        }
      }
    }
    data.frame(system_a = FOCAL, system_b = s, n_a = n_a, n_b = nb, n_both = n_both,
               method = "fitPagel", statistic = pg$lr, p_value = pg$p, p_adj = NA_real_,
               direction = sprintf("%s (OR=%s)", direction_of(or), fmt(or)),
               n = n, odds_ratio = or, or_ci_low = ft$conf.int[1], or_ci_high = ft$conf.int[2],
               fisher_p = ft$p.value, logL_independent = pg$logL_independent,
               logL_dependent = pg$logL_dependent, phyloglm_logodds = pl_est,
               phyloglm_se = pl_se, phyloglm_p = pl_p,
               note = clean_note(c(sprintf("tier set %s; n_excluded_tier=%d", paste(tiers, collapse = "+"), n_excl),
                                   pg$flags, pl_note)),
               stringsAsFactors = FALSE)
  }
  res <- if (opt$threads > 1) parallel::mclapply(seq_along(testable), one_system, mc.cores = opt$threads)
         else lapply(seq_along(testable), one_system)
  bad <- vapply(res, function(r) !is.data.frame(r), logical(1))
  if (any(bad)) {
    first <- res[[which(bad)[1]]]
    stop("per-system test failed for ", testable[which(bad)[1]], ": ",
         if (inherits(first, "try-error")) as.character(first) else "worker returned no result")
  }
  out <- do.call(rbind, res)
  out$p_adj <- stats::p.adjust(out$p_value, method = fdr)
  out <- out[order(out$system_b), OUT_COLS, drop = FALSE]
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  utils::write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  log_msg(sprintf("tested %d system(s) vs Prometheus (n=%d, n_a=%d; %s) -> %s",
                  nrow(out), n, n_a, fdr, opt$out))
  invisible(NULL)
}

# Run only under Rscript; source() exposes the helpers for unit tests.
if (sys.nframe() == 0L) main()
