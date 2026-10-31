#!/usr/bin/env Rscript
# =====================================================================
# Stage 5.5 -- ancestral gain/loss of Prometheus on the rooted archaeal species tree
# by stochastic character mapping, and the methanogen-gain enrichment test.
#
# Presence (Contract v2 item 6): for a tier set T and a genome's best validated tier
# prometheus_tier (A > B > C, "" = no validated locus): PRESENT if the tier is in T,
# ABSENT if "", EXCLUDED (tip dropped, counted) otherwise. --tiers gives the primary set;
# --tiers-sensitivity (optional) repeats everything for a second set (rows tagged by the
# additive tier_set column).
#
# Model (phytools; verified against 2.1.1 locally and the 2.3-0 source pinned in
# rstats.yaml -- same make.simmap/fitMk/sim.Mk interfaces):
#   * Mk model choice: fitMk ER vs ARD with the same root prior; ARD only when it lowers
#     AIC by >= 2 (parsimony rule), otherwise ER.
#   * make.simmap(tree, x, model, nsim, Q = --q-method ("empirical" = ML Q, the default;
#     "mcmc" = Q sampled from its posterior), pi = --root-prior ("fitzjohn" default:
#     FitzJohn et al. 2009 root weighting; "equal"; "estimated" = stationary)).
#     make.simmap internally reorders the tree "cladewise", so every per-edge quantity is
#     computed on the returned simmap trees, never on the input edge order.
#   * n_gains / n_losses: mean 0->1 / 1->0 transitions per map (phytools::countSimmap).
#   * n_independent_acquisitions: mean number of DISTINCT gain events (plus ancestral
#     presence at the root) that seed >= 1 extant carrier with no intervening loss.
#   * gains_on_methanogen_branches: mean gains on edges whose descendant tips are all
#     methanogens. The trait is master$methanogen exactly as in the association tests:
#     TRUE only for mcrA-defined methanogens; ANME, mcr 'other' (alkane oxidisers) and
#     ambiguous/unplaced McrA hosts are FALSE (non-methanogen), which is conservative for
#     the enrichment (a gain in such a clade never counts as a methanogen gain).
#   * the presence derived from prometheus_tier is cross-checked against the master's
#     prometheus_status_primary / _sensitivity columns when present (mismatch = error).
#   * structured null: the methanogen trait is simulated on the same tree under its own
#     fitted Mk model (ER/ARD by AIC, FitzJohn root prior), keeping only simulations whose
#     methanogen-tip count is within --prevalence-tol of the observed one (rejection; the
#     unconditional structured null is used, and labelled, if too few are accepted); the
#     FIXED gain maps are rescored on each simulated trait. One-sided p = (k + 1)/(n + 1),
#     k = #null >= observed. Free tip shuffling is never used.
# Inputs that are genuinely too small (trivial tree, < 5 matched tips, presence or trait
# invariant) give NA rows with an explanatory note; any other error exits non-zero.
# Output columns: metric, value, note, tier_set (gains_losses.tsv); kind, value, tier_set
# (null file). Plot: densityMap of presence (primary tier set).
# =====================================================================

suppressWarnings(suppressMessages({
  library(optparse)
  library(ape)
  library(phytools)
}))

.script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  ""
}
# plotting_theme.R sits next to this script (Rscript --file=...), or in the working
# directory when the file is source()d with chdir = TRUE (unit tests).
.theme <- Filter(file.exists, c(file.path(dirname(.script_path()), "plotting_theme.R"),
                                file.path(getwd(), "plotting_theme.R")))
if (!length(.theme)) stop("plotting_theme.R not found next to ancestral_states.R")
source(.theme[1])

OUT_COLS <- c("metric", "value", "note", "tier_set")
METRICS <- c("n_gains", "n_losses", "rate_gain", "rate_loss", "n_independent_acquisitions",
             "gains_on_methanogen_branches", "expected_gains_methanogen", "enrichment_p",
             "n_tips_analysed", "n_tips_excluded_tier", "n_tips_present")

as_flag <- function(x) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(s))
  out[s %in% c("true", "1", "yes", "t")] <- TRUE
  out[s %in% c("false", "0", "no", "f")] <- FALSE
  out
}

parse_tiers <- function(s) {
  t <- toupper(trimws(unlist(strsplit(as.character(s), "[,+;/ ]+"))))
  t <- t[nzchar(t)]
  if (!length(t) || !all(t %in% c("A", "B", "C")))
    stop("invalid tier set '", s, "' (use letters A,B,C)")
  sort(unique(t))
}

# Contract v2 item 6: 1 present, 0 absent, NA excluded.
tier_presence <- function(tier, tiers) {
  tier <- toupper(trimws(as.character(tier)))
  tier[is.na(tier)] <- ""
  bad <- !(tier %in% c("", "A", "B", "C"))
  if (any(bad)) stop("unexpected prometheus_tier value(s): ",
                     paste(unique(tier[bad]), collapse = ","))
  ifelse(tier == "", 0L, ifelse(tier %in% tiers, 1L, NA_integer_))
}

na_rows <- function(note, tier_lab, extra = list()) {
  d <- data.frame(metric = METRICS, value = NA_real_, note = note, tier_set = tier_lab,
                  stringsAsFactors = FALSE)
  for (k in names(extra)) d$value[d$metric == k] <- extra[[k]]
  d
}

# ---- tree helpers (on the simmap edge order) -------------------------------
descendant_all <- function(tree, tipvals) {
  # tipvals: logical matrix ntip x nsim (rows in tree$tip.label order). Returns a logical
  # matrix (Ntip + Nnode) x nsim: TRUE when every descendant tip is TRUE.
  n <- Ntip(tree); nn <- n + tree$Nnode
  out <- matrix(TRUE, nn, ncol(tipvals))
  out[seq_len(n), ] <- tipvals
  po <- ape::reorder.phylo(tree, "postorder", index.only = TRUE)
  for (e in po) {
    p <- tree$edge[e, 1]; c <- tree$edge[e, 2]
    out[p, ] <- out[p, ] & out[c, ]
  }
  out
}

edge_gain_means <- function(smaps) {
  # mean number of 0->1 transitions on each edge over the maps
  ne <- nrow(smaps[[1]]$edge)
  g <- numeric(ne)
  for (s in smaps) {
    for (i in seq_len(ne)) {
      st <- names(s$maps[[i]])
      if (length(st) > 1) g[i] <- g[i] + sum(st[-length(st)] == "0" & st[-1] == "1")
    }
  }
  g / length(smaps)
}

acquisitions_one <- function(s) {
  # distinct gain events (or root presence) seeding >= 1 extant carrier, no later loss
  tr <- s; n <- Ntip(tr)
  ne <- nrow(tr$edge)
  pe <- integer(n + tr$Nnode); pe[tr$edge[, 2]] <- seq_len(ne)   # parent edge of a node
  seed <- character(ne)
  for (e in ape::reorder.phylo(tr, "cladewise", index.only = TRUE)) {
    st <- names(tr$maps[[e]]); k <- length(st)
    if (st[k] != "1") { seed[e] <- NA; next }
    z <- which(st == "0")
    if (length(z)) {
      seed[e] <- paste0(e, ":", max(z) + 1)
    } else {
      par <- tr$edge[e, 1]
      seed[e] <- if (pe[par] == 0) "root" else seed[pe[par]]
    }
  }
  tip_e <- pe[seq_len(n)]
  st_tip <- vapply(tip_e, function(e) { v <- names(tr$maps[[e]]); v[length(v)] }, "")
  length(unique(seed[tip_e[st_tip == "1"]]))
}

# Two-state Mk simulation, vectorised over replicates (closed-form P(t)).
sim_binary <- function(tree, q01, q10, root_p1, nsim) {
  n <- Ntip(tree); nn <- n + tree$Nnode
  st <- matrix(FALSE, nn, nsim)
  st[n + 1, ] <- runif(nsim) < root_p1
  r <- q01 + q10
  for (e in ape::reorder.phylo(tree, "cladewise", index.only = TRUE)) {
    p <- tree$edge[e, 1]; c <- tree$edge[e, 2]; t <- tree$edge.length[e]
    decay <- if (r > 0) 1 - exp(-r * t) else 0
    p01 <- if (r > 0) q01 / r * decay else 0
    p10 <- if (r > 0) q10 / r * decay else 0
    u <- runif(nsim)
    st[c, ] <- ifelse(st[p, ], u >= p10, u < p01)
  }
  st[seq_len(n), , drop = FALSE]
}

fit_mk_aic <- function(tree, x, root_prior) {
  fits <- list(ER = fitMk(tree, x, model = "ER", pi = root_prior),
               ARD = fitMk(tree, x, model = "ARD", pi = root_prior))
  aic <- vapply(fits, AIC, numeric(1))
  best <- if (aic[["ARD"]] <= aic[["ER"]] - 2) "ARD" else "ER"
  list(fit = fits[[best]], model = best, aic = aic)
}

q_rates <- function(fit) {
  Q <- unclass(as.Qmatrix(fit))
  c(q01 = Q["0", "1"], q10 = Q["1", "0"])
}

root_p1 <- function(fit, q) {
  pi <- fit$pi
  if (is.numeric(pi) && all(c("0", "1") %in% names(pi)) && is.finite(pi[["1"]]))
    return(as.numeric(pi[["1"]]))
  if (sum(q) > 0) return(as.numeric(q[["q01"]] / sum(q)))   # stationary distribution
  0.5
}

analyse <- function(tree, master, tiers, opt, status_col = NULL) {
  tier_lab <- paste(tiers, collapse = ",")
  tip_ids <- ifelse(nzchar(master$tree_tip), master$tree_tip, master$genome_id)
  pres <- setNames(tier_presence(master$prometheus_tier, tiers), tip_ids)
  # The master table's own status column (scaffolds, same Contract v2 item 6 rule) must agree.
  if (!is.null(status_col) && status_col %in% names(master)) {
    st <- tolower(trimws(master[[status_col]]))
    ref <- ifelse(st == "present", 1L, ifelse(st == "absent", 0L, ifelse(st == "excluded", NA_integer_, -1L)))
    if (any(ref == -1L, na.rm = TRUE)) stop("unexpected values in master$", status_col)
    bad <- xor(is.na(ref), is.na(pres)) | (!is.na(ref) & !is.na(pres) & ref != pres)
    if (any(bad)) stop(sprintf("presence from prometheus_tier with tiers {%s} disagrees with master$%s for %d genome(s) (e.g. %s): the tier sets used by scaffolds and evolution differ",
                               tier_lab, status_col, sum(bad), paste(head(tip_ids[bad], 3), collapse = ",")))
  }
  meth <- setNames(as_flag(master$methanogen), tip_ids)
  on_tree <- intersect(tree$tip.label, names(pres))
  n_excl <- sum(is.na(pres[on_tree]))
  keep <- on_tree[!is.na(pres[on_tree])]
  base <- list(n_tips_excluded_tier = n_excl, n_tips_analysed = length(keep))
  message(sprintf("[%s] %d species-tree tips matched; %d excluded by the tier rule; %d analysed",
                  tier_lab, length(on_tree), n_excl, length(keep)))
  if (length(keep) < 5)
    return(list(rows = na_rows(sprintf("not computed: %d analysable tips (< 5)", length(keep)),
                               tier_lab, base)))
  tr <- keep.tip(tree, keep)
  x <- pres[tr$tip.label]
  base$n_tips_present <- sum(x)
  if (length(unique(x)) < 2)
    return(list(rows = na_rows("not computed: presence invariant after the tier rule",
                               tier_lab, base)))
  xf <- setNames(factor(x, levels = c(0, 1)), tr$tip.label)

  # ---- presence: model choice + stochastic maps ------------------------------
  mk <- fit_mk_aic(tr, xf, opt$root_prior)
  smaps <- make.simmap(tr, xf, model = mk$model, nsim = opt$nsim, Q = opt$q_method,
                       pi = opt$root_prior, message = FALSE)
  if (inherits(smaps, "simmap")) smaps <- structure(list(smaps), class = c("multiSimmap", "multiPhylo"))
  ref <- smaps[[1]]
  for (s in smaps) if (!identical(s$edge, ref$edge)) stop("simmap edge matrices differ between maps")
  cs <- countSimmap(smaps, states = c("0", "1"), message = FALSE)
  if (is.list(cs) && !is.null(cs$Tr)) cs <- cs$Tr      # single-map return shape
  cs <- as.matrix(cs)
  gains <- if ("0,1" %in% colnames(cs)) cs[, "0,1"] else rep(0, nrow(cs))
  losses <- if ("1,0" %in% colnames(cs)) cs[, "1,0"] else rep(0, nrow(cs))
  acq <- vapply(smaps, acquisitions_one, numeric(1))
  Qs <- lapply(smaps, function(s) if (!is.null(s$Q)) s$Q else attr(s, "Q"))  # phytools stores $Q
  if (any(vapply(Qs, is.null, logical(1)))) stop("simmap object carries no fitted Q")
  rg <- mean(vapply(Qs, function(Q) Q["0", "1"], numeric(1)))
  rl <- mean(vapply(Qs, function(Q) Q["1", "0"], numeric(1)))
  qnote <- sprintf("%s Mk (AIC ER=%.2f ARD=%.2f; ARD needs dAIC>=2), Q=%s, root prior=%s, %d maps",
                   mk$model, mk$aic[["ER"]], mk$aic[["ARD"]], opt$q_method, opt$root_prior, opt$nsim)

  # ---- methanogen enrichment -------------------------------------------------
  obs <- expv <- pval <- NA_real_
  enote <- "not computed"; null_vals <- numeric(0)
  mv <- meth[ref$tip.label]
  if (anyNA(mv)) {
    enote <- sprintf("not computed: methanogen status missing for %d tips", sum(is.na(mv)))
  } else if (length(unique(mv)) < 2) {
    enote <- "not computed: methanogen status invariant on the analysed tips"
  } else {
    gbar <- edge_gain_means(smaps)
    child <- ref$edge[, 2]
    obs_all <- descendant_all(ref, matrix(mv, ncol = 1))
    obs <- sum(gbar * obs_all[child, 1])
    mf <- setNames(factor(as.integer(mv), levels = c(0, 1)), ref$tip.label)
    mmk <- fit_mk_aic(ref, mf, opt$root_prior)
    q <- q_rates(mmk$fit); p1 <- root_p1(mmk$fit, q)
    n_obs <- sum(mv); tol <- max(1, round(opt$prevalence_tol * n_obs))
    acc <- list(); n_acc <- 0; n_tried <- 0
    while (n_acc < opt$nnull && n_tried < opt$max_null_sims) {
      b <- min(opt$batch, opt$max_null_sims - n_tried)
      sims <- sim_binary(ref, q[["q01"]], q[["q10"]], p1, b)
      n_tried <- n_tried + b
      ok <- abs(colSums(sims) - n_obs) <= tol & colSums(sims) > 0 & colSums(sims) < Ntip(ref)
      if (any(ok)) { acc[[length(acc) + 1]] <- sims[, ok, drop = FALSE]; n_acc <- n_acc + sum(ok) }
    }
    if (n_acc >= opt$min_accept) {
      S <- do.call(cbind, acc)[, seq_len(min(n_acc, opt$nnull)), drop = FALSE]
      cond <- sprintf("prevalence-matched (|n_meth - %d| <= %d; accepted %d of %d)", n_obs, tol,
                      n_acc, n_tried)
    } else {
      S <- sim_binary(ref, q[["q01"]], q[["q10"]], p1, opt$nnull)
      cond <- sprintf("UNCONDITIONAL (only %d of %d simulations matched the prevalence)", n_acc, n_tried)
    }
    D <- descendant_all(ref, S)
    null_vals <- as.numeric(crossprod(D[child, , drop = FALSE], gbar))
    expv <- mean(null_vals)
    pval <- (sum(null_vals >= obs - 1e-12) + 1) / (length(null_vals) + 1)
    enote <- sprintf("one-sided P(null >= observed) = (k+1)/(n+1), n=%d; methanogen trait simulated under its fitted %s Mk (q01=%.4g q10=%.4g, P(root=1)=%.3f), %s; gain maps held fixed",
                     length(null_vals), mmk$model, q[["q01"]], q[["q10"]], p1, cond)
  }
  rows <- data.frame(
    metric = METRICS,
    value = c(mean(gains), mean(losses), rg, rl, mean(acq), obs, expv, pval,
              base$n_tips_analysed, base$n_tips_excluded_tier, base$n_tips_present),
    note = c(paste("mean 0->1 transitions per map;", qnote),
             paste("mean 1->0 transitions per map;", qnote),
             "Q[0,1] gain rate (mean over maps)", "Q[1,0] loss rate (mean over maps)",
             sprintf("mean distinct gain events (or root presence) seeding >=1 extant carrier with no later loss (sd %.3f)", sd(acq)),
             "mean gains on edges whose descendant tips are all methanogens",
             "mean of the structured null", enote,
             "species-tree tips analysed", sprintf("tips excluded: best tier outside {%s}", tier_lab),
             sprintf("tips present (tier in {%s})", tier_lab)),
    tier_set = tier_lab, stringsAsFactors = FALSE)
  list(rows = rows, smaps = smaps, null = null_vals, obs = obs, n_acq = mean(acq))
}

main <- function() {
  option_list <- list(
    make_option("--tree", type = "character"),
    make_option("--master", type = "character"),
    make_option("--nsim", type = "integer", default = 1000),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--tiers", type = "character", default = "A"),
    make_option("--tiers-sensitivity", type = "character", dest = "tiers_sens", default = ""),
    make_option("--q-method", type = "character", dest = "q_method", default = "empirical"),
    make_option("--root-prior", type = "character", dest = "root_prior", default = "fitzjohn"),
    make_option("--nnull", type = "integer", default = 1000),
    make_option("--prevalence-tol", type = "double", dest = "prevalence_tol", default = 0.1),
    make_option("--max-null-sims", type = "integer", dest = "max_null_sims", default = 200000),
    make_option("--min-accept", type = "integer", dest = "min_accept", default = 200),
    make_option("--batch", type = "integer", default = 5000),
    make_option("--out", type = "character"),
    make_option("--out-null", type = "character", dest = "out_null", default = NULL),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  for (k in c("tree", "master", "out", "plot_prefix"))
    if (is.null(opt[[k]])) stop("--", gsub("_", "-", k), " is required")
  if (!opt$q_method %in% c("empirical", "mcmc")) stop("--q-method must be empirical or mcmc")
  if (!opt$root_prior %in% c("fitzjohn", "equal", "estimated")) stop("bad --root-prior")
  set.seed(opt$seed)
  dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
  tier_sets <- list(parse_tiers(opt$tiers))
  if (nzchar(opt$tiers_sens)) {
    ts <- parse_tiers(opt$tiers_sens)
    if (!identical(ts, tier_sets[[1]])) tier_sets[[2]] <- ts
  }

  master <- read.delim(opt$master, sep = "\t", colClasses = "character", na.strings = character(0),
                       check.names = FALSE)
  need <- c("genome_id", "prometheus_tier", "methanogen")
  miss <- setdiff(need, names(master))
  if (length(miss)) stop("master table lacks column(s): ", paste(miss, collapse = ", "))
  if (!"tree_tip" %in% names(master)) master$tree_tip <- ""
  master$genome_id <- trimws(master$genome_id); master$tree_tip <- trimws(master$tree_tip)

  tree_txt <- trimws(paste(readLines(opt$tree, warn = FALSE), collapse = ""))
  trivial <- !nzchar(tree_txt) || grepl("^\\(\\s*\\)\\s*;?$", tree_txt)
  res <- list()
  if (trivial) {
    for (t in tier_sets) res[[length(res) + 1]] <- list(rows = na_rows("not computed: trivial species tree",
                                                                        paste(t, collapse = ",")))
  } else {
    tree <- read.tree(opt$tree)
    if (!is.rooted(tree)) stop("species tree is not rooted; ancestral reconstruction needs the rooted ",
                               "tree written by scaffolds/root_tree.py")
    if (anyDuplicated(tree$tip.label)) stop("duplicated tip labels in the species tree")
    if (!is.binary(tree)) tree <- multi2di(tree, random = FALSE)
    if (is.null(tree$edge.length)) stop("species tree has no branch lengths")
    tree$edge.length[!is.finite(tree$edge.length) | tree$edge.length <= 0] <- 1e-6
    scol <- c("prometheus_status_primary", "prometheus_status_sensitivity")
    for (i in seq_along(tier_sets)) res[[i]] <- analyse(tree, master, tier_sets[[i]], opt, scol[i])
  }
  out <- do.call(rbind, lapply(res, `[[`, "rows"))
  write.table(out, opt$out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  if (!is.null(opt$out_null)) {
    nd <- do.call(rbind, lapply(seq_along(res), function(i) {
      r <- res[[i]]; lab <- paste(tier_sets[[i]], collapse = ",")
      d <- data.frame(kind = character(0), value = numeric(0), tier_set = character(0))
      if (length(r$null)) d <- rbind(d, data.frame(kind = "sim", value = r$null, tier_set = lab))
      if (!is.null(r$obs) && is.finite(r$obs)) d <- rbind(d, data.frame(kind = "observed", value = r$obs, tier_set = lab))
      d
    }))
    if (is.null(nd) || !nrow(nd)) nd <- data.frame(kind = character(0), value = numeric(0), tier_set = character(0))
    write.table(nd, opt$out_null, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  prim <- res[[1]]
  if (is.null(prim$smaps)) {
    placeholder_plot(opt$plot_prefix, paste("Ancestral states:", prim$rows$note[1]))
  } else {
    sm <- prim$smaps; nt <- Ntip(sm[[1]])
    fig <- function() {
      dm <- densityMap(sm, states = c("0", "1"), plot = FALSE, res = 200)
      dm <- setMap(dm, c("#D9D9D9", "#0072B2"))
      plot(dm, lwd = 2, ftype = if (nt <= 60) "i" else "off", fsize = 0.5, outline = FALSE,
           legend = FALSE)
      graphics::title(main = sprintf("P(Prometheus present), tiers {%s}; independent acquisitions = %.2f",
                                     paste(tier_sets[[1]], collapse = ","), prim$n_acq), cex.main = 0.7)
    }
    save_plot(fig, opt$plot_prefix, width = 5.5, height = max(3, min(0.12 * nt, 11)))
  }
  message("ancestral_states: wrote ", opt$out)
}

if (sys.nframe() == 0L) main()   # run as a script; source() only defines the functions
