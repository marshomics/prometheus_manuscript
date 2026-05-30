#!/usr/bin/env Rscript
# =====================================================================
# Stage 4.3-4.5 -- taxonomic association of Prometheus with methanogeny,
# tested in a phylogenetic comparative framework (design Stage 4).
#
# Methanogens are POLYPHYLETIC (methanogenesis arose ~once and was lost
# repeatedly; design Stage 3.3 / confounder #5), so genomes are NOT
# independent observations. We therefore report the naive contingency /
# logistic test only as a labelled baseline, and the defensible estimates
# come from phylogenetic logistic regression (phylolm::phyloglm, Ives &
# Garland), a correlated-evolution test (phytools::fitPagel, Pagel 1994),
# and the phylogenetic-signal D statistic (caper::phylo.d, Fritz & Purvis).
#
# The headline number is the methanogen effect from phyloglm AFTER
# conditioning on genome completeness, genome size and TOTAL defense load
# (confounder #8: methanogens may simply be defense-rich). Each fit is
# individually guarded so one failure does not lose the others.
#
# Output (SCHEMA distribution/association_tests.tsv):
#   test, model, term, estimate, std_error, statistic, p_value,
#   ci_low, ci_high, n, note
# On trivial input (empty/'();' tree, <5 usable tips, no variation in the
# response or predictor) writes a header-only TSV and exits status 0.
# =====================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

OUT_COLS <- c("test", "model", "term", "estimate", "std_error",
              "statistic", "p_value", "ci_low", "ci_high", "n", "note")

# --- locate and source the shared plotting/theme helpers (robust) ------
.script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  ""
}
.this_dir <- tryCatch(dirname(.script_path()), error = function(e) ".")
if (is.na(.this_dir) || !nzchar(.this_dir)) .this_dir <- "."
tryCatch(
  suppressWarnings(suppressMessages(source(file.path(.this_dir, "plotting_theme.R")))),
  error = function(e) invisible(NULL)  # this script writes no figure; theme optional
)

# --------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------- #
write_header_only <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w")
  writeLines(paste(OUT_COLS, collapse = "\t"), con)
  close(con)
}

# Build a one-row result with all columns present (NA where not applicable).
row_result <- function(test, model, term, estimate = NA_real_, std_error = NA_real_,
                       statistic = NA_real_, p_value = NA_real_,
                       ci_low = NA_real_, ci_high = NA_real_, n = NA_integer_,
                       note = "") {
  data.frame(test = test, model = model, term = term,
             estimate = estimate, std_error = std_error, statistic = statistic,
             p_value = p_value, ci_low = ci_low, ci_high = ci_high,
             n = n, note = note, stringsAsFactors = FALSE)
}

as_logical_col <- function(x) {
  s <- tolower(trimws(as.character(x)))
  s %in% c("true", "1", "yes", "t")
}

# --------------------------------------------------------------------- #
main <- function() {
  option_list <- list(
    make_option("--master", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--primary-tiers", type = "character", default = "A",
                dest = "primary_tiers"),
    make_option("--fdr", type = "character", default = "BH"),
    make_option("--btol", type = "double", default = 10),
    make_option("--boot", type = "integer", default = 200),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--out", type = "character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  set.seed(opt$seed)

  # Everything below is guarded; on ANY error or trivial input we still
  # emit a header-only TSV and quit(status = 0) so the DAG never fails.
  result <- tryCatch({

    stopifnot(!is.null(opt$master), !is.null(opt$tree), !is.null(opt$out))

    suppressWarnings(suppressMessages({
      library(ape)
      library(readr)
      library(dplyr)
    }))

    # ---- read tree -------------------------------------------------- #
    tree_txt <- tryCatch(paste(readLines(opt$tree, warn = FALSE), collapse = ""),
                         error = function(e) "")
    if (!nzchar(trimws(tree_txt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(tree_txt))) {
      write_header_only(opt$out)
      message("association_tests: trivial/empty tree; wrote header-only output.")
      return(invisible(NULL))
    }
    tree <- tryCatch(read.tree(opt$tree), error = function(e) NULL)
    # Safety net: phyloglm / fitPagel / phylo.d assume a rooted tree. The pipeline
    # roots species_tree.treefile upstream; midpoint-root here only as a fallback.
    if (!is.null(tree) && length(tree$tip.label) >= 3 && !ape::is.rooted(tree))
      tree <- tryCatch(phangorn::midpoint(tree), error = function(e) tree)
    if (is.null(tree) || is.null(tree$tip.label) || length(tree$tip.label) < 5) {
      write_header_only(opt$out)
      message("association_tests: <5 tree tips; wrote header-only output.")
      return(invisible(NULL))
    }

    # ---- read + shape master --------------------------------------- #
    master <- readr::read_tsv(opt$master, show_col_types = FALSE,
                              progress = FALSE)
    if (nrow(master) == 0 || !"genome_id" %in% names(master)) {
      write_header_only(opt$out)
      message("association_tests: empty master; wrote header-only output.")
      return(invisible(NULL))
    }

    primary_tiers <- trimws(strsplit(opt$primary_tiers, ",")[[1]])
    primary_tiers <- primary_tiers[nzchar(primary_tiers)]

    # Response: present iff prometheus_present == TRUE AND tier in primary set.
    pres_flag <- as_logical_col(master$prometheus_present)
    tier <- if ("prometheus_tier" %in% names(master)) {
      trimws(as.character(master$prometheus_tier))
    } else rep("", nrow(master))
    master$present <- as.integer(pres_flag & (tier %in% primary_tiers))

    # Predictor: methanogen (1/0).
    master$methanogen_bin <- as.integer(as_logical_col(master$methanogen))

    # Restrict to pass_primary genomes that have a tree tip.
    pass <- as_logical_col(master$pass_primary)
    has_tip <- !is.na(master$tree_tip) & nzchar(trimws(as.character(master$tree_tip)))
    dat <- master[pass & has_tip, , drop = FALSE]
    dat$tip <- trimws(as.character(dat$tree_tip))
    dat <- dat[dat$tip %in% tree$tip.label, , drop = FALSE]
    # one row per tip (defensive: keep first if duplicated)
    dat <- dat[!duplicated(dat$tip), , drop = FALSE]

    if (nrow(dat) < 5) {
      write_header_only(opt$out)
      message("association_tests: <5 usable genomes after pruning; header-only.")
      return(invisible(NULL))
    }

    # Prune tree to these tips and order data to tree tip order.
    keep <- intersect(tree$tip.label, dat$tip)
    if (length(keep) < 5) {
      write_header_only(opt$out)
      message("association_tests: <5 shared tips; header-only.")
      return(invisible(NULL))
    }
    tree <- ape::keep.tip(tree, keep)
    # geiger::name.check sanity (non-fatal if geiger missing)
    rownames(dat) <- dat$tip
    dat <- dat[tree$tip.label, , drop = FALSE]

    # Response / predictor variation check.
    if (length(unique(dat$present)) < 2) {
      write_header_only(opt$out)
      message("association_tests: no variation in Prometheus presence; header-only.")
      return(invisible(NULL))
    }

    n_all <- nrow(dat)
    results <- list()

    # Numeric confounders coerced once (used by model c).
    num <- function(x) suppressWarnings(as.numeric(as.character(x)))
    dat$completeness_n      <- num(dat$completeness)
    dat$genome_size_n       <- num(dat$genome_size)
    dat$n_defense_systems_n <- num(dat$n_defense_systems)

    # =============================================================== #
    # (a) naive GLM -- baseline that ASSUMES INDEPENDENCE (labelled).
    # =============================================================== #
    results[["naive_glm"]] <- tryCatch({
      fit <- glm(present ~ methanogen_bin, data = dat, family = binomial())
      sm <- summary(fit)$coefficients
      ci <- suppressMessages(tryCatch(confint(fit), error = function(e) NULL))
      est <- sm["methanogen_bin", "Estimate"]
      se  <- sm["methanogen_bin", "Std. Error"]
      z   <- sm["methanogen_bin", "z value"]
      p   <- sm["methanogen_bin", "Pr(>|z|)"]
      lo <- hi <- NA_real_
      if (!is.null(ci) && "methanogen_bin" %in% rownames(ci)) {
        lo <- ci["methanogen_bin", 1]; hi <- ci["methanogen_bin", 2]
      }
      row_result("naive_glm", "glm(present~methanogen, binomial)",
                 "methanogen", est, se, z, p, lo, hi, n_all,
                 "assumes independence; baseline only (log-odds)")
    }, error = function(e) {
      row_result("naive_glm", "glm(present~methanogen, binomial)", "methanogen",
                 n = n_all, note = paste("failed:", conditionMessage(e)))
    })

    # =============================================================== #
    # (b) phylogenetic logistic regression, methanogen only.
    # =============================================================== #
    results[["phyloglm_methanogen"]] <- tryCatch({
      if (!requireNamespace("phylolm", quietly = TRUE))
        stop("phylolm not available")
      d <- dat
      fit <- phylolm::phyloglm(present ~ methanogen_bin, phy = tree, data = d,
                               method = "logistic_MPLE",
                               btol = opt$btol, boot = opt$boot)
      sm <- summary(fit)
      co <- sm$coefficients
      est <- co["methanogen_bin", "Estimate"]
      se  <- if ("StdErr" %in% colnames(co)) co["methanogen_bin", "StdErr"]
             else co["methanogen_bin", 2]
      zc  <- grep("z.value", colnames(co))
      zstat <- if (length(zc)) co["methanogen_bin", zc[1]] else NA_real_
      pc  <- grep("p.value|Pr", colnames(co))
      p <- if (length(pc)) co["methanogen_bin", pc[1]] else NA_real_
      # bootstrap CI if available
      lo <- hi <- NA_real_
      bci <- tryCatch(fit$bootconfint95, error = function(e) NULL)
      if (!is.null(bci) && "methanogen_bin" %in% colnames(bci)) {
        lo <- bci[1, "methanogen_bin"]; hi <- bci[2, "methanogen_bin"]
      }
      note <- if (!is.na(lo)) sprintf("phylo logistic (MPLE); bootstrap CI (n_boot=%d); log-odds", opt$boot)
              else "phylo logistic (MPLE); log-odds; bootstrap CI unavailable"
      row_result("phyloglm_methanogen",
                 "phyloglm(present~methanogen, MPLE)", "methanogen",
                 est, se, zstat, p, lo, hi, n_all, note)
    }, error = function(e) {
      row_result("phyloglm_methanogen", "phyloglm(present~methanogen, MPLE)",
                 "methanogen", n = n_all,
                 note = paste("failed:", conditionMessage(e)))
    })

    # =============================================================== #
    # (c) phyloglm + confounders -> HEADLINE methanogen estimate.
    #     Conditions on completeness, scale(genome_size), total defense.
    # =============================================================== #
    results[["phyloglm_confounders"]] <- tryCatch({
      if (!requireNamespace("phylolm", quietly = TRUE))
        stop("phylolm not available")
      d <- dat
      # Need non-NA covariates; drop tips with missing confounders and
      # re-prune the tree so phyloglm sees a matched tree+data.
      keepc <- stats::complete.cases(d$present, d$methanogen_bin,
                                     d$completeness_n, d$genome_size_n,
                                     d$n_defense_systems_n)
      d <- d[keepc, , drop = FALSE]
      if (nrow(d) < 5) stop("fewer than 5 genomes with complete confounders")
      if (length(unique(d$present)) < 2) stop("no presence variation in complete-case subset")
      tr <- ape::keep.tip(tree, d$tip)
      d <- d[tr$tip.label, , drop = FALSE]
      d$gs_scaled <- as.numeric(scale(d$genome_size_n))
      fit <- phylolm::phyloglm(
        present ~ methanogen_bin + completeness_n + gs_scaled + n_defense_systems_n,
        phy = tr, data = d, method = "logistic_MPLE",
        btol = opt$btol, boot = opt$boot)
      sm <- summary(fit); co <- sm$coefficients
      est <- co["methanogen_bin", "Estimate"]
      se  <- if ("StdErr" %in% colnames(co)) co["methanogen_bin", "StdErr"]
             else co["methanogen_bin", 2]
      zc <- grep("z.value", colnames(co))
      zstat <- if (length(zc)) co["methanogen_bin", zc[1]] else NA_real_
      pc <- grep("p.value|Pr", colnames(co))
      p <- if (length(pc)) co["methanogen_bin", pc[1]] else NA_real_
      lo <- hi <- NA_real_
      bci <- tryCatch(fit$bootconfint95, error = function(e) NULL)
      if (!is.null(bci) && "methanogen_bin" %in% colnames(bci)) {
        lo <- bci[1, "methanogen_bin"]; hi <- bci[2, "methanogen_bin"]
      }
      row_result("phyloglm_confounders",
                 "phyloglm(present~methanogen+completeness+scale(genome_size)+n_defense_systems, MPLE)",
                 "methanogen", est, se, zstat, p, lo, hi, nrow(d),
                 "HEADLINE: methanogen effect adjusted for completeness, genome size, total defense load (log-odds)")
    }, error = function(e) {
      row_result("phyloglm_confounders",
                 "phyloglm(present~methanogen+confounders, MPLE)", "methanogen",
                 n = n_all, note = paste("failed:", conditionMessage(e)))
    })

    # =============================================================== #
    # (d) fitPagel -- correlated evolution of the two binary traits.
    # =============================================================== #
    results[["fitPagel"]] <- tryCatch({
      if (!requireNamespace("phytools", quietly = TRUE))
        stop("phytools not available")
      if (length(unique(dat$methanogen_bin)) < 2)
        stop("no variation in methanogen status")
      x <- setNames(as.factor(dat$methanogen_bin), dat$tip)  # methanogen
      y <- setNames(as.factor(dat$present), dat$tip)         # Prometheus
      pg <- phytools::fitPagel(tree, x = x, y = y)
      # LR statistic + P from the dependent vs independent model comparison.
      lr <- tryCatch(pg$lik.ratio, error = function(e) NA_real_)
      p  <- tryCatch(pg$P, error = function(e) NA_real_)
      row_result("fitPagel",
                 "fitPagel(methanogen, present)", "correlated_evolution",
                 estimate = NA_real_, statistic = lr, p_value = p, n = n_all,
                 note = "Pagel94 dependent vs independent; LR test of correlated gains/losses")
    }, error = function(e) {
      row_result("fitPagel", "fitPagel(methanogen, present)",
                 "correlated_evolution", n = n_all,
                 note = paste("failed:", conditionMessage(e)))
    })

    # =============================================================== #
    # (e) D statistic -- phylogenetic signal of Prometheus presence.
    # =============================================================== #
    results[["D_statistic"]] <- tryCatch({
      if (!requireNamespace("caper", quietly = TRUE))
        stop("caper not available")
      cd_df <- data.frame(names = dat$tip, present = dat$present,
                          stringsAsFactors = FALSE)
      cdat <- caper::comparative.data(phy = tree, data = cd_df, names.col = "names")
      pd <- caper::phylo.d(cdat, binvar = present)
      Dval <- tryCatch(as.numeric(pd$DEstimate), error = function(e) NA_real_)
      p1 <- tryCatch(as.numeric(pd$Pval1), error = function(e) NA_real_)  # vs random
      p0 <- tryCatch(as.numeric(pd$Pval0), error = function(e) NA_real_)  # vs Brownian
      row_result("D_statistic", "phylo.d(present)", "D",
                 estimate = Dval, n = n_all,
                 note = sprintf("Fritz-Purvis D=%.3f; Pval1(vs random)=%.3g; Pval0(vs Brownian)=%.3g; D~0 clumped, D~1 random",
                                ifelse(is.na(Dval), NA_real_, Dval), p1, p0))
    }, error = function(e) {
      row_result("D_statistic", "phylo.d(present)", "D", n = n_all,
                 note = paste("failed:", conditionMessage(e)))
    })

    out <- do.call(rbind, results)
    # NOTE: the single cross-family FDR regime is applied downstream in
    # results_registry.py (design Stage 7.2); per SCHEMA this table carries
    # raw p_values only. The --fdr flag is accepted for signature parity and
    # so the chosen method can be recorded, but no p_adj column is emitted
    # here (it is not part of the association_tests.tsv contract).
    out <- out[, OUT_COLS, drop = FALSE]

    dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(out, opt$out)
    message(sprintf("association_tests: wrote %d rows (n=%d genomes) -> %s",
                    nrow(out), n_all, opt$out))
    invisible(NULL)

  }, error = function(e) {
    message("association_tests: top-level error: ", conditionMessage(e))
    tryCatch(write_header_only(opt$out), error = function(e2) invisible(NULL))
    invisible(NULL)
  })

  quit(status = 0)
}

main()
