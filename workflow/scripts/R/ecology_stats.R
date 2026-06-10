#!/usr/bin/env Rscript
# =====================================================================
# Stage 6 -- ecology: does Prometheus presence track environment, host
# association, methanogenesis niche, or viral pressure BEYOND what the
# phylogeny explains (design Stage 6; confounder #11: biome and clade are
# deeply confounded in archaea)? Every predictor is tested in the SAME
# phylogenetic logistic framework as Stage 4 (phylolm::phyloglm), so a
# habitat signal cannot masquerade as a clade signal.
#
# Predictors = the SPECIES-LEVEL binary niche/group presence-absence
# columns carried in the master table (every column matching ^(group|niche)_,
# 0/1, autodetected — the legacy --predictors comma list is accepted but
# ignored) PLUS two continuous viral-pressure predictors appended from
# viral_summary: n_proviruses and frac_rloop_proviruses (the anti-viral /
# R-loop tie-back, Stage 6.3-6.4).
#
# Each binary group_/niche_ predictor is tested only when both states are
# present AND its prevalence (number of 1s) is >= 3; near-constant columns are
# skipped. Continuous viral predictors are tested when they have >= 5 non-NA
# values with variation. Every reported term's p-value is BH-adjusted across
# ALL terms.
#
# Output (SCHEMA ecology/ecology_associations.tsv):
#   predictor, level, model, estimate, std_error, p_value, p_adj, n, note
#   (one row per niche/group predictor [level = "1"] or viral predictor)
# Trivial input -> header-only TSV, quit(status=0).
# =====================================================================

suppressWarnings(suppressMessages({ library(optparse) }))

OUT_COLS <- c("predictor", "level", "model", "estimate", "std_error",
              "p_value", "p_adj", "n", "note")

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
  error = function(e) invisible(NULL)
)

write_header_only <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w"); writeLines(paste(OUT_COLS, collapse = "\t"), con); close(con)
}

as_logical_col <- function(x) {
  s <- tolower(trimws(as.character(x))); s %in% c("true", "1", "yes", "t")
}

erow <- function(predictor, level, model, estimate = NA_real_, std_error = NA_real_,
                 p_value = NA_real_, n = NA_integer_, note = "") {
  data.frame(predictor = predictor, level = level, model = model,
             estimate = estimate, std_error = std_error, p_value = p_value,
             p_adj = NA_real_, n = n, note = note, stringsAsFactors = FALSE)
}

main <- function() {
  option_list <- list(
    make_option("--master", type = "character"),
    make_option("--viral", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--predictors", type = "character",
                default = "biome_envo,oxygen,temp_class,salinity_class,host_associated"),
    make_option("--fdr", type = "character", default = "BH"),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--out", type = "character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  set.seed(opt$seed)

  tryCatch({
    stopifnot(!is.null(opt$master), !is.null(opt$tree), !is.null(opt$out))

    suppressWarnings(suppressMessages({
      library(ape); library(readr); library(dplyr)
    }))

    tree_txt <- tryCatch(paste(readLines(opt$tree, warn = FALSE), collapse = ""),
                         error = function(e) "")
    if (!nzchar(trimws(tree_txt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(tree_txt))) {
      write_header_only(opt$out)
      message("ecology_stats: trivial/empty tree; header-only."); return(invisible(NULL))
    }
    tree <- tryCatch(read.tree(opt$tree), error = function(e) NULL)
    # Safety net: phyloglm assumes a rooted tree. The pipeline roots
    # species_tree.treefile upstream; midpoint-root here only as a fallback.
    if (!is.null(tree) && length(tree$tip.label) >= 3 && !ape::is.rooted(tree))
      tree <- tryCatch(phangorn::midpoint(tree), error = function(e) tree)
    if (is.null(tree) || length(tree$tip.label) < 5) {
      write_header_only(opt$out)
      message("ecology_stats: <5 tree tips; header-only."); return(invisible(NULL))
    }

    master <- readr::read_tsv(opt$master, show_col_types = FALSE, progress = FALSE)
    if (nrow(master) == 0 || !"genome_id" %in% names(master)) {
      write_header_only(opt$out)
      message("ecology_stats: empty master; header-only."); return(invisible(NULL))
    }
    master$genome_id <- trimws(as.character(master$genome_id))

    # ---- merge viral summary by genome_id -------------------------- #
    if (!is.null(opt$viral) && file.exists(opt$viral)) {
      viral <- tryCatch(readr::read_tsv(opt$viral, show_col_types = FALSE,
                                        progress = FALSE),
                        error = function(e) NULL)
      if (!is.null(viral) && nrow(viral) > 0 && "genome_id" %in% names(viral)) {
        viral$genome_id <- trimws(as.character(viral$genome_id))
        keep_v <- intersect(c("genome_id", "n_proviruses", "frac_rloop_proviruses"),
                            names(viral))
        master <- dplyr::left_join(master, viral[, keep_v, drop = FALSE],
                                   by = "genome_id")
      }
    }
    if (!"n_proviruses" %in% names(master)) master$n_proviruses <- NA_real_
    if (!"frac_rloop_proviruses" %in% names(master)) master$frac_rloop_proviruses <- NA_real_

    # response
    master$present <- as.integer(as_logical_col(master$prometheus_present))

    # restrict to pass_primary genomes on the tree
    pass <- as_logical_col(master$pass_primary)
    tt <- if ("tree_tip" %in% names(master)) trimws(as.character(master$tree_tip)) else ""
    has_tip <- nzchar(tt)
    dat <- master[pass & has_tip, , drop = FALSE]
    dat$tip <- trimws(as.character(dat$tree_tip))
    dat <- dat[dat$tip %in% tree$tip.label, , drop = FALSE]
    dat <- dat[!duplicated(dat$tip), , drop = FALSE]
    if (nrow(dat) < 5 || length(unique(dat$present)) < 2) {
      write_header_only(opt$out)
      message("ecology_stats: <5 usable genomes or no presence variation; header-only.")
      return(invisible(NULL))
    }

    have_phylolm <- requireNamespace("phylolm", quietly = TRUE)
    if (!have_phylolm) {
      write_header_only(opt$out)
      message("ecology_stats: phylolm unavailable; header-only."); return(invisible(NULL))
    }

    # predictor list: AUTODETECT the binary species-level group_/niche_
    # presence-absence columns from the master table (the legacy --predictors
    # comma list is accepted but ignored), then append the continuous viral
    # predictors. group_ columns are reported before niche_, each block sorted.
    eco_cols <- grep("^(group|niche)_", names(dat), value = TRUE)
    bin_preds <- c(sort(grep("^group_", eco_cols, value = TRUE)),
                   sort(grep("^niche_", eco_cols, value = TRUE)))
    cont_preds <- c("n_proviruses", "frac_rloop_proviruses")
    all_preds <- unique(c(bin_preds, cont_preds))

    rows <- list()

    fit_one <- function(pname, is_continuous) {
      # Build a per-predictor data frame, drop NA, re-prune tree.
      if (!pname %in% names(dat)) {
        return(erow(pname, NA_character_, "phyloglm", note = "predictor not in master; skipped"))
      }
      d <- dat
      if (is_continuous) {
        d$pred <- suppressWarnings(as.numeric(as.character(d[[pname]])))
        # Provirus counts scale with genome completeness and size, so a raw
        # present~n_proviruses signal can be a genome-size/quality artefact. Adjust
        # for completeness and scaled genome size (as the headline methanogen model
        # does) before reporting any viral-pressure association.
        d$completeness_n <- suppressWarnings(as.numeric(as.character(d$completeness)))
        d$genome_size_n  <- suppressWarnings(as.numeric(as.character(d$genome_size)))
        d <- d[!is.na(d$pred) & !is.na(d$present), , drop = FALSE]
        if (nrow(d) < 5) return(erow(pname, NA_character_, "phyloglm", n = nrow(d),
                                     note = "fewer than 5 non-NA values; skipped"))
        if (length(unique(d$pred)) < 2) return(erow(pname, NA_character_, "phyloglm",
                                                    n = nrow(d), note = "no variation; skipped"))
        if (length(unique(d$present)) < 2) return(erow(pname, NA_character_, "phyloglm",
                                                       n = nrow(d), note = "no presence variation in subset; skipped"))
        covars <- character(0)
        if (all(!is.na(d$completeness_n)) && length(unique(d$completeness_n)) > 1)
          covars <- c(covars, "completeness_n")
        if (all(!is.na(d$genome_size_n)) && length(unique(d$genome_size_n)) > 1) {
          d$gs_scaled <- as.numeric(scale(d$genome_size_n)); covars <- c(covars, "gs_scaled")
        }
        if (nrow(d) < 5 || length(unique(d$present)) < 2 || length(unique(d$pred)) < 2)
          return(erow(pname, "continuous", "phyloglm", n = nrow(d),
                      note = "insufficient data after covariate adjustment; skipped"))
        tr <- ape::keep.tip(tree, d$tip)
        rownames(d) <- d$tip; d <- d[tr$tip.label, , drop = FALSE]
        rhs <- paste(c("pred", covars), collapse = " + ")
        fit <- phylolm::phyloglm(stats::as.formula(paste("present ~", rhs)), phy = tr,
                                 data = d, method = "logistic_MPLE", btol = 20)
        co <- summary(fit)$coefficients
        rn <- "pred"
        est <- co[rn, "Estimate"]
        se  <- if ("StdErr" %in% colnames(co)) co[rn, "StdErr"] else co[rn, 2]
        pc  <- grep("p.value|Pr", colnames(co))
        p   <- if (length(pc)) co[rn, pc[1]] else NA_real_
        model_str <- if (length(covars)) {
          sprintf("phyloglm(present~pred+%s, MPLE)", paste(covars, collapse = "+"))
        } else {
          "phyloglm(present~pred, MPLE)"
        }
        adj_note <- if (length(covars)) {
          paste("adjusted for", paste(covars, collapse = "+"))
        } else {
          "UNADJUSTED (completeness/genome_size constant or NA)"
        }
        return(erow(pname, "continuous", model_str, est, se, p, nrow(d),
                    paste("phylo logistic; continuous viral-pressure predictor;",
                          adj_note, "(log-odds per unit)")))
      } else {
        # binary 0/1 group_/niche_ presence-absence predictor.
        # robust coerce: 1/true/yes/t -> 1, everything else (incl NA/blank) -> 0
        d$pred <- as.integer(as_logical_col(d[[pname]]))
        d <- d[!is.na(d$pred) & !is.na(d$present), , drop = FALSE]
        if (nrow(d) < 5) return(erow(pname, NA_character_, "phyloglm", n = nrow(d),
                                     note = "fewer than 5 non-NA values; skipped"))
        prev <- sum(d$pred == 1L)
        # need both states present AND prevalence >= 3 (skip near-constant)
        if (length(unique(d$pred)) < 2)
          return(erow(pname, "1", "phyloglm", n = nrow(d),
                      note = "predictor constant in subset; skipped"))
        if (prev < 3 || (nrow(d) - prev) < 3)
          return(erow(pname, "1", "phyloglm", n = nrow(d),
                      note = sprintf("prevalence too low (n1=%d, n0=%d, min 3); skipped",
                                     prev, nrow(d) - prev)))
        if (length(unique(d$present)) < 2)
          return(erow(pname, "1", "phyloglm", n = nrow(d),
                      note = "no presence variation in subset; skipped"))
        tr <- ape::keep.tip(tree, d$tip)
        rownames(d) <- d$tip; d <- d[tr$tip.label, , drop = FALSE]
        fit <- phylolm::phyloglm(present ~ pred, phy = tr, data = d,
                                 method = "logistic_MPLE", btol = 20)
        co <- summary(fit)$coefficients
        rn <- "pred"
        est <- co[rn, "Estimate"]
        se  <- if ("StdErr" %in% colnames(co)) co[rn, "StdErr"] else co[rn, 2]
        pc  <- grep("p.value|Pr", colnames(co))
        p   <- if (length(pc)) co[rn, pc[1]] else NA_real_
        return(erow(pname, "1", "phyloglm(present~niche, MPLE)", est, se, p, nrow(d),
                    "phylo logistic; species in group/niche (1) vs not (log-odds)"))
      }
    }

    for (pn in all_preds) {
      is_cont <- pn %in% cont_preds
      r <- tryCatch(fit_one(pn, is_cont), error = function(e)
        erow(pn, "1", "phyloglm", note = paste("failed:", conditionMessage(e))))
      rows[[length(rows) + 1]] <- r
    }
    if (length(rows) == 0) {
      write_header_only(opt$out)
      message("ecology_stats: no usable predictor; header-only.")
      return(invisible(NULL))
    }

    out <- do.call(rbind, rows)
    # BH across all reported terms with a p-value.
    method <- if (toupper(opt$fdr) %in% c("BH", "FDR")) "BH" else opt$fdr
    pv <- which(!is.na(out$p_value))
    if (length(pv)) out$p_adj[pv] <- stats::p.adjust(out$p_value[pv], method = method)
    out <- out[, OUT_COLS, drop = FALSE]

    dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(out, opt$out)
    message(sprintf("ecology_stats: %d predictor-terms -> %s", nrow(out), opt$out))
    invisible(NULL)

  }, error = function(e) {
    message("ecology_stats: top-level error: ", conditionMessage(e))
    tryCatch(write_header_only(opt$out), error = function(e2) invisible(NULL))
    invisible(NULL)
  })

  quit(status = 0)
}

main()
