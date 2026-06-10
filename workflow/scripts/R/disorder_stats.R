#!/usr/bin/env Rscript
# =====================================================================
# Disorder as a comparative trait, phylogeny-controlled (exploratory).
#
# Per-genome consensus disorder (mean over that genome's intact Prometheus loci;
# LOW-COMPLEXITY-EXCLUDED fraction preferred) is regressed on methanogen status and
# on each environment/niche, using phylogenetic regression (phylolm, Pagel's lambda)
# on the species tree -- the SAME non-independence control as the rest of the
# pipeline, so a lineage/habitat disorder signal can't masquerade as clade structure.
# Because disordered-residue content tracks growth temperature, any configured
# temperature-proxy columns are added as covariates to the niche models.
#
# Output (disorder_associations.tsv): predictor, level, model, estimate, std_error,
#   p_value, p_adj, n, note. Sequence-based prediction -> EXPLORATORY only. Trivial
# input -> header-only, quit(status=0).
# =====================================================================

suppressWarnings(suppressMessages({ library(optparse) }))
OUT_COLS <- c("predictor", "level", "model", "estimate", "std_error",
              "p_value", "p_adj", "n", "note")

write_header_only <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w"); writeLines(paste(OUT_COLS, collapse = "\t"), con); close(con)
}
as_logical_col <- function(x) tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "t")
erow <- function(predictor, level, model, estimate = NA_real_, std_error = NA_real_,
                 p_value = NA_real_, n = NA_integer_, note = "") {
  data.frame(predictor = predictor, level = level, model = model, estimate = estimate,
             std_error = std_error, p_value = p_value, p_adj = NA_real_, n = n,
             note = note, stringsAsFactors = FALSE)
}

main <- function() {
  opt <- parse_args(OptionParser(option_list = list(
    make_option("--master", type = "character"),
    make_option("--disorder", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--temp-cols", type = "character", default = "", dest = "temp_cols"),
    make_option("--fdr", type = "character", default = "BH"),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--out", type = "character"))))
  set.seed(opt$seed)

  tryCatch({
    stopifnot(!is.null(opt$master), !is.null(opt$disorder), !is.null(opt$tree), !is.null(opt$out))
    suppressWarnings(suppressMessages({ library(ape); library(readr); library(dplyr) }))

    tt <- tryCatch(paste(readLines(opt$tree, warn = FALSE), collapse = ""), error = function(e) "")
    if (!nzchar(trimws(tt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(tt))) {
      write_header_only(opt$out); message("disorder_stats: trivial tree."); return(invisible(NULL))
    }
    tree <- tryCatch(read.tree(opt$tree), error = function(e) NULL)
    if (!is.null(tree) && length(tree$tip.label) >= 3 && !ape::is.rooted(tree))
      tree <- tryCatch(phangorn::midpoint(tree), error = function(e) tree)
    if (is.null(tree) || length(tree$tip.label) < 5) {
      write_header_only(opt$out); message("disorder_stats: <5 tips."); return(invisible(NULL))
    }

    dis <- tryCatch(readr::read_tsv(opt$disorder, show_col_types = FALSE, progress = FALSE),
                    error = function(e) NULL)
    if (is.null(dis) || nrow(dis) == 0 || !"genome_id" %in% names(dis)) {
      write_header_only(opt$out); message("disorder_stats: empty disorder table."); return(invisible(NULL))
    }
    # response: LCR-excluded disorder fraction (preferred), else raw fraction; mean per genome
    resp_col <- if ("frac_disordered_nonlcr" %in% names(dis)) "frac_disordered_nonlcr" else "frac_disordered"
    dis$resp <- suppressWarnings(as.numeric(dis[[resp_col]]))
    pg <- dis %>% group_by(genome_id) %>% summarise(disorder = mean(resp, na.rm = TRUE), .groups = "drop")

    master <- readr::read_tsv(opt$master, show_col_types = FALSE, progress = FALSE)
    if (nrow(master) == 0 || !"genome_id" %in% names(master)) {
      write_header_only(opt$out); return(invisible(NULL))
    }
    master$genome_id <- trimws(as.character(master$genome_id))
    pg$genome_id <- trimws(as.character(pg$genome_id))
    dat <- dplyr::left_join(master, pg, by = "genome_id")
    dat <- dat[is.finite(dat$disorder), , drop = FALSE]
    if ("pass_primary" %in% names(dat)) dat <- dat[as_logical_col(dat$pass_primary), , drop = FALSE]
    dat$tip <- if ("tree_tip" %in% names(dat)) trimws(as.character(dat$tree_tip)) else dat$genome_id
    dat <- dat[dat$tip %in% tree$tip.label, , drop = FALSE]
    dat <- dat[!duplicated(dat$tip), , drop = FALSE]
    if (nrow(dat) < 5 || length(unique(round(dat$disorder, 6))) < 2) {
      write_header_only(opt$out); message("disorder_stats: <5 genomes or no disorder variation."); return(invisible(NULL))
    }
    if (!requireNamespace("phylolm", quietly = TRUE)) {
      write_header_only(opt$out); message("disorder_stats: phylolm unavailable."); return(invisible(NULL))
    }

    temp_cols <- trimws(strsplit(opt$temp_cols, ",")[[1]]); temp_cols <- temp_cols[nzchar(temp_cols)]
    temp_cols <- intersect(temp_cols, names(dat))

    fit_one <- function(pname, level_lab, covars = character(0)) {
      d <- dat
      d$pred <- if (pname == "methanogen") as.integer(as_logical_col(d$methanogen)) else
                as.integer(as_logical_col(d[[pname]]))
      d <- d[!is.na(d$pred), , drop = FALSE]
      if (length(unique(d$pred)) < 2) return(erow(pname, level_lab, "phylolm", n = nrow(d), note = "predictor constant; skipped"))
      n1 <- sum(d$pred == 1L)
      if (pname != "methanogen" && (n1 < 3 || (nrow(d) - n1) < 3))
        return(erow(pname, level_lab, "phylolm", n = nrow(d), note = sprintf("prevalence too low (n1=%d); skipped", n1)))
      use_cov <- covars[vapply(covars, function(c) length(unique(d[[c]])) > 1, logical(1))]
      tr <- ape::keep.tip(tree, d$tip); rownames(d) <- d$tip; d <- d[tr$tip.label, , drop = FALSE]
      rhs <- paste(c("pred", use_cov), collapse = " + ")
      fit <- tryCatch(phylolm::phylolm(stats::as.formula(paste("disorder ~", rhs)),
                                       data = d, phy = tr, model = "lambda"),
                      error = function(e) NULL)
      if (is.null(fit)) return(erow(pname, level_lab, "phylolm(lambda)", n = nrow(d), note = "fit failed"))
      co <- summary(fit)$coefficients
      est <- co["pred", "Estimate"]
      se <- if ("StdErr" %in% colnames(co)) co["pred", "StdErr"] else co["pred", 2]
      pc <- grep("p.value|Pr", colnames(co)); p <- if (length(pc)) co["pred", pc[1]] else NA_real_
      note <- sprintf("phylo regression (Pagel lambda); response=%s%s; EXPLORATORY",
                      resp_col, if (length(use_cov)) paste0("; adj=", paste(use_cov, collapse = "+")) else "")
      erow(pname, level_lab, "phylolm(disorder~pred, lambda)", est, se, p, nrow(d), note)
    }

    rows <- list()
    if ("methanogen" %in% names(dat)) rows[["methanogen"]] <- fit_one("methanogen", "1")
    eco_cols <- grep("^(group|niche)_", names(dat), value = TRUE)
    eco_cols <- setdiff(eco_cols, temp_cols)
    for (pn in c(sort(grep("^group_", eco_cols, value = TRUE)), sort(grep("^niche_", eco_cols, value = TRUE)))) {
      rows[[pn]] <- tryCatch(fit_one(pn, "1", covars = temp_cols),
                             error = function(e) erow(pn, "1", "phylolm", note = paste("failed:", conditionMessage(e))))
    }
    if (length(rows) == 0) { write_header_only(opt$out); return(invisible(NULL)) }
    out <- do.call(rbind, rows)
    pv <- which(!is.na(out$p_value))
    if (length(pv)) out$p_adj[pv] <- stats::p.adjust(out$p_value[pv],
                                                     method = if (toupper(opt$fdr) %in% c("BH", "FDR")) "BH" else opt$fdr)
    out <- out[, OUT_COLS, drop = FALSE]
    dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(out, opt$out)
    message(sprintf("disorder_stats: %d predictor-terms -> %s", nrow(out), opt$out))
    invisible(NULL)
  }, error = function(e) {
    message("disorder_stats: top-level error: ", conditionMessage(e))
    tryCatch(write_header_only(opt$out), error = function(e2) invisible(NULL))
  })
  quit(status = 0)
}
main()
