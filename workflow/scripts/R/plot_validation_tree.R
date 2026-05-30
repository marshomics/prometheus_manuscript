#!/usr/bin/env Rscript
# =====================================================================
# Stage 1/7 figure -- the Prometheus validation tree (design Stage 1.4 /
# 7.1). A ggtree of the Stage-1 guide tree that places every candidate
# protein against the curated references. Tips are coloured by category:
#   * reference tips from --labels: Prometheus_ref / Paralog_ref
#   * candidate tips from --catalog: by confidence_tier A / B / C / Rejected
#   * anything else: grey ('other')
# This is the figure that shows candidates falling WITH Prometheus vs
# inside the housekeeping DNA2/NurA-HerA paralog clade -- the visual
# justification for the catalog's accept/reject calls.
#
# --labels is a HEADERLESS two-column TSV: <tip_name>\t<Prometheus_ref|Paralog_ref>
# --catalog is prometheus_catalog.tsv (protein_id/locus_id + confidence_tier);
#   catalog protein_id matches a gene-tree tip 'genome__orig' exactly.
# Saves {plot-prefix}.png/.svg via save_plot. Trivial tree -> placeholder.
# quit(status=0) always.
# =====================================================================

suppressWarnings(suppressMessages({ library(optparse) }))

.script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  ""
}
.this_dir <- tryCatch(dirname(.script_path()), error = function(e) ".")
if (is.na(.this_dir) || !nzchar(.this_dir)) .this_dir <- "."
.theme_ok <- tryCatch({
  suppressWarnings(suppressMessages(source(file.path(.this_dir, "plotting_theme.R")))); TRUE
}, error = function(e) FALSE)

.placeholder <- function(prefix, msg) {
  if (.theme_ok && exists("placeholder_plot"))
    tryCatch(placeholder_plot(prefix, msg), error = function(e) NULL)
  else for (dev in list(
    function() grDevices::png(paste0(prefix, ".png"), 4, 2.2, units = "in", res = 300),
    function() grDevices::svg(paste0(prefix, ".svg"), 4, 2.2)))
    tryCatch({ dev(); graphics::plot.new(); graphics::text(0.5, 0.5, msg)
               grDevices::dev.off() }, error = function(e) NULL)
}

is_trivial_tree <- function(path) {
  txt <- tryCatch(paste(readLines(path, warn = FALSE), collapse = ""),
                  error = function(e) "")
  !nzchar(trimws(txt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(txt))
}

main <- function() {
  option_list <- list(
    make_option("--tree", type = "character"),
    make_option("--labels", type = "character"),
    make_option("--catalog", type = "character"),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  tryCatch({
    stopifnot(!is.null(opt$tree), !is.null(opt$plot_prefix))

    suppressWarnings(suppressMessages({ library(ape) }))

    if (is_trivial_tree(opt$tree)) {
      .placeholder(opt$plot_prefix, "Validation tree: empty tree")
      message("plot_validation_tree: trivial tree; placeholder."); return(invisible(NULL))
    }
    tree <- tryCatch(read.tree(opt$tree), error = function(e) NULL)
    if (is.null(tree) || length(tree$tip.label) < 3) {
      .placeholder(opt$plot_prefix, "Validation tree: <3 tips")
      message("plot_validation_tree: <3 tips; placeholder."); return(invisible(NULL))
    }

    have_ggtree <- requireNamespace("ggtree", quietly = TRUE) &&
                   requireNamespace("ggplot2", quietly = TRUE)

    # ---- reference labels (HEADERLESS: tip<TAB>group) -------------- #
    ref <- setNames(character(0), character(0))
    if (!is.null(opt$labels) && file.exists(opt$labels) &&
        file.info(opt$labels)$size > 0) {
      lab <- tryCatch(utils::read.table(opt$labels, sep = "\t", header = FALSE,
                                        stringsAsFactors = FALSE, quote = "",
                                        comment.char = ""),
                      error = function(e) NULL)
      if (!is.null(lab) && ncol(lab) >= 2) {
        ref <- setNames(trimws(as.character(lab[[2]])), trimws(as.character(lab[[1]])))
      }
    }

    # ---- catalog tiers (protein_id/locus_id -> confidence_tier) ---- #
    cat_tier <- setNames(character(0), character(0))
    if (!is.null(opt$catalog) && file.exists(opt$catalog)) {
      catdf <- tryCatch(
        suppressWarnings(suppressMessages(
          readr::read_tsv(opt$catalog, show_col_types = FALSE, progress = FALSE))),
        error = function(e) NULL)
      if (!is.null(catdf) && nrow(catdf) > 0) {
        idcol <- intersect(c("protein_id", "locus_id"), names(catdf))
        if (length(idcol) && "confidence_tier" %in% names(catdf)) {
          ids <- trimws(as.character(catdf[[idcol[1]]]))
          cat_tier <- setNames(trimws(as.character(catdf$confidence_tier)), ids)
        }
      }
    }

    # ---- assign each tip a category ------------------------------- #
    tips <- tree$tip.label
    category <- vapply(tips, function(t) {
      if (t %in% names(ref)) return(ref[[t]])           # Prometheus_ref / Paralog_ref
      if (t %in% names(cat_tier)) {
        tv <- cat_tier[[t]]
        if (tv %in% c("A", "B", "C", "Rejected")) return(paste0("Candidate_", tv))
      }
      "other"
    }, character(1))

    cat_levels <- c("Prometheus_ref", "Paralog_ref",
                    "Candidate_A", "Candidate_B", "Candidate_C", "Candidate_Rejected",
                    "other")
    cat_colors <- c(
      Prometheus_ref     = "#000000",
      Paralog_ref        = "#CC79A7",
      Candidate_A        = "#0072B2",
      Candidate_B        = "#56B4E9",
      Candidate_C        = "#E69F00",
      Candidate_Rejected = "#D55E00",
      other              = "#BBBBBB")
    present_levels <- intersect(cat_levels, unique(category))

    n_tips <- length(tips)
    show_labels <- n_tips <= 90
    circular <- n_tips > 60
    h <- if (circular) max(6, min(0.10 * n_tips, 14)) else max(4, min(0.13 * n_tips, 16))

    drawn <- FALSE
    if (have_ggtree) {
      # ggtree is designed to be attached; the %<+% data-join operator has no
      # clean ggtree:: infix form, so attach the package (already confirmed
      # present via requireNamespace above).
      suppressWarnings(suppressMessages({ library(ggtree); library(ggplot2) }))
      df <- data.frame(label = tips,
                       category = factor(category, levels = present_levels),
                       stringsAsFactors = FALSE)
      p <- tryCatch({
        gt <- ggtree::ggtree(tree,
                             layout = if (circular) "circular" else "rectangular",
                             linewidth = 0.25)
        gt <- gt %<+% df
        gt <- gt + ggtree::geom_tippoint(ggplot2::aes(color = category), size = 1.1) +
          ggplot2::scale_color_manual(values = cat_colors[present_levels],
                                      name = "Tip category", drop = FALSE,
                                      na.value = "#BBBBBB")
        if (show_labels)
          gt <- gt + ggtree::geom_tiplab(size = 1.2, offset = 0.002)
        gt + ggplot2::ggtitle("Prometheus validation tree: candidates vs references") +
          ggplot2::theme(legend.position = "right",
                         legend.key.size = ggplot2::unit(0.6, "lines"),
                         legend.text = ggplot2::element_text(size = 6),
                         legend.title = ggplot2::element_text(size = 7),
                         plot.title = ggplot2::element_text(size = 8, face = "bold"))
      }, error = function(e) { message("plot_validation_tree: ggtree error: ",
                                       conditionMessage(e)); NULL })
      if (!is.null(p) && exists("save_plot")) {
        w <- if (circular) h else 7.5
        save_plot(p, opt$plot_prefix, width = w, height = h)
        drawn <- TRUE
      }
    }

    # ---- base-graphics fallback (ape) ----------------------------- #
    if (!drawn) {
      tipcols <- unname(cat_colors[category]); tipcols[is.na(tipcols)] <- "#BBBBBB"
      draw <- function() {
        ape::plot.phylo(tree, type = if (circular) "fan" else "phylogram",
                        show.tip.label = show_labels, cex = 0.4,
                        tip.color = tipcols, no.margin = FALSE,
                        edge.width = 0.5)
        ape::tiplabels(pch = 19, col = tipcols, cex = 0.5)
        graphics::legend("bottomleft", legend = present_levels,
                         col = cat_colors[present_levels], pch = 19,
                         cex = 0.5, bty = "n")
        graphics::title("Prometheus validation tree (candidates vs references)",
                        cex.main = 0.7)
      }
      if (exists("save_plot")) {
        w <- if (circular) h else 7.5
        save_plot(draw, opt$plot_prefix, width = w, height = h)
        drawn <- TRUE
      }
    }
    if (!drawn) .placeholder(opt$plot_prefix, "Validation tree (plotting unavailable)")

    message(sprintf("plot_validation_tree: %d tips (%d categories) -> %s.{png,svg}",
                    n_tips, length(present_levels), opt$plot_prefix))
    invisible(NULL)

  }, error = function(e) {
    message("plot_validation_tree: top-level error: ", conditionMessage(e))
    tryCatch(.placeholder(opt$plot_prefix, "Validation tree: error"),
             error = function(e2) invisible(NULL))
    invisible(NULL)
  })

  quit(status = 0)
}

main()
