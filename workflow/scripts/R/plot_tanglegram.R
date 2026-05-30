#!/usr/bin/env Rscript
# =====================================================================
# Stage 5/7 figure -- gene-tree vs species-tree tanglegram (design Stage
# 5.2/7.1). Each Prometheus gene-tree tip (genome__orig) is linked to its
# host genome in the species tree; links are coloured by the host's
# methanogen status. Visualises whether Prometheus's gene history tracks
# the archaeal species history (vertical) or crosses it (transfer).
#
# Uses phytools::cophylo(species_tree, gene_tree, assoc) when available
# (rotates tips to minimise crossings) and falls back to ape::cophyloplot.
# Saves {plot-prefix}.png/.svg via save_plot. Trivial trees ->
# placeholder_plot, quit(status=0).
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
  else {
    for (dev in list(
      function() grDevices::png(paste0(prefix, ".png"), 4, 2.2, units = "in", res = 300),
      function() grDevices::svg(paste0(prefix, ".svg"), 4, 2.2)))
      tryCatch({ dev(); graphics::plot.new(); graphics::text(0.5, 0.5, msg)
                 grDevices::dev.off() }, error = function(e) NULL)
  }
}

is_trivial_tree <- function(path) {
  txt <- tryCatch(paste(readLines(path, warn = FALSE), collapse = ""),
                  error = function(e) "")
  !nzchar(trimws(txt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(txt))
}

as_logical_col <- function(x) {
  s <- tolower(trimws(as.character(x))); s %in% c("true", "1", "yes", "t")
}

main <- function() {
  option_list <- list(
    make_option("--gene-tree", type = "character", dest = "gene_tree"),
    make_option("--species-tree", type = "character", dest = "species_tree"),
    make_option("--master", type = "character"),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  tryCatch({
    stopifnot(!is.null(opt$gene_tree), !is.null(opt$species_tree),
              !is.null(opt$plot_prefix))

    suppressWarnings(suppressMessages({ library(ape) }))

    if (is_trivial_tree(opt$gene_tree) || is_trivial_tree(opt$species_tree)) {
      .placeholder(opt$plot_prefix, "Tanglegram: empty gene or species tree")
      message("plot_tanglegram: trivial tree(s); placeholder."); return(invisible(NULL))
    }
    gtree <- tryCatch(read.tree(opt$gene_tree), error = function(e) NULL)
    stree <- tryCatch(read.tree(opt$species_tree), error = function(e) NULL)
    if (is.null(gtree) || is.null(stree) ||
        length(gtree$tip.label) < 3 || length(stree$tip.label) < 3) {
      .placeholder(opt$plot_prefix, "Tanglegram: <3 tips in a tree")
      message("plot_tanglegram: <3 tips; placeholder."); return(invisible(NULL))
    }

    # gene tip 'genome__orig' -> genome (species tip)
    g_geno <- sub("__.*", "", gtree$tip.label)
    # association: only gene tips whose genome is a species-tree tip
    on_species <- g_geno %in% stree$tip.label
    if (sum(on_species) < 2) {
      .placeholder(opt$plot_prefix, "Tanglegram: <2 gene tips map to species tree")
      message("plot_tanglegram: too few links; placeholder."); return(invisible(NULL))
    }
    assoc <- cbind(stree = g_geno[on_species],
                   gene  = gtree$tip.label[on_species])

    # methanogen colour by host genome (from master)
    link_col <- rep("#999999", nrow(assoc))  # grey default
    if (!is.null(opt$master) && file.exists(opt$master)) {
      master <- tryCatch(
        suppressWarnings(suppressMessages(
          readr::read_tsv(opt$master, show_col_types = FALSE, progress = FALSE))),
        error = function(e) NULL)
      if (!is.null(master) && all(c("genome_id", "methanogen") %in% names(master))) {
        meth <- setNames(as_logical_col(master$methanogen),
                         trimws(as.character(master$genome_id)))
        mv <- meth[assoc[, "stree"]]
        link_col <- ifelse(!is.na(mv) & mv, "#009E73", "#E69F00")  # green meth / orange non
      }
    }

    n_tips <- max(length(gtree$tip.label), length(stree$tip.label))
    h <- max(3.5, min(0.16 * n_tips, 13))
    show_labels <- n_tips <= 70

    have_cophylo <- requireNamespace("phytools", quietly = TRUE)

    draw <- function() {
      if (have_cophylo) {
        obj <- tryCatch(
          phytools::cophylo(stree, gtree, assoc = assoc, rotate = TRUE),
          error = function(e) NULL)
        if (!is.null(obj)) {
          # link colours must align to obj$assoc order; cophylo keeps assoc rows.
          lc <- link_col
          tryCatch(
            plot(obj, link.col = lc, link.lwd = 1.2,
                 fsize = if (show_labels) 0.4 else 1e-6,
                 link.type = "curved", pts = FALSE),
            error = function(e)
              plot(obj, fsize = if (show_labels) 0.4 else 1e-6))
          graphics::title(main = "Prometheus gene tree vs archaeal species tree\n(links: green=methanogen, orange=non-methanogen host)",
                          cex.main = 0.7)
          return(invisible(NULL))
        }
      }
      # fallback: ape::cophyloplot
      ape::cophyloplot(stree, gtree, assoc = assoc, length.line = 1,
                       space = 60, gap = 3, col = link_col,
                       show.tip.label = show_labels, font = 1)
      graphics::title(main = "Prometheus gene tree vs species tree (cophyloplot)",
                      cex.main = 0.7)
    }

    if (.theme_ok && exists("save_plot")) {
      save_plot(draw, opt$plot_prefix, width = 7.5, height = h)
    } else {
      .placeholder(opt$plot_prefix, "Tanglegram (plotting helper unavailable)")
    }
    message(sprintf("plot_tanglegram: %d links (gene<->species) -> %s.{png,svg}",
                    nrow(assoc), opt$plot_prefix))
    invisible(NULL)

  }, error = function(e) {
    message("plot_tanglegram: top-level error: ", conditionMessage(e))
    tryCatch(.placeholder(opt$plot_prefix, "Tanglegram: error"),
             error = function(e2) invisible(NULL))
    invisible(NULL)
  })

  quit(status = 0)
}

main()
