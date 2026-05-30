#!/usr/bin/env Rscript
# =====================================================================
# Stage 7 master figure -- the annotated archaeal species tree (design
# Stage 7.1). A ggtree of the species tree with concentric annotation
# layers added by gheatmap (ggnewscale resets the fill scale between
# layers so each variable keeps its own colour-blind-safe palette):
#   1. prometheus_tier  (A / B / C / none)
#   2. methanogen       (TRUE / FALSE)
#   3. methanogenesis_subtype
#   4. niche_richness   (rowSums of the species-level niche_* columns =
#                        number of niches the species occupies; continuous heat)
#      + a few high-level group_* presence rings (group_Human, group_Marine,
#        group_Environments) when those columns are present
#   5. proviral load    (n_proviruses, continuous heat)
# Title notes the total inferred gains from --gains (gains_losses.tsv).
# Saves {plot-prefix}.png/.svg via save_plot. Trivial / oversized-but-OK
# trees are guarded with placeholder_plot. quit(status=0) always.
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
as_logical_col <- function(x) {
  s <- tolower(trimws(as.character(x))); s %in% c("true", "1", "yes", "t")
}

# colour-blind-safe categorical palette generator (Okabe-Ito recycled)
ok_palette <- function(levels) {
  pal <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9",
           "#D55E00", "#F0E442", "#999999", "#000000")
  setNames(pal[(seq_along(levels) - 1) %% length(pal) + 1], levels)
}

main <- function() {
  option_list <- list(
    make_option("--tree", type = "character"),
    make_option("--master", type = "character"),
    make_option("--viral", type = "character"),
    make_option("--gains", type = "character"),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix")
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  tryCatch({
    stopifnot(!is.null(opt$tree), !is.null(opt$master), !is.null(opt$plot_prefix))

    suppressWarnings(suppressMessages({
      library(ape); library(readr); library(dplyr)
    }))

    if (is_trivial_tree(opt$tree)) {
      .placeholder(opt$plot_prefix, "Master tree: empty species tree")
      message("plot_tree_master: trivial tree; placeholder."); return(invisible(NULL))
    }
    tree <- tryCatch(read.tree(opt$tree), error = function(e) NULL)
    if (is.null(tree) || length(tree$tip.label) < 3) {
      .placeholder(opt$plot_prefix, "Master tree: <3 tips")
      message("plot_tree_master: <3 tips; placeholder."); return(invisible(NULL))
    }

    have_ggtree <- requireNamespace("ggtree", quietly = TRUE) &&
                   requireNamespace("ggplot2", quietly = TRUE)
    if (!have_ggtree) {
      .placeholder(opt$plot_prefix, "Master tree: ggtree/ggplot2 unavailable")
      message("plot_tree_master: ggtree missing; placeholder."); return(invisible(NULL))
    }
    # Attach so gheatmap + ggnewscale `+`-chaining dispatch reliably.
    suppressWarnings(suppressMessages({ library(ggtree); library(ggplot2) }))

    master <- readr::read_tsv(opt$master, show_col_types = FALSE, progress = FALSE)
    if (nrow(master) == 0 || !"genome_id" %in% names(master)) {
      .placeholder(opt$plot_prefix, "Master tree: empty master table")
      return(invisible(NULL))
    }
    master$genome_id <- trimws(as.character(master$genome_id))

    # merge viral n_proviruses
    if (!is.null(opt$viral) && file.exists(opt$viral)) {
      viral <- tryCatch(readr::read_tsv(opt$viral, show_col_types = FALSE,
                                        progress = FALSE), error = function(e) NULL)
      if (!is.null(viral) && "genome_id" %in% names(viral) &&
          "n_proviruses" %in% names(viral)) {
        viral$genome_id <- trimws(as.character(viral$genome_id))
        master <- dplyr::left_join(
          master, viral[, c("genome_id", "n_proviruses")], by = "genome_id")
      }
    }
    if (!"n_proviruses" %in% names(master)) master$n_proviruses <- NA_real_

    # total gains from gains_losses.tsv
    total_gains <- NA_real_
    if (!is.null(opt$gains) && file.exists(opt$gains)) {
      gl <- tryCatch(readr::read_tsv(opt$gains, show_col_types = FALSE,
                                     progress = FALSE), error = function(e) NULL)
      if (!is.null(gl) && all(c("metric", "value") %in% names(gl))) {
        hit <- gl$value[gl$metric == "n_independent_acquisitions"]
        if (length(hit)) total_gains <- suppressWarnings(as.numeric(hit[1]))
      }
    }

    # align annotation frame to tree tips (species tips = genome_id)
    ann <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
    m <- master[match(ann$label, master$genome_id), , drop = FALSE]
    tier <- trimws(as.character(m$prometheus_tier)); tier[is.na(tier) | tier == ""] <- "none"
    ann$tier <- factor(tier, levels = c("A", "B", "C", "none"))
    ann$methanogen <- ifelse(as_logical_col(m$methanogen), "yes", "no")
    sub <- trimws(as.character(m$methanogenesis_subtype))
    sub[is.na(sub) | sub == "" | sub == "NA"] <- "none"
    ann$subtype <- sub

    # ---- ecology rings from the species-level niche_*/group_* columns ----- #
    as_int01 <- function(x) {
      v <- as_logical_col(x)
      v2 <- suppressWarnings(as.numeric(as.character(x)))
      as.integer(v | (!is.na(v2) & v2 != 0))
    }
    # niche richness = how many niche_* the species occupies (rowSums of 0/1)
    niche_cols <- grep("^niche_", names(m), value = TRUE)
    if (length(niche_cols)) {
      nm <- as.data.frame(lapply(m[, niche_cols, drop = FALSE], as_int01))
      ann$niche_richness <- rowSums(nm, na.rm = TRUE)
    } else {
      ann$niche_richness <- NA_real_
    }
    # optional high-level group rings (only those actually present in the table)
    group_ring_cols <- intersect(c("group_Human", "group_Marine", "group_Environments"),
                                 names(m))
    for (gc in group_ring_cols) {
      ann[[gc]] <- ifelse(as_int01(m[[gc]]) == 1L, "yes", "no")
    }

    ann$n_proviruses <- suppressWarnings(as.numeric(as.character(m$n_proviruses)))

    n_tips <- length(tree$tip.label)
    show_labels <- n_tips <= 80
    circular <- n_tips > 40
    h <- if (circular) max(6, min(0.10 * n_tips, 14)) else max(4, min(0.14 * n_tips, 16))

    ggt <- ggtree::ggtree(tree,
                          layout = if (circular) "circular" else "rectangular",
                          linewidth = 0.25)
    if (show_labels) {
      ggt <- ggt + ggtree::geom_tiplab(size = 1.2,
                                       align = TRUE, linesize = 0.1,
                                       offset = 0.001)
    }

    rn <- ann; rownames(rn) <- rn$label

    # helper to add one categorical ring with its own fill scale
    add_ring <- function(p, df1col, palette, name, width = 0.06, offset = 0) {
      p2 <- p
      if (requireNamespace("ggnewscale", quietly = TRUE))
        p2 <- p2 + ggnewscale::new_scale_fill()
      ggtree::gheatmap(p2, df1col, width = width, offset = offset,
                       colnames = FALSE, color = NA) +
        ggplot2::scale_fill_manual(values = palette, name = name,
                                   na.value = "grey95")
    }

    # helper to add one continuous heat ring with its own fill scale
    add_heat_ring <- function(p, df1col, name, low, high, width = 0.06, offset = 0) {
      p2 <- p
      if (requireNamespace("ggnewscale", quietly = TRUE))
        p2 <- p2 + ggnewscale::new_scale_fill()
      ggtree::gheatmap(p2, df1col, width = width, offset = offset,
                       colnames = FALSE, color = NA) +
        ggplot2::scale_fill_gradient(low = low, high = high, name = name,
                                     na.value = "grey95")
    }

    # offsets scale with tree depth; use a fraction of max branch distance
    depth <- tryCatch(max(ape::node.depth.edgelength(tree)), error = function(e) 1)
    if (!is.finite(depth) || depth <= 0) depth <- 1
    step <- depth * 0.08

    p <- ggt
    # ring 1: tier
    p <- add_ring(p, data.frame(tier = rn$tier, row.names = rownames(rn)),
                  c(A = "#0072B2", B = "#56B4E9", C = "#999999", none = "grey92"),
                  "Prometheus tier", width = 0.06, offset = step * 0)
    # ring 2: methanogen
    p <- add_ring(p, data.frame(methanogen = rn$methanogen, row.names = rownames(rn)),
                  c(yes = "#009E73", no = "#F0E442"),
                  "Methanogen", width = 0.06, offset = step * 1)
    # ring 3: subtype
    p <- add_ring(p, data.frame(subtype = rn$subtype, row.names = rownames(rn)),
                  ok_palette(sort(unique(rn$subtype))),
                  "Methanogenesis\nsubtype", width = 0.06, offset = step * 2)
    # Running offset index so the ecology rings (variable in number) and the
    # proviral ring never overlap. Rings 1-3 used indices 0,1,2.
    ri <- 3
    # ring 4: niche richness (continuous heat = number of niches occupied)
    if ("niche_richness" %in% names(rn) && any(is.finite(rn$niche_richness))) {
      p <- add_heat_ring(p, data.frame(niche_richness = rn$niche_richness,
                                       row.names = rownames(rn)),
                         "Niche richness\n(n niches)", "#F7FBFF", "#08306B",
                         width = 0.06, offset = step * ri)
      ri <- ri + 1
    }
    # ring 4b: optional high-level group presence rings (those present in table)
    for (gc in group_ring_cols) {
      lbl <- sub("^group_", "", gc)
      p <- add_ring(p, setNames(data.frame(rn[[gc]], row.names = rownames(rn)), gc),
                    c(yes = "#CC79A7", no = "grey92"),
                    sprintf("In %s", lbl), width = 0.05, offset = step * ri)
      ri <- ri + 1
    }
    # final ring: proviral load (continuous)
    p <- add_heat_ring(p, data.frame(n_proviruses = rn$n_proviruses,
                                     row.names = rownames(rn)),
                       "Proviral load\n(n_proviruses)", "#FFF7EC", "#7F0000",
                       width = 0.06, offset = step * ri)

    ttl <- if (!is.na(total_gains))
      sprintf("Archaeal Prometheus distribution (inferred independent acquisitions = %.2f)", total_gains)
      else "Archaeal Prometheus distribution"
    p <- p + ggplot2::ggtitle(ttl) +
      ggplot2::theme(legend.position = "right",
                     legend.key.size = ggplot2::unit(0.5, "lines"),
                     legend.text = ggplot2::element_text(size = 5),
                     legend.title = ggplot2::element_text(size = 6),
                     plot.title = ggplot2::element_text(size = 8, face = "bold"))

    if (exists("save_plot")) {
      w <- if (circular) h else 7.5
      save_plot(p, opt$plot_prefix, width = w, height = h)
    } else {
      .placeholder(opt$plot_prefix, ttl)
    }
    message(sprintf("plot_tree_master: %d tips, %s layout -> %s.{png,svg}",
                    n_tips, if (circular) "circular" else "rect", opt$plot_prefix))
    invisible(NULL)

  }, error = function(e) {
    message("plot_tree_master: top-level error: ", conditionMessage(e))
    tryCatch(.placeholder(opt$plot_prefix, "Master tree: error"),
             error = function(e2) invisible(NULL))
    invisible(NULL)
  })

  quit(status = 0)
}

main()
