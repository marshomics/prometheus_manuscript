#!/usr/bin/env Rscript
# =====================================================================
# Stage 7.1 master figure -- the annotated archaeal species tree.
#
# ggtree of the ROOTED species tree (scaffolds/species_tree.treefile) with
# concentric gheatmap rings; ggnewscale::new_scale_fill() gives every ring
# its own fill scale and legend (inner -> outer):
#   1. Prometheus tier      A / B / C (best validated tier) / absent
#   2. mcr metabolism       methanogen / ANME (methanotroph) / neither
#   3. methanogenesis subtype (methanogens only; others "not methanogen")
#   4. niche richness       number of niche_* columns = 1 (if any niche_*; NA cells add nothing)
#   5. group rings          group_Human / group_Marine / group_Environments
#   6. proviral load        viral_summary n_proviruses
#   7. R-loop-prone fraction frac_rloop_proviruses (NA when no provirus)
# Tips are matched to master_table rows by tree_tip, then genome_id; a tip
# with no master / viral row is painted "no data" (white), NEVER as an
# absence. The title reports the mean number of independent acquisitions
# (gains_losses.tsv n_independent_acquisitions).
#
# Failure policy (Contract v2): missing packages, unreadable inputs, schema
# violations, no tip matching the master table, or a render error stop the
# script with a non-zero status. A labelled placeholder is written only for a
# genuinely empty input (an empty / '();' species tree or < 3 tips).
#
# Optional --annotation-out writes the per-tip ring values (figure source data).
# =====================================================================

suppressPackageStartupMessages(library(optparse))

.script_dir <- function() {
  m <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else "."
}

NO_DATA <- "no data"
NO_DATA_FILL <- "#FFFFFF"
ABSENT_FILL <- "#E8E8E8"
NA_TOKENS <- c("", "na", "nan", "none", "null")

as_bool <- function(x) {
  s <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(s))
  out[s %in% c("true", "t", "1", "yes", "1.0")] <- TRUE
  out[s %in% c("false", "f", "0", "no", "0.0")] <- FALSE
  out
}

is_missing_tok <- function(x) is.na(x) | tolower(trimws(as.character(x))) %in% NA_TOKENS

as_num <- function(x) {
  x <- as.character(x); x[is_missing_tok(x)] <- NA_character_
  v <- suppressWarnings(as.numeric(x))
  bad <- !is.na(x) & is.na(v)
  if (any(bad)) stop(sprintf("non-numeric value(s) %s", paste(head(unique(x[bad]), 3),
                                                            collapse = ", ")))
  v
}

read_tab <- function(path, what) {
  if (is.null(path) || !nzchar(path)) stop(sprintf("--%s is required", what))
  if (!file.exists(path)) stop(sprintf("%s not found: %s", what, path))
  if (file.info(path)$size == 0)
    stop(sprintf("%s is a 0-byte file (a header is expected): %s", what, path))
  utils::read.delim(path, sep = "\t", header = TRUE, colClasses = "character",
                    quote = "\"", na.strings = character(0), check.names = FALSE,
                    comment.char = "")
}

is_trivial_tree <- function(path) {
  txt <- trimws(paste(readLines(path, warn = FALSE), collapse = ""))
  !nzchar(txt) || grepl("^\\(\\s*\\)\\s*;?$", txt)
}

strip_quotes <- function(x) sub("^'(.*)'$", "\\1", sub("^\"(.*)\"$", "\\1", x))

ANN_COLS <- c("label", "matched", "tier", "metabolism", "subtype", "n_proviruses")
write_empty_annotation <- function(path) {
  if (!nzchar(path)) return(invisible(NULL))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(ANN_COLS, collapse = "\t"), path)
}

# Left-join `col` of `tab` (keyed by key_col) onto the tip labels; NA when absent.
lookup <- function(tab, key_col, col, labels) {
  tab[[col]][match(labels, trimws(tab[[key_col]]))]
}

main <- function() {
  opt <- parse_args(OptionParser(option_list = list(
    make_option("--tree", type = "character"),
    make_option("--master", type = "character"),
    make_option("--viral", type = "character"),
    make_option("--gains", type = "character"),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix"),
    make_option("--annotation-out", type = "character", dest = "annotation_out",
                default = "", help = "optional TSV of the per-tip ring values")
  )))
  for (a in c("tree", "master", "viral", "gains", "plot_prefix"))
    if (is.null(opt[[a]]) || !nzchar(opt[[a]])) stop(sprintf("--%s is required",
                                                               gsub("_", "-", a)))
  # presence check without loading (cheap); the namespaces are loaded before plotting
  for (p in c("ape", "ggplot2", "ggtree", "ggnewscale"))
    if (!nzchar(system.file(package = p))) stop(sprintf("R package '%s' is required", p))
  source(file.path(.script_dir(), "plotting_theme.R"))   # save_plot, placeholder_plot
  dir.create(dirname(opt$plot_prefix), showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(opt$tree)) stop("species tree not found: ", opt$tree)
  if (is_trivial_tree(opt$tree)) {
    # genuinely empty input (species_tree writes '();' for < 4 genomes)
    placeholder_plot(opt$plot_prefix, "Master tree: empty species tree")
    write_empty_annotation(opt$annotation_out)
    message("plot_tree_master: empty species tree; placeholder written")
    return(invisible(NULL))
  }
  tree <- ape::read.tree(opt$tree)
  if (is.null(tree) || is.null(tree$tip.label))
    stop("could not parse the species tree: ", opt$tree)
  tree$tip.label <- strip_quotes(trimws(tree$tip.label))
  if (anyDuplicated(tree$tip.label))
    stop("duplicated tip labels in the species tree: ",
         paste(head(unique(tree$tip.label[duplicated(tree$tip.label)]), 3), collapse = ", "))
  if (length(tree$tip.label) < 3) {
    placeholder_plot(opt$plot_prefix, "Master tree: fewer than 3 tips")
    write_empty_annotation(opt$annotation_out)
    message("plot_tree_master: < 3 tips; placeholder written")
    return(invisible(NULL))
  }
  tips <- tree$tip.label
  n_tips <- length(tips)
  if (!ape::is.rooted(tree))
    message("plot_tree_master: WARNING the species tree is unrooted; it is drawn from ",
            "its basal node (species_tree.treefile is expected to be rooted upstream)")

  # ---- master table -> tips (tree_tip first, genome_id second) ------------ #
  master <- read_tab(opt$master, "master")
  need <- c("genome_id", "prometheus_tier", "methanogen")
  if (!all(need %in% names(master)))
    stop("master table lacks column(s): ", paste(setdiff(need, names(master)), collapse = ", "))
  if (nrow(master) == 0) stop("master table is empty but the species tree has ", n_tips, " tips")
  gid <- trimws(master$genome_id)
  if (anyDuplicated(gid)) stop("duplicated genome_id in master table: ",
                               paste(head(unique(gid[duplicated(gid)]), 3), collapse = ", "))
  ttip <- if ("tree_tip" %in% names(master)) trimws(master$tree_tip) else rep("", nrow(master))
  ttip[is_missing_tok(ttip)] <- NA_character_
  idx <- match(tips, ttip)
  idx[is.na(idx)] <- match(tips[is.na(idx)], gid)
  matched <- !is.na(idx)
  if (!any(matched))
    stop(sprintf("no species-tree tip matches master_table tree_tip/genome_id (tips e.g. %s; master e.g. %s)",
                 paste(head(tips, 3), collapse = ", "), paste(head(gid, 3), collapse = ", ")))
  if (!all(matched))
    message(sprintf("plot_tree_master: WARNING %d/%d tip(s) have no master row (painted '%s'): %s",
                    sum(!matched), n_tips, NO_DATA, paste(head(tips[!matched], 5), collapse = ", ")))
  m <- master[idx, , drop = FALSE]          # NA rows for unmatched tips

  ann <- data.frame(label = tips, matched = matched, stringsAsFactors = FALSE)

  # ring 1: best validated tier (Contract v2 item 6: B/C are never shown as absent)
  tier <- trimws(m$prometheus_tier)
  bad <- matched & !(is_missing_tok(tier) | tier %in% c("A", "B", "C"))
  if (any(bad)) stop("unexpected prometheus_tier value(s): ",
                     paste(unique(tier[bad]), collapse = ", "), " (expected A/B/C or '')")
  ann$tier <- ifelse(!matched, NO_DATA, ifelse(is_missing_tok(tier), "absent", tier))

  # ring 2: mcr metabolism (ANME carry McrA but are methanotrophs: methanogen=False)
  meth <- as_bool(m$methanogen)
  anme <- if ("is_anme" %in% names(m)) as_bool(m$is_anme) else rep(FALSE, n_tips)
  anme[matched & is.na(anme)] <- FALSE
  ann$metabolism <- ifelse(!matched | is.na(meth), NO_DATA,
                           ifelse(anme, "ANME (methanotroph)",
                                  ifelse(meth, "methanogen", "neither")))

  # ring 3: methanogenesis subtype, defined for methanogens only
  sub <- if ("methanogenesis_subtype" %in% names(m)) trimws(m$methanogenesis_subtype)
         else rep(NA_character_, n_tips)
  sub[is_missing_tok(sub)] <- "unresolved"
  ann$subtype <- ifelse(ann$metabolism == NO_DATA, NO_DATA,
                        ifelse(ann$metabolism == "methanogen", sub, "not methanogen"))

  # ring 4: niche richness; ring 5: selected high-level groups
  niche_cols <- grep("^niche_", names(m), value = TRUE)
  if (length(niche_cols)) {
    nm <- vapply(niche_cols, function(cn) as_num(m[[cn]]), numeric(n_tips))
    nm <- matrix(nm, nrow = n_tips)
    # observed-niche count: an unknown (NA) cell adds nothing (unknown is not absence, and
    # not presence); NA ("no data") only when every niche of the species is unknown
    n_known <- rowSums(!is.na(nm))
    ann$niche_richness <- ifelse(matched & n_known > 0, rowSums(nm != 0, na.rm = TRUE), NA_real_)
  }
  group_cols <- intersect(c("group_Human", "group_Marine", "group_Environments"), names(m))
  for (gc in group_cols) {
    v <- as_num(m[[gc]])
    ann[[gc]] <- ifelse(!matched | is.na(v), NO_DATA, ifelse(v != 0, "yes", "no"))
  }

  # rings 6-7: proviral load and R-loop-prone fraction (viral_summary)
  viral <- read_tab(opt$viral, "viral")
  if (nrow(viral) == 0) {
    message("plot_tree_master: WARNING viral summary is header-only; viral rings are 'no data'")
    viral <- data.frame(genome_id = character(0), n_proviruses = character(0),
                        frac_rloop_proviruses = character(0))
  }
  if (!all(c("genome_id", "n_proviruses") %in% names(viral)))
    stop("viral summary lacks genome_id / n_proviruses")
  if (anyDuplicated(trimws(viral$genome_id))) stop("duplicated genome_id in viral summary")
  ann$n_proviruses <- as_num(lookup(viral, "genome_id", "n_proviruses", tips))
  if ("frac_rloop_proviruses" %in% names(viral)) {
    fr <- as_num(lookup(viral, "genome_id", "frac_rloop_proviruses", tips))
    fr[!is.na(ann$n_proviruses) & ann$n_proviruses == 0] <- NA_real_   # undefined
    ann$frac_rloop <- fr
  }

  # ---- title: independent acquisitions ----------------------------------- #
  gl <- read_tab(opt$gains, "gains")
  total_gains <- NA_real_
  if (nrow(gl) && all(c("metric", "value") %in% names(gl))) {
    hit <- gl$value[trimws(gl$metric) == "n_independent_acquisitions"]
    if (length(hit)) total_gains <- as_num(hit[1])
  } else if (nrow(gl)) stop("gains_losses table lacks metric/value columns")

  # ---- tree + rings ------------------------------------------------------- #
  # inputs are validated; now load the plotting stack (a broken install fails here)
  suppressPackageStartupMessages({ library(ggplot2); library(ggtree); library(ggnewscale) })
  circular <- n_tips > 40
  show_labels <- n_tips <= 80
  ggt <- ggtree::ggtree(tree, layout = if (circular) "circular" else "rectangular",
                        linewidth = 0.25)
  xr <- diff(range(ggt$data$x, na.rm = TRUE)); if (!is.finite(xr) || xr <= 0) xr <- 1
  xmax <- max(ggt$data$x[ggt$data$isTip], na.rm = TRUE)
  ring_w <- if (circular) 0.045 else 0.05     # fraction of the tree's x range
  gap <- 0.02 * xr; pitch <- (ring_w + 0.012) * xr
  ring_i <- 0L
  p <- ggt

  key_guide <- function(order) ggplot2::guide_legend(
    order = order, override.aes = list(colour = "grey55", linewidth = 0.2))
  add_ring <- function(p, column, values, title, kind = c("discrete", "continuous"),
                       low = NULL, high = NULL, limits = NULL) {
    kind <- match.arg(kind)
    d <- data.frame(ann[[column]], row.names = ann$label, check.names = FALSE)
    names(d) <- column
    off <- gap + ring_i * pitch - ring_w * xr / 2     # tile left edge at xmax + gap + i*pitch
    p <- p + ggnewscale::new_scale_fill()
    p <- suppressMessages(ggtree::gheatmap(p, d, offset = off, width = ring_w,
                                           colnames = FALSE, color = NA))
    # gheatmap adds a default fill scale; replace it (ggplot2 notes the replacement)
    if (kind == "discrete") {
      unknown <- setdiff(unique(d[[1]]), names(values))
      if (length(unknown)) stop(sprintf("ring '%s': no colour for value(s) %s", column,
                                        paste(unknown, collapse = ", ")))
      lv <- intersect(names(values), unique(d[[1]]))
      p <- suppressMessages(p + ggplot2::scale_fill_manual(
        values = values, breaks = lv, name = title, na.value = NO_DATA_FILL,
        guide = key_guide(ring_i + 1L)))
    } else {
      p <- suppressMessages(p + ggplot2::scale_fill_gradient(
        low = low, high = high, name = title, na.value = NO_DATA_FILL, limits = limits,
        guide = ggplot2::guide_colourbar(order = ring_i + 1L)))
    }
    ring_i <<- ring_i + 1L
    p
  }
  tc <- if (exists("TIER_COLORS")) TIER_COLORS else c(A = "#0072B2", B = "#56B4E9", C = "#999999")
  p <- add_ring(p, "tier", c(A = tc[["A"]], B = tc[["B"]], C = tc[["C"]],
                             absent = ABSENT_FILL, `no data` = NO_DATA_FILL),
                "Prometheus tier\n(best validated)")
  p <- add_ring(p, "metabolism", c(methanogen = "#009E73", `ANME (methanotroph)` = "#D55E00",
                                   neither = ABSENT_FILL, `no data` = NO_DATA_FILL),
                "mcr metabolism")
  # Subtype ring. methanogenesis_pathways.py writes a composite capacity label
  # ("CO2-reducing+acetoclastic+methyl-dismutating"); the legacy single-word labels are kept so
  # an older master table still plots. Colours are assigned DETERMINISTICALLY and never
  # recycled: two different subtypes must not share a colour. Composites beyond the reserved
  # palette are pooled into one honestly named level rather than silently re-coloured.
  sub_pal <- c(`CO2-reducing` = "#0072B2", hydrogenotrophic = "#0072B2",
               acetoclastic = "#56B4E9", methylotrophic = "#E69F00",
               `methyl-dismutating` = "#D55E00", `methyl-reducing` = "#CC79A7",
               methyl_reducing = "#CC79A7", none_detected = "#F0E442",
               unresolved = "#999999")
  COMBO <- c("#009E73", "#117733", "#44AA99", "#882255", "#AA4499", "#661100")
  extra <- sort(setdiff(unique(ann$subtype), c(names(sub_pal), "not methanogen", NO_DATA)))
  if (length(extra)) {
    freq <- sort(table(ann$subtype[ann$subtype %in% extra]), decreasing = TRUE)
    keep <- names(freq)[seq_len(min(length(freq), length(COMBO)))]
    sub_pal <- c(sub_pal, setNames(COMBO[seq_along(keep)], keep))
    pooled <- setdiff(extra, keep)
    if (length(pooled)) {
      message("plot_tree_master: pooling ", length(pooled),
              " rare subtype combination(s) into 'other combination': ",
              paste(pooled, collapse = ", "))
      ann$subtype[ann$subtype %in% pooled] <- "other combination"
      sub_pal <- c(sub_pal, `other combination` = "#000000")
    }
  }
  sub_pal <- sub_pal[!duplicated(names(sub_pal))]
  p <- add_ring(p, "subtype", c(sub_pal, `not methanogen` = ABSENT_FILL,
                                `no data` = NO_DATA_FILL), "Methanogenesis\nsubtype")
  if ("niche_richness" %in% names(ann))
    p <- add_ring(p, "niche_richness", NULL, "Niche richness\n(n niches)", "continuous",
                  low = "#DEEBF7", high = "#08306B")
  for (gc in group_cols)
    p <- add_ring(p, gc, c(yes = "#CC79A7", no = ABSENT_FILL, `no data` = NO_DATA_FILL),
                  sprintf("In %s", sub("^group_", "", gc)))
  p <- add_ring(p, "n_proviruses", NULL, "Proviral load\n(n proviruses)", "continuous",
                low = "#FEE6CE", high = "#7F2704")
  if ("frac_rloop" %in% names(ann))
    p <- add_ring(p, "frac_rloop", NULL, "R-loop-prone\nprovirus fraction", "continuous",
                  low = "#E5F5E0", high = "#00441B", limits = c(0, 1))
  ring_end <- xmax + gap + ring_i * pitch

  # tip labels outside the rings; dotted guides from each tip to the first ring
  lab_size <- if (n_tips <= 40) 1.6 else 1.2
  if (show_labels) {
    tipd <- ggt$data[ggt$data$isTip, c("x", "y")]
    p <- p + ggplot2::geom_segment(data = tipd, ggplot2::aes(x = x, xend = xmax + gap * 0.8,
                                                            y = y, yend = y),
                                   inherit.aes = FALSE, linetype = "dotted",
                                   linewidth = 0.15, colour = "grey60")
    p <- p + ggtree::geom_tiplab(size = lab_size, align = TRUE, linetype = NA,
                                 offset = ring_end - xmax + 0.01 * xr)
    char_frac <- lab_size / 25.4 * 0.62 * max(nchar(tips)) / (if (circular) 3.2 else 5.0)
    x_hi <- ring_end + 0.01 * xr + char_frac * (ring_end + 0.01 * xr) / max(0.2, 1 - char_frac)
  } else {
    x_hi <- ring_end
  }
  p <- p + ggplot2::geom_blank(data = data.frame(x = x_hi, y = 1), ggplot2::aes(x = x, y = y),
                               inherit.aes = FALSE)
  if (!circular) p <- suppressMessages(p + ggplot2::coord_cartesian(clip = "off"))

  ttl <- if (is.finite(total_gains))
    sprintf("Archaeal Prometheus distribution (%d tips; mean independent acquisitions = %.2f)",
            n_tips, total_gains) else sprintf("Archaeal Prometheus distribution (%d tips)", n_tips)
  p <- p + ggplot2::labs(title = ttl,
                         subtitle = "Rings, inner to outer: tier, mcr metabolism, subtype, ecology, proviral load, R-loop fraction; white = no data") +
    ggplot2::theme(legend.position = "right", legend.box = "vertical",
                   legend.key.size = ggplot2::unit(0.45, "lines"),
                   legend.text = ggplot2::element_text(size = 5),
                   legend.title = ggplot2::element_text(size = 5.5),
                   legend.spacing.y = ggplot2::unit(0.15, "lines"),
                   plot.title = ggplot2::element_text(size = 7.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 5.5, colour = "grey30"),
                   plot.margin = ggplot2::margin(4, 4, 4, 4))

  # ---- size, pre-render (errors surface here), save, verify -------------- #
  # ggplotGrob builds the whole plot, so any data/scale error stops the script here
  # (save_plot in plotting_theme.R would otherwise swallow it into a blank figure).
  # The legend stack is measured in a null device so the figure is tall enough.
  grDevices::pdf(NULL)
  g <- tryCatch(ggplot2::ggplotGrob(p), finally = grDevices::dev.off())
  grDevices::pdf(NULL)
  legend_h <- tryCatch({
    k <- grep("^guide-box", g$layout$name)
    max(c(0, vapply(k, function(j) grid::convertHeight(grid::grobHeight(g$grobs[[j]]),
                                                       "in", valueOnly = TRUE), numeric(1))))
  }, finally = grDevices::dev.off())
  if (circular) {
    h <- max(6, min(0.05 * n_tips, 10), legend_h + 0.8); w <- h + 1.8
  } else {
    h <- max(4, min(0.14 * n_tips, 16), legend_h + 0.8); w <- 7.5
  }
  render_err <- character(0)
  withCallingHandlers(
    save_plot(p, opt$plot_prefix, width = w, height = h),
    message = function(cnd) {
      if (grepl("render error", conditionMessage(cnd))) render_err <<- c(render_err, conditionMessage(cnd))
    })
  if (length(render_err)) stop("rendering failed: ", paste(render_err, collapse = "; "))
  for (ext in plot_formats()) {   # config plotting.formats (plotting_theme.R)
    f <- paste0(opt$plot_prefix, ".", ext)
    if (!file.exists(f) || file.info(f)$size == 0) stop("figure not written: ", f)
  }
  if (nzchar(opt$annotation_out)) {
    dir.create(dirname(opt$annotation_out), showWarnings = FALSE, recursive = TRUE)
    utils::write.table(ann, opt$annotation_out, sep = "\t", quote = FALSE,
                       row.names = FALSE, na = "")
  }
  message(sprintf("plot_tree_master: %d tips (%d matched), %s layout, %d rings, %.1fx%.1f in -> %s.{png,svg}",
                  n_tips, sum(matched), if (circular) "circular" else "rectangular",
                  ring_i, w, h, opt$plot_prefix))
  invisible(NULL)
}

main()
