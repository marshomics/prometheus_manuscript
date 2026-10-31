#!/usr/bin/env Rscript
# =====================================================================
# Stage 1/7 figure -- the Prometheus validation tree (design Stage 1.4 /
# 7.1). A ggtree of the Stage-1 guide tree that places every candidate
# protein against the curated references. Tips are coloured by category:
#   * reference tips from --labels: Prometheus_ref / Paralog_ref
#   * candidate tips from --catalog: by confidence_tier A / B / C / Rejected
#   * anything else: grey ('other')
# Internal nodes with SH-aLRT >= 80 and UFBoot >= 95 (IQ-TREE 'a/b' node labels; the
# thresholds validate_prometheus uses by default) carry a black dot. For display the
# tree is rooted on the paralog outgroup (smallest clade holding every Paralog_ref and no
# Prometheus_ref, as in validate_prometheus) with ape::root(..., edgelabel = TRUE) so the
# support labels stay on the right edges; otherwise it is drawn as inferred and the
# subtitle says so.
#
# --labels is a HEADERLESS two-column TSV: <tip_name>\t<Prometheus_ref|Paralog_ref>
# (ids sanitised identically in the tree; see combine_proteomes.py reference mode)
# --catalog is prometheus_catalog.tsv (protein_id/locus_id + confidence_tier).
# Saves {plot-prefix}.png/.svg via save_plot (plotting_theme.R).
# An empty / placeholder tree ("();", fewer than 3 tips) -> labelled placeholder figure.
# Any other failure (missing package, unreadable input, plotting error) exits non-zero:
# Contract v2 item 4 (no silent failure).
# =====================================================================

suppressWarnings(suppressMessages({ library(optparse) }))

.script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  ""
}
.this_dir <- dirname(.script_path())
if (is.na(.this_dir) || !nzchar(.this_dir)) .this_dir <- "."
source(file.path(.this_dir, "plotting_theme.R"))   # save_plot, placeholder_plot

is_trivial_tree <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  !nzchar(trimws(txt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(txt))
}

parse_support <- function(lbl) {
  # 'a/b' -> c(alrt = a, ufb = b); single value -> UFBoot-scale; else NA
  out <- c(alrt = NA_real_, ufb = NA_real_)
  if (is.null(lbl) || is.na(lbl) || !nzchar(lbl)) return(out)
  v <- suppressWarnings(as.numeric(strsplit(lbl, "/", fixed = TRUE)[[1]]))
  if (any(is.na(v))) return(out)
  if (length(v) == 1) out["ufb"] <- v
  if (length(v) == 2) out[] <- v
  if (length(v) == 3) out[] <- v[c(1, 3)]
  out
}

main <- function() {
  option_list <- list(
    make_option("--tree", type = "character"),
    make_option("--labels", type = "character"),
    make_option("--catalog", type = "character"),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$tree) || is.null(opt$plot_prefix))
    stop("--tree and --plot-prefix are required")
  for (p in c("ape", "ggtree", "ggplot2"))
    if (!requireNamespace(p, quietly = TRUE)) stop("required R package '", p, "' is not installed")

  if (is_trivial_tree(opt$tree)) {
    placeholder_plot(opt$plot_prefix, "Validation tree: no tree (no candidates or < 4 sequences)")
    message("plot_validation_tree: trivial tree; placeholder.")
    return(invisible(NULL))
  }
  tree <- ape::read.tree(opt$tree)
  if (is.null(tree)) stop("could not parse ", opt$tree)
  if (length(tree$tip.label) < 3) {
    placeholder_plot(opt$plot_prefix, "Validation tree: < 3 tips")
    return(invisible(NULL))
  }

  # ---- reference labels (HEADERLESS: tip<TAB>group) ---------------- #
  ref <- setNames(character(0), character(0))
  if (!is.null(opt$labels) && file.exists(opt$labels) && file.info(opt$labels)$size > 0) {
    lab <- utils::read.table(opt$labels, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                             quote = "", comment.char = "", colClasses = "character")
    if (ncol(lab) < 2) stop(opt$labels, " must have two columns (tip, group)")
    ref <- setNames(trimws(lab[[2]]), trimws(lab[[1]]))
  }

  # ---- catalog tiers (protein_id/locus_id -> confidence_tier) ------ #
  cat_tier <- setNames(character(0), character(0))
  if (!is.null(opt$catalog) && file.exists(opt$catalog)) {
    catdf <- utils::read.delim(opt$catalog, sep = "\t", quote = "", stringsAsFactors = FALSE,
                               colClasses = "character", check.names = FALSE)
    if (nrow(catdf) > 0) {
      idcol <- intersect(c("protein_id", "locus_id"), names(catdf))
      if (!length(idcol) || !"confidence_tier" %in% names(catdf))
        stop(opt$catalog, " lacks protein_id/locus_id or confidence_tier")
      cat_tier <- setNames(trimws(catdf$confidence_tier), trimws(catdf[[idcol[1]]]))
    }
  }

  # ---- display rooting on the paralog outgroup --------------------- #
  # Same convention as validate_prometheus: the outgroup is the smallest clade holding
  # every Paralog_ref (candidates nested among them included) and no Prometheus_ref.
  para <- intersect(names(ref)[ref == "Paralog_ref"], tree$tip.label)
  prom <- intersect(names(ref)[ref == "Prometheus_ref"], tree$tip.label)
  root_note <- "drawn as inferred (unrooted)"
  if (length(para) >= 1 && length(prom) >= 1) {
    t1 <- ape::root(tree, outgroup = prom[1], resolve.root = TRUE, edgelabel = TRUE)
    out_tips <- if (length(para) == 1) para else
      ape::extract.clade(t1, ape::getMRCA(t1, para))$tip.label
    if (!any(prom %in% out_tips)) {
      tree <- ape::root(t1, outgroup = out_tips, resolve.root = TRUE, edgelabel = TRUE)
      root_note <- sprintf("rooted on the paralog outgroup (%d tips)", length(out_tips))
    } else {
      root_note <- "drawn as inferred (Paralog_ref and Prometheus_ref tips not separable)"
    }
  }

  tips <- tree$tip.label
  category <- vapply(tips, function(t) {
    if (t %in% names(ref)) return(ref[[t]])
    if (t %in% names(cat_tier)) {
      tv <- cat_tier[[t]]
      if (tv %in% c("A", "B", "C", "Rejected")) return(paste0("Candidate_", tv))
    }
    "other"
  }, character(1))
  cat_levels <- c("Prometheus_ref", "Paralog_ref", "Candidate_A", "Candidate_B",
                  "Candidate_C", "Candidate_Rejected", "other")
  cat_colors <- c(Prometheus_ref = "#000000", Paralog_ref = "#CC79A7",
                  Candidate_A = "#0072B2", Candidate_B = "#56B4E9", Candidate_C = "#E69F00",
                  Candidate_Rejected = "#D55E00", other = "#BBBBBB")
  present_levels <- intersect(cat_levels, unique(category))

  # supported internal nodes (node labels index internal nodes Ntip+1 ...)
  nl <- if (is.null(tree$node.label)) rep("", tree$Nnode) else tree$node.label
  sup <- t(vapply(nl, parse_support, numeric(2)))
  supported <- which(!is.na(sup[, "ufb"]) & sup[, "ufb"] >= 95 &
                     (is.na(sup[, "alrt"]) | sup[, "alrt"] >= 80))
  sup_nodes <- length(tips) + supported

  n_tips <- length(tips)
  show_labels <- n_tips <= 90
  circular <- n_tips > 60
  h <- if (circular) max(6, min(0.10 * n_tips, 14)) else max(4, min(0.13 * n_tips, 16))

  suppressWarnings(suppressMessages({ library(ggtree); library(ggplot2) }))
  df <- data.frame(label = tips, category = factor(category, levels = present_levels),
                   stringsAsFactors = FALSE)
  gt <- ggtree::ggtree(tree, layout = if (circular) "circular" else "rectangular",
                       linewidth = 0.25)
  gt <- gt %<+% df
  gt <- gt + ggtree::geom_tippoint(ggplot2::aes(color = category), size = 1.1) +
    ggplot2::scale_color_manual(values = cat_colors[present_levels], name = "Tip category",
                                drop = FALSE, na.value = "#BBBBBB")
  if (length(sup_nodes))
    gt <- gt + ggtree::geom_point2(ggplot2::aes(subset = node %in% sup_nodes),
                                   size = 0.5, color = "black")
  if (show_labels) gt <- gt + ggtree::geom_tiplab(size = 1.2, offset = 0.002)
  gt <- gt + ggplot2::labs(
      title = "Prometheus validation tree: candidates vs references",
      subtitle = sprintf("%s; dots = SH-aLRT >= 80 & UFBoot >= 95 (%d nodes)",
                         root_note, length(sup_nodes))) +
    ggplot2::theme(legend.position = "right",
                   legend.key.size = ggplot2::unit(0.6, "lines"),
                   legend.text = ggplot2::element_text(size = 6),
                   legend.title = ggplot2::element_text(size = 7),
                   plot.title = ggplot2::element_text(size = 8, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 6))
  save_plot(gt, opt$plot_prefix, width = if (circular) h else 7.5, height = h)
  message(sprintf("plot_validation_tree: %d tips (%d categories, %d supported nodes; %s) -> %s.{png,svg}",
                  n_tips, length(present_levels), length(sup_nodes), root_note, opt$plot_prefix))
  invisible(NULL)
}

status <- tryCatch({ main(); 0L }, error = function(e) {
  message("plot_validation_tree: ERROR: ", conditionMessage(e)); 1L })
quit(status = status)
