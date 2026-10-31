#!/usr/bin/env Rscript
# =====================================================================
# Stage 5 figure -- Prometheus gene tree vs archaeal species tree (tanglegram).
# Each gene-tree tip (genome__gene) is linked to its host genome on the species tree;
# links are coloured by host status (methanogen / ANME / other). Shows whether the gene
# history tracks the species history (vertical) or crosses it (transfer).
#
# * The species tree is pruned to the genomes that carry a gene in the gene tree (the
#   carriers are what the comparison is about; thousands of empty tips are unreadable).
# * The IQ-TREE gene tree is unrooted; it is midpoint-rooted FOR DISPLAY ONLY (stated in
#   the title). GeneRax infers its own root for the reconciliation.
# * phytools::cophylo(rotate = TRUE) minimises crossings; link colours are taken from the
#   association matrix cophylo returns, so they stay matched to the links.
# Genuinely empty input (trivial tree, < 2 links) -> labelled placeholder; any other
# error exits non-zero.
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
.theme <- Filter(file.exists, c(file.path(dirname(.script_path()), "plotting_theme.R"),
                                file.path(getwd(), "plotting_theme.R")))
if (!length(.theme)) stop("plotting_theme.R not found next to plot_tanglegram.R")
source(.theme[1])

# other = non-methanogen host (no mcrA, or MCR-like alkane oxidiser); unknown = McrA host
# whose placement is ambiguous/unplaced, or a genome missing from the master table
COLS <- c(methanogen = "#009E73", ANME = "#CC79A7", other = "#E69F00", unknown = "#999999")

is_trivial_tree <- function(path) {
  txt <- trimws(paste(readLines(path, warn = FALSE), collapse = ""))
  !nzchar(txt) || grepl("^\\(\\s*\\)\\s*;?$", txt)
}

flag <- function(x) tolower(trimws(as.character(x))) %in% c("true", "1", "yes", "t")

host_class <- function(genomes, master_path) {
  out <- setNames(rep("unknown", length(genomes)), genomes)
  if (is.null(master_path) || !nzchar(master_path)) return(out)
  m <- read.delim(master_path, sep = "\t", colClasses = "character", check.names = FALSE)
  if (!all(c("genome_id", "methanogen") %in% names(m)))
    stop("master table lacks genome_id/methanogen")
  anme <- if ("is_anme" %in% names(m)) flag(m$is_anme) else rep(FALSE, nrow(m))
  grp <- if ("mcr_group" %in% names(m)) tolower(trimws(m$mcr_group)) else rep("", nrow(m))
  anme <- anme | grp == "anme"
  cls <- ifelse(flag(m$methanogen), "methanogen", ifelse(anme, "ANME",
         ifelse(grp %in% c("ambiguous", "unplaced", "not_placed"), "unknown", "other")))
  cls <- setNames(cls, trimws(m$genome_id))
  hit <- genomes %in% names(cls)
  out[hit] <- cls[genomes[hit]]
  out
}

main <- function() {
  option_list <- list(
    make_option("--gene-tree", type = "character", dest = "gene_tree"),
    make_option("--species-tree", type = "character", dest = "species_tree"),
    make_option("--master", type = "character", default = ""),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  for (k in c("gene_tree", "species_tree", "plot_prefix"))
    if (is.null(opt[[k]])) stop("--", gsub("_", "-", k), " is required")

  if (is_trivial_tree(opt$gene_tree) || is_trivial_tree(opt$species_tree)) {
    placeholder_plot(opt$plot_prefix, "Tanglegram: empty gene or species tree")
    message("plot_tanglegram: trivial tree(s); placeholder")
    return(invisible(NULL))
  }
  gtree <- read.tree(opt$gene_tree)
  stree <- read.tree(opt$species_tree)
  if (!is.rooted(stree)) stop("species tree is not rooted (scaffolds/root_tree.py must root it)")
  gtree$node.label <- NULL; stree$node.label <- NULL
  g_geno <- sub("__.*$", "", gtree$tip.label)
  on_sp <- g_geno %in% stree$tip.label
  if (sum(on_sp) < 2) {
    placeholder_plot(opt$plot_prefix, sprintf("Tanglegram: %d gene tip(s) map to the species tree", sum(on_sp)))
    message("plot_tanglegram: < 2 links; placeholder")
    return(invisible(NULL))
  }
  n_unmapped <- sum(!on_sp)
  if (n_unmapped) gtree <- drop.tip(gtree, gtree$tip.label[!on_sp])
  g_geno <- sub("__.*$", "", gtree$tip.label)
  carriers <- unique(g_geno)
  if (length(carriers) < 2) {
    placeholder_plot(opt$plot_prefix, "Tanglegram: genes come from < 2 genomes")
    return(invisible(NULL))
  }
  stree <- keep.tip(stree, carriers)
  if (Ntip(gtree) >= 3) gtree <- midpoint.root(gtree)
  if (!is.null(gtree$edge.length)) gtree$edge.length[gtree$edge.length <= 0] <- 1e-6
  assoc <- cbind(g_geno, gtree$tip.label)       # (species tip, gene tip)
  cls <- host_class(carriers, opt$master)

  obj <- cophylo(stree, gtree, assoc = assoc, rotate = TRUE, print = FALSE)
  a <- obj$assoc
  link_col <- COLS[cls[a[, 1]]]
  n_tips <- max(Ntip(gtree), Ntip(stree))
  show_labels <- n_tips <= 70
  draw <- function() {
    # outer margins hold the title (top) and the legend + omitted-gene note (bottom), so
    # neither is clipped nor drawn over the trees (cophylo uses the whole plot region)
    op <- graphics::par(oma = c(1.8, 0, 1.6, 0))
    on.exit(graphics::par(op))
    plot(obj, link.col = unname(link_col), link.lwd = 1.1, link.type = "curved",
         fsize = if (show_labels) 0.4 else 1e-6, pts = FALSE)
    graphics::mtext(sprintf(
      "Species tree (%d carrier genomes)  vs  Prometheus gene tree (%d genes, midpoint-rooted for display)",
      Ntip(stree), Ntip(gtree)), side = 3, outer = TRUE, line = 0.4, cex = 0.6, font = 2)
    present <- names(COLS)[names(COLS) %in% unique(cls)]
    graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    graphics::plot.new()
    graphics::legend("bottomleft", legend = paste("host:", present), col = COLS[present], lwd = 2,
                     bty = "n", cex = 0.5, horiz = TRUE, inset = c(0.005, 0.002))
    if (n_unmapped)
      graphics::legend("bottomright", legend = sprintf(
        "%d gene(s) whose genome is not on the species tree omitted", n_unmapped),
        bty = "n", cex = 0.5, inset = c(0.005, 0.002))
  }
  save_plot(draw, opt$plot_prefix, width = 7.5, height = max(3.5, min(0.16 * n_tips, 13)))
  message(sprintf("plot_tanglegram: %d links, %d carrier genomes, %d unmapped genes -> %s.{png,svg}",
                  nrow(a), Ntip(stree), n_unmapped, opt$plot_prefix))
}

if (sys.nframe() == 0L) main()
