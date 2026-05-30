#!/usr/bin/env Rscript
# =====================================================================
# Stage 4.6 -- co-occurrence / mutual exclusivity of Prometheus with the
# other DefenseFinder systems, with a PHYLOGENY-AWARE null (design Stage
# 4.6 / confounder #5). The naive cooccur test treats genomes as
# independent; instead we test correlated EVOLUTION of each system pair
# across the species tree with phytools::fitPagel (binary dependent vs
# independent model), which asks whether the two systems' gains and losses
# are coupled along branches rather than merely co-present at the tips.
#
# A Fisher exact odds ratio gives the direction (OR>1 co-occurrence, OR<1
# mutual exclusion). p-values are BH-adjusted across the tested systems.
#
# Output (SCHEMA distribution/cooccurrence.tsv):
#   system_a (='Prometheus'), system_b, n_a, n_b, n_both, method,
#   statistic, p_value, p_adj, direction
# Trivial input -> header-only TSV, quit(status=0).
# =====================================================================

suppressWarnings(suppressMessages({
  library(optparse)
}))

OUT_COLS <- c("system_a", "system_b", "n_a", "n_b", "n_both",
              "method", "statistic", "p_value", "p_adj", "direction")

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

# Find the Prometheus column in the defense matrix (case-insensitive).
find_prometheus_col <- function(cols) {
  hit <- cols[tolower(cols) == "prometheus"]
  if (length(hit)) return(hit[1])
  hit <- cols[grepl("prometheus", tolower(cols))]
  if (length(hit)) return(hit[1])
  NA_character_
}

main <- function() {
  option_list <- list(
    make_option("--matrix", type = "character"),
    make_option("--master", type = "character"),
    make_option("--tree", type = "character"),
    make_option("--fdr", type = "character", default = "BH"),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--out", type = "character")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  set.seed(opt$seed)

  tryCatch({
    stopifnot(!is.null(opt$matrix), !is.null(opt$tree), !is.null(opt$out))

    suppressWarnings(suppressMessages({
      library(ape); library(readr); library(dplyr)
    }))

    # ---- tree ------------------------------------------------------- #
    tree_txt <- tryCatch(paste(readLines(opt$tree, warn = FALSE), collapse = ""),
                         error = function(e) "")
    if (!nzchar(trimws(tree_txt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(tree_txt))) {
      write_header_only(opt$out)
      message("cooccurrence: trivial/empty tree; header-only."); return(invisible(NULL))
    }
    tree <- tryCatch(read.tree(opt$tree), error = function(e) NULL)
    # Safety net: fitPagel assumes a rooted tree. The pipeline roots
    # species_tree.treefile upstream; midpoint-root here only as a fallback.
    if (!is.null(tree) && length(tree$tip.label) >= 3 && !ape::is.rooted(tree))
      tree <- tryCatch(phangorn::midpoint(tree), error = function(e) tree)
    if (is.null(tree) || length(tree$tip.label) < 5) {
      write_header_only(opt$out)
      message("cooccurrence: <5 tree tips; header-only."); return(invisible(NULL))
    }

    # ---- defense matrix -------------------------------------------- #
    mat <- readr::read_tsv(opt$matrix, show_col_types = FALSE, progress = FALSE)
    if (nrow(mat) == 0 || !"genome_id" %in% names(mat)) {
      write_header_only(opt$out)
      message("cooccurrence: empty matrix; header-only."); return(invisible(NULL))
    }
    mat$genome_id <- trimws(as.character(mat$genome_id))

    # Map genome_id -> tree tip. The species-tree tips ARE genome_ids;
    # the master's tree_tip column confirms which genomes are on the tree.
    tip_for <- setNames(mat$genome_id, mat$genome_id)
    if (!is.null(opt$master) && file.exists(opt$master)) {
      master <- tryCatch(readr::read_tsv(opt$master, show_col_types = FALSE,
                                         progress = FALSE),
                         error = function(e) NULL)
      if (!is.null(master) && all(c("genome_id", "tree_tip") %in% names(master))) {
        tt <- trimws(as.character(master$tree_tip))
        gg <- trimws(as.character(master$genome_id))
        ok <- nzchar(tt)
        tip_for <- setNames(tt[ok], gg[ok])
      }
    }
    mat$tip <- unname(tip_for[mat$genome_id])
    mat <- mat[!is.na(mat$tip) & mat$tip %in% tree$tip.label, , drop = FALSE]
    mat <- mat[!duplicated(mat$tip), , drop = FALSE]
    if (nrow(mat) < 5) {
      write_header_only(opt$out)
      message("cooccurrence: <5 genomes on tree; header-only."); return(invisible(NULL))
    }

    prom_col <- find_prometheus_col(names(mat))
    if (is.na(prom_col)) {
      write_header_only(opt$out)
      message("cooccurrence: no Prometheus column; header-only."); return(invisible(NULL))
    }

    # Prune tree to matched tips; align matrix rows to tip order.
    tree <- ape::keep.tip(tree, mat$tip)
    rownames(mat) <- mat$tip
    mat <- mat[tree$tip.label, , drop = FALSE]
    n <- nrow(mat)

    # 0/1 coercion robust to "True"/"1"/1/TRUE; any NA/blank -> 0 so a stray
    # missing value never poisons a system's counts (defense_matrix is strict
    # 0/1 per SCHEMA, but we do not assume it).
    bin <- function(v) {
      num <- suppressWarnings(as.numeric(as.character(v)))
      b <- as_logical_col(v) | (!is.na(num) & num == 1)
      as.integer(ifelse(is.na(b), FALSE, b))
    }
    prom <- bin(mat[[prom_col]])
    if (length(unique(prom)) < 2) {
      write_header_only(opt$out)
      message("cooccurrence: Prometheus invariant on tree; header-only."); return(invisible(NULL))
    }
    names(prom) <- mat$tip
    n_a <- sum(prom)

    # Candidate systems: every other column that is 0/1 and not metadata.
    meta_cols <- c("genome_id", "tip", prom_col)
    sys_cols <- setdiff(names(mat), meta_cols)
    # keep only columns that look binary
    is_binary <- vapply(sys_cols, function(c) {
      v <- suppressWarnings(as.numeric(mat[[c]]))
      all(is.na(v) | v %in% c(0, 1)) || is.logical(mat[[c]])
    }, logical(1))
    sys_cols <- sys_cols[is_binary]

    have_phytools <- requireNamespace("phytools", quietly = TRUE)

    rows <- list()
    for (s in sys_cols) {
      sysv <- bin(mat[[s]]); names(sysv) <- mat$tip
      nb <- sum(sysv)
      # prevalence window [3, n-3] (design: avoid near-constant systems)
      if (nb < 3 || nb > (n - 3)) next
      n_both <- sum(prom == 1 & sysv == 1)

      # ---- Fisher exact for DIRECTION (OR) ---- #
      tab <- table(factor(prom, levels = c(0, 1)),
                   factor(sysv, levels = c(0, 1)))
      direction <- "ambiguous"
      orval <- NA_real_
      ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
      if (!is.null(ft)) {
        orval <- as.numeric(ft$estimate)
        if (is.finite(orval)) {
          direction <- if (orval > 1) "co-occurrence" else if (orval < 1) "exclusion" else "neutral"
        }
      }

      # ---- phytools::fitPagel for correlated EVOLUTION ---- #
      stat <- NA_real_; p <- NA_real_; method <- "fitPagel"
      if (have_phytools) {
        pg <- tryCatch(
          phytools::fitPagel(tree,
                             x = setNames(as.factor(prom), names(prom)),
                             y = setNames(as.factor(sysv), names(sysv))),
          error = function(e) NULL)
        if (!is.null(pg)) {
          stat <- tryCatch(pg$lik.ratio, error = function(e) NA_real_)
          p    <- tryCatch(pg$P, error = function(e) NA_real_)
        }
      }
      # Fall back to Fisher p only if fitPagel could not run (labelled).
      if (is.na(p) && !is.null(ft)) {
        p <- as.numeric(ft$p.value); stat <- orval; method <- "fisher(fallback)"
      }
      dir_lab <- sprintf("%s (OR=%.3g)", direction,
                         ifelse(is.na(orval), NA_real_, orval))

      rows[[length(rows) + 1]] <- data.frame(
        system_a = "Prometheus", system_b = s,
        n_a = n_a, n_b = nb, n_both = n_both,
        method = method, statistic = stat, p_value = p,
        p_adj = NA_real_, direction = dir_lab, stringsAsFactors = FALSE)
    }

    if (length(rows) == 0) {
      write_header_only(opt$out)
      message("cooccurrence: no systems in prevalence window [3, n-3]; header-only.")
      return(invisible(NULL))
    }

    out <- do.call(rbind, rows)
    # BH across the tested systems.
    method <- if (toupper(opt$fdr) %in% c("BH", "FDR")) "BH" else opt$fdr
    pv <- which(!is.na(out$p_value))
    if (length(pv)) out$p_adj[pv] <- stats::p.adjust(out$p_value[pv], method = method)
    out <- out[, OUT_COLS, drop = FALSE]

    dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(out, opt$out)
    message(sprintf("cooccurrence: tested %d systems vs Prometheus (n=%d) -> %s",
                    nrow(out), n, opt$out))
    invisible(NULL)

  }, error = function(e) {
    message("cooccurrence: top-level error: ", conditionMessage(e))
    tryCatch(write_header_only(opt$out), error = function(e2) invisible(NULL))
    invisible(NULL)
  })

  quit(status = 0)
}

main()
