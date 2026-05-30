#!/usr/bin/env Rscript
# =====================================================================
# Stage 5.5 -- ancestral gain/loss of Prometheus across the archaeal
# species tree by stochastic character mapping (design Stage 5.5;
# confounder #5/#9). This is the EVOLUTIONARY form of the association
# question and the one least confounded by tip-level clade structure.
#
# phytools::make.simmap(tree, x, model='ARD', nsim) draws nsim stochastic
# histories under ML-fitted asymmetric (all-rates-different) gain/loss
# rates. Across the simulations we summarise:
#   n_gains  (0->1 transitions, mean over sims)  = n_independent_acquisitions
#   n_losses (1->0 transitions, mean over sims)
#   rate_gain, rate_loss (from the fitted Q matrix)
#
# Methanogen enrichment of gains: for each mapped history we count gains
# that occur on branches leading into methanogen tips/clades (a branch is
# "methanogen" if its child subtree is entirely methanogen). We compare
# the observed mean gains-on-methanogen-branches to a null built by
# permuting the methanogen labels across tips (>=200 permutations,
# re-mapping each time) and report a one-sided enrichment p.
#
# Output (SCHEMA evolution/ancestral/gains_losses.tsv): rows
#   n_gains, n_losses, rate_gain, rate_loss, n_independent_acquisitions,
#   gains_on_methanogen_branches, expected_gains_methanogen, enrichment_p
# columns: metric, value, note.
# PLOT: {plot-prefix}.png/.svg -- densityMap (or simmap node pies) titled
# with n_independent_acquisitions. Trivial input -> header-only TSV +
# placeholder figure, quit(status=0).
# =====================================================================

suppressWarnings(suppressMessages({ library(optparse) }))

OUT_COLS <- c("metric", "value", "note")

.script_path <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) return(normalizePath(sub("^--file=", "", m[1])))
  ""
}
.this_dir <- tryCatch(dirname(.script_path()), error = function(e) ".")
if (is.na(.this_dir) || !nzchar(.this_dir)) .this_dir <- "."
.theme_ok <- tryCatch({
  suppressWarnings(suppressMessages(source(file.path(.this_dir, "plotting_theme.R"))))
  TRUE
}, error = function(e) FALSE)

# Fallback placeholder if the shared helper failed to load.
.placeholder <- function(prefix, msg) {
  if (.theme_ok && exists("placeholder_plot")) {
    tryCatch(placeholder_plot(prefix, msg), error = function(e) NULL)
  } else {
    tryCatch({
      grDevices::png(paste0(prefix, ".png"), width = 4, height = 2.2,
                     units = "in", res = 300)
      graphics::plot.new(); graphics::text(0.5, 0.5, msg); grDevices::dev.off()
    }, error = function(e) NULL)
    tryCatch({
      grDevices::svg(paste0(prefix, ".svg"), width = 4, height = 2.2)
      graphics::plot.new(); graphics::text(0.5, 0.5, msg); grDevices::dev.off()
    }, error = function(e) NULL)
  }
}

write_header_only <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  con <- file(path, "w"); writeLines(paste(OUT_COLS, collapse = "\t"), con); close(con)
}

as_logical_col <- function(x) {
  s <- tolower(trimws(as.character(x))); s %in% c("true", "1", "yes", "t")
}

mrow <- function(metric, value, note = "") {
  data.frame(metric = metric, value = value, note = note, stringsAsFactors = FALSE)
}

# Count 0->1 and 1->0 changes in one simmap tree from its $maps.
# Returns list(gains, losses) and per-edge child-state (the state at the
# tipward end of each edge) for the methanogen-branch test.
count_changes <- function(smap) {
  gains <- 0L; losses <- 0L
  child_state <- character(nrow(smap$edge))
  for (i in seq_along(smap$maps)) {
    seg <- smap$maps[[i]]
    st <- names(seg)
    child_state[i] <- st[length(st)]
    if (length(st) >= 2) {
      for (j in 2:length(st)) {
        a <- st[j - 1]; b <- st[j]
        if (a == "0" && b == "1") gains <- gains + 1L
        if (a == "1" && b == "0") losses <- losses + 1L
      }
    }
  }
  list(gains = gains, losses = losses, child_state = child_state)
}

# For each edge, is the descendant subtree entirely methanogen?
# meth is a named 0/1 vector over tips.
edge_into_methanogen <- function(tree, meth) {
  ntip <- length(tree$tip.label)
  # descendant tips per node via ape::prop.part-like; use phangorn-free approach
  # Build children lists.
  desc_tips <- vector("list", max(tree$edge))
  # postorder: compute tip sets
  # Initialise tips.
  for (t in seq_len(ntip)) desc_tips[[t]] <- tree$tip.label[t]
  # Children-before-parents ordering so each parent's tip set is complete
  # by the time we reach it. reorder.phylo shares the tree's edge indexing,
  # which is also the indexing of make.simmap's $maps (same tree object).
  ord <- ape::reorder.phylo(tree, "postorder")
  for (e in seq_len(nrow(ord$edge))) {
    parent <- ord$edge[e, 1]; child <- ord$edge[e, 2]
    if (is.null(desc_tips[[parent]])) desc_tips[[parent]] <- character(0)
    if (child <= ntip) {
      desc_tips[[parent]] <- c(desc_tips[[parent]], ord$tip.label[child])
    } else {
      desc_tips[[parent]] <- c(desc_tips[[parent]], desc_tips[[child]])
    }
  }
  # Now evaluate per edge of the ORIGINAL tree ordering.
  res <- logical(nrow(tree$edge))
  for (e in seq_len(nrow(tree$edge))) {
    child <- tree$edge[e, 2]
    tips_c <- if (child <= ntip) tree$tip.label[child] else desc_tips[[child]]
    if (length(tips_c) == 0) { res[e] <- FALSE; next }
    mv <- meth[tips_c]
    res[e] <- all(!is.na(mv)) && all(mv == 1)
  }
  res
}

# Count gains on methanogen-leading branches for one simmap.
gains_on_meth <- function(smap, meth_edges) {
  cnt <- 0L
  for (i in seq_along(smap$maps)) {
    if (!meth_edges[i]) next
    st <- names(smap$maps[[i]])
    if (length(st) >= 2) for (j in 2:length(st))
      if (st[j - 1] == "0" && st[j] == "1") cnt <- cnt + 1L
  }
  cnt
}

main <- function() {
  option_list <- list(
    make_option("--tree", type = "character"),
    make_option("--master", type = "character"),
    make_option("--nsim", type = "integer", default = 1000),
    make_option("--seed", type = "integer", default = 12345),
    make_option("--out", type = "character"),
    make_option("--plot-prefix", type = "character", dest = "plot_prefix")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  set.seed(opt$seed)

  tryCatch({
    stopifnot(!is.null(opt$tree), !is.null(opt$master), !is.null(opt$out),
              !is.null(opt$plot_prefix))

    suppressWarnings(suppressMessages({
      library(ape); library(readr); library(dplyr)
    }))

    tree_txt <- tryCatch(paste(readLines(opt$tree, warn = FALSE), collapse = ""),
                         error = function(e) "")
    if (!nzchar(trimws(tree_txt)) || grepl("^\\(\\s*\\)\\s*;?$", trimws(tree_txt))) {
      write_header_only(opt$out)
      .placeholder(opt$plot_prefix, "Ancestral states: empty tree")
      message("ancestral_states: trivial tree; header-only + placeholder.")
      return(invisible(NULL))
    }
    tree <- tryCatch(read.tree(opt$tree), error = function(e) NULL)
    # Safety net: make.simmap requires a rooted tree. The pipeline roots
    # species_tree.treefile upstream; if an unrooted tree is ever passed in, apply a
    # midpoint root rather than mapping gains/losses on an arbitrary root.
    if (!is.null(tree) && length(tree$tip.label) >= 3 && !ape::is.rooted(tree))
      tree <- tryCatch(phangorn::midpoint(tree), error = function(e) tree)
    if (is.null(tree) || length(tree$tip.label) < 5) {
      write_header_only(opt$out)
      .placeholder(opt$plot_prefix, "Ancestral states: <5 tips")
      message("ancestral_states: <5 tips; header-only + placeholder.")
      return(invisible(NULL))
    }

    master <- readr::read_tsv(opt$master, show_col_types = FALSE, progress = FALSE)
    if (nrow(master) == 0 || !"genome_id" %in% names(master)) {
      write_header_only(opt$out)
      .placeholder(opt$plot_prefix, "Ancestral states: empty master")
      return(invisible(NULL))
    }

    # tip-keyed presence (0/1). tree_tip == genome_id for the species tree.
    tt <- if ("tree_tip" %in% names(master)) trimws(as.character(master$tree_tip)) else ""
    gid <- trimws(as.character(master$genome_id))
    use_tip <- ifelse(nzchar(tt), tt, gid)
    present <- as.integer(as_logical_col(master$prometheus_present))
    meth <- as.integer(as_logical_col(master$methanogen))
    names(present) <- use_tip; names(meth) <- use_tip

    keep <- intersect(tree$tip.label, names(present))
    keep <- keep[!duplicated(keep)]
    if (length(keep) < 5) {
      write_header_only(opt$out)
      .placeholder(opt$plot_prefix, "Ancestral states: <5 matched tips")
      return(invisible(NULL))
    }
    tree <- ape::keep.tip(tree, keep)
    x <- present[tree$tip.label]
    mvec <- meth[tree$tip.label]

    if (length(unique(x)) < 2) {
      write_header_only(opt$out)
      .placeholder(opt$plot_prefix, "Ancestral states: presence invariant")
      message("ancestral_states: no presence variation; header-only + placeholder.")
      return(invisible(NULL))
    }
    if (!requireNamespace("phytools", quietly = TRUE)) {
      write_header_only(opt$out)
      .placeholder(opt$plot_prefix, "Ancestral states: phytools unavailable")
      message("ancestral_states: phytools missing; header-only + placeholder.")
      return(invisible(NULL))
    }

    # Ensure branch lengths for simmap; if absent, set unit lengths.
    if (is.null(tree$edge.length) || all(is.na(tree$edge.length))) {
      tree <- ape::compute.brlen(tree, 1)
    }

    nsim <- max(1L, as.integer(opt$nsim))
    xf <- setNames(as.factor(x), tree$tip.label)

    # ---- stochastic character mapping (ARD) ------------------------ #
    smaps <- phytools::make.simmap(tree, xf, model = "ARD", nsim = nsim,
                                   message = FALSE)
    # Normalise to a multiSimmap list regardless of nsim. nsim==1 returns a
    # single "simmap"; wrap it. Keep BOTH classes so densityMap() (which
    # dispatches on multiPhylo) still accepts the object.
    if (inherits(smaps, "simmap") && !inherits(smaps, "multiSimmap")) {
      smaps <- structure(list(smaps), class = c("multiSimmap", "multiPhylo"))
    } else if (!inherits(smaps, "multiPhylo")) {
      class(smaps) <- c("multiSimmap", "multiPhylo")
    }

    # gains/losses per sim
    gl <- lapply(smaps, count_changes)
    gains_vec  <- vapply(gl, function(z) z$gains, numeric(1))
    losses_vec <- vapply(gl, function(z) z$losses, numeric(1))
    mean_gains  <- mean(gains_vec)
    mean_losses <- mean(losses_vec)

    # rates from fitted Q (average over sims' attached Q if present).
    rate_gain <- NA_real_; rate_loss <- NA_real_
    Qs <- lapply(smaps, function(s) attr(s, "Q"))
    Qs <- Qs[!vapply(Qs, is.null, logical(1))]
    if (length(Qs)) {
      g <- vapply(Qs, function(Q) {
        if (all(c("0", "1") %in% rownames(Q))) Q["0", "1"] else NA_real_
      }, numeric(1))
      l <- vapply(Qs, function(Q) {
        if (all(c("0", "1") %in% rownames(Q))) Q["1", "0"] else NA_real_
      }, numeric(1))
      rate_gain <- mean(g, na.rm = TRUE); rate_loss <- mean(l, na.rm = TRUE)
    }

    # ---- methanogen enrichment of gains --------------------------- #
    enrichment_p <- NA_real_
    obs_meth_gains <- NA_real_; exp_meth_gains <- NA_real_
    if (length(unique(mvec)) >= 2 && !anyNA(mvec)) {
      meth_named <- setNames(as.integer(mvec), tree$tip.label)
      meth_edges <- tryCatch(edge_into_methanogen(tree, meth_named),
                             error = function(e) rep(FALSE, nrow(tree$edge)))
      obs_per_sim <- vapply(smaps, function(s) gains_on_meth(s, meth_edges), numeric(1))
      obs_meth_gains <- mean(obs_per_sim)

      # Null: permute methanogen labels across tips, recompute meth_edges,
      # and score the SAME mapped histories (so only the label assignment
      # is randomised, holding the inferred gain placements fixed). >=200 perms.
      nperm <- max(200L, nsim)
      null_means <- numeric(nperm)
      for (p in seq_len(nperm)) {
        permv <- sample(meth_named)
        names(permv) <- tree$tip.label
        me <- tryCatch(edge_into_methanogen(tree, permv),
                       error = function(e) rep(FALSE, nrow(tree$edge)))
        # average over a small subset of sims for speed if nsim large
        idx <- if (length(smaps) > 50) sample(seq_along(smaps), 50) else seq_along(smaps)
        null_means[p] <- mean(vapply(smaps[idx], function(s) gains_on_meth(s, me), numeric(1)))
      }
      exp_meth_gains <- mean(null_means)
      # one-sided: P(null >= observed)
      enrichment_p <- (sum(null_means >= obs_meth_gains) + 1) / (nperm + 1)
    }

    n_indep <- mean_gains  # independent acquisitions = mean inferred gains

    out <- rbind(
      mrow("n_gains", round(mean_gains, 4),
           sprintf("mean 0->1 over %d simmaps (ARD)", nsim)),
      mrow("n_losses", round(mean_losses, 4),
           sprintf("mean 1->0 over %d simmaps (ARD)", nsim)),
      mrow("rate_gain", ifelse(is.na(rate_gain), NA_real_, round(rate_gain, 6)),
           "Q[0,1] gain rate (mean over sims)"),
      mrow("rate_loss", ifelse(is.na(rate_loss), NA_real_, round(rate_loss, 6)),
           "Q[1,0] loss rate (mean over sims)"),
      mrow("n_independent_acquisitions", round(n_indep, 4),
           "mean independent gains across the archaeal tree"),
      mrow("gains_on_methanogen_branches",
           ifelse(is.na(obs_meth_gains), NA_real_, round(obs_meth_gains, 4)),
           "observed mean gains on branches leading into methanogen clades"),
      mrow("expected_gains_methanogen",
           ifelse(is.na(exp_meth_gains), NA_real_, round(exp_meth_gains, 4)),
           "null mean from permuting methanogen labels across tips"),
      mrow("enrichment_p",
           ifelse(is.na(enrichment_p), NA_real_, signif(enrichment_p, 4)),
           "one-sided P(null >= observed); >=200 label permutations")
    )

    dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
    readr::write_tsv(out, opt$out)

    # ---- figure: densityMap of presence, titled by n_indep -------- #
    title <- sprintf("Prometheus gains/losses (n independent acquisitions = %.2f)",
                     n_indep)
    drawn <- FALSE
    fig <- function() {
      dm <- tryCatch(
        phytools::densityMap(smaps, states = c("0", "1"), plot = FALSE,
                             res = 200),
        error = function(e) NULL)
      if (!is.null(dm)) {
        dm <- tryCatch(phytools::setMap(dm, c("#D9D9D9", "#0072B2")),
                       error = function(e) dm)
        plot(dm, lwd = 2, ftype = if (length(tree$tip.label) <= 60) "i" else "off",
             fsize = 0.5, outline = FALSE, legend = FALSE)
      } else {
        # fallback: plot one simmap with node states
        phytools::plotSimmap(smaps[[1]], colors = c("0" = "#D9D9D9", "1" = "#0072B2"),
                             ftype = if (length(tree$tip.label) <= 60) "i" else "off",
                             fsize = 0.5, lwd = 2)
      }
      graphics::title(main = title, cex.main = 0.8)
    }
    if (.theme_ok && exists("save_plot")) {
      tryCatch({
        h <- max(3, min(0.12 * length(tree$tip.label), 11))
        save_plot(fig, opt$plot_prefix, width = 5.5, height = h)
        drawn <- TRUE
      }, error = function(e) message("ancestral_states: plot error: ", conditionMessage(e)))
    }
    if (!drawn) .placeholder(opt$plot_prefix, title)

    message(sprintf("ancestral_states: n_indep=%.2f gains=%.2f losses=%.2f enrich_p=%s -> %s",
                    n_indep, mean_gains, mean_losses,
                    ifelse(is.na(enrichment_p), "NA", signif(enrichment_p, 3)), opt$out))
    invisible(NULL)

  }, error = function(e) {
    message("ancestral_states: top-level error: ", conditionMessage(e))
    tryCatch(write_header_only(opt$out), error = function(e2) invisible(NULL))
    tryCatch(.placeholder(opt$plot_prefix, "Ancestral states: error"),
             error = function(e2) invisible(NULL))
    invisible(NULL)
  })

  quit(status = 0)
}

main()
