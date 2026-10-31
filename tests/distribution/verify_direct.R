#!/usr/bin/env Rscript
# Independent reference computations for run_tests.py (does NOT source the pipeline
# scripts): rebuilds the analysis set with its own implementation of the Contract v2
# item 6 tier rule and calls the underlying tools directly.
#   Rscript verify_direct.R assoc <master> <tree> <tiers> <seed> <d_perm> <out.json>
#   Rscript verify_direct.R cooc  <master> <tree> <tiers> <matrix> <systems,csv> <out.json>
suppressPackageStartupMessages({
  library(ape); library(phylolm); library(phytools); library(caper); library(jsonlite)
})
a <- commandArgs(trailingOnly = TRUE)
mode <- a[1]; master_p <- a[2]; tree_p <- a[3]; tiers <- strsplit(a[4], ",")[[1]]
m <- read.delim(master_p, colClasses = "character", na.strings = character(0), quote = "")
tr <- read.tree(text = sub("^\\[&R\\]\\s*", "", paste(readLines(tree_p), collapse = "")))
tr$node.label <- NULL
tr$tip.label <- gsub("'", "", tr$tip.label)
tier <- trimws(m$prometheus_tier)
y <- ifelse(tier %in% tiers, 1L, ifelse(tier == "", 0L, NA_integer_))
cand <- m$pass_primary == "True" & m$tree_tip != "" & m$tree_tip %in% tr$tip.label
use <- cand & !is.na(y)
d <- data.frame(tip = m$tree_tip[use], present = y[use],
                methanogen = as.integer(m$methanogen[use] == "True"), stringsAsFactors = FALSE)
t2 <- keep.tip(tr, d$tip)
d <- d[match(t2$tip.label, d$tip), ]
rownames(d) <- d$tip
res <- list(n = nrow(d), n_present = sum(d$present), n_excluded = sum(cand & is.na(y)))

if (mode == "assoc") {
  seed <- as.integer(a[5]); dperm <- as.integer(a[6]); out <- a[7]
  g <- glm(present ~ methanogen, family = binomial, data = d)
  res$glm_est <- unname(coef(g)["methanogen"])
  res$glm_lrt_p <- pchisq(g$null.deviance - g$deviance, 1, lower.tail = FALSE)
  ft <- fisher.test(table(factor(d$methanogen, 0:1), factor(d$present, 0:1)))
  res$fisher_logor <- log(unname(ft$estimate)); res$fisher_p <- ft$p.value
  fp <- fitPagel(t2, setNames(factor(d$methanogen, 0:1), d$tip), setNames(factor(d$present, 0:1), d$tip))
  res$pagel_lr <- as.numeric(fp$lik.ratio); res$pagel_p <- as.numeric(fp$P)
  res$pagel_indep_logL <- as.numeric(fp$independent.logL)
  sx <- fitMk(t2, setNames(factor(d$methanogen, 0:1), d$tip), model = "ARD")
  sy <- fitMk(t2, setNames(factor(d$present, 0:1), d$tip), model = "ARD")
  res$sep_logL_sum <- as.numeric(sx$logLik) + as.numeric(sy$logLik)
  f0 <- suppressWarnings(phyloglm(present ~ methanogen, data = d, phy = t2, method = "logistic_MPLE"))
  res$pglm_default_est <- unname(coef(f0)["methanogen"]); res$pglm_default_pen <- f0$penlogLik
  res$pglm_default_conv <- f0$convergence
  cd <- comparative.data(t2, data.frame(tip = d$tip, present = d$present), tip)
  set.seed(seed)
  pd <- phylo.d(cd, binvar = present, permut = dperm)
  res$D <- as.numeric(pd$DEstimate)
  res$D_p_random <- (sum(pd$Permutations$random <= pd$Parameters$Observed) + 1) / (dperm + 1)
  writeLines(toJSON(res, auto_unbox = TRUE, digits = NA), out)
} else if (mode == "cooc") {
  mat <- read.delim(a[5], colClasses = "character", quote = "")
  systems <- strsplit(a[6], ",")[[1]]; out <- a[7]
  mat <- mat[match(m$genome_id[use], mat$genome_id), ]
  mat <- mat[match(t2$tip.label, m$tree_tip[use]), ]
  res$systems <- list()
  for (s in systems) {
    sv <- setNames(factor(as.integer(mat[[s]]), 0:1), t2$tip.label)
    pv <- setNames(factor(d$present, 0:1), d$tip)
    fp <- fitPagel(t2, pv, sv, lik.func = "pruning")
    ft <- fisher.test(table(factor(d$present, 0:1), factor(as.integer(mat[[s]]), 0:1)))
    res$systems[[s]] <- list(lr = as.numeric(fp$lik.ratio), p = as.numeric(fp$P),
                             n_b = sum(as.integer(mat[[s]])), or = unname(ft$estimate))
  }
  writeLines(toJSON(res, auto_unbox = TRUE, digits = NA), out)
}
