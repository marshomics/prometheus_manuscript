#!/usr/bin/env Rscript
# Synthetic Stage-4 inputs with a KNOWN generative signal (used by run_tests.py).
#   Rscript sim_data.R <outdir> <scenario> <n_tips> <seed>
# scenario = signal  : methanogen has many independent origins (fast Mk); Prometheus
#                      presence ~ phylogenetic logistic model (phylolm::rbinTrait) with
#                      beta = (-1, 1.5) on (intercept, methanogen), alpha = 3/T.
#            null    : same, but beta_methanogen = 0.
#            clumped : presence = one large clade (strong phylogenetic signal, D << 1).
# Writes tree.nwk (rooted, '[&R] ' prefix, IQ-TREE-style duplicated quoted '100/100'
# support labels), master.tsv (SCHEMA master_table columns), matrix.tsv (defense
# matrix with a co-occurring, an excluded, a random and a rare system).
# Genome ids are zero-padded digits ('000001') so numeric-looking ids are exercised.
# Tier design (Contract v2 item 6): a few present genomes are B-only or C-only, three
# genomes fail pass_primary and two are not on the tree.
suppressPackageStartupMessages({ library(ape); library(phylolm); library(phytools) })
a <- commandArgs(trailingOnly = TRUE)
outdir <- a[1]; scen <- a[2]; n <- as.integer(a[3]); seed <- as.integer(a[4])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)
ids <- sprintf("%06d", seq_len(n))
tr <- rcoal(n, tip.label = ids)
tr$edge.length <- tr$edge.length / max(node.depth.edgelength(tr))   # height 1
Q <- matrix(c(-3, 3, 3, -3), 2, 2, dimnames = list(c("0", "1"), c("0", "1")))
meth <- as.integer(as.character(sim.Mk(tr, Q)[ids]))
names(meth) <- ids
if (scen == "clumped") {
  big <- extract.clade(tr, n + 2)$tip.label            # one daughter clade of the root
  if (length(big) < n / 4 || length(big) > 3 * n / 4)
    big <- extract.clade(tr, n + 3)$tip.label
  y <- setNames(as.integer(ids %in% big), ids)
} else {
  b1 <- if (scen == "signal") 1.5 else 0
  y <- rbinTrait(n = 1, phy = tr, beta = c(-1, b1), alpha = 3, X = cbind(1, meth))
  y <- setNames(as.integer(y), ids)
}
tier <- ifelse(y == 1, "A", "")
pres_idx <- which(y == 1)
bset <- pres_idx[seq(2, length(pres_idx), by = 7)]          # B-only genomes
cset <- setdiff(pres_idx[seq(5, length(pres_idx), by = 11)], bset)  # C-only genomes
tier[bset] <- "B"; tier[cset] <- "C"
pass <- rep("True", n); pass[c(3, 17, 29)] <- "False"
tip <- ids; tip[c(8, 41)] <- ""
gs <- round(exp(rnorm(n, log(2.5e6), 0.3)))
comp <- round(runif(n, 70, 100), 2)
dex <- rpois(n, 5)
m <- data.frame(genome_id = ids, representative = "True", tree_tip = tip,
  gtdb_domain = "Archaea",
  gtdb_phylum = ifelse(seq_len(n) %% 3 == 0, "Halobacteriota", "Methanobacteriota"),
  gtdb_class = "c__x", gtdb_order = "o__x",
  gtdb_family = paste0("Fam", (seq_len(n) %% 4)), gtdb_genus = paste0("Gen", seq_len(n) %% 9),
  gtdb_species = paste0("sp", seq_len(n)), completeness = comp, contamination = 1,
  quality_tier = "A", pass_primary = pass, genome_size = gs,
  mcrA = ifelse(meth == 1, "True", "False"), n_mcr = meth,
  methanogen = ifelse(meth == 1, "True", "False"), is_anme = "False",
  methanogenesis_subtype = ifelse(meth == 1, "unresolved", "NA"),
  prometheus_present = ifelse(tier != "", "True", "False"), prometheus_tier = tier,
  prometheus_n_loci = as.integer(tier != ""), detect_defensefinder = as.integer(tier != ""),
  n_defense_systems = dex + as.integer(tier != ""), n_defense_systems_excl_prometheus = dex,
  stringsAsFactors = FALSE)
write.table(m, file.path(outdir, "master.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# defense matrix: raw DefenseFinder-style 0/1 (Prometheus column = any validated tier)
prom_df <- as.integer(tier != "")
mat <- data.frame(genome_id = ids, Prometheus = prom_df,
  CoSys = ifelse(y == 1, 1L, rbinom(n, 1, 0.08)),
  ExSys = ifelse(y == 1, rbinom(n, 1, 0.03), rbinom(n, 1, 0.6)),
  RandSys = rbinom(n, 1, 0.4),
  RareSys = as.integer(seq_len(n) %in% c(1, 2)), stringsAsFactors = FALSE)
write.table(mat, file.path(outdir, "matrix.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
nl <- tr; nl$node.label <- rep("100/100", nl$Nnode)
txt <- write.tree(nl)
txt <- gsub(")100/100", ")'100/100'", txt, fixed = TRUE)
writeLines(paste0("[&R] ", txt), file.path(outdir, "tree.nwk"))
cat(sprintf("scenario=%s n=%d present=%d methanogen=%d B=%d C=%d\n", scen, n, sum(y), sum(meth),
            length(bset), length(cset)))
