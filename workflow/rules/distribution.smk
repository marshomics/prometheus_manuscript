# =====================================================================
# Stage 4 — distribution and taxonomic association (phylogeny-aware).
# Prometheus presence follows SCHEMA Contract v2 item 6: for a tier set T a genome is
# present if its best validated tier is in T, absent if it has none, and EXCLUDED
# (dropped and counted) when its best tier lies outside T. Every rule logs to LOGS and
# fails on any tool error; header-only outputs arise only from genuinely empty input
# (decided and logged inside the scripts).
# =====================================================================

rule defense_matrix:
    """Genome x DefenseFinder system-type presence/absence (Antidefense rows excluded)."""
    input:
        master=rules.build_master_table.output.master,
        df_systems=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_systems.tsv", genome=GENOMES),
    output:
        matrix=f"{RESULTS}/distribution/defense_matrix.tsv",
    params:
        df_dir=f"{RESULTS}/defense/defensefinder",
    log:
        f"{LOGS}/distribution/defense_matrix.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        mkdir -p $(dirname {output.matrix})
        python workflow/scripts/py/build_defense_matrix.py \
            --master {input.master} --df-dir {params.df_dir} --require-all \
            --out {output.matrix} > {log} 2>&1
        """


rule compute_prevalence:
    """Prevalence by GTDB rank (Wilson 95% CI) at the primary and the A+B sensitivity tier set."""
    input:
        master=rules.build_master_table.output.master,
    output:
        prevalence=f"{RESULTS}/distribution/prevalence_by_taxon.tsv",
    params:
        tiers=",".join(cfg("stats/confidence_tiers_primary", ["A"])),
        sens_tiers=",".join(cfg("stats/confidence_tiers_sensitivity", ["A", "B"])),
    log:
        f"{LOGS}/distribution/compute_prevalence.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        mkdir -p $(dirname {output.prevalence})
        python workflow/scripts/py/compute_prevalence.py \
            --master {input.master} --primary-tiers "{params.tiers}" \
            --sensitivity-tiers "{params.sens_tiers}" --out {output.prevalence} > {log} 2>&1
        """


rule association_tests:
    """Methanogen <-> Prometheus: naive baselines, multi-start phyloglm (unadjusted and
    confounder-adjusted HEADLINE), Pagel correlated evolution, Fritz-Purvis D."""
    input:
        master=rules.build_master_table.output.master,
        tree=rules.species_tree.output.tree,
    output:
        tests=f"{RESULTS}/distribution/association_tests.tsv",
    params:
        primary_tiers=",".join(cfg("stats/confidence_tiers_primary", ["A"])),
        fdr=cfg("stats/fdr_method", "BH"),
        btol=cfg("stats/phyloglm_btol", 10),
        boot=cfg("stats/phyloglm_boot", 200),
        d_perm=cfg("stats/phylo_d_permutations", 1000),
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/distribution/association_tests.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=240,
    shell:
        r"""
        mkdir -p $(dirname {output.tests})
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
        Rscript workflow/scripts/R/association_tests.R \
            --master {input.master} --tree {input.tree} \
            --primary-tiers "{params.primary_tiers}" --fdr {params.fdr} \
            --btol {params.btol} --boot {params.boot} --d-permutations {params.d_perm} \
            --seed {params.seed} --threads {threads} \
            --out {output.tests} > {log} 2>&1
        """


rule cooccurrence:
    """Prometheus vs every other system type: fitPagel correlated evolution (BH), Fisher
    direction, genome-size-adjusted phyloglm effect."""
    input:
        matrix=rules.defense_matrix.output.matrix,
        master=rules.build_master_table.output.master,
        tree=rules.species_tree.output.tree,
    output:
        cooc=f"{RESULTS}/distribution/cooccurrence.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
        primary_tiers=",".join(cfg("stats/confidence_tiers_primary", ["A"])),
        min_count=cfg("stats/cooccurrence_min_count", 5),
        btol=cfg("stats/phyloglm_btol", 10),
    log:
        f"{LOGS}/distribution/cooccurrence.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=240,
    shell:
        r"""
        mkdir -p $(dirname {output.cooc})
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
        Rscript workflow/scripts/R/cooccurrence.R \
            --matrix {input.matrix} --master {input.master} --tree {input.tree} \
            --fdr {params.fdr} --seed {params.seed} --primary-tiers "{params.primary_tiers}" \
            --min-count {params.min_count} --btol {params.btol} --threads {threads} \
            --out {output.cooc} > {log} 2>&1
        """


rule plot_distribution:
    """Prevalence bars, association forest and co-occurrence heatmap."""
    input:
        prevalence=rules.compute_prevalence.output.prevalence,
        tests=rules.association_tests.output.tests,
        cooc=rules.cooccurrence.output.cooc,
    output:
        expand(f"{RESULTS}/distribution/plots/prevalence_by_taxon.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/distribution/plots/association_forest.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/distribution/plots/cooccurrence_heatmap.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/distribution/plots",
        formats=",".join(PLOT_FORMATS),
    log:
        f"{LOGS}/distribution/plot_distribution.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_distribution.py \
            --prevalence {input.prevalence} --tests {input.tests} --cooccurrence {input.cooc} \
            --outdir {params.outdir} --formats "{params.formats}" > {log} 2>&1
        """
