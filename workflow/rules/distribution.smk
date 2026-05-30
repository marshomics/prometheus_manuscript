# =====================================================================
# Stage 4 — distribution and taxonomic association (phylogeny-aware).
# =====================================================================

rule defense_matrix:
    """Per-genome presence/absence of every defense system (for co-occurrence + total load)."""
    input:
        master=rules.build_master_table.output.master,
        df_systems=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_systems.tsv", genome=GENOMES),
    output:
        matrix=f"{RESULTS}/distribution/defense_matrix.tsv",
    params:
        df_dir=f"{RESULTS}/defense/defensefinder",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        mkdir -p $(dirname {output.matrix})
        python workflow/scripts/py/build_defense_matrix.py \
            --master {input.master} --df-dir {params.df_dir} --out {output.matrix}
        """


rule compute_prevalence:
    """Prevalence by GTDB rank with binomial CIs, reported alongside sampling + completeness."""
    input:
        master=rules.build_master_table.output.master,
    output:
        prevalence=f"{RESULTS}/distribution/prevalence_by_taxon.tsv",
    params:
        tiers=",".join(cfg("stats/confidence_tiers_primary", ["A"])),
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        mkdir -p $(dirname {output.prevalence})
        python workflow/scripts/py/compute_prevalence.py \
            --master {input.master} --primary-tiers {params.tiers} --out {output.prevalence}
        """


rule association_tests:
    """Phylogenetic logistic regression + correlated-evolution + D, with confounder conditioning."""
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
        Rscript workflow/scripts/R/association_tests.R \
            --master {input.master} --tree {input.tree} \
            --primary-tiers {params.primary_tiers} --fdr {params.fdr} \
            --btol {params.btol} --boot {params.boot} --seed {params.seed} \
            --out {output.tests} > {log} 2>&1
        """


rule cooccurrence:
    """Prometheus vs other systems: phylogenetically-aware co-occurrence / exclusivity."""
    input:
        matrix=rules.defense_matrix.output.matrix,
        master=rules.build_master_table.output.master,
        tree=rules.species_tree.output.tree,
    output:
        cooc=f"{RESULTS}/distribution/cooccurrence.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/distribution/cooccurrence.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=240,
    shell:
        r"""
        Rscript workflow/scripts/R/cooccurrence.R \
            --matrix {input.matrix} --master {input.master} --tree {input.tree} \
            --fdr {params.fdr} --seed {params.seed} --out {output.cooc} > {log} 2>&1
        """


rule plot_distribution:
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
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_distribution.py \
            --prevalence {input.prevalence} --tests {input.tests} --cooccurrence {input.cooc} \
            --outdir {params.outdir}
        """
