# =====================================================================
# Stage 7 — synthesis: integrated figures, one FDR regime, reproducibility.
# =====================================================================

rule plot_validation_tree:
    """Prometheus validation tree: candidates vs Prometheus / housekeeping-paralog references."""
    input:
        tree=rules.prometheus_validation_tree.output.tree,
        labels=rules.prometheus_validation_tree.output.labels,
        catalog=rules.validate_prometheus.output.catalog,
    output:
        expand(f"{RESULTS}/prometheus/plots/prometheus_validation_tree.{{ext}}", ext=PLOT_FORMATS),
    params:
        plotdir=f"{RESULTS}/prometheus/plots",
    log:
        f"{LOGS}/synthesis/validation_tree.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        mkdir -p {params.plotdir}
        Rscript workflow/scripts/R/plot_validation_tree.R \
            --tree {input.tree} --labels {input.labels} --catalog {input.catalog} \
            --plot-prefix {params.plotdir}/prometheus_validation_tree > {log} 2>&1
        """


rule master_tree_figure:
    """Annotated archaeal species tree: presence/tier, methanogen, subtype, ecology, proviral load."""
    input:
        tree=rules.species_tree.output.tree,
        master=rules.build_master_table.output.master,
        viral=rules.viral_summary.output.summary,
        gl=rules.ancestral_states.output.gl,
    output:
        expand(f"{RESULTS}/synthesis/figures/master_tree.{{ext}}", ext=PLOT_FORMATS),
    params:
        figdir=f"{RESULTS}/synthesis/figures",
    log:
        f"{LOGS}/synthesis/master_tree.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_SMALL
    resources:
        mem_mb=16000, runtime=120,
    shell:
        r"""
        mkdir -p {params.figdir}
        Rscript workflow/scripts/R/plot_tree_master.R \
            --tree {input.tree} --master {input.master} --viral {input.viral} --gains {input.gl} \
            --plot-prefix {params.figdir}/master_tree > {log} 2>&1
        """


rule results_registry:
    """Collect every hypothesis test; apply one FDR regime per family; mark confirmatory vs exploratory."""
    input:
        association=rules.association_tests.output.tests,
        cooccurrence=rules.cooccurrence.output.cooc,
        ecology=rules.ecology_stats.output.assoc,
        selection=rules.selection.output.summary,
        prevalence=rules.compute_prevalence.output.prevalence,
    output:
        registry=f"{RESULTS}/synthesis/results_registry.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        mkdir -p $(dirname {output.registry})
        python workflow/scripts/py/results_registry.py \
            --association {input.association} --cooccurrence {input.cooccurrence} \
            --ecology {input.ecology} --selection {input.selection} --prevalence {input.prevalence} \
            --fdr {params.fdr} --out {output.registry}
        """


rule plot_summary_forest:
    input:
        registry=rules.results_registry.output.registry,
    output:
        expand(f"{RESULTS}/synthesis/figures/summary_forest.{{ext}}", ext=PLOT_FORMATS),
    params:
        figdir=f"{RESULTS}/synthesis/figures",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_summary_forest.py \
            --registry {input.registry} --outdir {params.figdir}
        """
