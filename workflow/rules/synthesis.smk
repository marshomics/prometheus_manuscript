# =====================================================================
# Stage 7 — synthesis: integrated figures, one FDR regime, reproducibility.
#
# Config read here (defaults keep the rules valid before a key exists):
#   stats.fdr_method               BH (also BY, holm, bonferroni; anything else fails)
#   stats.confidence_tiers_primary the primary tier set; registry rows of other tier
#                                  sets are sensitivity rows outside the FDR family
#   stats.confirmatory             optional list of "family:test:term" keys replacing the
#                                  registry's built-in confirmatory set
#   stats.alpha                    p_adj threshold used to colour the summary forest
#   plotting.formats               figure formats (PLOT_FORMATS)
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
        prefix=lambda wc, output: os.path.splitext(output[0])[0],
    log:
        f"{LOGS}/synthesis/validation_tree.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        Rscript workflow/scripts/R/plot_validation_tree.R \
            --tree {input.tree} --labels {input.labels} --catalog {input.catalog} \
            --plot-prefix {params.prefix} > {log} 2>&1
        """


rule master_tree_figure:
    """Annotated ROOTED species tree: tier, mcr metabolism (methanogen / ANME), subtype,
    ecology, proviral load and R-loop fraction rings. Tips without a master/viral row are
    drawn as 'no data', never as absences. Also writes the per-tip ring values."""
    input:
        tree=rules.species_tree.output.tree,
        master=rules.build_master_table.output.master,
        viral=rules.viral_summary.output.summary,
        gl=rules.ancestral_states.output.gl,
    output:
        figs=expand(f"{RESULTS}/synthesis/figures/master_tree.{{ext}}", ext=PLOT_FORMATS),
        annotation=f"{RESULTS}/synthesis/figures/master_tree_annotation.tsv",
    params:
        prefix=lambda wc, output: os.path.splitext(output.figs[0])[0],
    log:
        f"{LOGS}/synthesis/master_tree.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_SMALL
    resources:
        mem_mb=16000, runtime=120,
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        Rscript workflow/scripts/R/plot_tree_master.R \
            --tree {input.tree} --master {input.master} --viral {input.viral} --gains {input.gl} \
            --plot-prefix {params.prefix} \
            --annotation-out {output.annotation} > {log} 2>&1
        """


rule results_registry:
    """Collect EVERY hypothesis test (association, co-occurrence, ecology, the
    within-methanogen pathway tests of methanogenesis_stats, selection, reconciliation,
    ancestral; disorder when enabled); apply one multiple-testing regime (BH) within each
    family over the rows with a p-value (NA p-values do not count in m; the naive
    independence baselines and non-primary tier sets are reported but kept outside the
    family; the methanogenesis family, always exploratory, counts its primary
    habitat-adjusted rows only); tag confirmatory vs exploratory by exact (family, test,
    term) keys. A non-empty family that lacks its confirmatory key fails the rule
    (--strict-confirmatory)."""
    input:
        association=rules.association_tests.output.tests,
        cooccurrence=rules.cooccurrence.output.cooc,
        ecology=rules.ecology_stats.output.assoc,
        methanogenesis=([rules.methanogenesis_stats.output.assoc]
                       if MGEN_ENABLED else []),   # opt-in (methanogenesis.enabled)
        selection=rules.selection.output.summary,
        reconciliation=rules.reconciliation.output.events,
        ancestral=rules.ancestral_states.output.gl,
        # disorder.smk is included after this file, so reference its output by path
        disorder=lambda wc: ([f"{RESULTS}/disorder/disorder_associations.tsv"]
                             if DISORDER_ENABLED else []),
    output:
        registry=f"{RESULTS}/synthesis/results_registry.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
        primary_tiers=",".join(cfg("stats/confidence_tiers_primary", ["A"]) or ["A"]),
        confirmatory=";".join(cfg("stats/confirmatory", []) or []),
    log:
        f"{LOGS}/synthesis/results_registry.log",
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
            --ecology {input.ecology} --methanogenesis "{input.methanogenesis}" \
            --selection {input.selection} \
            --reconciliation {input.reconciliation} --ancestral {input.ancestral} \
            --disorder "{input.disorder}" \
            --fdr {params.fdr} --primary-tiers "{params.primary_tiers}" \
            --confirmatory "{params.confirmatory}" --strict-confirmatory \
            --out {output.registry} > {log} 2>&1
        """


rule plot_summary_forest:
    """One panel per confirmatory registry row, each on its own scale and null value."""
    input:
        registry=rules.results_registry.output.registry,
    output:
        expand(f"{RESULTS}/synthesis/figures/summary_forest.{{ext}}", ext=PLOT_FORMATS),
    params:
        figdir=lambda wc, output: os.path.dirname(output[0]),
        alpha=cfg("stats/alpha", 0.05),
        formats=",".join(PLOT_FORMATS),
    log:
        f"{LOGS}/synthesis/summary_forest.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_summary_forest.py \
            --registry {input.registry} --outdir {params.figdir} \
            --alpha {params.alpha} --formats {params.formats} > {log} 2>&1
        """
