# =====================================================================
# Extra figures (added ALONGSIDE the existing plots; nothing is replaced).
# Each reuses plotting_style.py (600-dpi PNG + Illustrator-editable SVG) and
# degrades to a labelled placeholder on empty input, like the existing plots.
# Foreground: distribution/architecture, ecology, evolution.
# =====================================================================

rule plot_qc_extra:
    """Covariate diagnostic: genome size / defense load / completeness by presence."""
    input:
        master=rules.build_master_table.output.master,
    output:
        expand(f"{RESULTS}/qc/plots/prometheus_covariates.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/qc/plots",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_qc_extra.py --master {input.master} --outdir {params.outdir}
        """


rule plot_distribution_extra:
    """Prevalence caterpillar (Wilson CI) + directional co-occurrence odds ratios."""
    input:
        prevalence=rules.compute_prevalence.output.prevalence,
        cooccurrence=rules.cooccurrence.output.cooc,
        master=rules.build_master_table.output.master,
    output:
        expand(f"{RESULTS}/distribution/plots/prevalence_caterpillar.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/distribution/plots/cooccurrence_oddsratio.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/distribution/plots",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_distribution_extra.py \
            --prevalence {input.prevalence} --cooccurrence {input.cooccurrence} \
            --master {input.master} --outdir {params.outdir}
        """


rule plot_diversity_extra:
    """Domain-combination UpSet + architecture-by-phylum + subfamily-vs-identity."""
    input:
        typology=rules.architecture_typology.output.typology,
        catalog=rules.validate_prometheus.output.catalog,
        master=rules.build_master_table.output.master,
        subfamilies=rules.mmseqs_cluster.output.subfam,
    output:
        expand(f"{RESULTS}/diversity/plots/architecture_combinations.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/diversity/plots/architecture_by_phylum.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/diversity/plots/subfamily_identity_curve.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/diversity/plots",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_diversity_extra.py \
            --typology {input.typology} --catalog {input.catalog} --master {input.master} \
            --subfamilies {input.subfamilies} --outdir {params.outdir}
        """


rule plot_ecology_extra:
    """Viral/R-loop pressure vs presence + niche effect forest + niche x phylum grid."""
    input:
        master=rules.build_master_table.output.master,
        viral=rules.viral_summary.output.summary,
        associations=rules.ecology_stats.output.assoc,
    output:
        expand(f"{RESULTS}/ecology/plots/viral_pressure_prometheus.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/ecology/plots/niche_effect_forest.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/ecology/plots/prevalence_niche_phylum.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/ecology/plots",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_ecology_extra.py \
            --master {input.master} --viral {input.viral} --associations {input.associations} \
            --outdir {params.outdir}
        """


rule plot_evolution_extra:
    """Per-site dN/dS along the protein + reconciliation events + gain-enrichment null."""
    input:
        sites=rules.selection.output.sites,
        reconciliation=rules.reconciliation.output.events,
        gains=rules.ancestral_states.output.gl,
        gains_null=rules.ancestral_states.output.gl_null,
    output:
        expand(f"{RESULTS}/evolution/plots/dnds_along_protein.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/evolution/plots/reconciliation_events.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/evolution/plots/gain_enrichment_null.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/evolution/plots",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_evolution_extra.py \
            --sites {input.sites} --reconciliation {input.reconciliation} \
            --gains {input.gains} --gains-null {input.gains_null} --outdir {params.outdir}
        """
