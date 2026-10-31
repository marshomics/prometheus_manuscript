# =====================================================================
# Stage 2 — diversity of archaeal Prometheus (sequence, architecture, context).
# =====================================================================

rule mmseqs_cluster:
    """Cluster validated Prometheus proteins into subfamilies at graded identities
    (MMseqs2 easy-cluster, bidirectional coverage, exact identity, cluster reassignment).
    Columns are cluster_<identity %.2f> (cluster_0.30 ... as in SCHEMA)."""
    input:
        faa=rules.validate_prometheus.output.proteins,
    output:
        subfam=f"{RESULTS}/diversity/subfamilies.tsv",
        summary=f"{RESULTS}/diversity/subfamily_summary.tsv",
    params:
        identities=",".join(f"{float(x):.2f}" for x in cfg("diversity/mmseqs_identities", [0.30, 0.50, 0.70, 0.90])),
        coverage=cfg("diversity/mmseqs_coverage", 0.60),
        cov_mode=cfg("diversity/mmseqs_cov_mode", 0),
        sensitivity=cfg("diversity/mmseqs_sensitivity", 7.5),
        workdir=f"{RESULTS}/diversity/mmseqs",
    log:
        f"{LOGS}/diversity/mmseqs.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p {params.workdir}
        python workflow/scripts/py/cluster_diversity.py \
            --faa {input.faa} --identities {params.identities} --coverage {params.coverage} \
            --cov-mode {params.cov_mode} --sensitivity {params.sensitivity} \
            --threads {threads} --workdir {params.workdir} \
            --out {output.subfam} --summary {output.summary} > {log} 2>&1
        """


rule architecture_typology:
    """Enumerate domain architectures (canonical, optional-domain loss, core loss,
    rearrangements, fusions) through the shared domain spec."""
    input:
        domtbl=rules.prometheus_domain_scan.output.domtbl,
        catalog=rules.validate_prometheus.output.catalog,
    output:
        typology=f"{RESULTS}/diversity/architecture_typology.tsv",
    params:
        domains_json=lambda wc: json_domains(),
    log:
        f"{LOGS}/diversity/architecture_typology.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/architecture_typology.py \
            --domtbl {input.domtbl} --catalog {input.catalog} \
            --domains '{params.domains_json}' --out {output.typology} > {log} 2>&1
        """


rule plot_diversity:
    """Subfamily rarefaction (exact, incidence-based, over all analysis genomes) and
    per-locus architecture by phylum."""
    input:
        subfam=rules.mmseqs_cluster.output.subfam,
        typology=rules.architecture_typology.output.typology,
        catalog=rules.validate_prometheus.output.catalog,
        master=f"{RESULTS}/scaffolds/master_table.tsv",
    output:
        expand(f"{RESULTS}/diversity/plots/subfamily_rarefaction.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/diversity/plots/architecture_by_taxon.{{ext}}", ext=PLOT_FORMATS),
        curve=f"{RESULTS}/diversity/subfamily_rarefaction.tsv",
    params:
        outdir=f"{RESULTS}/diversity/plots",
        seed=cfg("seed", 12345),
        permutations=cfg("diversity/rarefaction_permutations", 200),
    log:
        f"{LOGS}/diversity/plot_diversity.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_diversity.py \
            --subfamilies {input.subfam} --typology {input.typology} --master {input.master} \
            --catalog {input.catalog} --seed {params.seed} --permutations {params.permutations} \
            --curve {output.curve} --outdir {params.outdir} > {log} 2>&1
        """
