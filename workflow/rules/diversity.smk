# =====================================================================
# Stage 2 — diversity of archaeal Prometheus (sequence, architecture, context).
# =====================================================================

rule mmseqs_cluster:
    """Cluster validated Prometheus proteins into subfamilies at graded identities."""
    input:
        faa=rules.validate_prometheus.output.proteins,
    output:
        subfam=f"{RESULTS}/diversity/subfamilies.tsv",
    params:
        identities=",".join(str(x) for x in cfg("diversity/mmseqs_identities", [0.30, 0.50, 0.70, 0.90])),
        coverage=cfg("diversity/mmseqs_coverage", 0.60),
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
            --threads {threads} --workdir {params.workdir} --out {output.subfam} > {log} 2>&1
        """


rule architecture_typology:
    """Enumerate domain architectures (canonical 4-domain, losses, insertions, reorderings)."""
    input:
        domtbl=rules.prometheus_domain_scan.output.domtbl,
        catalog=rules.validate_prometheus.output.catalog,
    output:
        typology=f"{RESULTS}/diversity/architecture_typology.tsv",
    params:
        domains_json=lambda wc: json_domains(),
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/architecture_typology.py \
            --domtbl {input.domtbl} --catalog {input.catalog} \
            --domains '{params.domains_json}' --out {output.typology}
        """


rule plot_diversity:
    input:
        subfam=rules.mmseqs_cluster.output.subfam,
        typology=rules.architecture_typology.output.typology,
        master=f"{RESULTS}/scaffolds/master_table.tsv",
    output:
        expand(f"{RESULTS}/diversity/plots/subfamily_rarefaction.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/diversity/plots/architecture_by_taxon.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/diversity/plots",
        seed=cfg("seed", 12345),
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_diversity.py \
            --subfamilies {input.subfam} --typology {input.typology} --master {input.master} \
            --seed {params.seed} --outdir {params.outdir}
        """
