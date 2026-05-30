# =====================================================================
# Stage 6 — ecology: environment association, viral pressure proxy,
# R-loop-targeting hypothesis. All tests are phylogenetically controlled.
# =====================================================================

rule genomad:
    """Predict proviruses / viral contigs (viral-pressure proxy)."""
    input:
        genome=genome_fasta,
    output:
        virus_fna=f"{RESULTS}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus.fna",
        virus_summary=f"{RESULTS}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv",
    params:
        db=cfg("db/genomad", ""),
        outdir=f"{RESULTS}/ecology/viral/{{genome}}",
    log:
        f"{LOGS}/genomad/{{genome}}.log",
    conda:
        "../envs/viral.yaml"
    threads: T_MED
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=12000, per_attempt=8000),
        runtime=240,
    shell:
        r"""
        mkdir -p {params.outdir}
        genomad end-to-end --cleanup --threads {threads} {input.genome} {params.outdir} {params.db} > {log} 2>&1 || true
        mkdir -p $(dirname {output.virus_fna})
        [ -s {output.virus_fna} ] || : > {output.virus_fna}
        [ -s {output.virus_summary} ] || printf "seq_name\tlength\ttopology\tn_genes\tvirus_score\n" > {output.virus_summary}
        """


rule checkv:
    """Assess provirus completeness/quality."""
    input:
        virus_fna=rules.genomad.output.virus_fna,
    output:
        quality=f"{RESULTS}/ecology/viral/{{genome}}/checkv/quality_summary.tsv",
    params:
        db=cfg("db/checkv", ""),
        outdir=f"{RESULTS}/ecology/viral/{{genome}}/checkv",
    log:
        f"{LOGS}/checkv/{{genome}}.log",
    conda:
        "../envs/viral.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p {params.outdir}
        if [ -s {input.virus_fna} ]; then
            checkv end_to_end {input.virus_fna} {params.outdir} -d {params.db} -t {threads} > {log} 2>&1 || true
        fi
        [ -s {output.quality} ] || printf "contig_id\tcontig_length\tcompleteness\tcheckv_quality\n" > {output.quality}
        """


rule rloop_predict:
    """Predict R-loop-forming sequences (RLFS) in proviruses (RIZ/REZ G-skew model)."""
    input:
        virus_fna=rules.genomad.output.virus_fna,
    output:
        rloop=f"{RESULTS}/ecology/viral/{{genome}}/rloop.tsv",
    params:
        riz_min_g=cfg("ecology/rloop/riz_min_g", 3),
        window=cfg("ecology/rloop/window", 100),
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/rloop_predict.py \
            --fasta {input.virus_fna} --genome {wildcards.genome} \
            --riz-min-g {params.riz_min_g} --window {params.window} --out {output.rloop}
        """


rule viral_summary:
    """Aggregate proviral load, CheckV quality, and predicted R-loop potential per genome."""
    input:
        virus=expand(f"{RESULTS}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv", genome=GENOMES),
        checkv=expand(f"{RESULTS}/ecology/viral/{{genome}}/checkv/quality_summary.tsv", genome=GENOMES),
        rloop=expand(f"{RESULTS}/ecology/viral/{{genome}}/rloop.tsv", genome=GENOMES),
        samples=config["samples"],
    output:
        summary=f"{RESULTS}/ecology/viral/viral_summary.tsv",
    params:
        viral_dir=f"{RESULTS}/ecology/viral",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/aggregate_viral.py \
            --viral-dir {params.viral_dir} --samples {input.samples} --out {output.summary}
        """


rule ecology_stats:
    """Phylogenetically-controlled tests: environment, niche, viral pressure, R-loop."""
    input:
        master=rules.build_master_table.output.master,
        viral=rules.viral_summary.output.summary,
        tree=rules.species_tree.output.tree,
    output:
        assoc=f"{RESULTS}/ecology/ecology_associations.tsv",
    params:
        predictors=",".join(cfg("ecology/predictors", ["biome_envo", "oxygen", "temp_class", "salinity_class", "host_associated"])),
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/ecology/ecology_stats.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=240,
    shell:
        r"""
        mkdir -p $(dirname {output.assoc})
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --viral {input.viral} --tree {input.tree} \
            --predictors {params.predictors} --fdr {params.fdr} --seed {params.seed} \
            --out {output.assoc} > {log} 2>&1
        """


rule plot_ecology:
    input:
        assoc=rules.ecology_stats.output.assoc,
        viral=rules.viral_summary.output.summary,
        master=rules.build_master_table.output.master,
    output:
        expand(f"{RESULTS}/ecology/plots/ecology_enrichment.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/ecology/plots/viral_rloop_association.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/ecology/plots",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_ecology.py \
            --associations {input.assoc} --viral {input.viral} --master {input.master} \
            --outdir {params.outdir}
        """
