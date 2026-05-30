# =====================================================================
# Stage 0 — standardized gene calls + per-genome QC table.
# Quality, taxonomy, assembly stats, species dereplication, and the niche
# union are all done in Phase 1 (rules/preprocess.smk); the analysis runs on
# the species representatives that preprocessing staged. No dRep / extra QC here.
# =====================================================================

rule prodigal_genes:
    """Consistent gene calls (one coordinate space for all downstream context work)."""
    input:
        genome=genome_fasta,
    output:
        faa=f"{RESULTS}/genes/{{genome}}/{{genome}}.faa",
        fna=f"{RESULTS}/genes/{{genome}}/{{genome}}.fna",
        gff=f"{RESULTS}/genes/{{genome}}/{{genome}}.gff",
    params:
        mode=cfg("qc/prodigal_mode", "single"),
    log:
        f"{LOGS}/prodigal/{{genome}}.log",
    conda:
        "../envs/qc.yaml"
    threads: T_SMALL
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=2000, per_attempt=2000),
        runtime=60,
    shell:
        r"""
        mkdir -p $(dirname {output.faa})
        # Fall back to 'meta' automatically for very small/fragmented inputs.
        bp=$(grep -v '^>' {input.genome} | tr -d '\n' | wc -c)
        mode={params.mode}
        if [ "$bp" -lt 20000 ]; then mode=meta; fi
        prodigal -i {input.genome} -a {output.faa} -d {output.fna} \
                 -f gff -o {output.gff} -p $mode -q 2> {log} || true
        # Guarantee outputs exist even for pathological inputs (keeps DAG valid).
        touch {output.faa} {output.fna} {output.gff}
        # Rename headers to genome-underscore-number, identical for the matching
        # protein and CDS so .faa and .fna stay consistent gene-by-gene.
        python workflow/scripts/py/rename_prodigal.py \
            --faa {output.faa} --fna {output.fna} --genome {wildcards.genome} >> {log} 2>&1
        """


rule genome_metadata:
    """Per-genome QC table for the analysis, taken from preprocessing (completeness,
    contamination, BBMap assembly stats, quality tier). No QC is recomputed here:
    every analysis genome already passed the Phase-1 quality filter."""
    input:
        meta=rules.filter_collapse.output.genome_meta,
    output:
        meta=f"{RESULTS}/qc/genome_metadata.tsv",
    threads: T_SMALL
    resources:
        mem_mb=2000, runtime=10,
    shell:
        r"""
        mkdir -p $(dirname {output.meta})
        cp {input.meta} {output.meta}
        """


rule plot_qc:
    input:
        meta=rules.genome_metadata.output.meta,
    output:
        expand(f"{RESULTS}/qc/plots/qc_completeness_contamination.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/qc/plots/qc_assembly_stats.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/qc/plots",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_qc.py --meta {input.meta} --outdir {params.outdir}
        """
