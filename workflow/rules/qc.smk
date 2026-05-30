# =====================================================================
# Stage 0 — genome set, quality, provenance, standardized gene calls.
# =====================================================================

rule prodigal_genes:
    """Consistent gene calls (one coordinate space for all downstream context work)."""
    input:
        genome=genome_fasta,
    output:
        faa=f"{RESULTS}/genes/{{genome}}/{{genome}}.faa",
        ffn=f"{RESULTS}/genes/{{genome}}/{{genome}}.ffn",
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
        prodigal -i {input.genome} -a {output.faa} -d {output.ffn} \
                 -f gff -o {output.gff} -p $mode -q 2> {log} || true
        # Guarantee outputs exist even for pathological inputs (keeps DAG valid).
        touch {output.faa} {output.ffn} {output.gff}
        """


rule checkm2:
    """Completeness/contamination (ML-based; better for reduced archaeal genomes)."""
    input:
        genomes=expand(os.path.join(GENOME_DIR, "{genome}." + EXT), genome=GENOMES),
    output:
        report=f"{RESULTS}/qc/checkm2/quality_report.tsv",
    params:
        db=cfg("db/checkm2", ""),
        outdir=f"{RESULTS}/qc/checkm2",
        ext=EXT,
    log:
        f"{LOGS}/checkm2/checkm2.log",
    conda:
        "../envs/qc.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=32000, per_attempt=16000),
        runtime=720,
    shell:
        r"""
        mkdir -p {params.outdir}
        checkm2 predict --threads {threads} --force \
            --database_path {params.db} \
            --input {input.genomes} -x {params.ext} \
            --output-directory {params.outdir} > {log} 2>&1
        test -f {output.report}
        """


rule seqkit_stats:
    """Assembly statistics (size, N50, contigs, GC) independent of CheckM2."""
    input:
        genomes=expand(os.path.join(GENOME_DIR, "{genome}." + EXT), genome=GENOMES),
    output:
        f"{RESULTS}/qc/seqkit_stats.tsv",
    conda:
        "../envs/qc.yaml"
    threads: T_MED
    resources:
        mem_mb=4000, runtime=60,
    shell:
        r"""
        mkdir -p $(dirname {output})
        seqkit stats -a -T -j {threads} {input.genomes} > {output}
        """


rule barrnap_rrna:
    input:
        genome=genome_fasta,
    output:
        gff=f"{RESULTS}/qc/rrna/{{genome}}.rrna.gff",
    conda:
        "../envs/qc.yaml"
    threads: T_SMALL
    resources:
        mem_mb=2000, runtime=30,
    shell:
        r"""
        mkdir -p $(dirname {output.gff})
        barrnap --kingdom arc --threads {threads} {input.genome} > {output.gff} 2>/dev/null || touch {output.gff}
        """


rule trnascan_trna:
    input:
        genome=genome_fasta,
    output:
        txt=f"{RESULTS}/qc/trna/{{genome}}.trna.txt",
    conda:
        "../envs/qc.yaml"
    threads: T_SMALL
    resources:
        mem_mb=2000, runtime=60,
    shell:
        r"""
        mkdir -p $(dirname {output.txt})
        tRNAscan-SE -A -q -o {output.txt} {input.genome} 2>/dev/null || touch {output.txt}
        """


rule genome_metadata:
    """Aggregate QC into one row per genome, assign MIMAG-style tiers + analysis flags."""
    input:
        checkm2=rules.checkm2.output.report,
        stats=rules.seqkit_stats.output,
        rrna=expand(f"{RESULTS}/qc/rrna/{{genome}}.rrna.gff", genome=GENOMES),
        trna=expand(f"{RESULTS}/qc/trna/{{genome}}.trna.txt", genome=GENOMES),
        samples=config["samples"],
    output:
        meta=f"{RESULTS}/qc/genome_metadata.tsv",
        genomeinfo=f"{RESULTS}/qc/genomeInfo.csv",
    params:
        comp_floor=cfg("qc/completeness_floor", 50.0),
        cont_ceil=cfg("qc/contamination_ceiling", 10.0),
        comp_high=cfg("qc/tier_high_completeness", 90.0),
        cont_high=cfg("qc/tier_high_contamination", 5.0),
        ext=EXT,
        rrna_dir=f"{RESULTS}/qc/rrna",
        trna_dir=f"{RESULTS}/qc/trna",
    conda:
        "../envs/qc.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/aggregate_qc.py \
            --checkm2 {input.checkm2} --stats {input.stats} \
            --rrna-dir {params.rrna_dir} --trna-dir {params.trna_dir} \
            --samples {input.samples} --ext {params.ext} \
            --comp-floor {params.comp_floor} --cont-ceil {params.cont_ceil} \
            --comp-high {params.comp_high} --cont-high {params.cont_high} \
            --out {output.meta} --genomeinfo {output.genomeinfo}
        """


rule dereplicate:
    """Species-level dereplication so prevalence is not inflated by redundant sequencing."""
    input:
        genomes=expand(os.path.join(GENOME_DIR, "{genome}." + EXT), genome=GENOMES),
        ginfo=rules.genome_metadata.output.genomeinfo,
    output:
        reps=f"{RESULTS}/qc/representatives.txt",
    params:
        do=str(DEREPLICATE),
        ani=cfg("drep_ani", 0.95),
        outdir=f"{RESULTS}/qc/drep",
        ids=" ".join(GENOMES),
        ext=EXT,
    log:
        f"{LOGS}/drep/drep.log",
    conda:
        "../envs/qc.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=16000, per_attempt=8000),
        runtime=480,
    shell:
        r"""
        mkdir -p $(dirname {output.reps})
        if [ "{params.do}" = "True" ]; then
            dRep dereplicate {params.outdir} -p {threads} -g {input.genomes} \
                 --genomeInfo {input.ginfo} -sa {params.ani} --S_algorithm fastANI \
                 > {log} 2>&1 || true
            if [ -d {params.outdir}/dereplicated_genomes ]; then
                ls {params.outdir}/dereplicated_genomes/ | sed 's/\.{params.ext}$//' > {output.reps}
            else
                printf "%s\n" {params.ids} > {output.reps}   # fallback: keep all
            fi
        else
            printf "%s\n" {params.ids} > {output.reps}
        fi
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
