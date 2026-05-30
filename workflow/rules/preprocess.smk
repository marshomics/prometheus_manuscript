# =====================================================================
# Phase 1 — preprocessing of the raw input genome set (e.g. 92k genomes).
#   batch split  -> GTDB-Tk (batched) + CheckM2 (batched) + BBMap stats (batched)
#   -> merge/annotate + quality score -> quality filter + collapse to one
#      representative per GTDB species (union niche/group presence-absence)
#   -> species-level table + staged representative genomes for the analysis.
# Run with:  snakemake preprocess --use-conda ...
# No dRep / other QC runs before the analysis; this replaces it.
# =====================================================================

wildcard_constraints:
    b = r"\d+"

PRE = f"{RESULTS}/preprocess"


rule batch_genomes:
    """Split the input genomes into deterministic batches; stage per-batch symlinks
    named {genome_id}.fna so GTDB-Tk/CheckM2/BBMap report stable ids."""
    input:
        metadata=INPUT_METADATA,
    output:
        ready=expand(f"{PRE}/batches/batch_{{b}}/.ready", b=BATCHES),
    params:
        genomes_dir=INPUT_GENOMES_DIR,
        id_col=ID_COL, file_col=FILE_COL,
        input_ext=INPUT_GENOME_EXT,
        batch_size=BATCH_SIZE,
        outdir=f"{PRE}/batches",
    log:
        f"{LOGS}/preprocess/batch_genomes.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        python workflow/scripts/py/batch_genomes.py \
            --metadata {input.metadata} --genomes-dir "{params.genomes_dir}" \
            --id-col {params.id_col} --file-col {params.file_col} --input-ext {params.input_ext} \
            --batch-size {params.batch_size} --outdir {params.outdir} > {log} 2>&1
        """


rule gtdbtk_batch:
    """GTDB-Tk classify_wf on one batch (~20k genomes)."""
    input:
        ready=f"{PRE}/batches/batch_{{b}}/.ready",
    output:
        summary=f"{PRE}/batches/batch_{{b}}/gtdbtk_summary.tsv",
    params:
        gdir=f"{PRE}/batches/batch_{{b}}/genomes",
        outdir=f"{PRE}/batches/batch_{{b}}/gtdbtk",
        data=cfg("db/gtdbtk_data", ""),
        extra=cfg("tree/gtdbtk_classify_extra", "--skip_ani_screen"),
    log:
        f"{LOGS}/preprocess/gtdbtk_batch_{{b}}.log",
    conda:
        "../envs/gtdbtk.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=128000, per_attempt=64000),
        runtime=4320,
    shell:
        r"""
        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.data}
        gtdbtk classify_wf --genome_dir {params.gdir} -x fna \
            --out_dir {params.outdir} --cpus {threads} {params.extra} > {log} 2>&1 || true
        # Merge ar53 (+ bac120 if any) summaries, keeping a single header.
        first=1; : > {output.summary}
        for s in $(find {params.outdir} -name '*summary.tsv' 2>/dev/null); do
            if [ "$first" = "1" ]; then cat "$s" > {output.summary}; first=0;
            else tail -n +2 "$s" >> {output.summary}; fi
        done
        [ -s {output.summary} ] || printf "user_genome\tclassification\n" > {output.summary}
        """


rule checkm2_batch:
    """CheckM2 completeness/contamination on one batch."""
    input:
        ready=f"{PRE}/batches/batch_{{b}}/.ready",
    output:
        report=f"{PRE}/batches/batch_{{b}}/checkm2_quality.tsv",
    params:
        gdir=f"{PRE}/batches/batch_{{b}}/genomes",
        outdir=f"{PRE}/batches/batch_{{b}}/checkm2",
        db=cfg("db/checkm2", ""),
    log:
        f"{LOGS}/preprocess/checkm2_batch_{{b}}.log",
    conda:
        "../envs/qc.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=48000, per_attempt=24000),
        runtime=2880,
    shell:
        r"""
        mkdir -p {params.outdir}
        checkm2 predict --threads {threads} --force --database_path {params.db} \
            --input {params.gdir} -x fna --output-directory {params.outdir} > {log} 2>&1 || true
        if [ -s {params.outdir}/quality_report.tsv ]; then
            cp {params.outdir}/quality_report.tsv {output.report}
        else
            printf "Name\tCompleteness\tContamination\n" > {output.report}
        fi
        """


rule bbmap_batch:
    """BBMap assembly statistics (statswrapper.sh format=6) on one batch."""
    input:
        ready=f"{PRE}/batches/batch_{{b}}/.ready",
    output:
        stats=f"{PRE}/batches/batch_{{b}}/bbmap_stats.tsv",
    params:
        gdir=f"{PRE}/batches/batch_{{b}}/genomes",
        chunk=2000,
    log:
        f"{LOGS}/preprocess/bbmap_batch_{{b}}.log",
    conda:
        "../envs/bbmap.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=720,
    shell:
        r"""
        python workflow/scripts/py/run_bbmap_stats.py \
            --genomes-dir {params.gdir} --ext fna --chunk {params.chunk} \
            --out {output.stats} > {log} 2>&1
        """


rule merge_annotations:
    """Join GTDB-Tk + CheckM2 + BBMap onto the input metadata; add quality_score."""
    input:
        metadata=INPUT_METADATA,
        gtdbtk=expand(f"{PRE}/batches/batch_{{b}}/gtdbtk_summary.tsv", b=BATCHES),
        checkm2=expand(f"{PRE}/batches/batch_{{b}}/checkm2_quality.tsv", b=BATCHES),
        bbmap=expand(f"{PRE}/batches/batch_{{b}}/bbmap_stats.tsv", b=BATCHES),
    output:
        annotated=f"{PRE}/annotated_metadata.tsv",
    params:
        gc_col=cfg("preprocess/bbmap/gc_col", "gc_avg"),
        length_col=cfg("preprocess/bbmap/length_col", "scaf_bp"),
        contigs_col=cfg("preprocess/bbmap/contigs_col", "n_contigs"),
        n50_col=cfg("preprocess/bbmap/n50_col", "ctg_L50"),
        id_col=ID_COL,
        batch_dir=f"{PRE}/batches",
    log:
        f"{LOGS}/preprocess/merge_annotations.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=16000, per_attempt=8000),
        runtime=120,
    shell:
        r"""
        python workflow/scripts/py/merge_annotations.py \
            --metadata {input.metadata} --batch-dir {params.batch_dir} --id-col {params.id_col} \
            --gc-col {params.gc_col} --length-col {params.length_col} \
            --contigs-col {params.contigs_col} --n50-col {params.n50_col} \
            --out {output.annotated} > {log} 2>&1
        """


rule filter_collapse:
    """Quality-filter, then collapse to one representative per GTDB species,
    unioning niche/group presence-absence across each species' genomes."""
    input:
        annotated=rules.merge_annotations.output.annotated,
    output:
        species_meta=f"{PRE}/species_metadata.tsv",
        reps_list=f"{PRE}/representatives.txt",
        genome_meta=f"{PRE}/genome_metadata.tsv",
        taxonomy=f"{PRE}/taxonomy.tsv",
        staged=f"{PRE}/genomes/.staged",
    params:
        genomes_dir=INPUT_GENOMES_DIR,
        input_ext=INPUT_GENOME_EXT,
        file_col=FILE_COL, src_col=SRC_COL,
        min_comp=cfg("preprocess/quality/min_completeness", 90.0),
        max_cont=cfg("preprocess/quality/max_contamination", 5.0),
        max_contigs=cfg("preprocess/quality/max_contigs", 1000),
        min_n50=cfg("preprocess/quality/min_n50", 5000),
        min_qs=cfg("preprocess/quality/min_quality_score", 50.0),
        niche_agg=cfg("preprocess/niche_aggregation", "all"),
        unclassified=cfg("preprocess/unclassified_species", "drop"),
        keep_domain=cfg("preprocess/keep_domain", "Archaea"),
        stage_dir=f"{PRE}/genomes",
    log:
        f"{LOGS}/preprocess/filter_collapse.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=16000, per_attempt=8000),
        runtime=120,
    shell:
        r"""
        python workflow/scripts/py/filter_collapse.py \
            --annotated {input.annotated} --genomes-dir "{params.genomes_dir}" \
            --input-ext {params.input_ext} \
            --file-col {params.file_col} --src-col {params.src_col} \
            --min-completeness {params.min_comp} --max-contamination {params.max_cont} \
            --max-contigs {params.max_contigs} --min-n50 {params.min_n50} \
            --min-quality-score {params.min_qs} \
            --niche-aggregation {params.niche_agg} --unclassified {params.unclassified} \
            --keep-domain "{params.keep_domain}" \
            --stage-dir {params.stage_dir} \
            --out-species {output.species_meta} \
            --out-reps-list {output.reps_list} --out-genome-meta {output.genome_meta} \
            --out-taxonomy {output.taxonomy} --out-staged {output.staged} > {log} 2>&1
        """


rule plot_preprocess:
    """Diagnostic: quality scatter + filter pass/fail + per-species collapse."""
    input:
        annotated=rules.merge_annotations.output.annotated,
        species=rules.filter_collapse.output.species_meta,
    output:
        expand(f"{PRE}/plots/preprocess_quality.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{PRE}/plots",
        min_comp=cfg("preprocess/quality/min_completeness", 90.0),
        max_cont=cfg("preprocess/quality/max_contamination", 5.0),
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/plot_preprocess.py \
            --annotated {input.annotated} --species {input.species} \
            --min-completeness {params.min_comp} --max-contamination {params.max_cont} \
            --outdir {params.outdir}
        """
