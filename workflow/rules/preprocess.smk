# =====================================================================
# Phase 1 — preprocessing of the raw input genome set (e.g. 92k genomes).
#   batch split  -> GTDB-Tk (batched) + CheckM2 (batched) + BBMap stats (batched)
#   -> merge/annotate + quality score -> quality filter + collapse to one
#      representative per GTDB species (union niche/group presence-absence)
#   -> species-level table + staged representative genomes for the analysis.
# Run with:  snakemake preprocess --use-conda ...
# No dRep / other QC runs before the analysis; this replaces it.
#
# Failure policy (Contract v2): a crashed tool or a missing database fails its rule.
# Genomes a tool ran on but did not report are flagged per genome (gtdbtk_status /
# checkm2_status / bbmap_status in annotated_metadata.tsv) and listed with the reason in
# preprocess/filter_report.tsv; nothing disappears silently.
# =====================================================================

wildcard_constraints:
    b = r"\d+"

PRE = f"{RESULTS}/preprocess"
# GTDB-Tk 2.4.0 refuses user genomes whose id equals a GTDB reference accession
# (gtdbtk/main.py _genomes_to_process); batch_genomes.py gives every genome this prefix
# in the GTDB-Tk batchfile and gtdbtk_batch strips it from the merged summary.
GTDBTK_ID_PREFIX = "user_"


rule batch_genomes:
    """Split the input genomes into deterministic batches (sorted by id); stage per-batch
    symlinks named {genome_id}.fna so CheckM2/BBMap report the genome_id, write the
    GTDB-Tk batchfile per batch, and a manifest of every genome. Fails on an invalid id or
    a missing/empty/compressed FASTA."""
    input:
        metadata=INPUT_METADATA,
    output:
        ready=expand(f"{PRE}/batches/batch_{{b}}/.ready", b=BATCHES),
        batchfiles=expand(f"{PRE}/batches/batch_{{b}}/gtdbtk_batchfile.tsv", b=BATCHES),
        manifest=f"{PRE}/batches/manifest.tsv",
    params:
        genomes_dir=INPUT_GENOMES_DIR,
        id_col=ID_COL, file_col=FILE_COL,
        input_ext=INPUT_GENOME_EXT,
        batch_size=BATCH_SIZE,
        expected_n=N_INPUT,
        gtdbtk_prefix=GTDBTK_ID_PREFIX,
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
            --batch-size {params.batch_size} --expected-n {params.expected_n} \
            --gtdbtk-id-prefix {params.gtdbtk_prefix} \
            --outdir {params.outdir} --manifest {output.manifest} > {log} 2>&1
        """


rule gtdbtk_batch:
    """GTDB-Tk classify_wf on one batch (~20k genomes). Output = the batch's ar53 + bac120
    classification summaries (classify/gtdbtk.{ar53,bac120}.summary.tsv, GTDB-Tk 2.4.0
    paths; never the root-level symlinks or identify/*markers_summary.tsv), one header,
    plus a marker_set column, with user_genome = genome_id. Genomes GTDB-Tk could not
    process are reported by GTDB-Tk itself as 'Unclassified*' rows and stay visible."""
    input:
        ready=f"{PRE}/batches/batch_{{b}}/.ready",
        batchfile=f"{PRE}/batches/batch_{{b}}/gtdbtk_batchfile.tsv",
    output:
        summary=f"{PRE}/batches/batch_{{b}}/gtdbtk_summary.tsv",
    params:
        outdir=f"{PRE}/batches/batch_{{b}}/gtdbtk",
        data=cfg("db/gtdbtk_data", ""),
        extra=cfg("tree/gtdbtk_classify_extra", "--skip_ani_screen"),
        prefix=GTDBTK_ID_PREFIX,
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
        if [ -n "{params.data}" ]; then
            if [ ! -d "{params.data}" ]; then
                echo "ERROR: GTDB-Tk reference data '{params.data}' (config db.gtdbtk_data) is not a directory" > {log}
                exit 1
            fi
            export GTDBTK_DATA_PATH="$(cd "{params.data}" && pwd)"
        fi
        if [ -z "${{GTDBTK_DATA_PATH:-}}" ] || [ ! -d "${{GTDBTK_DATA_PATH}}" ]; then
            echo "ERROR: GTDBTK_DATA_PATH is unset or not a directory; set config db.gtdbtk_data" > {log}
            exit 1
        fi
        # A fresh output directory: GTDB-Tk 2.4 reuses a completed ANI screen found in
        # out_dir/gtdbtk.json, which would be stale if the batch composition changed.
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        python -c "import gtdbtk; print('gtdbtk', gtdbtk.__version__)" > {log} 2>&1
        echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH" >> {log}
        gtdbtk classify_wf --batchfile {input.batchfile} --out_dir {params.outdir} \
            --cpus {threads} {params.extra} >> {log} 2>&1
        files=""
        for s in {params.outdir}/classify/gtdbtk.ar53.summary.tsv {params.outdir}/classify/gtdbtk.bac120.summary.tsv; do
            if [ -f "$s" ]; then files="$files $s"; fi
        done
        if [ -z "$files" ]; then
            echo "ERROR: GTDB-Tk exited 0 but wrote no classify/gtdbtk.{{ar53,bac120}}.summary.tsv" >> {log}
            exit 1
        fi
        LC_ALL=C awk -v pre="{params.prefix}" '
            BEGIN {{ FS = OFS = "\t" }}
            FNR == 1 {{
                ms = (FILENAME ~ /\.ar53\.summary\.tsv$/) ? "ar53" : "bac120"
                if (NR == 1) {{
                    hdr = $0
                    if ($1 != "user_genome") {{ print "ERROR: unexpected GTDB-Tk summary header: " $0 > "/dev/stderr"; bad = 1; exit 1 }}
                    print $0, "marker_set"; next
                }}
                if ($0 != hdr) {{ print "ERROR: ar53/bac120 summary headers differ" > "/dev/stderr"; bad = 1; exit 1 }}
                next
            }}
            {{
                if (index($1, pre) != 1) {{ print "ERROR: user_genome without the batch prefix: " $1 > "/dev/stderr"; bad = 1; exit 1 }}
                $1 = substr($1, length(pre) + 1)
                print $0, ms
            }}
            END {{ if (bad) exit 1 }}' $files > {output.summary}.tmp 2>> {log}
        mv {output.summary}.tmp {output.summary}
        echo "merged summaries:$files -> {output.summary} ($(($(wc -l < {output.summary}) - 1)) genomes)" >> {log}
        """


rule checkm2_batch:
    """CheckM2 completeness/contamination on one batch. quality_report.tsv `Name` is the
    file stem (checkm2/prodigal.py: os.path.splitext(basename)), i.e. the genome_id."""
    input:
        ready=f"{PRE}/batches/batch_{{b}}/.ready",
    output:
        report=f"{PRE}/batches/batch_{{b}}/checkm2_quality.tsv",
    params:
        gdir=f"{PRE}/batches/batch_{{b}}/genomes",
        outdir=f"{PRE}/batches/batch_{{b}}/checkm2",
        db=cfg("db/checkm2", ""),
        extra=cfg("preprocess/checkm2_extra", "--remove_intermediates"),
    log:
        f"{LOGS}/preprocess/checkm2_batch_{{b}}.log",
    conda:
        "../envs/checkm2.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=48000, per_attempt=24000),
        runtime=2880,
    shell:
        r"""
        if [ ! -s "{params.db}" ]; then
            echo "ERROR: CheckM2 DIAMOND database '{params.db}' (config db.checkm2) is missing or empty" > {log}
            exit 1
        fi
        rm -rf {params.outdir}
        echo "checkm2 $(checkm2 --version 2>&1)" > {log}
        checkm2 predict --threads {threads} --force --database_path {params.db} \
            --input {params.gdir} -x fna --output-directory {params.outdir} {params.extra} >> {log} 2>&1
        cp {params.outdir}/quality_report.tsv {output.report}
        """


rule bbmap_batch:
    """BBMap assembly statistics (statswrapper.sh format=6) on one batch; fails unless every
    genome of the batch gets exactly one stats row."""
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
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=720,
    shell:
        r"""
        python workflow/scripts/py/run_bbmap_stats.py \
            --genomes-dir {params.gdir} --ext fna --chunk {params.chunk} \
            --out {output.stats} > {log} 2>&1
        """


rule merge_annotations:
    """Join GTDB-Tk + CheckM2 + BBMap onto the input metadata; add quality_score and the
    per-tool status flags."""
    input:
        metadata=INPUT_METADATA,
        manifest=rules.batch_genomes.output.manifest,
        gtdbtk=expand(f"{PRE}/batches/batch_{{b}}/gtdbtk_summary.tsv", b=BATCHES),
        checkm2=expand(f"{PRE}/batches/batch_{{b}}/checkm2_quality.tsv", b=BATCHES),
        bbmap=expand(f"{PRE}/batches/batch_{{b}}/bbmap_stats.tsv", b=BATCHES),
    output:
        annotated=f"{PRE}/annotated_metadata.tsv",
    params:
        gc_col=cfg("preprocess/bbmap/gc_col", "gc_avg"),
        length_col=cfg("preprocess/bbmap/length_col", "contig_bp"),
        contigs_col=cfg("preprocess/bbmap/contigs_col", "n_contigs"),
        n50_col=cfg("preprocess/bbmap/n50_col", "ctg_L50"),
        id_col=ID_COL,
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
            --metadata {input.metadata} --manifest {input.manifest} --id-col {params.id_col} \
            --gtdbtk {input.gtdbtk} --checkm2 {input.checkm2} --bbmap {input.bbmap} \
            --gc-col {params.gc_col} --length-col {params.length_col} \
            --contigs-col {params.contigs_col} --n50-col {params.n50_col} \
            --out {output.annotated} > {log} 2>&1
        """


rule filter_collapse:
    """Quality-filter, then collapse to one representative per GTDB species, unioning
    niche/group presence-absence across each species' genomes. Stages the representatives
    where Phase 2 reads them (config genomes_dir / genome_extension)."""
    input:
        annotated=rules.merge_annotations.output.annotated,
    output:
        species_meta=f"{PRE}/species_metadata.tsv",
        reps_list=f"{PRE}/representatives.txt",
        genome_meta=f"{PRE}/genome_metadata.tsv",
        taxonomy=f"{PRE}/taxonomy.tsv",
        staged=f"{PRE}/genomes/.staged",
        report=f"{PRE}/filter_report.tsv",
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
        tier_hi_comp=cfg("qc/tier_high_completeness", 90.0),
        tier_hi_cont=cfg("qc/tier_high_contamination", 5.0),
        comp_floor=cfg("qc/completeness_floor", 50.0),
        cont_ceil=cfg("qc/contamination_ceiling", 10.0),
        stage_dir=GENOME_DIR,
        stage_ext=EXT,
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
            --tier-high-completeness {params.tier_hi_comp} --tier-high-contamination {params.tier_hi_cont} \
            --completeness-floor {params.comp_floor} --contamination-ceiling {params.cont_ceil} \
            --stage-dir {params.stage_dir} --stage-ext {params.stage_ext} \
            --out-species {output.species_meta} \
            --out-reps-list {output.reps_list} --out-genome-meta {output.genome_meta} \
            --out-taxonomy {output.taxonomy} --out-staged {output.staged} \
            --report {output.report} > {log} 2>&1
        """


rule plot_preprocess:
    """Diagnostic: quality scatter (all five filter thresholds applied via the filter
    report) + the input -> domain -> species -> passing -> representative funnel."""
    input:
        annotated=rules.merge_annotations.output.annotated,
        species=rules.filter_collapse.output.species_meta,
        report=rules.filter_collapse.output.report,
    output:
        expand(f"{PRE}/plots/preprocess_quality.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{PRE}/plots",
        min_comp=cfg("preprocess/quality/min_completeness", 90.0),
        max_cont=cfg("preprocess/quality/max_contamination", 5.0),
        keep_domain=cfg("preprocess/keep_domain", "Archaea"),
        formats=" ".join(PLOT_FORMATS),
    log:
        f"{LOGS}/preprocess/plot_preprocess.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/plot_preprocess.py \
            --annotated {input.annotated} --species {input.species} --report {input.report} \
            --min-completeness {params.min_comp} --max-contamination {params.max_cont} \
            --keep-domain "{params.keep_domain}" --formats {params.formats} \
            --outdir {params.outdir} > {log} 2>&1
        """
