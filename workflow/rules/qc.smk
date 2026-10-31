# =====================================================================
# Stage 0 — standardized gene calls + per-genome QC table.
# Quality, taxonomy, assembly stats, species dereplication, and the niche
# union are all done in Phase 1 (rules/preprocess.smk); the analysis runs on
# the species representatives that preprocessing staged. No dRep / extra QC here.
# =====================================================================

# Prodigal 2.6.3 refuses single-genome training below this many bases
# (main.c: MIN_SINGLE_GENOME 20000, exit code 10).
PRODIGAL_MIN_SINGLE = 20000


rule prodigal_genes:
    """Consistent gene calls (one coordinate space for all downstream context work).
    Contract v2: gene ids {genome}_{contig#}_{gene#} (identical in .faa/.fna/.gff) and the
    per-gene coordinate table genes.tsv. A crashed Prodigal, zero predicted genes or
    inconsistent Prodigal outputs fail the rule (an empty proteome would read as
    'Prometheus absent')."""
    input:
        genome=genome_fasta,
    output:
        faa=f"{RESULTS}/genes/{{genome}}/{{genome}}.faa",
        fna=f"{RESULTS}/genes/{{genome}}/{{genome}}.fna",
        gff=f"{RESULTS}/genes/{{genome}}/{{genome}}.gff",
        genes=f"{RESULTS}/genes/{{genome}}/{{genome}}.genes.tsv",
    params:
        mode=cfg("qc/prodigal_mode", "single"),
        table=cfg("qc/prodigal_translation_table", 11),
        min_single=PRODIGAL_MIN_SINGLE,
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
        mode="{params.mode}"
        case "$mode" in
            single|meta) ;;
            *) echo "ERROR: config qc.prodigal_mode must be 'single' or 'meta' (got '$mode')" > {log}; exit 1 ;;
        esac
        # Training length as Prodigal counts it (sequence.c read_seq_training): sequence
        # letters plus a 12-bp spacer per record when there is more than one record.
        # (Prodigal also counts the non-letters [\]^_` between 'Z' and 'a'; they never occur
        # in valid FASTA, and ignoring them can only choose meta, never Prodigal's exit 10.)
        bp=$(LC_ALL=C awk '/^>/ {{ n++; next }} {{ gsub(/[^A-Za-z]/, ""); L += length($0) }} END {{ if (n > 1) L += 12 * n; print L + 0 }}' {input.genome})
        if [ "$mode" = "single" ] && [ "$bp" -lt {params.min_single} ]; then
            mode=meta   # single mode cannot train on < {params.min_single} bp (Prodigal exits 10)
        fi
        echo "prodigal_genes {wildcards.genome}: training_length=$bp mode=$mode translation_table={params.table}" > {log}
        prodigal -v >> {log} 2>&1
        prodigal -i {input.genome} -a {output.faa}.raw -d {output.fna}.raw \
                 -f gff -o {output.gff}.raw -p $mode -g {params.table} -q 2>> {log}
        python workflow/scripts/py/rename_prodigal.py --genome {wildcards.genome} \
            --in-faa {output.faa}.raw --in-fna {output.fna}.raw --in-gff {output.gff}.raw \
            --faa {output.faa} --fna {output.fna} --gff {output.gff} \
            --genes-tsv {output.genes} >> {log} 2>&1
        rm -f {output.faa}.raw {output.fna}.raw {output.gff}.raw
        """


rule genome_metadata:
    """Per-genome QC table for the analysis, taken from preprocessing (completeness,
    contamination, BBMap assembly stats, quality tier). No QC is recomputed here:
    every analysis genome already passed the Phase-1 quality filter."""
    input:
        meta=rules.filter_collapse.output.genome_meta,
    output:
        meta=f"{RESULTS}/qc/genome_metadata.tsv",
    log:
        f"{LOGS}/qc/genome_metadata.log",
    threads: T_SMALL
    resources:
        mem_mb=2000, runtime=10,
    shell:
        r"""
        mkdir -p $(dirname {output.meta})
        cp {input.meta} {output.meta} 2> {log}
        echo "copied $(($(wc -l < {output.meta}) - 1)) genome row(s) from {input.meta}" >> {log}
        """


rule plot_qc:
    """Completeness-contamination scatter (by quality tier, with the configured qc
    thresholds) + assembly-statistic histograms of the analysis genomes."""
    input:
        meta=rules.genome_metadata.output.meta,
    output:
        expand(f"{RESULTS}/qc/plots/qc_completeness_contamination.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/qc/plots/qc_assembly_stats.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/qc/plots",
        comp_floor=cfg("qc/completeness_floor", 50.0),
        cont_ceil=cfg("qc/contamination_ceiling", 10.0),
        hi_comp=cfg("qc/tier_high_completeness", 90.0),
        hi_cont=cfg("qc/tier_high_contamination", 5.0),
        formats=" ".join(PLOT_FORMATS),
    log:
        f"{LOGS}/qc/plot_qc.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_qc.py --meta {input.meta} --outdir {params.outdir} \
            --completeness-floor {params.comp_floor} --contamination-ceiling {params.cont_ceil} \
            --tier-high-completeness {params.hi_comp} --tier-high-contamination {params.hi_cont} \
            --formats {params.formats} > {log} 2>&1
        """
