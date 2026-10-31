# =====================================================================
# Stage 1 — Prometheus identification, validation, curation.
# Detection comes ONLY from DefenseFinder's Prometheus calls (no PADLOC,
# no custom HMM sweep) -> domain architecture + catalytic motifs ->
# phylogenetic placement against Prometheus vs DNA2/NurA-HerA references ->
# +/-N-ORF defense-island context -> confidence tiers.
# Contract v2: contig-aware gene ids ({genome}_{contig#}_{gene#}), DefenseFinder with
# --db-type gembase, genes/{g}/{g}.genes.tsv, one exact-matching domain spec, and no
# silent failure (a crashed tool fails its rule; placeholders only for empty input).
# =====================================================================

DF_SYSTEM = cfg("prometheus/defensefinder_system", "Prometheus")
# Pinned detector (must equal workflow/envs/defensefinder.yaml) and models.
DF_SOFTWARE_PIN = cfg("prometheus/defensefinder_version", "3.0.0")
DF_MODELS_PIN = cfg("prometheus/defensefinder_models_pin", "3.1.0")
DF_MODELS_DIR = cfg("prometheus/defensefinder_models_dir", "")
# CasFinder 3.1.1 definitions (grammar vers 2.1) do not parse with macsyfinder 2.1.4.
DF_CASFINDER_PIN = cfg("prometheus/defensefinder_casfinder_pin", "3.1.0")
# DefenseFinder's own header lines for a genuinely empty proteome
# (defense_finder_posttreat/df_systems.py; macsypy TsvSystemSerializer + type/subtype/activity).
DF_SYSTEMS_HEADER = "sys_id\ttype\tsubtype\tactivity\tsys_beg\tsys_end\tprotein_in_syst\tgenes_count\tname_of_profiles_in_sys"
DF_GENES_HEADER = ("replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\tsys_loci\tlocus_num\t"
                   "sys_wholeness\tsys_score\tsys_occ\thit_gene_ref\thit_status\thit_seq_len\t"
                   "hit_i_eval\thit_score\thit_profile_cov\thit_seq_cov\thit_begin_match\t"
                   "hit_end_match\tcounterpart\tused_in\ttype\tsubtype\tactivity")


def _df_models_opt():
    return f"--models-dir {DF_MODELS_DIR}" if DF_MODELS_DIR else ""


rule defensefinder_check:
    """Preflight: fail fast unless the pinned detector can actually call Prometheus.

    Verifies the defense-finder software version (exact pin), the installed
    defense-finder-models (exact pin, >= 2.0.0) and CasFinder packages, the Prometheus
    model definition (mandatory gene Prometheus__ProA) and its profile, and -- as a
    positive control -- that DefenseFinder calls Prometheus on the consensus of its own
    ProA profile. Writes a provenance sentinel (versions, md5s)."""
    output:
        ok=f"{RESULTS}/defense/.defensefinder_ok",
    params:
        run=str(RUN_DEFENSEFINDER),
        system=DF_SYSTEM,
        min_sw=cfg("prometheus/defensefinder_min_software", "2.0.0"),
        models_version=cfg("prometheus/defensefinder_models_version", "2.0.0"),
        software_pin=DF_SOFTWARE_PIN,
        models_pin=DF_MODELS_PIN,
        casfinder_pin=DF_CASFINDER_PIN,
        models_dir=DF_MODELS_DIR,
        positive=str(cfg("prometheus/defensefinder_positive_control", True)),
        provided=cfg("prometheus/defensefinder_provided_dir", "resources/defensefinder"),
        workdir=f"{RESULTS}/defense/.defensefinder_positive_control",
    log:
        f"{LOGS}/defensefinder/preflight.log",
    conda:
        "../envs/defensefinder.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=60,
    shell:
        r"""
        mkdir -p $(dirname {output.ok})
        python workflow/scripts/py/check_defensefinder.py \
            --run {params.run} --system {params.system} \
            --min-software {params.min_sw} --models-version {params.models_version} \
            --software-version "{params.software_pin}" --models-pin "{params.models_pin}" \
            --casfinder-pin "{params.casfinder_pin}" \
            --models-dir "{params.models_dir}" --positive-control {params.positive} \
            --workdir {params.workdir} --threads {threads} \
            --provided-dir {params.provided} --out {output.ok} > {log} 2>&1
        """


rule defensefinder:
    """Run (or ingest) DefenseFinder per genome under the pinned version, then verify the
    tables (columns, ids present in this genome's proteome, v2+ model paths, gembase
    replicons). --db-type gembase: replicon = contig, so no system spans a contig junction.
    --skip-model-version-check: the preflight pins/verifies the models; without it
    DefenseFinder queries GitHub (fails on offline nodes)."""
    input:
        faa=f"{RESULTS}/genes/{{genome}}/{{genome}}.faa",
        models_ok=rules.defensefinder_check.output.ok,
    output:
        systems=f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_systems.tsv",
        genes=f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_genes.tsv",
    params:
        run=str(RUN_DEFENSEFINDER),
        provided=cfg("prometheus/defensefinder_provided_dir", "resources/defensefinder"),
        outdir=f"{RESULTS}/defense/defensefinder/{{genome}}",
        models_opt=_df_models_opt(),
        sys_header=DF_SYSTEMS_HEADER,
        genes_header=DF_GENES_HEADER,
    log:
        f"{LOGS}/defensefinder/{{genome}}.log",
    conda:
        "../envs/defensefinder.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p {params.outdir}
        if [ "{params.run}" = "True" ]; then
            if [ -s {input.faa} ]; then
                defense-finder run --db-type gembase --skip-model-version-check {params.models_opt} \
                    --index-dir {params.outdir}/index -w {threads} -o {params.outdir} \
                    {input.faa} > {log} 2>&1
                gembase="--expect-gembase"
            else
                # empty proteome (no genes): DefenseFinder's own empty-table headers
                printf '%s\n' "{params.sys_header}" > {output.systems}  # nonfatal: empty proteome, no genes to search
                printf '%s\n' "{params.genes_header}" > {output.genes}  # nonfatal: empty proteome, no genes to search
                echo "empty proteome {input.faa}: DefenseFinder not run" > {log}
                gembase=""
            fi
        else
            python workflow/scripts/py/check_defensefinder.py --ingest-copy \
                --provided-dir {params.provided} --genome {wildcards.genome} \
                --systems {output.systems} --genes {output.genes} > {log} 2>&1
            gembase=""
        fi
        python workflow/scripts/py/check_defensefinder.py --check-outputs $gembase \
            --systems {output.systems} --genes {output.genes} --faa {input.faa} \
            --genome {wildcards.genome} >> {log} 2>&1
        """


rule combined_proteome:
    """Concatenate proteomes with genome-prefixed IDs (genome__orig_id) for candidate
    extraction + mcr search. A missing per-genome proteome is an error."""
    input:
        faa=expand(f"{RESULTS}/genes/{{genome}}/{{genome}}.faa", genome=GENOMES),
        samples=config["samples"],
    output:
        faa=f"{RESULTS}/prometheus/all_proteins.faa",
        map=f"{RESULTS}/prometheus/protein_map.tsv",
    params:
        genes_dir=f"{RESULTS}/genes",
    log:
        f"{LOGS}/prometheus/combined_proteome.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/combine_proteomes.py \
            --samples {input.samples} --genes-dir {params.genes_dir} \
            --out-faa {output.faa} --out-map {output.map} > {log} 2>&1
        """


rule assemble_candidates:
    """Proteins of Prometheus__ProA genes in DefenseFinder systems whose type/subtype is
    exactly the configured system; genomes from the sample table only."""
    input:
        df_systems=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_systems.tsv", genome=GENOMES),
        df_genes=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_genes.tsv", genome=GENOMES),
        proteome=rules.combined_proteome.output.faa,
        pmap=rules.combined_proteome.output.map,
        samples=config["samples"],
    output:
        candidates=f"{RESULTS}/prometheus/candidates.faa",
        sources=f"{RESULTS}/prometheus/candidate_sources.tsv",
        raw=f"{RESULTS}/prometheus/raw_hits.tsv",
    params:
        system=DF_SYSTEM,
        gene=cfg("prometheus/defensefinder_gene", f"{DF_SYSTEM}__ProA"),
        df_dir=f"{RESULTS}/defense/defensefinder",
        strict=str(cfg("prometheus/strict_detection", True)),
    log:
        f"{LOGS}/prometheus/assemble_candidates.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/parse_defensefinder.py \
            --df-dir {params.df_dir} --proteome {input.proteome} --pmap {input.pmap} \
            --samples {input.samples} --system {params.system} --gene "{params.gene}" \
            --strict {params.strict} \
            --out-candidates {output.candidates} --out-sources {output.sources} \
            --out-raw {output.raw} > {log} 2>&1
        """


rule prometheus_domain_scan:
    """Domain architecture by hmmscan --cut_ga vs pressed Pfam-A. The main output (-o,
    --notextw) keeps the per-domain alignments used for residue-level catalytic motifs.
    An unpressed/missing Pfam or a crashed hmmscan fails the rule."""
    input:
        candidates=rules.assemble_candidates.output.candidates,
        pfam=cfg("db/pfam_hmm", "resources/db/pfam/Pfam-A.hmm"),
    output:
        domtbl=f"{RESULTS}/prometheus/domains.domtbl",
        hmmout=f"{RESULTS}/prometheus/domains.hmmscan.txt",
    log:
        f"{LOGS}/prometheus/domains.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=240,
    shell:
        r"""
        mkdir -p $(dirname {output.domtbl})
        if [ -s {input.candidates} ]; then
            for ext in h3m h3i h3f h3p; do
                [ -s {input.pfam}.$ext ] || {{ echo "ERROR: {input.pfam}.$ext missing: run 'hmmpress {input.pfam}'" > {log}; exit 1; }}
            done
            hmmscan --cpu {threads} --cut_ga --notextw -o {output.hmmout} \
                --domtblout {output.domtbl} {input.pfam} {input.candidates} > {log} 2>&1
        else
            # zero candidates (genuinely empty input): comment-only domtbl, empty alignments
            printf '# no candidates\n' > {output.domtbl}  # nonfatal: zero candidates, nothing to scan
            : > {output.hmmout}  # nonfatal: zero candidates, nothing to scan
            echo "no candidates: hmmscan not run" > {log}
        fi
        """


rule prometheus_validation_tree:
    """Place candidates against Prometheus + housekeeping DNA2/NurA-HerA references.
    Reference ids are sanitised once (combine_proteomes.py reference mode) so the
    alignment, tree tips and reference_labels.tsv agree. MAFFT -> ClipKIT
    (kpic-smart-gap) -> IQ-TREE with UFBoot2 + SH-aLRT (+ --bnni) and a fixed seed; the
    supports drive validate_prometheus' supported-clade placement."""
    input:
        candidates=rules.assemble_candidates.output.candidates,
        prom_ref=cfg("refs/prometheus_faa", "resources/refs/prometheus_reference.faa"),
        para_ref=cfg("refs/paralog_faa", "resources/refs/paralog_reference.faa"),
    output:
        tree=f"{RESULTS}/prometheus/validation_tree/validation.treefile",
        labels=f"{RESULTS}/prometheus/validation_tree/reference_labels.tsv",
    params:
        outdir=f"{RESULTS}/prometheus/validation_tree",
        model=cfg("prometheus/validation_tree_model", "MFP"),
        ufboot=cfg("tree/ufboot", 1000),
        alrt=cfg("tree/alrt", 1000),
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/prometheus/validation_tree.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=1440,
    shell:
        r"""
        mkdir -p {params.outdir}
        python workflow/scripts/py/combine_proteomes.py --candidates {input.candidates} \
            --ref Prometheus_ref={input.prom_ref} --ref Paralog_ref={input.para_ref} \
            --out-faa {params.outdir}/all.faa --out-labels {output.labels} \
            --out-map {params.outdir}/reference_map.tsv > {log} 2>&1
        n_cand=$(awk '/^>/ {{ n++ }} END {{ print n + 0 }}' {input.candidates})
        n_all=$(awk '/^>/ {{ n++ }} END {{ print n + 0 }}' {params.outdir}/all.faa)
        if [ "$n_cand" -eq 0 ] || [ "$n_all" -lt 4 ]; then
            # no candidates, or < 4 sequences (no testable internal edge): empty-input tree
            echo "();" > {output.tree}  # nonfatal: no candidates or fewer than 4 sequences, no testable edge
            echo "placeholder tree: $n_cand candidate(s), $n_all sequence(s)" >> {log}
            exit 0
        fi
        # --threadit 0: MAFFT's multithreaded iterative refinement differs run to run (MAFFT multithreading docs); this keeps it reproducible
        mafft --auto --thread {threads} --threadit 0 {params.outdir}/all.faa > {params.outdir}/all.aln 2>> {log}
        [ -s {params.outdir}/all.aln ] || {{ echo "ERROR: MAFFT wrote no alignment" >> {log}; exit 1; }}
        rm -f {params.outdir}/all.trim.aln
        clipkit {params.outdir}/all.aln -m kpic-smart-gap -o {params.outdir}/all.trim.aln >> {log} 2>&1
        [ -s {params.outdir}/all.trim.aln ] || {{ echo "ERROR: ClipKIT wrote no alignment (it exits 0 on errors)" >> {log}; exit 1; }}
        iqtree2 -s {params.outdir}/all.trim.aln -m {params.model} -B {params.ufboot} \
            --alrt {params.alrt} --bnni -T {threads} --seed {params.seed} \
            --prefix {params.outdir}/validation -redo >> {log} 2>&1
        [ -s {output.tree} ] || {{ echo "ERROR: IQ-TREE wrote no tree" >> {log}; exit 1; }}
        """


rule validate_prometheus:
    """Integrate all evidence into a confidence-tiered Prometheus catalog."""
    input:
        sources=rules.assemble_candidates.output.sources,
        raw=rules.assemble_candidates.output.raw,
        domtbl=rules.prometheus_domain_scan.output.domtbl,
        hmmout=rules.prometheus_domain_scan.output.hmmout,
        tree=rules.prometheus_validation_tree.output.tree,
        ref_labels=rules.prometheus_validation_tree.output.labels,
        candidates=rules.assemble_candidates.output.candidates,
        meta=rules.genome_metadata.output.meta,
        reps=rules.filter_collapse.output.reps_list,
        df_genes=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_genes.tsv", genome=GENOMES),
        genes_tsv=expand(f"{RESULTS}/genes/{{genome}}/{{genome}}.genes.tsv", genome=GENOMES),
    output:
        catalog=f"{RESULTS}/prometheus/prometheus_catalog.tsv",
        proteins=f"{RESULTS}/prometheus/prometheus_proteins.faa",
        sensitivity=f"{RESULTS}/prometheus/detection_sensitivity.tsv",
    params:
        domains_json=lambda wc: json_domains(),
        min_len=cfg("prometheus/min_len_aa", 800),
        max_len=cfg("prometheus/max_len_aa", 3000),
        df_dir=f"{RESULTS}/defense/defensefinder",
        genes_dir=f"{RESULTS}/genes",
        context_window=cfg("prometheus/context_window", 10),
        min_alrt=cfg("prometheus/placement_min_alrt", 80),
        min_ufboot=cfg("prometheus/placement_min_ufboot", 95),
        system=DF_SYSTEM,
    log:
        f"{LOGS}/prometheus/validate_prometheus.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/validate_prometheus.py \
            --sources {input.sources} --raw {input.raw} --domtbl {input.domtbl} \
            --hmmscan-out {input.hmmout} \
            --tree {input.tree} --ref-labels {input.ref_labels} \
            --candidates {input.candidates} --meta {input.meta} --reps {input.reps} \
            --df-dir {params.df_dir} --genes-dir {params.genes_dir} \
            --context-window {params.context_window} --system {params.system} \
            --min-alrt {params.min_alrt} --min-ufboot {params.min_ufboot} \
            --domains '{params.domains_json}' --min-len {params.min_len} --max-len {params.max_len} \
            --out-catalog {output.catalog} --out-proteins {output.proteins} \
            --out-sensitivity {output.sensitivity} > {log} 2>&1
        """


rule plot_prometheus:
    input:
        catalog=rules.validate_prometheus.output.catalog,
        sensitivity=rules.validate_prometheus.output.sensitivity,
    output:
        expand(f"{RESULTS}/prometheus/plots/prometheus_detection_summary.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{RESULTS}/prometheus/plots/prometheus_architecture.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/prometheus/plots",
    log:
        f"{LOGS}/prometheus/plot_prometheus.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_prometheus_validation.py \
            --catalog {input.catalog} --sensitivity {input.sensitivity} \
            --outdir {params.outdir} > {log} 2>&1
        """
