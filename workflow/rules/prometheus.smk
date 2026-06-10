# =====================================================================
# Stage 1 — Prometheus identification, validation, curation.
# Detection comes ONLY from DefenseFinder's Prometheus calls (no PADLOC,
# no custom HMM sweep) -> domain architecture -> phylogenetic separation
# from cellular paralogs -> defense-island context -> confidence tiers.
# =====================================================================

rule defensefinder_check:
    """Preflight: fail fast unless the detector can actually call Prometheus.

    Prometheus only exists in defense-finder-models v2.0.0+ (needs software >= 2.0.0);
    an older binary or model set detects zero Prometheus loci silently. This gate must
    pass before any genome is searched."""
    output:
        ok=f"{RESULTS}/defense/.defensefinder_ok",
    params:
        run=str(RUN_DEFENSEFINDER),
        system=cfg("prometheus/defensefinder_system", "Prometheus"),
        min_sw=cfg("prometheus/defensefinder_min_software", "2.0.0"),
        models_version=cfg("prometheus/defensefinder_models_version", "2.0.0"),
        provided=cfg("prometheus/defensefinder_provided_dir", "resources/defensefinder"),
    log:
        f"{LOGS}/defensefinder/preflight.log",
    conda:
        "../envs/defensefinder.yaml"
    threads: T_SMALL
    resources:
        mem_mb=2000, runtime=20,
    shell:
        r"""
        mkdir -p $(dirname {output.ok})
        python workflow/scripts/py/check_defensefinder.py \
            --run {params.run} --system {params.system} \
            --min-software {params.min_sw} --models-version {params.models_version} \
            --provided-dir {params.provided} --out {output.ok} > {log} 2>&1
        """


rule defensefinder:
    """Run (or ingest) DefenseFinder per genome under a pinned version."""
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
            defense-finder run -o {params.outdir} -w {threads} {input.faa} > {log} 2>&1 || true
        else
            cp {params.provided}/{wildcards.genome}/*_defense_finder_systems.tsv {output.systems} 2>/dev/null || true
            cp {params.provided}/{wildcards.genome}/*_defense_finder_genes.tsv   {output.genes}   2>/dev/null || true
        fi
        # Guarantee headers exist so empty genomes don't break aggregation.
        [ -s {output.systems} ] || printf "sys_id\ttype\tsubtype\tsys_beg\tsys_end\tprotein_in_syst\tname_of_profiles_in_sys\n" > {output.systems}
        [ -s {output.genes} ]   || printf "replicon\thit_id\tgene_name\thit_pos\tmodel_fqn\tsys_id\thit_i_eval\thit_score\thit_profile_cov\n" > {output.genes}
        """


rule combined_proteome:
    """Concatenate proteomes with genome-prefixed IDs (genome__orig_id) for candidate extraction + mcr search."""
    input:
        faa=expand(f"{RESULTS}/genes/{{genome}}/{{genome}}.faa", genome=GENOMES),
        samples=config["samples"],
    output:
        faa=f"{RESULTS}/prometheus/all_proteins.faa",
        map=f"{RESULTS}/prometheus/protein_map.tsv",
    params:
        genes_dir=f"{RESULTS}/genes",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/combine_proteomes.py \
            --samples {input.samples} --genes-dir {params.genes_dir} \
            --out-faa {output.faa} --out-map {output.map}
        """


rule assemble_candidates:
    """Extract Prometheus candidates from DefenseFinder's Prometheus system calls."""
    input:
        df_systems=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_systems.tsv", genome=GENOMES),
        df_genes=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_genes.tsv", genome=GENOMES),
        proteome=rules.combined_proteome.output.faa,
        pmap=rules.combined_proteome.output.map,
    output:
        candidates=f"{RESULTS}/prometheus/candidates.faa",
        sources=f"{RESULTS}/prometheus/candidate_sources.tsv",
        raw=f"{RESULTS}/prometheus/raw_hits.tsv",
    params:
        system=cfg("prometheus/defensefinder_system", "Prometheus"),
        df_dir=f"{RESULTS}/defense/defensefinder",
        strict=str(cfg("prometheus/strict_detection", True)),
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/parse_defensefinder.py \
            --df-dir {params.df_dir} --proteome {input.proteome} --pmap {input.pmap} \
            --system {params.system} --strict {params.strict} \
            --out-candidates {output.candidates} --out-sources {output.sources} --out-raw {output.raw}
        """


rule prometheus_domain_scan:
    """Domain architecture by HMM scan vs Pfam-A (names matched downstream; release-robust)."""
    input:
        candidates=rules.assemble_candidates.output.candidates,
        pfam=cfg("db/pfam_hmm", "resources/db/pfam/Pfam-A.hmm"),
    output:
        domtbl=f"{RESULTS}/prometheus/domains.domtbl",
    log:
        f"{LOGS}/prometheus/domains.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p $(dirname {output.domtbl})
        if [ -s {input.candidates} ] && [ -s {input.pfam} ]; then
            hmmscan --cpu {threads} --domtblout {output.domtbl} --cut_ga \
                    {input.pfam} {input.candidates} > {log} 2>&1 || true
        fi
        [ -s {output.domtbl} ] || printf "# no domains\n" > {output.domtbl}
        """


rule prometheus_validation_tree:
    """Place candidates against Prometheus + housekeeping DNA2/NurA/HerA references."""
    input:
        candidates=rules.assemble_candidates.output.candidates,
        prom_ref=cfg("refs/prometheus_faa", "resources/refs/prometheus_reference.faa"),
        para_ref=cfg("refs/paralog_faa", "resources/refs/paralog_reference.faa"),
    output:
        tree=f"{RESULTS}/prometheus/validation_tree/validation.treefile",
        labels=f"{RESULTS}/prometheus/validation_tree/reference_labels.tsv",
    params:
        outdir=f"{RESULTS}/prometheus/validation_tree",
    log:
        f"{LOGS}/prometheus/validation_tree.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=240,
    shell:
        r"""
        mkdir -p {params.outdir}
        # Label references by origin so the validator can assign each candidate to a clade.
        awk '/^>/{{print substr($1,2)"\tPrometheus_ref"}}' {input.prom_ref} > {output.labels}
        awk '/^>/{{print substr($1,2)"\tParalog_ref"}}'    {input.para_ref} >> {output.labels}
        cat {input.candidates} {input.prom_ref} {input.para_ref} > {params.outdir}/all.faa
        n=$(grep -c '>' {params.outdir}/all.faa || echo 0)
        if [ "$n" -ge 4 ]; then
            mafft --auto --thread {threads} {params.outdir}/all.faa > {params.outdir}/all.aln 2> {log}
            clipkit {params.outdir}/all.aln -m kpic-smart-gap -o {params.outdir}/all.trim.aln >> {log} 2>&1 || cp {params.outdir}/all.aln {params.outdir}/all.trim.aln
            iqtree2 -s {params.outdir}/all.trim.aln -m LG+G -fast -T {threads} \
                    --prefix {params.outdir}/validation -redo >> {log} 2>&1 || true
        fi
        [ -s {output.tree} ] || echo "();" > {output.tree}
        """


rule validate_prometheus:
    """Integrate all evidence into a confidence-tiered Prometheus catalog."""
    input:
        sources=rules.assemble_candidates.output.sources,
        raw=rules.assemble_candidates.output.raw,
        domtbl=rules.prometheus_domain_scan.output.domtbl,
        tree=rules.prometheus_validation_tree.output.tree,
        ref_labels=rules.prometheus_validation_tree.output.labels,
        candidates=rules.assemble_candidates.output.candidates,
        meta=rules.genome_metadata.output.meta,
        reps=rules.filter_collapse.output.reps_list,
        df_genes=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_genes.tsv", genome=GENOMES),
    output:
        catalog=f"{RESULTS}/prometheus/prometheus_catalog.tsv",
        proteins=f"{RESULTS}/prometheus/prometheus_proteins.faa",
        sensitivity=f"{RESULTS}/prometheus/detection_sensitivity.tsv",
    params:
        domains_json=lambda wc: json_domains(),
        min_len=cfg("prometheus/min_len_aa", 800),
        max_len=cfg("prometheus/max_len_aa", 3000),
        df_dir=f"{RESULTS}/defense/defensefinder",
        context_window=cfg("prometheus/context_window", 10),
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        python workflow/scripts/py/validate_prometheus.py \
            --sources {input.sources} --raw {input.raw} --domtbl {input.domtbl} \
            --tree {input.tree} --ref-labels {input.ref_labels} \
            --candidates {input.candidates} --meta {input.meta} --reps {input.reps} \
            --df-dir {params.df_dir} --context-window {params.context_window} \
            --domains '{params.domains_json}' --min-len {params.min_len} --max-len {params.max_len} \
            --out-catalog {output.catalog} --out-proteins {output.proteins} --out-sensitivity {output.sensitivity}
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
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_prometheus_validation.py \
            --catalog {input.catalog} --sensitivity {input.sensitivity} --outdir {params.outdir}
        """
