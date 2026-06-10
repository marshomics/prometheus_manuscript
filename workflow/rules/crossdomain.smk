# =====================================================================
# Bacterial / cross-domain branch (OPT-IN: bacteria.enabled: true).
#
# Ingests a PRE-COMPUTED bacterial DefenseFinder set (combined systems/genes
# tables + per-genome faa/cds), validates bacterial Prometheus with the SAME
# scripts as archaea, then adds the cross-domain analyses: the combined
# phylogeny (origin + bacteria<->archaea transfer direction), matched prevalence
# and architecture, cross-domain dN/dS, and co-occurrence. Built with the
# `crossdomain` target; the archaeal `all`/`preprocess` runs are untouched and
# are consumed here only as read-only inputs (archaeal validated proteins + CDS).
# =====================================================================

BD = f"{RESULTS}/bacteria"
CD = f"{RESULTS}/crossdomain"

def _bact(key, default=""):
    return cfg(f"bacteria/{key}", default)


# --------------------------------------------------------------------- #
# 1. Ingest bacterial DefenseFinder -> Prometheus candidates (+ split systems)
# --------------------------------------------------------------------- #
rule bact_extract:
    input:
        systems=_bact("df_systems"),
        genes=_bact("df_genes"),
    output:
        candidates=f"{BD}/prometheus/candidates.faa",
        sources=f"{BD}/prometheus/candidate_sources.tsv",
        raw=f"{BD}/prometheus/raw_hits.tsv",
        cds=f"{BD}/prometheus/bacterial_prometheus_cds.ffn",
        split_done=f"{BD}/defense/.split.done",
    params:
        faa_dir=_bact("faa_dir"),
        cds_dir=_bact("cds_dir"),
        sys_gcol=_bact("systems_genome_col", "genome_id"),
        gene_gcol=_bact("genes_genome_col", "genome_id"),
        system=_bact("system", "Prometheus"),
        defense_dir=f"{BD}/defense",
    log: f"{LOGS}/bacteria/extract.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        python workflow/scripts/py/extract_bacterial_prometheus.py \
            --df-systems {input.systems} --df-genes {input.genes} \
            --faa-dir {params.faa_dir} --cds-dir "{params.cds_dir}" \
            --systems-genome-col {params.sys_gcol} --genes-genome-col {params.gene_gcol} \
            --system {params.system} --out-defense-dir {params.defense_dir} \
            --out-candidates {output.candidates} --out-sources {output.sources} \
            --out-raw {output.raw} --out-cds {output.cds} > {log} 2>&1
        [ -f {output.cds} ] || : > {output.cds}
        touch {output.split_done}
        """


# --------------------------------------------------------------------- #
# 2. Bacterial Prometheus validation (REUSES validate_prometheus.py)
# --------------------------------------------------------------------- #
rule bact_domain_scan:
    input:
        candidates=rules.bact_extract.output.candidates,
        pfam=cfg("db/pfam_hmm", "resources/db/pfam/Pfam-A.hmm"),
    output:
        domtbl=f"{BD}/prometheus/domains.domtbl",
    log: f"{LOGS}/bacteria/domains.log"
    conda: "../envs/annotation.yaml"
    threads: T_MED
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        mkdir -p $(dirname {output.domtbl})
        if [ -s {input.candidates} ] && [ -s {input.pfam} ]; then
            hmmscan --cpu {threads} --domtblout {output.domtbl} --cut_ga \
                    {input.pfam} {input.candidates} > {log} 2>&1 || true
        fi
        [ -s {output.domtbl} ] || printf "# no domains\n" > {output.domtbl}
        """


rule bact_validation_tree:
    input:
        candidates=rules.bact_extract.output.candidates,
        prom_ref=cfg("refs/prometheus_faa", "resources/refs/prometheus_reference.faa"),
        para_ref=cfg("refs/paralog_faa", "resources/refs/paralog_reference.faa"),
    output:
        tree=f"{BD}/prometheus/validation_tree/validation.treefile",
        labels=f"{BD}/prometheus/validation_tree/reference_labels.tsv",
    params:
        outdir=f"{BD}/prometheus/validation_tree",
    log: f"{LOGS}/bacteria/validation_tree.log"
    conda: "../envs/phylo.yaml"
    threads: T_MED
    resources: mem_mb=8000, runtime=240
    shell:
        r"""
        mkdir -p {params.outdir}
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


rule bact_validate:
    input:
        sources=rules.bact_extract.output.sources,
        raw=rules.bact_extract.output.raw,
        domtbl=rules.bact_domain_scan.output.domtbl,
        tree=rules.bact_validation_tree.output.tree,
        ref_labels=rules.bact_validation_tree.output.labels,
        candidates=rules.bact_extract.output.candidates,
        meta=_bact("qc"),
        split_done=rules.bact_extract.output.split_done,
    output:
        catalog=f"{BD}/prometheus/bacterial_prometheus_catalog.tsv",
        proteins=f"{BD}/prometheus/bacterial_prometheus_proteins.faa",
        sensitivity=f"{BD}/prometheus/bacterial_detection_sensitivity.tsv",
    params:
        domains_json=lambda wc: json_domains(),
        min_len=cfg("prometheus/min_len_aa", 800),
        max_len=cfg("prometheus/max_len_aa", 3000),
        df_dir=f"{BD}/defense",
        ctx=cfg("prometheus/context_window", 10),
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        python workflow/scripts/py/validate_prometheus.py \
            --sources {input.sources} --raw {input.raw} --domtbl {input.domtbl} \
            --tree {input.tree} --ref-labels {input.ref_labels} \
            --candidates {input.candidates} --meta {input.meta} --reps {input.meta} \
            --df-dir {params.df_dir} --context-window {params.ctx} \
            --domains '{params.domains_json}' --min-len {params.min_len} --max-len {params.max_len} \
            --out-catalog {output.catalog} --out-proteins {output.proteins} \
            --out-sensitivity {output.sensitivity}
        """


# --------------------------------------------------------------------- #
# 3. Bacterial table + matched prevalence + architecture (REUSED scripts)
# --------------------------------------------------------------------- #
rule bact_table:
    input:
        taxonomy=_bact("taxonomy"),
        qc=_bact("qc"),
        catalog=rules.bact_validate.output.catalog,
        split_done=rules.bact_extract.output.split_done,
    output:
        master=f"{BD}/bacterial_master_table.tsv",
    params:
        defense_dir=f"{BD}/defense",
        min_comp=_bact("min_completeness", 90.0),
        genomes_dir=_bact("genomes_dir", ""),
        genome_ext=_bact("genome_extension", "fasta"),
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/build_bacterial_table.py \
            --taxonomy {input.taxonomy} --qc {input.qc} --catalog {input.catalog} \
            --defense-dir {params.defense_dir} --min-completeness {params.min_comp} \
            --genomes-dir "{params.genomes_dir}" --genome-ext {params.genome_ext} \
            --out {output.master}
        """


rule bact_prevalence:
    input:
        master=rules.bact_table.output.master,
    output:
        prevalence=f"{BD}/distribution/bacterial_prevalence_by_taxon.tsv",
    params:
        tiers=",".join(cfg("stats/confidence_tiers_primary", ["A"])),
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/compute_prevalence.py \
            --master {input.master} --primary-tiers "{params.tiers}" --out {output.prevalence}
        """


rule bact_architecture:
    input:
        domtbl=rules.bact_domain_scan.output.domtbl,
        catalog=rules.bact_validate.output.catalog,
    output:
        typology=f"{BD}/diversity/bacterial_architecture_typology.tsv",
    params:
        domains_json=lambda wc: json_domains(),
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/architecture_typology.py \
            --domtbl {input.domtbl} --catalog {input.catalog} \
            --domains '{params.domains_json}' --out {output.typology}
        """


# --------------------------------------------------------------------- #
# 4. Bacterial species tree (prune GTDB bac120) + defense matrix + co-occurrence
# --------------------------------------------------------------------- #
rule bact_species_tree:
    input:
        qc=_bact("qc"),
    output:
        tree=f"{BD}/scaffolds/bacterial_species_tree.treefile",
    params:
        ref_tree=_bact("gtdb_bac_tree"),
    log: f"{LOGS}/bacteria/species_tree.log"
    conda: "../envs/phylo.yaml"
    threads: T_SMALL
    resources: mem_mb=16000, runtime=120
    shell:
        r"""
        python workflow/scripts/py/prune_tree_to_taxa.py \
            --in-tree "{params.ref_tree}" --keep {input.qc} --out-tree {output.tree} > {log} 2>&1
        """


rule bact_defense_matrix:
    input:
        master=rules.bact_table.output.master,
        split_done=rules.bact_extract.output.split_done,
    output:
        matrix=f"{BD}/distribution/bacterial_defense_matrix.tsv",
    params:
        df_dir=f"{BD}/defense",
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/build_defense_matrix.py \
            --master {input.master} --df-dir {params.df_dir} --out {output.matrix}
        """


rule bact_cooccurrence:
    input:
        matrix=rules.bact_defense_matrix.output.matrix,
        master=rules.bact_table.output.master,
        tree=rules.bact_species_tree.output.tree,
    output:
        cooc=f"{BD}/distribution/bacterial_cooccurrence.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
    log: f"{LOGS}/bacteria/cooccurrence.log"
    conda: "../envs/rstats.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        Rscript workflow/scripts/R/cooccurrence.R \
            --matrix {input.matrix} --master {input.master} --tree {input.tree} \
            --fdr {params.fdr} --seed {params.seed} --out {output.cooc} > {log} 2>&1
        """


# --------------------------------------------------------------------- #
# 5. Cross-domain phylogeny (archaeal + bacterial Prometheus + references)
# --------------------------------------------------------------------- #
rule crossdomain_align:
    input:
        arch=rules.validate_prometheus.output.proteins,     # archaeal (read-only)
        bact=rules.bact_validate.output.proteins,
        prom_ref=cfg("refs/prometheus_faa", "resources/refs/prometheus_reference.faa"),
        para_ref=cfg("refs/paralog_faa", "resources/refs/paralog_reference.faa"),
    output:
        aln=f"{CD}/align/crossdomain.aln",
        trim=f"{CD}/align/crossdomain.trim.aln",
        labels=f"{CD}/align/reference_labels.tsv",
    log: f"{LOGS}/crossdomain/align.log"
    conda: "../envs/phylo.yaml"
    threads: T_LARGE
    resources: mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=16000, per_attempt=16000), runtime=720
    shell:
        r"""
        mkdir -p $(dirname {output.aln})
        awk '/^>/{{print substr($1,2)"\tPrometheus_ref"}}' {input.prom_ref} > {output.labels}
        awk '/^>/{{print substr($1,2)"\tParalog_ref"}}'    {input.para_ref} >> {output.labels}
        cat {input.arch} {input.bact} {input.prom_ref} {input.para_ref} > {output.aln}.in.faa
        n=$(grep -c '>' {output.aln}.in.faa || echo 0)
        if [ "$n" -ge 4 ]; then
            mafft --auto --thread {threads} {output.aln}.in.faa > {output.aln} 2> {log}
            clipkit {output.aln} -m kpic-smart-gap -o {output.trim} >> {log} 2>&1 || cp {output.aln} {output.trim}
        else
            cp {output.aln}.in.faa {output.aln}; cp {output.aln}.in.faa {output.trim}
        fi
        """


rule crossdomain_tree:
    input:
        trim=rules.crossdomain_align.output.trim,
    output:
        tree=f"{CD}/tree/crossdomain.treefile",
    params:
        workdir=f"{CD}/tree",
        model=cfg("evolution/gene_tree_model", "MFP"),
        ufboot=cfg("tree/ufboot", 1000),
        alrt=cfg("tree/alrt", 1000),
        seed=cfg("seed", 12345),
    log: f"{LOGS}/crossdomain/tree.log"
    conda: "../envs/phylo.yaml"
    threads: T_LARGE
    resources: mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=16000, per_attempt=16000), runtime=1440
    shell:
        r"""
        mkdir -p {params.workdir}
        n=$(grep -c '>' {input.trim} || echo 0)
        if [ "$n" -ge 4 ]; then
            iqtree2 -s {input.trim} -m {params.model} -B {params.ufboot} -alrt {params.alrt} \
                    -T {threads} --seed {params.seed} --prefix {params.workdir}/crossdomain -redo >> {log} 2>&1 || true
            cp {params.workdir}/crossdomain.treefile {output.tree} 2>/dev/null || echo "();" > {output.tree}
        else
            echo "();" > {output.tree}
        fi
        """


rule crossdomain_analysis:
    input:
        tree=rules.crossdomain_tree.output.tree,
        arch=rules.validate_prometheus.output.proteins,
        bact=rules.bact_validate.output.proteins,
        labels=rules.crossdomain_align.output.labels,
    output:
        summary=f"{CD}/crossdomain_summary.tsv",
        tips=f"{CD}/crossdomain_tip_domains.tsv",
    conda: "../envs/phylo.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        python workflow/scripts/py/crossdomain_analysis.py \
            --tree {input.tree} --arch-proteins {input.arch} --bact-proteins {input.bact} \
            --ref-labels {input.labels} --out {output.summary} --out-tips {output.tips}
        """


rule plot_crossdomain:
    input:
        tree=rules.crossdomain_tree.output.tree,
        tips=rules.crossdomain_analysis.output.tips,
        summary=rules.crossdomain_analysis.output.summary,
    output:
        expand(f"{CD}/figures/crossdomain_tree.{{ext}}", ext=PLOT_FORMATS),
    params:
        prefix=f"{CD}/figures/crossdomain_tree",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        python workflow/scripts/py/plot_crossdomain_tree.py \
            --tree {input.tree} --tips {input.tips} --summary {input.summary} \
            --out-prefix {params.prefix}
        """


# --------------------------------------------------------------------- #
# 6. Cross-domain dN/dS (codon alignment of archaeal+bacterial Prometheus)
# --------------------------------------------------------------------- #
rule crossdomain_codon_align:
    input:
        prot_aln=rules.crossdomain_align.output.aln,
        arch_cds=rules.combined_cds.output.ffn,            # archaeal CDS (read-only)
        bact_cds=rules.bact_extract.output.cds,
        proteins=rules.bact_validate.output.proteins,
    output:
        codon=f"{CD}/selection/crossdomain.codon.fasta",
        cds=f"{CD}/selection/crossdomain_all_cds.ffn",
    log: f"{LOGS}/crossdomain/codon_align.log"
    conda: "../envs/selection.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        mkdir -p $(dirname {output.codon})
        cat {input.arch_cds} {input.bact_cds} > {output.cds}
        python workflow/scripts/py/codon_align.py \
            --prot-aln {input.prot_aln} --cds {output.cds} --proteins {input.proteins} \
            --out {output.codon} > {log} 2>&1
        """


rule crossdomain_selection:
    input:
        codon=rules.crossdomain_codon_align.output.codon,
        gene_tree=rules.crossdomain_tree.output.tree,
        arch_catalog=rules.validate_prometheus.output.catalog,
        bact_catalog=rules.bact_validate.output.catalog,
        master=rules.build_master_table.output.master,     # archaeal master (RELAX off here)
        arch_domtbl=rules.prometheus_domain_scan.output.domtbl,
        bact_domtbl=rules.bact_domain_scan.output.domtbl,
    output:
        summary=f"{CD}/selection/crossdomain_selection_summary.tsv",
        sites=f"{CD}/selection/crossdomain_site_table.tsv",
    params:
        run=str(RUN_SELECTION),
        methods="FEL,MEME,BUSTED,aBSREL",          # cross-domain: no methanogen RELAX partition
        intact_only=str(cfg("evolution/selection/intact_only", True)),
        site_pvalue=cfg("evolution/selection/site_pvalue", 0.1),
        domains_json=lambda wc: json_domains(),
        workdir=f"{CD}/selection",
    log: f"{LOGS}/crossdomain/selection.log"
    conda: "../envs/selection.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=1440
    shell:
        r"""
        mkdir -p {params.workdir}
        cat {input.arch_catalog} <(tail -n +2 {input.bact_catalog}) > {params.workdir}/combined_catalog.tsv
        cat {input.arch_domtbl} {input.bact_domtbl} > {params.workdir}/combined.domtbl
        python workflow/scripts/py/run_selection.py \
            --codon {input.codon} --gene-tree {input.gene_tree} \
            --catalog {params.workdir}/combined_catalog.tsv --master {input.master} \
            --domtbl {params.workdir}/combined.domtbl --domains '{params.domains_json}' \
            --methods {params.methods} --relax-partition none \
            --intact-only {params.intact_only} --site-pvalue {params.site_pvalue} --run {params.run} \
            --threads {threads} --workdir {params.workdir} \
            --out-summary {output.summary} --out-sites {output.sites} > {log} 2>&1
        """


rule plot_crossdomain_extra:
    """Bacteria vs archaea: prevalence (Wilson CI) + domain-architecture fractions."""
    input:
        arch_master=rules.build_master_table.output.master,
        bact_master=rules.bact_table.output.master,
        arch_typology=rules.architecture_typology.output.typology,
        bact_typology=rules.bact_architecture.output.typology,
    output:
        expand(f"{CD}/figures/crossdomain_comparison.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{CD}/figures",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_crossdomain_extra.py \
            --arch-master {input.arch_master} --bact-master {input.bact_master} \
            --arch-typology {input.arch_typology} --bact-typology {input.bact_typology} \
            --outdir {params.outdir}
        """


# --------------------------------------------------------------------- #
# 7. OPTIONAL: anaerobe / environment specificity control
#    Activates only when bacteria.niche_metadata is set (gated in crossdomain_targets).
#    Tests whether bacterial Prometheus tracks an oxygen/biome signal in the SAME
#    phyloglm framework as archaea -- a clean control because bacteria can't be
#    methanogens, so a niche signal here means the archaeal methanogen association may
#    be partly a niche effect. Reuses ecology_stats.R unchanged.
# --------------------------------------------------------------------- #
rule bact_niche:
    input:
        master=rules.bact_table.output.master,
        tree=rules.bact_species_tree.output.tree,
    output:
        assoc=f"{BD}/ecology/bacterial_niche_associations.tsv",
    params:
        niche=_bact("niche_metadata"),
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
        niche_master=f"{BD}/ecology/bacterial_niche_master.tsv",
    log: f"{LOGS}/bacteria/niche.log"
    conda: "../envs/rstats.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=180
    shell:
        r"""
        mkdir -p $(dirname {output.assoc})
        python workflow/scripts/py/build_bacterial_niche_master.py \
            --master {input.master} --niche "{params.niche}" --out {params.niche_master}
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {params.niche_master} --tree {input.tree} \
            --fdr {params.fdr} --seed {params.seed} --out {output.assoc} > {log} 2>&1
        """


# --------------------------------------------------------------------- #
# 8. Bacterial provirus / R-loop mechanism test
#    Activates only when bacteria.genomes_dir (assemblies) is set. The R-loop
#    targeting was demonstrated on Pseudomonas (a bacterium), so replicating the
#    proviral-R-loop-vs-Prometheus signal in bacteria validates the proxy the archaeal
#    side relies on. Reuses the archaeal geNomad/CheckV/rloop_predict/aggregate_viral/
#    ecology_stats scripts UNCHANGED, just on the bacterial assemblies + species tree.
# --------------------------------------------------------------------- #
rule bact_genomad:
    """Predict bacterial proviruses / viral contigs from the provided assemblies."""
    input:
        genome=bact_genome_fasta,
    output:
        virus_fna=f"{BD}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus.fna",
        virus_summary=f"{BD}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv",
    params:
        db=cfg("db/genomad", ""),
        outdir=f"{BD}/ecology/viral/{{genome}}",
    log: f"{LOGS}/bacteria/genomad/{{genome}}.log"
    conda: "../envs/viral.yaml"
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


rule bact_checkv:
    """Bacterial provirus completeness/quality."""
    input:
        virus_fna=rules.bact_genomad.output.virus_fna,
    output:
        quality=f"{BD}/ecology/viral/{{genome}}/checkv/quality_summary.tsv",
    params:
        db=cfg("db/checkv", ""),
        outdir=f"{BD}/ecology/viral/{{genome}}/checkv",
    log: f"{LOGS}/bacteria/checkv/{{genome}}.log"
    conda: "../envs/viral.yaml"
    threads: T_MED
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        mkdir -p {params.outdir}
        if [ -s {input.virus_fna} ]; then
            checkv end_to_end {input.virus_fna} {params.outdir} -d {params.db} -t {threads} > {log} 2>&1 || true
        fi
        [ -s {output.quality} ] || printf "contig_id\tcontig_length\tcompleteness\tcheckv_quality\n" > {output.quality}
        """


rule bact_rloop:
    """R-loop-forming-sequence prediction in bacterial proviruses (QmRLFS RIZ+REZ)."""
    input:
        virus_fna=rules.bact_genomad.output.virus_fna,
    output:
        rloop=f"{BD}/ecology/viral/{{genome}}/rloop.tsv",
    params:
        riz_min_g=cfg("ecology/rloop/riz_min_g", 3),
        window=cfg("ecology/rloop/window", 100),
        require_rez=str(cfg("ecology/rloop/require_rez", True)),
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/rloop_predict.py \
            --fasta {input.virus_fna} --genome {wildcards.genome} \
            --riz-min-g {params.riz_min_g} --window {params.window} \
            --require-rez {params.require_rez} --out {output.rloop}
        """


rule bact_viral_summary:
    """Aggregate bacterial proviral load + CheckV + R-loop per genome."""
    input:
        virus=expand(f"{BD}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv", genome=BACT_GENOMES),
        checkv=expand(f"{BD}/ecology/viral/{{genome}}/checkv/quality_summary.tsv", genome=BACT_GENOMES),
        rloop=expand(f"{BD}/ecology/viral/{{genome}}/rloop.tsv", genome=BACT_GENOMES),
        samples=_bact("qc"),
    output:
        summary=f"{BD}/ecology/viral/bacterial_viral_summary.tsv",
    params:
        viral_dir=f"{BD}/ecology/viral",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/aggregate_viral.py \
            --viral-dir {params.viral_dir} --samples {input.samples} --out {output.summary}
        """


rule bact_ecology_stats:
    """Bacterial Prometheus vs viral pressure / R-loop, phylogeny-controlled (the
    mechanism replication). ecology_stats.R adjusts the continuous viral predictors for
    completeness + genome size; no niche columns in the bacterial master, so only the
    viral-pressure terms are tested."""
    input:
        master=rules.bact_table.output.master,
        viral=rules.bact_viral_summary.output.summary,
        tree=rules.bact_species_tree.output.tree,
    output:
        assoc=f"{BD}/ecology/bacterial_viral_associations.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
    log: f"{LOGS}/bacteria/viral_stats.log"
    conda: "../envs/rstats.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        mkdir -p $(dirname {output.assoc})
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --viral {input.viral} --tree {input.tree} \
            --fdr {params.fdr} --seed {params.seed} --out {output.assoc} > {log} 2>&1
        """


rule plot_viral_crossdomain:
    """Mechanism figure: R-loop-prone proviral fraction + provirus load vs Prometheus
    presence in BOTH domains, and Prometheus prevalence vs R-loop pressure."""
    input:
        arch_master=rules.build_master_table.output.master,
        arch_viral=rules.viral_summary.output.summary,
        bact_master=rules.bact_table.output.master,
        bact_viral=rules.bact_viral_summary.output.summary,
        bact_assoc=rules.bact_ecology_stats.output.assoc,
    output:
        expand(f"{CD}/figures/viral_rloop_crossdomain.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{CD}/figures",
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_viral_crossdomain.py \
            --arch-master {input.arch_master} --arch-viral {input.arch_viral} \
            --bact-master {input.bact_master} --bact-viral {input.bact_viral} \
            --bact-assoc {input.bact_assoc} --outdir {params.outdir}
        """
