# =====================================================================
# Stage 3 — reference scaffolds: GTDB taxonomy, archaeal species tree,
# gene-based methanogen classification, ecology harmonization, master table.
# =====================================================================

rule gtdbtk_classify:
    """Standardized GTDB taxonomy for every genome."""
    input:
        genomes=expand(os.path.join(GENOME_DIR, "{genome}." + EXT), genome=GENOMES),
    output:
        summary=f"{RESULTS}/scaffolds/taxonomy_gtdbtk.tsv",
    params:
        gdir=GENOME_DIR,
        ext=EXT,
        data=cfg("db/gtdbtk_data", ""),
        outdir=f"{RESULTS}/scaffolds/gtdbtk",
        # Version-specific ANI-screen handling. 2.4+ uses skani; `--skip_ani_screen`
        # avoids needing the screen DB. `--mash_db` was removed in 2.4 (do NOT pass it).
        # Override here if your GTDB-Tk version differs.
        extra=cfg("tree/gtdbtk_classify_extra", "--skip_ani_screen"),
    log:
        f"{LOGS}/gtdbtk/classify.log",
    conda:
        "../envs/gtdbtk.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=64000, per_attempt=32000),
        runtime=1440,
    shell:
        r"""
        mkdir -p {params.outdir} $(dirname {output.summary})
        export GTDBTK_DATA_PATH={params.data}
        gtdbtk classify_wf --genome_dir {params.gdir} -x {params.ext} \
            --out_dir {params.outdir} --cpus {threads} {params.extra} > {log} 2>&1 || true
        # The summary lands in out_dir/ or out_dir/classify/ depending on version;
        # locate it wherever it is and normalise to a stable path.
        s=$(find {params.outdir} -name '*ar53.summary.tsv' 2>/dev/null | head -n1)
        if [ -n "$s" ]; then cp "$s" {output.summary}; \
        else printf "user_genome\tclassification\n" > {output.summary}; fi
        """


rule gtdbtk_align:
    """Concatenated ar53 marker MSA for the de novo focal tree."""
    input:
        genomes=expand(os.path.join(GENOME_DIR, "{genome}." + EXT), genome=GENOMES),
    output:
        msa=f"{RESULTS}/scaffolds/gtdbtk_align/user_msa.fasta.gz",
    params:
        gdir=GENOME_DIR,
        ext=EXT,
        data=cfg("db/gtdbtk_data", ""),
        outdir=f"{RESULTS}/scaffolds/gtdbtk_align",
        scope=cfg("tree/scope", "focal"),
    log:
        f"{LOGS}/gtdbtk/align.log",
    conda:
        "../envs/gtdbtk.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=64000, per_attempt=32000),
        runtime=1440,
    shell:
        r"""
        mkdir -p {params.outdir}
        export GTDBTK_DATA_PATH={params.data}
        skip=""
        [ "{params.scope}" = "focal" ] && skip="--skip_gtdb_refs"
        gtdbtk identify --genome_dir {params.gdir} -x {params.ext} \
            --out_dir {params.outdir}/identify --cpus {threads} > {log} 2>&1 || true
        gtdbtk align --identify_dir {params.outdir}/identify \
            --out_dir {params.outdir} --cpus {threads} $skip >> {log} 2>&1 || true
        # Locate the user MSA (align/ subdir, name varies by version) and normalise.
        m=$(find {params.outdir} -name '*user_msa.fasta.gz' 2>/dev/null | head -n1)
        if [ -n "$m" ]; then cp "$m" {output.msa}; else : | gzip > {output.msa}; fi
        """


rule species_tree:
    """ML archaeal species tree (IQ-TREE2, two-step PMSF for deep-node robustness)."""
    input:
        msa=rules.gtdbtk_align.output.msa,
        taxonomy=rules.gtdbtk_classify.output.summary,
    output:
        tree=f"{RESULTS}/scaffolds/species_tree.treefile",
    params:
        workdir=f"{RESULTS}/scaffolds/species_tree",
        pmsf=str(cfg("tree/pmsf", True)),
        pmsf_model=cfg("tree/pmsf_model", "LG+C20+F+G"),
        model=cfg("tree/iqtree_model", "MFP"),
        ufboot=cfg("tree/ufboot", 1000),
        alrt=cfg("tree/alrt", 1000),
        seed=cfg("seed", 12345),
        rooting_method=cfg("tree/rooting_method", "outgroup"),
        outgroup=",".join(cfg("tree/outgroup", []) or []),
        outgroup_taxon=cfg("tree/outgroup_taxon", ""),
    log:
        f"{LOGS}/scaffolds/species_tree.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=32000, per_attempt=16000),
        runtime=2880,
    shell:
        r"""
        mkdir -p {params.workdir}
        msa={params.workdir}/msa.faa
        gunzip -c {input.msa} > $msa 2>/dev/null || cp {input.msa} $msa
        n=$(grep -c '>' $msa || echo 0)
        unrooted={params.workdir}/species.unrooted.treefile
        if [ "$n" -ge 4 ]; then
            if [ "{params.pmsf}" = "True" ]; then
                iqtree2 -s $msa -m LG+G -T {threads} -fast --seed {params.seed} \
                        --prefix {params.workdir}/guide -redo >> {log} 2>&1 || true
                iqtree2 -s $msa -m "{params.pmsf_model}" -ft {params.workdir}/guide.treefile \
                        -B {params.ufboot} -alrt {params.alrt} -T {threads} --seed {params.seed} \
                        --prefix {params.workdir}/species -redo >> {log} 2>&1 || true
            else
                iqtree2 -s $msa -m {params.model} -B {params.ufboot} -alrt {params.alrt} \
                        -T {threads} --seed {params.seed} --prefix {params.workdir}/species -redo >> {log} 2>&1 || true
            fi
            cp {params.workdir}/species.treefile $unrooted 2>/dev/null || echo "();" > $unrooted
        else
            echo "();" > $unrooted
        fi
        # Root at a defensible position (outgroup preferred; midpoint fallback) before
        # any comparative / ancestral use. The canonical species_tree.treefile is ROOTED.
        python workflow/scripts/py/root_tree.py \
            --in-tree $unrooted --taxonomy {input.taxonomy} \
            --method {params.rooting_method} --outgroup "{params.outgroup}" \
            --outgroup-taxon "{params.outgroup_taxon}" \
            --out-tree {output.tree} --log {params.workdir}/rooting_method.txt >> {log} 2>&1
        """


rule mcr_search:
    """Detect mcrABG to define methanogens from genes (not taxonomy)."""
    input:
        proteome=rules.combined_proteome.output.faa,
        pmap=rules.combined_proteome.output.map,
    output:
        hits=f"{RESULTS}/scaffolds/metabolism/mcr_hits.tsv",
    params:
        hmmdir=cfg("refs/mcr_hmm_dir", "resources/refs/mcr_hmm"),
    log:
        f"{LOGS}/scaffolds/mcr_search.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p $(dirname {output.hits})
        python workflow/scripts/py/mcr_search.py \
            --proteome {input.proteome} --pmap {input.pmap} --hmmdir {params.hmmdir} \
            --threads {threads} --out {output.hits} > {log} 2>&1
        """


rule methanogen_classify:
    """Methanogen vs ANME via mcrA presence + McrA phylogenetic placement."""
    input:
        hits=rules.mcr_search.output.hits,
        proteome=rules.combined_proteome.output.faa,
        ref_aln=cfg("refs/mcrA_ref_aln", "resources/refs/mcrA_reference.afa"),
        ref_labels=cfg("refs/mcrA_ref_labels", "resources/refs/mcrA_reference_labels.tsv"),
    output:
        calls=f"{RESULTS}/scaffolds/metabolism/methanogen_calls.tsv",
    params:
        required=",".join(cfg("methanogen/mcr_required", ["K00399"])),
        anme_tree=str(cfg("methanogen/anme_via_mcrA_tree", True)),
        workdir=f"{RESULTS}/scaffolds/metabolism/mcrA_placement",
    log:
        f"{LOGS}/scaffolds/methanogen_classify.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=180,
    shell:
        r"""
        mkdir -p {params.workdir}
        python workflow/scripts/py/classify_methanogen.py \
            --hits {input.hits} --proteome {input.proteome} \
            --ref-aln {input.ref_aln} --ref-labels {input.ref_labels} \
            --required {params.required} --anme-tree {params.anme_tree} \
            --threads {threads} --workdir {params.workdir} --out {output.calls} > {log} 2>&1
        """


rule build_master_table:
    """Join QC, taxonomy, tree tip, methanogen status, ecology, Prometheus, defense load."""
    input:
        meta=rules.genome_metadata.output.meta,
        taxonomy=rules.gtdbtk_classify.output.summary,
        tree=rules.species_tree.output.tree,
        methanogen=rules.methanogen_classify.output.calls,
        ecology=config["metadata"],
        catalog=rules.validate_prometheus.output.catalog,
        sensitivity=rules.validate_prometheus.output.sensitivity,
        reps=rules.dereplicate.output.reps,
        df_systems=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_systems.tsv", genome=GENOMES),
    output:
        master=f"{RESULTS}/scaffolds/master_table.tsv",
        ecology=f"{RESULTS}/scaffolds/ecology_harmonized.tsv",
    params:
        df_dir=f"{RESULTS}/defense/defensefinder",
        envo_map=cfg("ecology/envo_map", "config/envo_map.tsv"),
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/build_master_table.py \
            --meta {input.meta} --taxonomy {input.taxonomy} --tree {input.tree} \
            --methanogen {input.methanogen} --ecology {input.ecology} --envo-map {params.envo_map} \
            --catalog {input.catalog} --sensitivity {input.sensitivity} --reps {input.reps} \
            --df-dir {params.df_dir} \
            --out-master {output.master} --out-ecology {output.ecology}
        """
