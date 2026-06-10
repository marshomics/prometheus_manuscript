# =====================================================================
# Stage 5 — molecular evolution: family + per-domain phylogenetics,
# reconciliation (DTL), ancestral gain/loss, selection (dN/dS).
# =====================================================================

rule family_align:
    """Align validated Prometheus proteins; trim with ClipKIT."""
    input:
        faa=rules.validate_prometheus.output.proteins,
    output:
        aln=f"{RESULTS}/evolution/align/family.aln",
        trim=f"{RESULTS}/evolution/align/family.trim.aln",
    log:
        f"{LOGS}/evolution/family_align.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p $(dirname {output.aln})
        n=$(grep -c '>' {input.faa} || echo 0)
        if [ "$n" -ge 3 ]; then
            mafft --auto --thread {threads} {input.faa} > {output.aln} 2> {log}
            clipkit {output.aln} -m kpic-smart-gap -o {output.trim} >> {log} 2>&1 || cp {output.aln} {output.trim}
        else
            cp {input.faa} {output.aln}; cp {input.faa} {output.trim}
        fi
        """


rule family_tree:
    """ML Prometheus gene tree (IQ-TREE2, ModelFinder, UFBoot2 + SH-aLRT)."""
    input:
        trim=rules.family_align.output.trim,
    output:
        tree=f"{RESULTS}/evolution/tree/family.treefile",
        ufboot=f"{RESULTS}/evolution/tree/family.ufboot",
    params:
        workdir=f"{RESULTS}/evolution/tree",
        model=cfg("evolution/gene_tree_model", "MFP"),
        ufboot=cfg("tree/ufboot", 1000),
        alrt=cfg("tree/alrt", 1000),
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/evolution/family_tree.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=720,
    shell:
        r"""
        mkdir -p {params.workdir}
        n=$(grep -c '>' {input.trim} || echo 0)
        if [ "$n" -ge 4 ]; then
            iqtree2 -s {input.trim} -m {params.model} -B {params.ufboot} -alrt {params.alrt} \
                    --wbtl -T {threads} --seed {params.seed} --prefix {params.workdir}/family -redo >> {log} 2>&1 || true
            cp {params.workdir}/family.treefile {output.tree} 2>/dev/null || echo "();" > {output.tree}
            cp {params.workdir}/family.ufboot {output.ufboot} 2>/dev/null || echo "();" > {output.ufboot}
        else
            echo "();" > {output.tree}; echo "();" > {output.ufboot}
        fi
        """


rule domain_trees:
    """Per-domain alignments + trees; Robinson–Foulds vs the family tree (domain shuffling)."""
    input:
        proteins=rules.validate_prometheus.output.proteins,
        domtbl=rules.prometheus_domain_scan.output.domtbl,
        family_tree=rules.family_tree.output.tree,
    output:
        rf=f"{RESULTS}/evolution/domains/rf_distances.tsv",
    params:
        domains_json=lambda wc: json_domains(),
        workdir=f"{RESULTS}/evolution/domains",
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/evolution/domain_trees.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=360,
    shell:
        r"""
        mkdir -p {params.workdir}
        python workflow/scripts/py/domain_trees.py \
            --proteins {input.proteins} --domtbl {input.domtbl} --family-tree {input.family_tree} \
            --domains '{params.domains_json}' --threads {threads} --seed {params.seed} \
            --workdir {params.workdir} --out {output.rf} > {log} 2>&1
        """


rule combined_cds:
    """Concatenate CDS (ffn) with genome-prefixed IDs (for codon alignment)."""
    input:
        fna=expand(f"{RESULTS}/genes/{{genome}}/{{genome}}.fna", genome=GENOMES),
        samples=config["samples"],
    output:
        ffn=f"{RESULTS}/evolution/all_cds.ffn",
        map=f"{RESULTS}/evolution/cds_map.tsv",
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
            --samples {input.samples} --genes-dir {params.genes_dir} --suffix fna \
            --out-faa {output.ffn} --out-map {output.map}
        """


rule codon_align:
    """Codon (in-frame) alignment from the protein alignment + CDS (PAL2NAL)."""
    input:
        prot_aln=rules.family_align.output.aln,
        cds=rules.combined_cds.output.ffn,
        proteins=rules.validate_prometheus.output.proteins,
    output:
        codon=f"{RESULTS}/evolution/align/family.codon.fasta",
    log:
        f"{LOGS}/evolution/codon_align.log",
    conda:
        "../envs/selection.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p $(dirname {output.codon})
        python workflow/scripts/py/codon_align.py \
            --prot-aln {input.prot_aln} --cds {input.cds} --proteins {input.proteins} \
            --out {output.codon} > {log} 2>&1
        """


rule reconciliation:
    """Species-tree-aware DTL reconciliation (GeneRax; optional ALE)."""
    input:
        aln=rules.family_align.output.trim,
        gene_tree=rules.family_tree.output.tree,
        ufboot=rules.family_tree.output.ufboot,
        species_tree=rules.species_tree.output.tree,
    output:
        events=f"{RESULTS}/evolution/reconciliation/reconciliation_events.tsv",
    params:
        workdir=f"{RESULTS}/evolution/reconciliation",
        model=cfg("evolution/reconciliation/generax_model", "UndatedDTL"),
        run_generax=str(RUN_GENERAX),
        run_ale=str(RUN_ALE),
    log:
        f"{LOGS}/evolution/reconciliation.log",
    conda:
        "../envs/reconciliation.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=720,
    shell:
        r"""
        mkdir -p {params.workdir}
        python workflow/scripts/py/run_reconciliation.py \
            --aln {input.aln} --gene-tree {input.gene_tree} --ufboot {input.ufboot} \
            --species-tree {input.species_tree} --model {params.model} \
            --run-generax {params.run_generax} --run-ale {params.run_ale} \
            --threads {threads} --workdir {params.workdir} --out {output.events} > {log} 2>&1
        """


rule ancestral_states:
    """Ancestral gain/loss of Prometheus across the archaeal tree (phytools stochastic mapping)."""
    input:
        tree=rules.species_tree.output.tree,
        master=rules.build_master_table.output.master,
    output:
        gl=f"{RESULTS}/evolution/ancestral/gains_losses.tsv",
        plot_png=f"{RESULTS}/evolution/plots/ancestral_gain_loss.png",
        plot_svg=f"{RESULTS}/evolution/plots/ancestral_gain_loss.svg",
    params:
        nsim=cfg("evolution/ancestral/nsim", 1000),
        seed=cfg("seed", 12345),
        outdir=f"{RESULTS}/evolution/ancestral",
        plotdir=f"{RESULTS}/evolution/plots",
    log:
        f"{LOGS}/evolution/ancestral_states.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=360,
    shell:
        r"""
        mkdir -p {params.outdir} {params.plotdir}
        Rscript workflow/scripts/R/ancestral_states.R \
            --tree {input.tree} --master {input.master} --nsim {params.nsim} --seed {params.seed} \
            --out {output.gl} --plot-prefix {params.plotdir}/ancestral_gain_loss > {log} 2>&1
        """


rule selection:
    """dN/dS on catalytically-intact loci (HyPhy FEL/MEME/BUSTED/aBSREL/RELAX)."""
    input:
        codon=rules.codon_align.output.codon,
        gene_tree=rules.family_tree.output.tree,
        catalog=rules.validate_prometheus.output.catalog,
        master=rules.build_master_table.output.master,
        domtbl=rules.prometheus_domain_scan.output.domtbl,
    output:
        summary=f"{RESULTS}/evolution/selection/selection_summary.tsv",
        sites=f"{RESULTS}/evolution/selection/site_table.tsv",
    params:
        run=str(RUN_SELECTION),
        methods=",".join(SELECTION_METHODS),
        relax_partition=cfg("evolution/selection/relax_partition", "methanogen"),
        intact_only=str(cfg("evolution/selection/intact_only", True)),
        site_pvalue=cfg("evolution/selection/site_pvalue", 0.1),
        domains_json=lambda wc: json_domains(),
        workdir=f"{RESULTS}/evolution/selection",
    log:
        f"{LOGS}/evolution/selection.log",
    conda:
        "../envs/selection.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=1440,
    shell:
        r"""
        mkdir -p {params.workdir}
        python workflow/scripts/py/run_selection.py \
            --codon {input.codon} --gene-tree {input.gene_tree} --catalog {input.catalog} \
            --master {input.master} --domtbl {input.domtbl} --domains '{params.domains_json}' \
            --methods {params.methods} --relax-partition {params.relax_partition} \
            --intact-only {params.intact_only} --site-pvalue {params.site_pvalue} --run {params.run} \
            --threads {threads} --workdir {params.workdir} \
            --out-summary {output.summary} --out-sites {output.sites} > {log} 2>&1
        """


rule plot_tanglegram:
    input:
        gene_tree=rules.family_tree.output.tree,
        species_tree=rules.species_tree.output.tree,
        master=rules.build_master_table.output.master,
    output:
        expand(f"{RESULTS}/evolution/plots/genetree_speciestree_tanglegram.{{ext}}", ext=PLOT_FORMATS),
    params:
        plotdir=f"{RESULTS}/evolution/plots",
    log:
        f"{LOGS}/evolution/tanglegram.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        mkdir -p {params.plotdir}
        Rscript workflow/scripts/R/plot_tanglegram.R \
            --gene-tree {input.gene_tree} --species-tree {input.species_tree} --master {input.master} \
            --plot-prefix {params.plotdir}/genetree_speciestree_tanglegram > {log} 2>&1
        """


rule plot_selection:
    input:
        summary=rules.selection.output.summary,
        sites=rules.selection.output.sites,
    output:
        expand(f"{RESULTS}/evolution/plots/selection_by_domain.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/evolution/plots",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_selection.py \
            --summary {input.summary} --sites {input.sites} --outdir {params.outdir}
        """
