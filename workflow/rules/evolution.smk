# =====================================================================
# Stage 5 — molecular evolution: family + per-domain phylogenetics,
# reconciliation (DTL), ancestral gain/loss, selection (dN/dS).
#
# Config keys read here (cfg(path, default); defaults apply until the key exists):
#   evolution.align_tool          mafft | famsa          (mafft: L-INS-i up to linsi_max_seqs,
#   evolution.linsi_max_seqs      200                     FAMSA above it; famsa: FAMSA always)
#   evolution.trim_tool           clipkit | none ;  evolution.clipkit_mode  kpic-smart-gap
#   evolution.gene_tree_model     MFP ; evolution.domain_tree_model (default = gene_tree_model)
#   evolution.domain_min_taxa     6
#   evolution.codon_table         11 (NCBI; archaeal/bacterial initiators GTG/TTG)
#   evolution.max_cds_drop_frac   0.05 (codon_align fails above it)
#   evolution.reconciliation.{generax_model UndatedDTL, generax_strategy SPR,
#       subst_model LG+G4, ale_samples 100}
#   evolution.ancestral.{nsim 1000, q_method empirical, root_prior fitzjohn, nnull 1000,
#       prevalence_tol 0.1}; stats.confidence_tiers_primary / _sensitivity
#   evolution.selection.{relax_labelling clade, gard false, gard_mode Faster}
# =====================================================================

_EVO_LINSI_MAX = int(cfg("evolution/linsi_max_seqs", 200))
_EVO_TIERS = ",".join(cfg("stats/confidence_tiers_primary", ["A"]) or ["A"])
_EVO_TIERS_SENS = ",".join(cfg("stats/confidence_tiers_sensitivity", ["A", "B"]) or [])


rule family_align:
    """Align validated Prometheus proteins (MAFFT L-INS-i, or FAMSA for large sets); ClipKIT."""
    input:
        faa=rules.validate_prometheus.output.proteins,
    output:
        aln=f"{RESULTS}/evolution/align/family.aln",
        trim=f"{RESULTS}/evolution/align/family.trim.aln",
        method=f"{RESULTS}/evolution/align/family.align_method.txt",
    params:
        tool=cfg("evolution/align_tool", "mafft"),
        linsi_max=_EVO_LINSI_MAX,
        trim=cfg("evolution/trim_tool", "clipkit"),
        clipkit_mode=cfg("evolution/clipkit_mode", "kpic-smart-gap"),
        clean=f"{RESULTS}/evolution/align/family.input.faa",
    log:
        f"{LOGS}/evolution/family_align.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=720,
    shell:
        r"""
        mkdir -p $(dirname {output.aln}) $(dirname {log})
        # MAFFT needs an ABSOLUTE TMPDIR (it pushd's into its mktemp dir): keep an absolute
        # TMPDIR (e.g. the scheduler's), replace a relative or unset one by /tmp.
        case "${{TMPDIR:-}}" in /*) ;; *) export TMPDIR=/tmp ;; esac
        # Prodigal proteins end in '*': MAFFT silently deletes '*', so strip the terminal
        # stop explicitly (internal '*' -> X) and align exactly what codon_align expects.
        seqkit replace -s -p '\*$' -r '' {input.faa} 2> {log} | seqkit replace -s -p '\*' -r 'X' > {params.clean} 2>> {log}
        n=$(awk '/^>/{{n++}} END{{print n+0}}' {params.clean})
        if [ "$n" -eq 0 ]; then
            : > {output.aln}; : > {output.trim}  # nonfatal: no validated Prometheus proteins
            echo "none: no validated Prometheus proteins" > {output.method}
        elif [ "$n" -eq 1 ]; then
            cp {params.clean} {output.aln}; cp {params.clean} {output.trim}
            echo "none: a single sequence is its own alignment" > {output.method}
        else
            if [ "{params.tool}" = "famsa" ] || {{ [ "{params.tool}" = "mafft" ] && [ "$n" -gt {params.linsi_max} ]; }}; then
                famsa -t {threads} {params.clean} {output.aln} >> {log} 2>&1
                echo -e "FAMSA (famsa -t {threads})\tn=$n\tevolution.align_tool={params.tool}\tlinsi_max_seqs={params.linsi_max}" > {output.method}
            elif [ "{params.tool}" = "mafft" ]; then
                # --threadit 0: MAFFT's multithreaded iterative refinement differs run to run (MAFFT multithreading docs); this keeps it reproducible
                mafft --amino --localpair --maxiterate 1000 --thread {threads} --threadit 0 {params.clean} > {output.aln} 2>> {log}
                echo -e "MAFFT $(mafft --version 2>&1 | head -n1) L-INS-i (--localpair --maxiterate 1000 --threadit 0)\tn=$n\tlinsi_max_seqs={params.linsi_max}" > {output.method}
            else
                echo "unknown evolution.align_tool '{params.tool}' (mafft|famsa)" >&2; exit 1
            fi
            if [ "{params.trim}" = "clipkit" ]; then
                clipkit {output.aln} -m {params.clipkit_mode} -o {output.trim} >> {log} 2>&1
                echo -e "trim\tClipKIT -m {params.clipkit_mode}" >> {output.method}
            elif [ "{params.trim}" = "none" ]; then
                cp {output.aln} {output.trim}
                echo -e "trim\tnone" >> {output.method}
            else
                echo "unknown evolution.trim_tool '{params.trim}' (clipkit|none)" >&2; exit 1
            fi
        fi
        """


rule family_tree:
    """ML Prometheus gene tree (IQ-TREE2, ModelFinder, UFBoot2 + SH-aLRT; UFBoot trees for ALE)."""
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
        mkdir -p {params.workdir} $(dirname {log})
        n=$(awk '/^>/{{n++}} END{{print n+0}}' {input.trim})
        if [ "$n" -ge 4 ]; then
            # IQ-TREE writes <prefix>.treefile/.ufboot itself; these ARE the rule outputs, so
            # they are copied only if the paths differ (cp onto itself exits 1). Identical
            # sequences are re-added to the ML tree and to every UFBoot tree (IQ-TREE 2.3.6).
            prefix={params.workdir}/family
            iqtree2 -s {input.trim} -m {params.model} -B {params.ufboot} -alrt {params.alrt} \
                    --wbtl -T {threads} --seed {params.seed} --prefix "$prefix" -redo > {log} 2>&1
            for ext in treefile ufboot; do
                if [ ! -s "$prefix.$ext" ]; then echo "IQ-TREE wrote no $prefix.$ext" >&2; exit 1; fi
            done
            if [ ! "$prefix.treefile" -ef {output.tree} ]; then cp "$prefix.treefile" {output.tree}; fi
            if [ ! "$prefix.ufboot" -ef {output.ufboot} ]; then cp "$prefix.ufboot" {output.ufboot}; fi
            ntip=$(awk '{{n += gsub(/,/, ",")}} END{{print n + 1}}' {output.tree})
            if [ "$ntip" -ne "$n" ]; then
                echo "family tree {output.tree} has $ntip tips but the alignment has $n sequences" >&2; exit 1
            fi
        else
            echo "only $n sequence(s): no family tree (IQ-TREE/UFBoot need >= 4)" > {log}
            echo "();" > {output.tree}; echo "();" > {output.ufboot}  # nonfatal: < 4 sequences, no tree is defined
        fi
        """


rule domain_trees:
    """Per-domain trees vs the family tree: RF, supported conflicts, AU test (domain shuffling)."""
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
        model=cfg("evolution/domain_tree_model", cfg("evolution/gene_tree_model", "MFP")),
        ufboot=cfg("tree/ufboot", 1000),
        min_taxa=cfg("evolution/domain_min_taxa", 6),
        tool=cfg("evolution/align_tool", "mafft"),
        linsi_max=_EVO_LINSI_MAX,
        trim=cfg("evolution/trim_tool", "clipkit"),
        clipkit_mode=cfg("evolution/clipkit_mode", "kpic-smart-gap"),
    log:
        f"{LOGS}/evolution/domain_trees.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=720,
    shell:
        r"""
        mkdir -p {params.workdir} $(dirname {log})
        python workflow/scripts/py/domain_trees.py \
            --proteins {input.proteins} --domtbl {input.domtbl} --family-tree {input.family_tree} \
            --domains '{params.domains_json}' --threads {threads} --seed {params.seed} \
            --model {params.model} --ufboot {params.ufboot} --min-taxa {params.min_taxa} \
            --align-tool {params.tool} --linsi-max {params.linsi_max} --trim {params.trim} \
            --clipkit-mode {params.clipkit_mode} --workdir {params.workdir} --out {output.rf} > {log} 2>&1
        """


rule combined_cds:
    """Concatenate CDS (fna) with genome-prefixed ids identical to the protein ids."""
    input:
        fna=expand(f"{RESULTS}/genes/{{genome}}/{{genome}}.fna", genome=GENOMES),
        samples=config["samples"],
    output:
        ffn=f"{RESULTS}/evolution/all_cds.ffn",
        map=f"{RESULTS}/evolution/cds_map.tsv",
    params:
        genes_dir=f"{RESULTS}/genes",
    log:
        f"{LOGS}/evolution/combined_cds.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=60,
    shell:
        r"""
        mkdir -p $(dirname {output.ffn}) $(dirname {log})
        python workflow/scripts/py/combine_proteomes.py \
            --samples {input.samples} --genes-dir {params.genes_dir} --suffix fna \
            --out-faa {output.ffn} --out-map {output.map} > {log} 2>&1
        """


rule codon_align:
    """Codon (in-frame) alignment: pre-checked pep/CDS pairs -> PAL2NAL -codontable 11 -> verified."""
    input:
        prot_aln=rules.family_align.output.aln,
        cds=rules.combined_cds.output.ffn,
        proteins=rules.validate_prometheus.output.proteins,
    output:
        codon=f"{RESULTS}/evolution/align/family.codon.fasta",
        report=f"{RESULTS}/evolution/align/family.codon_report.tsv",
    params:
        table=cfg("evolution/codon_table", 11),
        max_drop=cfg("evolution/max_cds_drop_frac", 0.05),
    log:
        f"{LOGS}/evolution/codon_align.log",
    conda:
        "../envs/selection.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        mkdir -p $(dirname {output.codon}) $(dirname {log})
        python workflow/scripts/py/codon_align.py \
            --prot-aln {input.prot_aln} --cds {input.cds} --proteins {input.proteins} \
            --codon-table {params.table} --max-drop-frac {params.max_drop} \
            --report {output.report} --out {output.codon} > {log} 2>&1
        """


rule reconciliation:
    """Species-tree-aware DTL reconciliation (GeneRax UndatedDTL; optional ALE on UFBoot trees)."""
    input:
        aln=rules.family_align.output.trim,
        gene_tree=rules.family_tree.output.tree,
        ufboot=rules.family_tree.output.ufboot,
        species_tree=rules.species_tree.output.tree,
        master=rules.build_master_table.output.master,
    output:
        events=f"{RESULTS}/evolution/reconciliation/reconciliation_events.tsv",
        transfers=f"{RESULTS}/evolution/reconciliation/reconciliation_transfers.tsv",
    params:
        workdir=f"{RESULTS}/evolution/reconciliation",
        model=cfg("evolution/reconciliation/generax_model", "UndatedDTL"),
        strategy=cfg("evolution/reconciliation/generax_strategy", "SPR"),
        subst_model=cfg("evolution/reconciliation/subst_model", "LG+G4"),
        ale_samples=cfg("evolution/reconciliation/ale_samples", 100),
        seed=cfg("seed", 12345),
        run_generax=str(RUN_GENERAX),
        run_ale=str(RUN_ALE),
    log:
        f"{LOGS}/evolution/reconciliation.log",
    conda:
        "../envs/reconciliation.yaml"
    # GeneRax parallelises over families via MPI; one family = one process.
    threads: T_SMALL
    resources:
        mem_mb=16000, runtime=1440,
    shell:
        r"""
        mkdir -p {params.workdir} $(dirname {log})
        python workflow/scripts/py/run_reconciliation.py \
            --aln {input.aln} --gene-tree {input.gene_tree} --ufboot {input.ufboot} \
            --species-tree {input.species_tree} --master {input.master} --model {params.model} \
            --strategy {params.strategy} --subst-model '{params.subst_model}' --seed {params.seed} \
            --run-generax {params.run_generax} --run-ale {params.run_ale} --ale-samples {params.ale_samples} \
            --threads {threads} --workdir {params.workdir} \
            --out {output.events} --out-transfers {output.transfers} > {log} 2>&1
        """


rule ancestral_states:
    """Gain/loss by stochastic mapping (phytools) + structured methanogen-gain enrichment null."""
    input:
        tree=rules.species_tree.output.tree,
        master=rules.build_master_table.output.master,
    output:
        gl=f"{RESULTS}/evolution/ancestral/gains_losses.tsv",
        gl_null=f"{RESULTS}/evolution/ancestral/gains_losses_null.tsv",
        plots=expand(f"{RESULTS}/evolution/plots/ancestral_gain_loss.{{ext}}", ext=PLOT_FORMATS),
    params:
        nsim=cfg("evolution/ancestral/nsim", 1000),
        nnull=cfg("evolution/ancestral/nnull", 1000),
        q_method=cfg("evolution/ancestral/q_method", "empirical"),
        root_prior=cfg("evolution/ancestral/root_prior", "fitzjohn"),
        prevalence_tol=cfg("evolution/ancestral/prevalence_tol", 0.1),
        tiers=_EVO_TIERS,
        tiers_sens=_EVO_TIERS_SENS,
        seed=cfg("seed", 12345),
        outdir=f"{RESULTS}/evolution/ancestral",
        plotdir=f"{RESULTS}/evolution/plots",
    log:
        f"{LOGS}/evolution/ancestral_states.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_SMALL
    resources:
        mem_mb=16000, runtime=720,
    shell:
        r"""
        mkdir -p {params.outdir} {params.plotdir} $(dirname {log})
        Rscript workflow/scripts/R/ancestral_states.R \
            --tree {input.tree} --master {input.master} --nsim {params.nsim} --nnull {params.nnull} \
            --q-method {params.q_method} --root-prior {params.root_prior} \
            --prevalence-tol {params.prevalence_tol} --tiers '{params.tiers}' \
            --tiers-sensitivity '{params.tiers_sens}' --seed {params.seed} \
            --out {output.gl} --out-null {output.gl_null} \
            --plot-prefix {params.plotdir}/ancestral_gain_loss > {log} 2>&1
        """


rule selection:
    """dN/dS on catalytically-intact loci (HyPhy FEL/MEME/BUSTED/aBSREL/RELAX; optional GARD)."""
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
        relax_labelling=cfg("evolution/selection/relax_labelling", "clade"),
        intact_only=str(cfg("evolution/selection/intact_only", True)),
        site_pvalue=cfg("evolution/selection/site_pvalue", 0.1),
        gard=str(cfg("evolution/selection/gard", False)),
        gard_mode=cfg("evolution/selection/gard_mode", "Faster"),
        seed=cfg("seed", 12345),
        domains_json=lambda wc: json_domains(),
        workdir=f"{RESULTS}/evolution/selection",
    log:
        f"{LOGS}/evolution/selection.log",
    conda:
        "../envs/selection.yaml"
    threads: T_MED
    resources:
        mem_mb=16000, runtime=2880,
    shell:
        r"""
        mkdir -p {params.workdir} $(dirname {log})
        python workflow/scripts/py/run_selection.py \
            --codon {input.codon} --gene-tree {input.gene_tree} --catalog {input.catalog} \
            --master {input.master} --domtbl {input.domtbl} --domains '{params.domains_json}' \
            --methods {params.methods} --relax-partition {params.relax_partition} \
            --relax-labelling {params.relax_labelling} --gard {params.gard} --gard-mode {params.gard_mode} \
            --intact-only {params.intact_only} --site-pvalue {params.site_pvalue} --run {params.run} \
            --threads {threads} --seed {params.seed} --workdir {params.workdir} \
            --out-summary {output.summary} --out-sites {output.sites} > {log} 2>&1
        """


rule plot_tanglegram:
    """Gene tree vs species tree (carrier genomes), links coloured by host methanogen/ANME status."""
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
        mkdir -p {params.plotdir} $(dirname {log})
        Rscript workflow/scripts/R/plot_tanglegram.R \
            --gene-tree {input.gene_tree} --species-tree {input.species_tree} --master {input.master} \
            --plot-prefix {params.plotdir}/genetree_speciestree_tanglegram > {log} 2>&1
        """


rule plot_selection:
    """Fraction of selected sites per domain (FEL/MEME, Wilson CI) + per-site beta/alpha."""
    input:
        summary=rules.selection.output.summary,
        sites=rules.selection.output.sites,
    output:
        expand(f"{RESULTS}/evolution/plots/selection_by_domain.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{RESULTS}/evolution/plots",
    log:
        f"{LOGS}/evolution/plot_selection.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        mkdir -p {params.outdir} $(dirname {log})
        python workflow/scripts/py/plot_selection.py \
            --summary {input.summary} --sites {input.sites} --outdir {params.outdir} > {log} 2>&1
        """
