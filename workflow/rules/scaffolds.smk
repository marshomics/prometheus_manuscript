# =====================================================================
# Stage 3 — reference scaffolds: GTDB-Tk marker MSA, rooted species tree, gene-based
# methanogen classification, master table. GTDB taxonomy + genome QC come from
# Phase-1 preprocessing (rules/preprocess.smk); only the species-tree marker MSA is
# computed here, on the species representatives.
#
# No silent failure (Contract v2 item 4): a failing tool fails its rule. The only
# placeholder (a '();' tree when the MSA holds < 4 genomes) is marked '# nonfatal:'.
# Defaults of config keys that config.yaml does not (yet) define are given in cfg().
# =====================================================================

def _scaf_bool(x):
    return str(x).strip().lower() in ("true", "1", "yes")


_SCAF_SCOPE = str(cfg("tree/scope", "focal")).strip().lower()
_SCAF_MARKERS = str(cfg("tree/marker_set", "ar53")).strip()
if _SCAF_SCOPE not in ("focal", "domain"):
    raise ValueError(f"[scaffolds] tree.scope must be 'focal' or 'domain' (got {_SCAF_SCOPE!r})")
if _SCAF_MARKERS not in ("ar53", "bac120"):
    raise ValueError(f"[scaffolds] tree.marker_set must be 'ar53' or 'bac120' (got {_SCAF_MARKERS!r})")
_SCAF_DOMAIN = _SCAF_SCOPE == "domain"
_SCAF_ANME = _scaf_bool(cfg("methanogen/anme_via_mcrA_tree", True))
_SCAF_PARTIAL = cfg("methanogen/partial_score_fraction", 0.5)
# Methanogenesis pathway capacity (design §3.3): opt-in, because it needs the KOfam profile
# subset staged by workflow/scripts/py/fetch_kofam_subset.py (resources/README.md).
_MG_ENABLED = MGEN_ENABLED          # common.smk (one definition for every consumer)
_MG_EXPECT = str(cfg("methanogenesis/taxonomy_expectation", "") or "")


rule gtdbtk_align:
    """GTDB-Tk 2.4 identify + align on exactly the analysed genomes (batchfile).
    tree.scope=focal: the user MSA (--skip_gtdb_refs). tree.scope=domain: the MSA that
    also holds the GTDB reference genomes (references are pruned again after rooting).
    Genomes that GTDB-Tk leaves out of the MSA are listed, with GTDB-Tk's own reason,
    in gtdbtk_dropped.tsv."""
    input:
        genomes=expand(os.path.join(GENOME_DIR, "{genome}." + EXT), genome=GENOMES),
        reps=rules.filter_collapse.output.reps_list,
    output:
        msa=f"{RESULTS}/scaffolds/gtdbtk_align/user_msa.fasta.gz",
        dropped=f"{RESULTS}/scaffolds/gtdbtk_align/gtdbtk_dropped.tsv",
    params:
        gdir=GENOME_DIR,
        ext=EXT,
        data=cfg("db/gtdbtk_data", ""),
        outdir=f"{RESULTS}/scaffolds/gtdbtk_align",
        refs="" if _SCAF_DOMAIN else "--skip_gtdb_refs",
        msa_name=(f"gtdbtk.{_SCAF_MARKERS}.msa.fasta.gz" if _SCAF_DOMAIN
                  else f"gtdbtk.{_SCAF_MARKERS}.user_msa.fasta.gz"),
        markers=_SCAF_MARKERS,
        other="bac120" if _SCAF_MARKERS == "ar53" else "ar53",
        min_perc_aa=cfg("tree/min_perc_aa", 10),
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
        out={params.outdir}
        mkdir -p "$out"
        if [ -z "{params.data}" ] || [ ! -d "{params.data}" ]; then
            echo "ERROR: db.gtdbtk_data='{params.data}' is not a GTDB-Tk reference data directory" > {log}
            exit 1
        fi
        export GTDBTK_DATA_PATH="$(cd "{params.data}" && pwd)"
        gabs="$(cd "{params.gdir}" && pwd)"
        # Exactly the analysed genomes (absolute path <TAB> GTDB-Tk id), not the directory.
        # GTDB-Tk 2.4.0 aborts the whole run when a user id equals a GTDB reference id, and
        # those include the bare NCBI accession of every GTDB genome (gtdbtk/tools.py
        # get_reference_ids; main.py _genomes_to_process), so a genome named by its
        # accession (GCA_...) would stop the species tree. GTDB-Tk therefore sees
        # 'user_<genome_id>'; the prefix is stripped from every output below.
        awk -v d="$gabs" -v x="{params.ext}" 'NF && $1 != "genome_id" {{print d "/" $1 "." x "\tuser_" $1}}' \
            {input.reps} > "$out/batchfile.tsv"
        # Fresh GTDB-Tk directories: a stale identify dir would leak earlier genomes into align.
        rm -rf "$out/identify" "$out/align"
        gtdbtk identify --batchfile "$out/batchfile.tsv" --out_dir "$out" \
            --cpus {threads} > {log} 2>&1
        gtdbtk align --identify_dir "$out" --out_dir "$out" --cpus {threads} \
            --min_perc_aa {params.min_perc_aa} {params.refs} >> {log} 2>&1
        # Genomes GTDB-Tk left out of the MSA, each with GTDB-Tk's own reason.
        {{
          f="$out/align/gtdbtk.{params.markers}.filtered.tsv"
          if [ -s "$f" ]; then awk -F'\t' -v m={params.markers} '{{print $1 "\t" m " --min_perc_aa filter: " $2}}' "$f"; fi
          f="$out/align/gtdbtk.{params.other}.user_msa.fasta.gz"
          if [ -s "$f" ]; then gunzip -c "$f" | awk -v o={params.other} '/^>/ {{sub(/^>/, ""); print $1 "\tassigned to the " o " marker set by GTDB-Tk"}}'; fi
          f="$out/align/gtdbtk.{params.other}.filtered.tsv"
          if [ -s "$f" ]; then awk -F'\t' -v o={params.other} '{{print $1 "\tassigned to the " o " marker set by GTDB-Tk; " $2}}' "$f"; fi
          f="$out/align/intermediate_results/gtdbtk.align.failed.tsv"
          if [ -s "$f" ]; then awk -F'\t' '{{print $1 "\t" $2}}' "$f"; fi
          f="$out/identify/gtdbtk.failed_genomes.tsv"
          if [ -s "$f" ]; then awk -F'\t' '{{print $1 "\tGTDB-Tk identify failed: " $2}}' "$f"; fi
        }} > "$out/dropped.raw.tsv"
        awk -F'\t' 'BEGIN {{OFS = "\t"}} {{sub(/^user_/, "", $1); print}}' "$out/dropped.raw.tsv" > {output.dropped}
        msa="$out/align/{params.msa_name}"
        if [ ! -s "$msa" ]; then
            echo "ERROR: GTDB-Tk wrote no $msa (no analysed genome reached the {params.markers} MSA; see {output.dropped})" >> {log}
            exit 1
        fi
        # Back to the pipeline's genome ids (GTDB reference headers, tree.scope=domain, start
        # with RS_/GB_ and are untouched).
        gunzip -c "$msa" | awk '/^>user_/ {{sub(/^>user_/, ">")}} {{print}}' | gzip -c > {output.msa}
        echo "MSA sequences: $(gunzip -c {output.msa} | awk '/^>/{{c++}} END{{print c+0}}'); genomes dropped by GTDB-Tk: $(awk 'END{{print NR}}' {output.dropped})" >> {log}
        """


rule species_tree:
    """ML species tree (IQ-TREE 2: two-step PMSF, or ModelFinder), rooted by root_tree.py
    (outgroup preferred; labelled midpoint fallback) and pruned to the analysed genomes.
    rooting_method.txt records the root actually applied; species_tree_membership.tsv
    says why any analysed genome is not a tip."""
    input:
        msa=rules.gtdbtk_align.output.msa,
        dropped=rules.gtdbtk_align.output.dropped,
        taxonomy=rules.filter_collapse.output.taxonomy,
        reps=rules.filter_collapse.output.reps_list,
    output:
        tree=f"{RESULTS}/scaffolds/species_tree.treefile",
        rooting=f"{RESULTS}/scaffolds/species_tree/rooting_method.txt",
        membership=f"{RESULTS}/scaffolds/species_tree/species_tree_membership.tsv",
        full_tree=f"{RESULTS}/scaffolds/species_tree/species_tree.full.treefile",
    params:
        workdir=f"{RESULTS}/scaffolds/species_tree",
        pmsf=str(_scaf_bool(cfg("tree/pmsf", True))),
        guide_model=cfg("tree/pmsf_guide_model", "LG+F+G"),
        pmsf_model=cfg("tree/pmsf_model", "LG+C20+F+G"),
        model=cfg("tree/iqtree_model", "MFP"),
        ufboot=cfg("tree/ufboot", 1000),
        alrt=cfg("tree/alrt", 1000),
        seed=cfg("seed", 12345),
        rooting_method=cfg("tree/rooting_method", "outgroup"),
        outgroup=",".join(str(g) for g in (cfg("tree/outgroup", []) or [])),
        outgroup_taxon=cfg("tree/outgroup_taxon", "") or "",
        on_failure=cfg("tree/outgroup_on_failure", "fail"),
    log:
        f"{LOGS}/scaffolds/species_tree.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_LARGE
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt,
                                          base=128000 if _SCAF_DOMAIN else 32000,
                                          per_attempt=64000 if _SCAF_DOMAIN else 16000),
        runtime=10080 if _SCAF_DOMAIN else 2880,
    shell:
        r"""
        wd={params.workdir}
        mkdir -p "$wd"
        msa="$wd/msa.faa"
        # Sequence id only on header lines (GTDB reference headers carry the taxonomy).
        gunzip -c {input.msa} | awk '/^>/ {{print $1; next}} {{print}}' > "$msa"
        # GTDB reference taxonomy from those headers (tree.scope=domain; empty for focal).
        gunzip -c {input.msa} | awk '/^>/ && index($0, " ") {{id = substr($1, 2); tx = $0; sub(/^>[^ ]+ /, "", tx); print id "\t" tx}}' > "$wd/ref_taxonomy.tsv"
        n=$(awk '/^>/ {{c++}} END {{print c+0}}' "$msa")
        unrooted="$wd/species.unrooted.treefile"
        rm -f "$wd"/guide.* "$wd"/species.*
        echo "species_tree: $n sequences in the marker MSA" > {log}
        if [ "$n" -ge 4 ]; then
            if [ "{params.pmsf}" = "True" ]; then
                # PMSF: guide tree, then posterior mean site frequencies under the mixture model.
                iqtree2 -s "$msa" -m "{params.guide_model}" -T {threads} --seed {params.seed} \
                    --prefix "$wd/guide" -redo >> {log} 2>&1
                iqtree2 -s "$msa" -m "{params.pmsf_model}" -ft "$wd/guide.treefile" \
                    -B {params.ufboot} -alrt {params.alrt} -T {threads} --seed {params.seed} \
                    --prefix "$wd/species" -redo >> {log} 2>&1
            else
                iqtree2 -s "$msa" -m "{params.model}" -B {params.ufboot} -alrt {params.alrt} \
                    -T {threads} --seed {params.seed} --prefix "$wd/species" -redo >> {log} 2>&1
            fi
            cp "$wd/species.treefile" "$unrooted"
        else
            echo "();" > "$unrooted"  # nonfatal: fewer than 4 genomes in the MSA, no resolvable topology; recorded in membership + rooting files
        fi
        python workflow/scripts/py/root_tree.py \
            --in-tree "$unrooted" --taxonomy {input.taxonomy} --ref-taxonomy "$wd/ref_taxonomy.tsv" \
            --method {params.rooting_method} --outgroup "{params.outgroup}" \
            --outgroup-taxon "{params.outgroup_taxon}" --on-outgroup-failure {params.on_failure} \
            --keep-tips {input.reps} --msa "$msa" --gtdbtk-dropped {input.dropped} \
            --out-tree {output.tree} --out-full-tree {output.full_tree} \
            --out-membership {output.membership} --log {output.rooting} >> {log} 2>&1
        """


rule mcr_search:
    """Detect mcrA/B/G (methanogens are defined from genes, not taxonomy): hmmsearch with
    the curated KOfam / GA thresholds applied exactly as KofamScan does, plus contig-edge
    fragments (Prodigal partial flags)."""
    input:
        proteome=rules.combined_proteome.output.faa,
        pmap=rules.combined_proteome.output.map,
        genes=expand(f"{RESULTS}/genes/{{genome}}/{{genome}}.genes.tsv", genome=GENOMES),
    output:
        hits=f"{RESULTS}/scaffolds/metabolism/mcr_hits.tsv",
        all_hits=f"{RESULTS}/scaffolds/metabolism/mcr_hits_all.tsv",
    params:
        hmmdir=cfg("refs/mcr_hmm_dir", "resources/refs/mcr_hmm"),
        score_fraction=cfg("methanogen/kofam_score_fraction", 1.0),
        partial=(f"--partial-score-fraction {_SCAF_PARTIAL}"
                 if _SCAF_PARTIAL not in (None, "", False) else ""),
        ko_list=cfg("db/kofam_ko_list", ""),
        kos=",".join(str(k) for k in (cfg("methanogen/mcr_kos", ["K00399", "K00401", "K00402"]) or [])),
        genes_dir=f"{RESULTS}/genes",
    log:
        f"{LOGS}/scaffolds/mcr_search.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=120,
    shell:
        r"""
        python workflow/scripts/py/mcr_search.py \
            --proteome {input.proteome} --pmap {input.pmap} --hmmdir {params.hmmdir} \
            --score-fraction {params.score_fraction} {params.partial} \
            --genes-dir {params.genes_dir} --ko-list "{params.ko_list}" --kos "{params.kos}" \
            --threads {threads} --out {output.hits} --out-all {output.all_hits} > {log} 2>&1
        """


rule methanogen_classify:
    """Methanogen vs ANME vs other by credible Mcr + McrA nearest-labelled-clade placement
    (MAFFT --add --keeplength, IQ-TREE). methanogenesis_subtype is never asserted from mcr."""
    input:
        hits=rules.mcr_search.output.hits,
        proteome=rules.combined_proteome.output.faa,
        pmap=rules.combined_proteome.output.map,
        ref_aln=[cfg("refs/mcrA_ref_aln", "resources/refs/mcrA_reference.afa")] if _SCAF_ANME else [],
        ref_labels=[cfg("refs/mcrA_ref_labels", "resources/refs/mcrA_reference_labels.tsv")] if _SCAF_ANME else [],
    output:
        calls=f"{RESULTS}/scaffolds/metabolism/methanogen_calls.tsv",
        placement=f"{RESULTS}/scaffolds/metabolism/mcrA_placement.tsv",
    params:
        required=",".join(str(k) for k in (cfg("methanogen/mcr_required", ["K00399"]) or [])),
        anme_tree=str(_SCAF_ANME),
        workdir=f"{RESULTS}/scaffolds/metabolism/mcrA_placement",
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/scaffolds/methanogen_classify.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=180,
    shell:
        r"""
        python workflow/scripts/py/classify_methanogen.py \
            --hits {input.hits} --proteome {input.proteome} --pmap {input.pmap} \
            --ref-aln "{input.ref_aln}" --ref-labels "{input.ref_labels}" \
            --required {params.required} --anme-tree {params.anme_tree} \
            --threads {threads} --seed {params.seed} --workdir {params.workdir} \
            --out {output.calls} --out-placement {output.placement} > {log} 2>&1
        """


rule methanogenesis_search:
    """Detect the methanogenesis pathway KOs in the methanogens with the curated KOfam
    thresholds (identical acceptance to mcr_search; kofam.py), rescuing contig-edge fragments
    and genes Prodigal truncated at an in-frame UAG (pyrrolysine) or UGA (selenocysteine)."""
    input:
        proteome=rules.combined_proteome.output.faa,
        pmap=rules.combined_proteome.output.map,
        calls=rules.methanogen_classify.output.calls,
        rules=cfg("methanogenesis/rules",
                  "workflow/resources/methanogenesis/rules/methanogenesis_rules.yaml"),
        genes=expand(f"{RESULTS}/genes/{{genome}}/{{genome}}.genes.tsv", genome=GENOMES),
    output:
        hits=f"{RESULTS}/scaffolds/metabolism/methanogenesis_ko_hits.tsv",
        all_hits=f"{RESULTS}/scaffolds/metabolism/methanogenesis_ko_hits_all.tsv",
        readthrough=f"{RESULTS}/scaffolds/metabolism/methanogenesis_readthrough.faa",
    params:
        hmmdir=cfg("db/kofam_profiles", "resources/db/kofam/methanogenesis_profiles"),
        ko_list=cfg("db/kofam_ko_list", ""),
        genes_dir=f"{RESULTS}/genes",
        genome_dir=GENOME_DIR,
        genome_ext=EXT,
        score_fraction=cfg("methanogen/kofam_score_fraction", 1.0),
        partial=(f"--partial-score-fraction {_SCAF_PARTIAL}"
                 if _SCAF_PARTIAL not in (None, "", False) else ""),
        readthrough=str(_scaf_bool(cfg("methanogenesis/readthrough", True))),
    log:
        f"{LOGS}/scaffolds/methanogenesis_search.log",
    conda:
        "../envs/phylo.yaml"
    threads: T_MED
    resources:
        mem_mb=8000, runtime=240,
    shell:
        r"""
        python workflow/scripts/py/methanogenesis_search.py \
            --proteome {input.proteome} --pmap {input.pmap} --genes-dir {params.genes_dir} \
            --hmmdir {params.hmmdir} --rules {input.rules} --ko-list "{params.ko_list}" \
            --methanogen-calls {input.calls} --genome-dir {params.genome_dir} \
            --genome-ext {params.genome_ext} --score-fraction {params.score_fraction} \
            {params.partial} --readthrough {params.readthrough} --threads {threads} \
            --out {output.hits} --out-all {output.all_hits} \
            --out-readthrough {output.readthrough} > {log} 2>&1
        """


rule methanogenesis_pathways:
    """Pathway capacity per methanogen: CO2-reducing, acetoclastic, methylotrophic,
    methyl-dismutating, methyl-reducing, each present / absent / unresolved with a
    completeness-aware bound on a false absence, plus substrate flags, a composite subtype
    label and the GTDB-taxon reconciliation (design §3.3)."""
    input:
        hits=rules.methanogenesis_search.output.hits,
        calls=rules.methanogen_classify.output.calls,
        meta=rules.genome_metadata.output.meta,
        taxonomy=rules.filter_collapse.output.taxonomy,
        rules=rules.methanogenesis_search.input.rules,
        expectation=[_MG_EXPECT] if _MG_EXPECT else [],
    output:
        pathways=f"{RESULTS}/scaffolds/metabolism/methanogenesis_pathways.tsv",
        reconciliation=f"{RESULTS}/scaffolds/metabolism/methanogenesis_reconciliation.tsv",
    params:
        alpha=cfg("methanogenesis/absence_alpha", 0.05),
        expectation=(f"--expectation {_MG_EXPECT}" if _MG_EXPECT else ""),
        taxonomy=lambda wc, input: f"--taxonomy {input.taxonomy}",
    log:
        f"{LOGS}/scaffolds/methanogenesis_pathways.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/methanogenesis_pathways.py \
            --hits {input.hits} --methanogen {input.calls} --meta {input.meta} \
            --rules {input.rules} {params.taxonomy} {params.expectation} \
            --absence-alpha {params.alpha} --out {output.pathways} \
            --out-reconciliation {output.reconciliation} > {log} 2>&1
        """


rule build_master_table:
    """Join QC, taxonomy, tree tip, methanogen status, ecology, Prometheus, defense load."""
    input:
        meta=rules.genome_metadata.output.meta,
        taxonomy=rules.filter_collapse.output.taxonomy,
        tree=rules.species_tree.output.tree,
        membership=rules.species_tree.output.membership,
        methanogen=rules.methanogen_classify.output.calls,
        pathways=([rules.methanogenesis_pathways.output.pathways] if _MG_ENABLED else []),
        ecology=config["metadata"],
        catalog=rules.validate_prometheus.output.catalog,
        sensitivity=rules.validate_prometheus.output.sensitivity,
        reps=rules.filter_collapse.output.reps_list,
        df_systems=expand(f"{RESULTS}/defense/defensefinder/{{genome}}/{{genome}}_defense_finder_systems.tsv", genome=GENOMES),
    output:
        master=f"{RESULTS}/scaffolds/master_table.tsv",
        ecology=f"{RESULTS}/scaffolds/ecology_harmonized.tsv",
        missingness=f"{RESULTS}/scaffolds/ecology_missingness.tsv",
    params:
        df_dir=f"{RESULTS}/defense/defensefinder",
        envo_map=cfg("ecology/envo_map", "") or "",
        system=cfg("prometheus/defensefinder_system", "Prometheus"),
        primary=",".join(str(t) for t in (cfg("stats/confidence_tiers_primary", ["A"]) or ["A"])),
        pathways=lambda wc, input: (f"--pathways {input.pathways}" if _MG_ENABLED else ""),
    log:
        f"{LOGS}/scaffolds/build_master_table.log",
    conda:
        "../envs/annotation.yaml"
    threads: T_SMALL
    resources:
        mem_mb=8000, runtime=30,
    shell:
        r"""
        python workflow/scripts/py/build_master_table.py \
            --meta {input.meta} --taxonomy {input.taxonomy} --tree {input.tree} \
            --tree-membership {input.membership} \
            --methanogen {input.methanogen} --ecology {input.ecology} --envo-map "{params.envo_map}" \
            --catalog {input.catalog} --sensitivity {input.sensitivity} --reps {input.reps} \
            --df-dir {params.df_dir} --system "{params.system}" --primary-tiers "{params.primary}" \
            {params.pathways} \
            --out-master {output.master} --out-ecology {output.ecology} \
            --out-missingness {output.missingness} > {log} 2>&1
        """
