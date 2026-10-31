# =====================================================================
# Bacterial / cross-domain branch (OPT-IN: bacteria.enabled: true).
#
# Ingests a PRE-COMPUTED bacterial DefenseFinder set (combined systems/genes
# tables + per-genome faa/cds), validates bacterial Prometheus with the SAME
# scripts as archaea, then adds the cross-domain analyses: the combined
# phylogeny (origin + bacteria<->archaea transfer direction), matched prevalence
# and architecture, cross-domain dN/dS (+ a synonymous-saturation check), and
# co-occurrence. Built with the `crossdomain` target; the archaeal `all`/`preprocess`
# runs are untouched and are consumed here only as read-only inputs (archaeal
# validated proteins + CDS, master table, typology, domtbl, viral summary).
#
# Contract v2: every tool failure fails its rule; placeholders ('();' trees,
# header-only tables) are written only for a genuinely empty input and are marked
# `# nonfatal:`. Every rule logs. Sequence ids entering MAFFT/IQ-TREE are sanitised to
# [A-Za-z0-9_.-] (IQ-TREE rewrites ':' ',' etc., which would silently break the
# reference-label lookups), and duplicated ids are an error.
# =====================================================================

import re

BD = f"{RESULTS}/bacteria"
CD = f"{RESULTS}/crossdomain"


def _bact(key, default=""):
    return cfg(f"bacteria/{key}", default)


def _opt_input(key):
    """Optional input file: [] (no dependency) when the config key is unset."""
    p = _bact(key, "")
    return [p] if p else []


def _cd(key, default):
    return cfg(f"bacteria/crossdomain/{key}", default)


def _csv(value, default):
    """'a,b' from a YAML list or a comma/space separated string (a bare ",".join() would
    turn the string "A,B" into "A,,,B")."""
    items = value if isinstance(value, (list, tuple)) else re.split(r"[,\s]+", str(value))
    return ",".join(str(x).strip() for x in items if str(x).strip()) or default


def _tiers(key, default):
    return _csv(cfg(key, default), default)


# Cross-domain dN/dS needs bacterial CDS (bacteria.cds_dir). Without them the codon
# alignment would hold archaea only -- not a cross-domain test -- and codon_align.py
# (fails when > --max-drop-frac of the expected loci lack a CDS) would abort the run.
BACT_CDS = bool(_bact("cds_dir", ""))
# validate_prometheus.py --system takes ONE focal type (excluded from the defense
# context); bacteria.system may list several exact names -- the first is the type.
BACT_FOCAL = str(_bact("system", "Prometheus")).split(",")[0].strip() or "Prometheus"


# Sanitise FASTA ids / label files to the newick- and IQ-TREE-safe alphabet.
_AWK_FASTA = r"""awk '/^>/{{id=substr($1,2); gsub(/[^A-Za-z0-9_.-]/,"_",id); print ">"id; next}} {{print}}'"""
_AWK_DUPS = r"""awk '/^>/{{c[$1]++}} END{{for (k in c) if (c[k] > 1) printf "%s ", k}}'"""
_AWK_COUNT = r"""awk '/^>/{{c++}} END{{print c+0}}'"""


def _awk_labels(group):
    """awk printing '<sanitised id>\t<group>' for every FASTA header."""
    return ("awk -v grp=" + group + r""" '/^>/{{id=substr($1,2); gsub(/[^A-Za-z0-9_.-]/,"_",id); print id "\t" grp}}'""")


CHECKV_HEADER = ("contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\t"
                 "viral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\t"
                 "completeness_method\\tcontamination\\tkmer_freq\\twarnings")


def _write_bact_samples(src, table, reps, logf):
    """Write the bacterial analysis universe (= BACT_GENOMES, the ids the per-genome
    viral rules are expanded over) as a genome_id table (+ completeness/contamination/
    genome_size copied from the QC table when present) and a one-id-per-line list in
    the format of the archaeal representatives.txt."""
    os.makedirs(os.path.dirname(os.path.abspath(table)), exist_ok=True)
    with open(logf, "w") as lg:
        if not BACT_GENOMES:
            raise ValueError(f"[crossdomain] no bacterial genome ids could be read from {src}")
        df = pd.read_csv(src, sep="\t", dtype=str, comment="#")
        idc = "genome_id" if "genome_id" in df.columns else (
            "user_genome" if "user_genome" in df.columns else df.columns[0])
        df[idc] = df[idc].astype(str).str.strip()
        low = {c.lower(): c for c in df.columns}
        keep = [low[c] for c in ("completeness", "contamination", "genome_size") if c in low]
        sub = df.drop_duplicates(subset=[idc]).set_index(idc)
        out = pd.DataFrame({"genome_id": BACT_GENOMES})
        for c in keep:
            out[c.lower()] = [sub.at[g, c] if g in sub.index else "" for g in BACT_GENOMES]
        out.to_csv(table, sep="\t", index=False)
        with open(reps, "w") as fh:
            fh.write("\n".join(BACT_GENOMES) + "\n")
        lg.write(f"{len(BACT_GENOMES)} bacterial genomes from {src} (columns: {keep})\n")


# --------------------------------------------------------------------- #
# 0. Bacterial analysis universe
# --------------------------------------------------------------------- #
rule bact_samples:
    """genome_id table + representatives list of the bacterial genome universe."""
    input:
        src=_bact("qc") or _bact("taxonomy"),
    output:
        table=f"{BD}/bacterial_samples.tsv",
        reps=f"{BD}/bacterial_representatives.txt",
    log: f"{LOGS}/bacteria/samples.log"
    threads: 1
    resources: mem_mb=2000, runtime=10
    run:
        _write_bact_samples(input.src, output.table, output.reps, log[0])


# --------------------------------------------------------------------- #
# 1. Ingest bacterial DefenseFinder -> Prometheus candidates (+ split tables)
# --------------------------------------------------------------------- #
rule bact_extract:
    input:
        systems=_bact("df_systems"),
        genes=_bact("df_genes"),
        samples=rules.bact_samples.output.table,
    output:
        candidates=f"{BD}/prometheus/candidates.faa",
        sources=f"{BD}/prometheus/candidate_sources.tsv",
        raw=f"{BD}/prometheus/raw_hits.tsv",
        cds=f"{BD}/prometheus/bacterial_prometheus_cds.ffn",
        report=f"{BD}/prometheus/extraction_report.tsv",
        split_done=f"{BD}/defense/.split.done",
    params:
        faa_dir=_bact("faa_dir"),
        cds_dir=_bact("cds_dir"),
        sys_gcol=_bact("systems_genome_col", "genome_id"),
        gene_gcol=_bact("genes_genome_col", "genome_id"),
        system=_bact("system", "Prometheus"),
        allow_missing=str(_bact("allow_missing_proteins", False)),
        defense_dir=f"{BD}/defense",
        genes_dir=f"{BD}/genes",
    log: f"{LOGS}/bacteria/extract.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        python workflow/scripts/py/extract_bacterial_prometheus.py \
            --df-systems {input.systems} --df-genes {input.genes} --genomes {input.samples} \
            --faa-dir "{params.faa_dir}" --cds-dir "{params.cds_dir}" \
            --systems-genome-col {params.sys_gcol} --genes-genome-col {params.gene_gcol} \
            --system "{params.system}" --allow-missing {params.allow_missing} \
            --out-defense-dir {params.defense_dir} --out-genes-dir {params.genes_dir} \
            --out-candidates {output.candidates} --out-sources {output.sources} \
            --out-raw {output.raw} --out-cds {output.cds} --out-report {output.report} > {log} 2>&1
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
        hmmout=f"{BD}/prometheus/domains.hmmscan.txt",
    log: f"{LOGS}/bacteria/domains.log"
    conda: "../envs/annotation.yaml"
    threads: T_MED
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.domtbl})
        if [ -s {input.candidates} ]; then
            for ext in h3m h3i h3f h3p; do
                [ -s {input.pfam}.$ext ] || {{ echo "ERROR: {input.pfam}.$ext missing: run 'hmmpress {input.pfam}'"; exit 1; }}
            done
            # -o --notextw: per-domain alignments, read by validate_prometheus.py
            # --hmmscan-out for the Walker A / motif II / nuclease-core checks (as archaea)
            hmmscan --cpu {threads} --cut_ga --notextw -o {output.hmmout} \
                    --domtblout {output.domtbl} {input.pfam} {input.candidates}
        else
            echo "no bacterial candidates: hmmscan not run"
            printf "# no candidates\n" > {output.domtbl}  # nonfatal: zero bacterial Prometheus candidates
            : > {output.hmmout}  # nonfatal: zero bacterial Prometheus candidates
        fi
        """


rule bact_validation_tree:
    """Place the bacterial candidates against the SAME references with the SAME settings
    as the archaeal prometheus_validation_tree (combine_proteomes.py reference mode ->
    MAFFT -> ClipKIT kpic-smart-gap -> IQ-TREE model + UFBoot2 + SH-aLRT + --bnni, fixed
    seed). validate_prometheus.py places a candidate only in a clade whose supports reach
    prometheus.placement_min_alrt/_ufboot, so a support-less (-fast) bacterial tree would
    make every bacterial locus 'unsupported' and tier A unreachable in bacteria only."""
    input:
        candidates=rules.bact_extract.output.candidates,
        prom_ref=cfg("refs/prometheus_faa", "resources/refs/prometheus_reference.faa"),
        para_ref=cfg("refs/paralog_faa", "resources/refs/paralog_reference.faa"),
    output:
        tree=f"{BD}/prometheus/validation_tree/validation.treefile",
        labels=f"{BD}/prometheus/validation_tree/reference_labels.tsv",
    params:
        outdir=f"{BD}/prometheus/validation_tree",
        model=cfg("prometheus/validation_tree_model", "MFP"),
        ufboot=cfg("tree/ufboot", 1000),
        alrt=cfg("tree/alrt", 1000),
        seed=cfg("seed", 12345),
    log: f"{LOGS}/bacteria/validation_tree.log"
    conda: "../envs/phylo.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=1440
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        python workflow/scripts/py/combine_proteomes.py --candidates {input.candidates} \
            --ref Prometheus_ref={input.prom_ref} --ref Paralog_ref={input.para_ref} \
            --out-faa {params.outdir}/all.faa --out-labels {output.labels} \
            --out-map {params.outdir}/reference_map.tsv
        n_cand=$(""" + _AWK_COUNT + r""" {input.candidates})
        n_all=$(""" + _AWK_COUNT + r""" {params.outdir}/all.faa)
        if [ "$n_cand" -eq 0 ] || [ "$n_all" -lt 4 ]; then
            echo "placeholder tree: $n_cand candidate(s), $n_all sequence(s)"
            echo "();" > {output.tree}  # nonfatal: no candidates or < 4 sequences, no testable edge
            exit 0
        fi
        # --threadit 0: MAFFT's multithreaded iterative refinement differs run to run (MAFFT multithreading docs); this keeps it reproducible
        mafft --auto --thread {threads} --threadit 0 {params.outdir}/all.faa > {params.outdir}/all.aln
        [ -s {params.outdir}/all.aln ] || {{ echo "ERROR: MAFFT wrote no alignment"; exit 1; }}
        rm -f {params.outdir}/all.trim.aln
        # ClipKIT 2.14 exits 0 (writing nothing) on a missing/unreadable input
        clipkit {params.outdir}/all.aln -m kpic-smart-gap -o {params.outdir}/all.trim.aln
        [ -s {params.outdir}/all.trim.aln ] || {{ echo "ERROR: ClipKIT wrote no alignment"; exit 1; }}
        iqtree2 -s {params.outdir}/all.trim.aln -m {params.model} -B {params.ufboot} \
            --alrt {params.alrt} --bnni -T {threads} --seed {params.seed} \
            --prefix {params.outdir}/validation -redo
        [ -s {output.tree} ] || {{ echo "ERROR: IQ-TREE wrote no tree"; exit 1; }}
        """


rule bact_validate:
    input:
        sources=rules.bact_extract.output.sources,
        raw=rules.bact_extract.output.raw,
        domtbl=rules.bact_domain_scan.output.domtbl,
        hmmout=rules.bact_domain_scan.output.hmmout,
        tree=rules.bact_validation_tree.output.tree,
        ref_labels=rules.bact_validation_tree.output.labels,
        candidates=rules.bact_extract.output.candidates,
        meta=rules.bact_samples.output.table,
        reps=rules.bact_samples.output.reps,
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
        genes_dir=f"{BD}/genes",
        ctx=cfg("prometheus/context_window", 10),
        min_alrt=cfg("prometheus/placement_min_alrt", 80),
        min_ufboot=cfg("prometheus/placement_min_ufboot", 95),
        system=BACT_FOCAL,
    log: f"{LOGS}/bacteria/validate.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        python workflow/scripts/py/validate_prometheus.py \
            --sources {input.sources} --raw {input.raw} --domtbl {input.domtbl} \
            --hmmscan-out {input.hmmout} \
            --tree {input.tree} --ref-labels {input.ref_labels} \
            --candidates {input.candidates} --meta {input.meta} --reps {input.reps} \
            --df-dir {params.df_dir} --genes-dir {params.genes_dir} \
            --context-window {params.ctx} --system "{params.system}" \
            --min-alrt {params.min_alrt} --min-ufboot {params.min_ufboot} \
            --domains '{params.domains_json}' --min-len {params.min_len} --max-len {params.max_len} \
            --out-catalog {output.catalog} --out-proteins {output.proteins} \
            --out-sensitivity {output.sensitivity} > {log} 2>&1
        """


# --------------------------------------------------------------------- #
# 3. Bacterial table + matched prevalence + architecture (REUSED scripts)
# --------------------------------------------------------------------- #
rule bact_table:
    input:
        taxonomy=_bact("taxonomy"),
        qc=_bact("qc"),
        samples=rules.bact_samples.output.table,
        catalog=rules.bact_validate.output.catalog,
        split_done=rules.bact_extract.output.split_done,
    output:
        master=f"{BD}/bacterial_master_table.tsv",
    params:
        defense_dir=f"{BD}/defense",
        min_comp=_bact("min_completeness", 90.0),
        # matched to the archaeal analysis set (preprocess keeps contamination <= 5)
        max_contam=_bact("max_contamination", cfg("preprocess/quality/max_contamination", 5.0)),
        genomes_dir=_bact("genomes_dir", ""),
        genome_ext=_bact("genome_extension", "fasta"),
        system=_bact("system", "Prometheus"),
        tiers=_tiers("stats/confidence_tiers_primary", "A"),
        sens_tiers=_tiers("stats/confidence_tiers_sensitivity", "A,B"),
    log: f"{LOGS}/bacteria/table.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/build_bacterial_table.py \
            --taxonomy {input.taxonomy} --qc {input.qc} --catalog {input.catalog} \
            --samples {input.samples} --defense-dir {params.defense_dir} \
            --min-completeness {params.min_comp} --max-contamination {params.max_contam} \
            --genomes-dir "{params.genomes_dir}" --genome-ext {params.genome_ext} \
            --system "{params.system}" --primary-tiers "{params.tiers}" \
            --sensitivity-tiers "{params.sens_tiers}" --out {output.master} > {log} 2>&1
        """


rule bact_prevalence:
    input:
        master=rules.bact_table.output.master,
    output:
        prevalence=f"{BD}/distribution/bacterial_prevalence_by_taxon.tsv",
    params:
        tiers=_tiers("stats/confidence_tiers_primary", "A"),
        sens_tiers=_tiers("stats/confidence_tiers_sensitivity", "A,B"),
    log: f"{LOGS}/bacteria/prevalence.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/compute_prevalence.py \
            --master {input.master} --primary-tiers "{params.tiers}" \
            --sensitivity-tiers "{params.sens_tiers}" --out {output.prevalence} > {log} 2>&1
        """


rule bact_architecture:
    input:
        domtbl=rules.bact_domain_scan.output.domtbl,
        catalog=rules.bact_validate.output.catalog,
    output:
        typology=f"{BD}/diversity/bacterial_architecture_typology.tsv",
    params:
        domains_json=lambda wc: json_domains(),
    log: f"{LOGS}/bacteria/architecture.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/architecture_typology.py \
            --domtbl {input.domtbl} --catalog {input.catalog} \
            --domains '{params.domains_json}' --out {output.typology} > {log} 2>&1
        """


# --------------------------------------------------------------------- #
# 4. Bacterial species tree (prune GTDB bac120) + defense matrix + co-occurrence
# --------------------------------------------------------------------- #
rule bact_species_tree:
    """GTDB bac120 pruned to the bacterial genomes, tips relabelled to genome_id.
    Without bacteria.gtdb_bac_tree the phylogeny-controlled bacterial tests are skipped
    (documented optional input); a configured tree that matches too few genomes fails."""
    input:
        samples=rules.bact_samples.output.table,
        qc=_bact("qc"),
        taxonomy=_bact("taxonomy"),
        ref_tree=lambda wc: _opt_input("gtdb_bac_tree"),
        ref_taxonomy=lambda wc: _opt_input("gtdb_bac_taxonomy"),
    output:
        tree=f"{BD}/scaffolds/bacterial_species_tree.treefile",
        report=f"{BD}/scaffolds/bacterial_species_tree_tips.tsv",
    params:
        min_frac=_bact("min_tree_match_fraction", 0.9),
    log: f"{LOGS}/bacteria/species_tree.log"
    conda: "../envs/phylo.yaml"
    threads: T_SMALL
    resources: mem_mb=16000, runtime=120
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.tree})
        if [ -n "{input.ref_tree}" ]; then
            python workflow/scripts/py/prune_tree_to_taxa.py \
                --in-tree {input.ref_tree} --keep {input.samples} --quality {input.qc} \
                --taxonomy {input.taxonomy} --ref-taxonomy "{input.ref_taxonomy}" \
                --min-match-fraction {params.min_frac} \
                --out-tree {output.tree} --out-report {output.report}
        else
            echo "bacteria.gtdb_bac_tree not set: no bacterial species tree"
            echo "();" > {output.tree}  # nonfatal: optional bacteria.gtdb_bac_tree not configured
            printf "genome_id\ttree_tip\tmatch\tnote\n" > {output.report}  # nonfatal: no reference tree
        fi
        """


rule bact_defense_matrix:
    input:
        master=rules.bact_table.output.master,
        split_done=rules.bact_extract.output.split_done,
    output:
        matrix=f"{BD}/distribution/bacterial_defense_matrix.tsv",
    params:
        df_dir=f"{BD}/defense",
    log: f"{LOGS}/bacteria/defense_matrix.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/build_defense_matrix.py \
            --master {input.master} --df-dir {params.df_dir} --require-all \
            --out {output.matrix} > {log} 2>&1
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
        primary_tiers=_tiers("stats/confidence_tiers_primary", "A"),
        min_count=cfg("stats/cooccurrence_min_count", 5),
        btol=cfg("stats/phyloglm_btol", 10),
    log: f"{LOGS}/bacteria/cooccurrence.log"
    conda: "../envs/rstats.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=240
    shell:
        r"""
        mkdir -p $(dirname {output.cooc})
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
        Rscript workflow/scripts/R/cooccurrence.R \
            --matrix {input.matrix} --master {input.master} --tree {input.tree} \
            --fdr {params.fdr} --seed {params.seed} --primary-tiers "{params.primary_tiers}" \
            --min-count {params.min_count} --btol {params.btol} --threads {threads} \
            --out {output.cooc} > {log} 2>&1
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
        fasta=f"{CD}/align/crossdomain.in.faa",
        aln=f"{CD}/align/crossdomain.aln",
        trim=f"{CD}/align/crossdomain.trim.aln",
        labels=f"{CD}/align/reference_labels.tsv",
    log: f"{LOGS}/crossdomain/align.log"
    conda: "../envs/phylo.yaml"
    threads: T_LARGE
    resources: mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=16000, per_attempt=16000), runtime=720
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.aln})
        """ + _awk_labels("Prometheus_ref") + r""" {input.prom_ref} > {output.labels}
        """ + _awk_labels("Paralog_ref") + r""" {input.para_ref} >> {output.labels}
        """ + _AWK_FASTA + r""" {input.arch} {input.bact} {input.prom_ref} {input.para_ref} > {output.fasta}
        dups=$(""" + _AWK_DUPS + r""" {output.fasta})
        if [ -n "$dups" ]; then echo "duplicated sequence ids after sanitising: $dups" >&2; exit 1; fi
        n=$(""" + _AWK_COUNT + r""" {output.fasta})
        na=$(""" + _AWK_COUNT + r""" {input.arch})
        nb=$(""" + _AWK_COUNT + r""" {input.bact})
        if [ "$n" -ge 4 ] && [ "$na" -ge 1 ] && [ "$nb" -ge 1 ]; then
            # --threadit 0: MAFFT's multithreaded iterative refinement differs run to run (MAFFT multithreading docs); this keeps it reproducible
            mafft --auto --thread {threads} --threadit 0 {output.fasta} > {output.aln}
            [ -s {output.aln} ] || {{ echo "ERROR: MAFFT wrote no alignment"; exit 1; }}
            # ClipKIT 2.14 exits 0 (writing nothing) on a missing/unreadable input
            clipkit {output.aln} -m kpic-smart-gap -o {output.trim}
            [ -s {output.trim} ] || {{ echo "ERROR: ClipKIT wrote no alignment"; exit 1; }}
        else
            echo "archaea=$na bacteria=$nb total=$n: no cross-domain alignment"
            cp {output.fasta} {output.aln}  # nonfatal: one domain has no validated Prometheus or < 4 sequences
            cp {output.fasta} {output.trim}  # nonfatal: as above (the tree rule writes '();')
        fi
        """


rule crossdomain_tree:
    input:
        trim=rules.crossdomain_align.output.trim,
        arch=rules.validate_prometheus.output.proteins,
        bact=rules.bact_validate.output.proteins,
    output:
        tree=f"{CD}/tree/crossdomain.treefile",
        ufboot=f"{CD}/tree/crossdomain.ufboot",
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
        exec > {log} 2>&1
        mkdir -p {params.workdir}
        n=$(""" + _AWK_COUNT + r""" {input.trim})
        na=$(""" + _AWK_COUNT + r""" {input.arch})
        nb=$(""" + _AWK_COUNT + r""" {input.bact})
        if [ "$n" -ge 4 ] && [ "$na" -ge 1 ] && [ "$nb" -ge 1 ]; then
            # IQ-TREE writes {output.tree} and {output.ufboot} itself (prefix = workdir/crossdomain)
            iqtree2 -s {input.trim} -m {params.model} -B {params.ufboot} -alrt {params.alrt} \
                    --wbtl -T {threads} --seed {params.seed} --prefix {params.workdir}/crossdomain -redo
            [ -s {output.tree} ] && [ -s {output.ufboot} ] || {{ echo "ERROR: IQ-TREE wrote no tree/UFBoot trees"; exit 1; }}
        else
            echo "archaea=$na bacteria=$nb total=$n: no cross-domain tree"
            echo "();" > {output.tree}  # nonfatal: one domain has no validated Prometheus or < 4 sequences
            : > {output.ufboot}  # nonfatal: no tree, so no UFBoot replicates
        fi
        """


rule crossdomain_codon_align:
    """Codon alignment of the archaeal + bacterial validated Prometheus. codon_align.py
    gets the COMBINED validated protein list (archaeal + bacterial), so neither domain is
    dropped, and a CDS file restricted to those ids (not the whole archaeal CDS set)."""
    input:
        prot_aln=rules.crossdomain_align.output.aln,
        arch=rules.validate_prometheus.output.proteins,    # archaeal (read-only)
        bact=rules.bact_validate.output.proteins,
        arch_cds=rules.combined_cds.output.ffn,            # archaeal CDS (read-only)
        bact_cds=rules.bact_extract.output.cds,
    output:
        codon=f"{CD}/selection/crossdomain.codon.fasta",
        cds=f"{CD}/selection/crossdomain_all_cds.ffn",
        proteins=f"{CD}/selection/crossdomain_proteins.faa",
        report=f"{CD}/selection/crossdomain_codon_report.tsv",
    params:
        bact_cds=str(BACT_CDS),
        max_drop=_cd("max_codon_drop_frac", cfg("evolution/max_cds_drop_frac", 0.05)),
        table=cfg("evolution/codon_table", 11),
    log: f"{LOGS}/crossdomain/codon_align.log"
    conda: "../envs/selection.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.codon})
        cat {input.arch} {input.bact} > {output.proteins}
        awk 'NR==FNR {{ if (/^>/) keep[substr($1,2)] = 1; next }}
             /^>/ {{ p = (substr($1,2) in keep) }} p' \
            {output.proteins} {input.arch_cds} {input.bact_cds} > {output.cds}
        if [ "{params.bact_cds}" = "True" ]; then
            python workflow/scripts/py/codon_align.py \
                --prot-aln {input.prot_aln} --cds {output.cds} --proteins {output.proteins} \
                --codon-table {params.table} --max-drop-frac {params.max_drop} \
                --report {output.report} --out {output.codon}
        else
            echo "bacteria.cds_dir is not set: no bacterial CDS, so no cross-domain codon alignment"
            : > {output.codon}  # nonfatal: cross-domain dN/dS not applicable without bacterial CDS
            printf "id\tstatus\tin_proteins\tn_aa\tn_codons\tnote\n" > {output.report}  # nonfatal: codon alignment not built (no bacterial CDS)
        fi
        """


rule crossdomain_analysis:
    input:
        tree=rules.crossdomain_tree.output.tree,
        ufboot=rules.crossdomain_tree.output.ufboot,
        arch=rules.validate_prometheus.output.proteins,
        bact=rules.bact_validate.output.proteins,
        labels=rules.crossdomain_align.output.labels,
        codon=rules.crossdomain_codon_align.output.codon,
    output:
        summary=f"{CD}/crossdomain_summary.tsv",
        tips=f"{CD}/crossdomain_tip_domains.tsv",
        clades=f"{CD}/crossdomain_archaeal_clades.tsv",
        rooted=f"{CD}/tree/crossdomain_rooted.treefile",
        saturation=f"{CD}/selection/crossdomain_saturation.tsv",
    params:
        min_ufboot=_cd("min_ufboot", 95),
        min_alrt=_cd("min_alrt", 80),
        n_subsample=_cd("n_subsample", 200),
        subsample_min=_cd("subsample_min", 5),
        robust=_cd("robust_fraction", 0.95),
        ref_domain=_cd("ref_domain", "Bacteria"),
        n_pairs=_cd("saturation_pairs", 200),
        sat_ds=_cd("saturation_ds", 3.0),
        seed=cfg("seed", 12345),
    log: f"{LOGS}/crossdomain/analysis.log"
    conda: "../envs/phylo.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        python workflow/scripts/py/crossdomain_analysis.py \
            --tree {input.tree} --ufboot {input.ufboot} \
            --arch-proteins {input.arch} --bact-proteins {input.bact} --ref-labels {input.labels} \
            --min-ufboot {params.min_ufboot} --min-alrt {params.min_alrt} \
            --n-subsample {params.n_subsample} --subsample-min {params.subsample_min} \
            --robust-fraction {params.robust} \
            --ref-domain {params.ref_domain} --seed {params.seed} \
            --codon {input.codon} --n-pairs {params.n_pairs} --saturation-ds {params.sat_ds} \
            --out {output.summary} --out-tips {output.tips} --out-clades {output.clades} \
            --out-rooted-tree {output.rooted} --out-saturation {output.saturation} > {log} 2>&1
        """


rule plot_crossdomain:
    input:
        tree=rules.crossdomain_analysis.output.rooted,
        tips=rules.crossdomain_analysis.output.tips,
        summary=rules.crossdomain_analysis.output.summary,
    output:
        expand(f"{CD}/figures/crossdomain_tree.{{ext}}", ext=PLOT_FORMATS),
    params:
        prefix=f"{CD}/figures/crossdomain_tree",
        formats=",".join(PLOT_FORMATS),
        min_ufboot=_cd("min_ufboot", 95),
        min_alrt=_cd("min_alrt", 80),
    log: f"{LOGS}/crossdomain/plot_tree.log"
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        python workflow/scripts/py/plot_crossdomain_tree.py \
            --tree {input.tree} --tips {input.tips} --summary {input.summary} \
            --min-ufboot {params.min_ufboot} --min-alrt {params.min_alrt} \
            --formats {params.formats} --out-prefix {params.prefix} > {log} 2>&1
        """


# --------------------------------------------------------------------- #
# 6. Cross-domain dN/dS (codon alignment of archaeal+bacterial Prometheus)
# --------------------------------------------------------------------- #
rule crossdomain_selection_tree:
    """Cross-domain gene tree pruned to exactly the codon-alignment taxa (the reference
    and paralog tips have no CDS), so HyPhy always gets a matching tree."""
    input:
        tree=rules.crossdomain_tree.output.tree,
        codon=rules.crossdomain_codon_align.output.codon,
    output:
        tree=f"{CD}/selection/crossdomain_selection.nwk",
    log: f"{LOGS}/crossdomain/selection_tree.log"
    conda: "../envs/phylo.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        exec > {log} 2>&1
        nc=$(""" + _AWK_COUNT + r""" {input.codon})
        ncomma=$(awk '{{n += gsub(/,/, ",")}} END{{print n+0}}' {input.tree})
        if [ "$nc" -ge 3 ] && [ "$ncomma" -ge 2 ]; then
            python workflow/scripts/py/prune_tree_to_taxa.py --in-tree {input.tree} \
                --keep-fasta {input.codon} --match exact --no-relabel --min-match-fraction 1.0 \
                --out-tree {output.tree}
        else
            echo "codon records=$nc: nothing to test"
            echo "();" > {output.tree}  # nonfatal: < 3 codon-aligned sequences or no cross-domain tree
        fi
        """


rule crossdomain_selection:
    input:
        codon=rules.crossdomain_codon_align.output.codon,
        gene_tree=rules.crossdomain_selection_tree.output.tree,
        arch_catalog=rules.validate_prometheus.output.catalog,
        bact_catalog=rules.bact_validate.output.catalog,
        master=rules.build_master_table.output.master,     # archaeal master (RELAX off here)
        arch_domtbl=rules.prometheus_domain_scan.output.domtbl,
        bact_domtbl=rules.bact_domain_scan.output.domtbl,
    output:
        summary=f"{CD}/selection/crossdomain_selection_summary.tsv",
        sites=f"{CD}/selection/crossdomain_site_table.tsv",
    params:
        run=str(bool(RUN_SELECTION) and BACT_CDS),
        no_cds=str(bool(RUN_SELECTION) and not BACT_CDS),
        methods="FEL,MEME,BUSTED,aBSREL",          # cross-domain: no methanogen RELAX partition
        intact_only=str(cfg("evolution/selection/intact_only", True)),
        site_pvalue=cfg("evolution/selection/site_pvalue", 0.1),
        gard=str(cfg("evolution/selection/gard", False)),
        gard_mode=cfg("evolution/selection/gard_mode", "Faster"),
        seed=cfg("seed", 12345),
        domains_json=lambda wc: json_domains(),
        workdir=f"{CD}/selection",
    log: f"{LOGS}/crossdomain/selection.log"
    conda: "../envs/selection.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=1440
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.workdir}
        awk 'FNR == 1 && NR != 1 {{ next }} {{ print }}' {input.arch_catalog} {input.bact_catalog} \
            > {params.workdir}/combined_catalog.tsv
        cat {input.arch_domtbl} {input.bact_domtbl} > {params.workdir}/combined.domtbl
        python workflow/scripts/py/run_selection.py \
            --codon {input.codon} --gene-tree {input.gene_tree} \
            --catalog {params.workdir}/combined_catalog.tsv --master {input.master} \
            --domtbl {params.workdir}/combined.domtbl --domains '{params.domains_json}' \
            --methods {params.methods} --relax-partition none \
            --gard {params.gard} --gard-mode {params.gard_mode} \
            --intact-only {params.intact_only} --site-pvalue {params.site_pvalue} --run {params.run} \
            --threads {threads} --seed {params.seed} --workdir {params.workdir} \
            --out-summary {output.summary} --out-sites {output.sites}
        if [ "{params.no_cds}" = "True" ]; then
            # run_selection.py --run False labels its rows 'selection disabled'; the actual
            # reason here is the missing bacterial CDS, so the note column says so
            awk -F'\t' -v OFS='\t' -v why="not run: bacteria.cds_dir not set (cross-domain dN/dS needs bacterial CDS)" \
                'NR == 1 {{ for (i = 1; i <= NF; i++) if ($i == "note") c = i; if (!c) exit 3; print; next }}
                 {{ $c = why; print }}' {output.summary} > {output.summary}.tmp
            mv {output.summary}.tmp {output.summary}
        fi
        """


rule plot_crossdomain_extra:
    """Bacteria vs archaea: prevalence (Wilson CI, Contract-v2 presence rule) +
    domain-architecture fractions."""
    input:
        arch_master=rules.build_master_table.output.master,
        bact_master=rules.bact_table.output.master,
        arch_typology=rules.architecture_typology.output.typology,
        bact_typology=rules.bact_architecture.output.typology,
    output:
        expand(f"{CD}/figures/crossdomain_comparison.{{ext}}", ext=PLOT_FORMATS),
        table=f"{CD}/figures/crossdomain_comparison.tsv",
    params:
        outdir=f"{CD}/figures",
        domains_json=lambda wc: json_domains(),
        tiers=_tiers("stats/confidence_tiers_primary", "A"),
        formats=",".join(PLOT_FORMATS),
    log: f"{LOGS}/crossdomain/plot_comparison.log"
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_crossdomain_extra.py \
            --arch-master {input.arch_master} --bact-master {input.bact_master} \
            --arch-typology {input.arch_typology} --bact-typology {input.bact_typology} \
            --domains '{params.domains_json}' --primary-tiers "{params.tiers}" \
            --formats {params.formats} --out-table {output.table} \
            --outdir {params.outdir} > {log} 2>&1
        """


# --------------------------------------------------------------------- #
# 7. OPTIONAL: anaerobe / environment specificity control
#    Activates only when bacteria.niche_metadata is set (gated in crossdomain_targets).
#    Tests whether bacterial Prometheus tracks an oxygen/biome signal in the SAME
#    phyloglm framework as archaea -- a clean control because bacteria can't be
#    methanogens, so a niche signal here means the archaeal methanogen association may
#    be partly a niche effect. Reuses ecology_stats.R unchanged (no --viral: only the
#    niche_* predictors are informative; the viral terms are reported as skipped).
# --------------------------------------------------------------------- #
rule bact_niche_master:
    """Niche labels -> niche_* indicators on the annotated bacterial genomes (python env;
    the R env of bact_niche has no python/pandas)."""
    input:
        master=rules.bact_table.output.master,
        niche=lambda wc: _opt_input("niche_metadata"),
    output:
        master=f"{BD}/ecology/bacterial_niche_master.tsv",
    log: f"{LOGS}/bacteria/niche_master.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/build_bacterial_niche_master.py \
            --master {input.master} --niche "{input.niche}" --out {output.master} > {log} 2>&1
        """


rule bact_niche:
    input:
        master=rules.bact_niche_master.output.master,
        tree=rules.bact_species_tree.output.tree,
    output:
        assoc=f"{BD}/ecology/bacterial_niche_associations.tsv",
        assoc_sens=f"{BD}/ecology/bacterial_niche_associations_sensitivity.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
        primary_tiers=_tiers("stats/confidence_tiers_primary", "A"),
        sens_tiers=_tiers("stats/confidence_tiers_sensitivity", "A,B"),
        btol=cfg("stats/phyloglm_btol", 10),
        covariates=_csv(_bact("covariates", cfg("ecology/covariates", ["completeness", "genome_size"])),
                        "completeness,genome_size"),
    log: f"{LOGS}/bacteria/niche.log"
    conda: "../envs/rstats.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=180
    shell:
        r"""
        mkdir -p $(dirname {output.assoc})
        echo "== primary tiers {params.primary_tiers}" > {log}
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --tree {input.tree} \
            --primary-tiers "{params.primary_tiers}" --covariates '{params.covariates}' \
            --fdr {params.fdr} --btol {params.btol} --seed {params.seed} \
            --out {output.assoc} >> {log} 2>&1
        echo "== sensitivity tiers {params.sens_tiers}" >> {log}
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --tree {input.tree} \
            --primary-tiers "{params.sens_tiers}" --covariates '{params.covariates}' \
            --fdr {params.fdr} --btol {params.btol} --seed {params.seed} \
            --out {output.assoc_sens} >> {log} 2>&1
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
    """Predict bacterial proviruses / viral contigs from the provided assemblies.
    geNomad 1.8 names its outputs after the input FASTA stem (= genome_id) and always
    writes {stem}_virus.fna + {stem}_virus_summary.tsv (header-only when nothing is
    found), so no placeholder is needed: a missing output means geNomad failed."""
    input:
        genome=bact_genome_fasta,
    output:
        virus_fna=f"{BD}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus.fna",
        virus_summary=f"{BD}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv",
    params:
        db=cfg("db/genomad", ""),
        outdir=f"{BD}/ecology/viral/{{genome}}",
        splits=int(cfg("ecology/genomad/splits", 0)),
    log: f"{LOGS}/bacteria/genomad/{{genome}}.log"
    conda: "../envs/viral.yaml"
    threads: T_MED
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=12000, per_attempt=8000),
        runtime=240,
    shell:
        r"""
        exec > {log} 2>&1
        if [ ! -f "{params.db}/version.txt" ]; then
            echo "geNomad database not found at '{params.db}' (config db.genomad = the genomad_db directory made by 'genomad download-database', with version.txt)" >&2
            exit 1
        fi
        mkdir -p {params.outdir}
        # --restart: geNomad 1.8.0 otherwise skips every step whose output exists, i.e. a
        # re-run after a failure silently resumes stale intermediates
        genomad end-to-end --cleanup --restart --threads {threads} --splits {params.splits} \
            {input.genome} {params.outdir} "{params.db}"
        """


rule bact_checkv:
    """Bacterial provirus completeness/quality. CheckV 1.0.3 exits 0 WITHOUT writing
    quality_summary.tsv on an empty FASTA, so a genome without viral sequences gets the
    CheckV header (genuinely empty input); a non-empty input must produce the table."""
    input:
        virus_fna=rules.bact_genomad.output.virus_fna,
    output:
        quality=f"{BD}/ecology/viral/{{genome}}/checkv/quality_summary.tsv",
    params:
        db=cfg("db/checkv", ""),
        outdir=f"{BD}/ecology/viral/{{genome}}/checkv",
        header=CHECKV_HEADER,
    log: f"{LOGS}/bacteria/checkv/{{genome}}.log"
    conda: "../envs/viral.yaml"
    threads: T_MED
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.outdir}
        if [ ! -s {input.virus_fna} ]; then
            echo "no viral sequences: CheckV not run"
            printf "{params.header}\n" > {output.quality}  # nonfatal: geNomad found no viral sequence in this genome
            exit 0
        fi
        if [ ! -d "{params.db}/genome_db" ] || [ ! -d "{params.db}/hmm_db" ]; then
            echo "CheckV database not found at '{params.db}' (config db.checkv = the checkv-db-v* directory made by 'checkv download_database')" >&2
            exit 1
        fi
        # --restart: CheckV 1.0.3 otherwise resumes from a stale tmp/; --remove_tmp drops it
        checkv end_to_end {input.virus_fna} {params.outdir} -d "{params.db}" -t {threads} \
            --restart --remove_tmp
        # CheckV exits 0 WITHOUT results on some fatal errors: require the 1.0.3 table
        if [ ! -s {output.quality} ] || [ "$(head -n 1 {output.quality})" != "$(printf '{params.header}')" ]; then
            echo "CheckV exited 0 but wrote no quality_summary.tsv with the CheckV 1.0.3 header (see the messages above)" >&2
            exit 1
        fi
        """


rule bact_rloop:
    """R-loop-forming-sequence prediction in bacterial proviruses (QmRLFS RIZ+REZ)."""
    input:
        virus_fna=rules.bact_genomad.output.virus_fna,
    output:
        rloop=f"{BD}/ecology/viral/{{genome}}/rloop.tsv",
        calls=f"{BD}/ecology/viral/{{genome}}/rloop_rlfs_calls.tsv",
    params:
        riz_min_g=cfg("ecology/rloop/riz_min_g", 3),
        window=cfg("ecology/rloop/window", 100),
        require_rez=str(cfg("ecology/rloop/require_rez", True)),
        min_rlfs=cfg("ecology/rloop/min_rlfs_score", 1),
        models=_tiers("ecology/rloop/models", "m1,m2"),
    log: f"{LOGS}/bacteria/rloop/{{genome}}.log"
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/rloop_predict.py \
            --fasta {input.virus_fna} --genome {wildcards.genome} \
            --riz-min-g {params.riz_min_g} --window {params.window} \
            --require-rez {params.require_rez} --min-rlfs {params.min_rlfs} \
            --models {params.models} --calls {output.calls} --out {output.rloop} > {log} 2>&1
        """


rule bact_viral_summary:
    """Aggregate bacterial proviral load + CheckV + R-loop per genome, over exactly the
    genomes the per-genome viral rules ran on (bacterial_samples.tsv = BACT_GENOMES)."""
    input:
        virus=expand(f"{BD}/ecology/viral/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv", genome=BACT_GENOMES),
        checkv=expand(f"{BD}/ecology/viral/{{genome}}/checkv/quality_summary.tsv", genome=BACT_GENOMES),
        rloop=expand(f"{BD}/ecology/viral/{{genome}}/rloop.tsv", genome=BACT_GENOMES),
        samples=rules.bact_samples.output.table,
    output:
        summary=f"{BD}/ecology/viral/bacterial_viral_summary.tsv",
    params:
        viral_dir=f"{BD}/ecology/viral",
    log: f"{LOGS}/bacteria/viral_summary.log"
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/aggregate_viral.py \
            --viral-dir {params.viral_dir} --samples {input.samples} \
            --out {output.summary} > {log} 2>&1
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
        assoc_sens=f"{BD}/ecology/bacterial_viral_associations_sensitivity.tsv",
    params:
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
        primary_tiers=_tiers("stats/confidence_tiers_primary", "A"),
        sens_tiers=_tiers("stats/confidence_tiers_sensitivity", "A,B"),
        btol=cfg("stats/phyloglm_btol", 10),
        covariates=_csv(_bact("covariates", cfg("ecology/covariates", ["completeness", "genome_size"])),
                        "completeness,genome_size"),
    log: f"{LOGS}/bacteria/viral_stats.log"
    conda: "../envs/rstats.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        mkdir -p $(dirname {output.assoc})
        echo "== primary tiers {params.primary_tiers}" > {log}
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --viral {input.viral} --tree {input.tree} \
            --primary-tiers "{params.primary_tiers}" --covariates '{params.covariates}' \
            --fdr {params.fdr} --btol {params.btol} --seed {params.seed} \
            --out {output.assoc} >> {log} 2>&1
        echo "== sensitivity tiers {params.sens_tiers}" >> {log}
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --viral {input.viral} --tree {input.tree} \
            --primary-tiers "{params.sens_tiers}" --covariates '{params.covariates}' \
            --fdr {params.fdr} --btol {params.btol} --seed {params.seed} \
            --out {output.assoc_sens} >> {log} 2>&1
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
        table=f"{CD}/figures/viral_rloop_crossdomain.tsv",
    params:
        outdir=f"{CD}/figures",
        tiers=_tiers("stats/confidence_tiers_primary", "A"),
        formats=",".join(PLOT_FORMATS),
    log: f"{LOGS}/crossdomain/plot_viral.log"
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_viral_crossdomain.py \
            --arch-master {input.arch_master} --arch-viral {input.arch_viral} \
            --bact-master {input.bact_master} --bact-viral {input.bact_viral} \
            --bact-assoc {input.bact_assoc} --primary-tiers "{params.tiers}" \
            --formats {params.formats} --out-table {output.table} \
            --outdir {params.outdir} > {log} 2>&1
        """
