# =====================================================================
# Stage 6 -- ecology: viral pressure (geNomad 1.8.0 + CheckV 1.0.3), R-loop-forming
# sequences in the viral sequences (QmRLFS, rloop_predict.py) and phylogenetically
# controlled niche / viral-pressure associations (phyloglm, ecology_stats.R).
# Every rule fails when its tool fails (SCHEMA Contract v2 item 4); the only placeholder
# is CheckV's header-only table when geNomad found no viral sequence (a genuinely empty
# input), marked `# nonfatal:`. Presence follows Contract v2 item 6 (tier sets below).
# =====================================================================
import sys as _eco_sys

ECO = f"{RESULTS}/ecology"
ECO_VIRAL = f"{ECO}/viral"

# CheckV 1.0.3 quality_summary.tsv header (checkv/modules/quality_summary.py)
CHECKV_QS_HEADER = "\t".join([
    "contig_id", "contig_length", "provirus", "proviral_length", "gene_count",
    "viral_genes", "host_genes", "checkv_quality", "miuvig_quality", "completeness",
    "completeness_method", "contamination", "kmer_freq", "warnings"])


def _eco_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        value = value.split(",")
    return [str(v).strip() for v in value if str(v).strip()]


def _eco_predictors():
    """--predictors for ecology_stats.R. The documented predictor set (SCHEMA, README) is
    every species-level group_*/niche_* column of the master table ("auto"); config
    ecology.predictors may restrict it to named group_/niche_ columns or globs. Legacy
    ENVO field names (biome_envo, oxygen, temp_class, ...) are never master-table columns,
    so they are dropped with a warning instead of leaving the analysis without any niche
    predictor."""
    toks = _eco_list(cfg("ecology/predictors", "auto"))
    keep = [t for t in toks if t.lower() == "auto" or t.startswith(("group_", "niche_"))]
    dropped = [t for t in toks if t not in keep]
    if dropped:
        print(f"[ecology] config ecology.predictors {dropped} are not group_/niche_ master "
              f"columns; testing {keep or ['auto']} (all group_/niche_ columns when 'auto')",
              file=_eco_sys.stderr)
    return ",".join(keep) if keep else "auto"


ECO_TIERS = ",".join(_eco_list(cfg("stats/confidence_tiers_primary", ["A"]))) or "A"
ECO_TIERS_SENS = ",".join(_eco_list(cfg("stats/confidence_tiers_sensitivity", ["A", "B"]))) or "A,B"


rule genomad:
    """Proviruses and viral contigs with geNomad 1.8.0 end-to-end (defaults: --min-score
    0.7, no score calibration). geNomad names its outputs after the input file stem, i.e.
    {genome} for the staged {genome}.<ext> assembly, and always writes both outputs
    (header-only / empty when nothing is found). --restart: never resume stale
    intermediates of an earlier failed run. ecology.genomad.splits (default 0 = geNomad's
    default) splits the MMseqs2 marker search to bound memory."""
    input:
        genome=genome_fasta,
    output:
        virus_fna=f"{ECO_VIRAL}/{{genome}}/{{genome}}_summary/{{genome}}_virus.fna",
        virus_summary=f"{ECO_VIRAL}/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv",
    params:
        db=cfg("db/genomad", ""),
        outdir=f"{ECO_VIRAL}/{{genome}}",
        splits=int(cfg("ecology/genomad/splits", 0)),
    log:
        f"{LOGS}/genomad/{{genome}}.log",
    conda:
        "../envs/viral.yaml"
    threads: T_MED
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=12000, per_attempt=8000),
        runtime=240,
    shell:
        r"""
        if [ ! -f "{params.db}/version.txt" ]; then
            echo "geNomad database not found at '{params.db}' (config db.genomad = the genomad_db directory made by 'genomad download-database', with version.txt)" > {log}
            exit 1
        fi
        mkdir -p {params.outdir}
        genomad end-to-end --cleanup --restart --threads {threads} --splits {params.splits} \
            {input.genome} {params.outdir} "{params.db}" > {log} 2>&1
        """


rule checkv:
    """CheckV 1.0.3 end_to_end on the geNomad viral sequences. CheckV exits 0 WITHOUT
    results on some fatal errors (empty or duplicate-id FASTA, missing executables), so the
    output is checked explicitly; --restart stops a re-run resuming a stale tmp/ (CheckV
    resumes by default), --remove_tmp drops the intermediates."""
    input:
        virus_fna=rules.genomad.output.virus_fna,
    output:
        quality=f"{ECO_VIRAL}/{{genome}}/checkv/quality_summary.tsv",
    params:
        db=cfg("db/checkv", ""),
        outdir=f"{ECO_VIRAL}/{{genome}}/checkv",
        header=CHECKV_QS_HEADER,
    log:
        f"{LOGS}/checkv/{{genome}}.log",
    conda:
        "../envs/viral.yaml"
    threads: T_MED
    resources:
        mem_mb=8000,
        runtime=120,
    shell:
        r"""
        mkdir -p {params.outdir}
        if [ ! -s {input.virus_fna} ]; then
            echo "geNomad reported no viral sequence for {wildcards.genome}; CheckV not run" > {log}
            printf '%s\n' "{params.header}" > {output.quality}  # nonfatal: genuinely empty input (geNomad found no viral sequence)
            exit 0
        fi
        if [ ! -d "{params.db}/genome_db" ] || [ ! -d "{params.db}/hmm_db" ]; then
            echo "CheckV database not found at '{params.db}' (config db.checkv = the checkv-db-v* directory made by 'checkv download_database')" > {log}
            exit 1
        fi
        checkv end_to_end {input.virus_fna} {params.outdir} -d "{params.db}" -t {threads} \
            --restart --remove_tmp > {log} 2>&1
        if [ ! -s {output.quality} ] || ! head -n 1 {output.quality} | grep -q '^contig_id'; then
            echo "CheckV exited 0 but wrote no quality_summary.tsv (see the messages above)" >> {log}
            exit 1
        fi
        """


rule rloop_predict:
    """QmRLFS R-loop-forming sequences in every geNomad viral sequence: rloop_predict.py
    reproduces QmRLFS-finder v2.0 call-for-call (m1+m2 RIZ >= 50 % G, linker 0-49 nt,
    REZ 100-2000 nt >= 40 % G, both strands). config ecology.rloop: riz_min_g (m1 G-run,
    3 = QmRLFS), window (minimum REZ length, 100 = QmRLFS), require_rez (true = QmRLFS),
    min_rlfs_score (RLFS loci needed to call a sequence R-loop-prone, 1), models."""
    input:
        virus_fna=rules.genomad.output.virus_fna,
    output:
        rloop=f"{ECO_VIRAL}/{{genome}}/rloop.tsv",
        calls=f"{ECO_VIRAL}/{{genome}}/rloop_rlfs_calls.tsv",
    params:
        riz_min_g=cfg("ecology/rloop/riz_min_g", 3),
        window=cfg("ecology/rloop/window", 100),
        require_rez=str(cfg("ecology/rloop/require_rez", True)),
        min_rlfs=cfg("ecology/rloop/min_rlfs_score", 1),
        models=",".join(_eco_list(cfg("ecology/rloop/models", ["m1", "m2"]))) or "m1,m2",
    log:
        f"{LOGS}/rloop/{{genome}}.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000,
        runtime=60,
    shell:
        r"""
        python workflow/scripts/py/rloop_predict.py \
            --fasta {input.virus_fna} --genome {wildcards.genome} \
            --riz-min-g {params.riz_min_g} --window {params.window} \
            --require-rez {params.require_rez} --min-rlfs {params.min_rlfs} \
            --models {params.models} --calls {output.calls} --out {output.rloop} > {log} 2>&1
        """


rule viral_summary:
    """Per-genome proviral load (integrated proviruses only), whole viral contigs, CheckV
    quality and the R-loop-prone fraction; fails on missing or inconsistent inputs."""
    input:
        virus=expand(f"{ECO_VIRAL}/{{genome}}/{{genome}}_summary/{{genome}}_virus_summary.tsv", genome=GENOMES),
        checkv=expand(f"{ECO_VIRAL}/{{genome}}/checkv/quality_summary.tsv", genome=GENOMES),
        rloop=expand(f"{ECO_VIRAL}/{{genome}}/rloop.tsv", genome=GENOMES),
        samples=config["samples"],
    output:
        summary=f"{ECO_VIRAL}/viral_summary.tsv",
    params:
        viral_dir=ECO_VIRAL,
    log:
        f"{LOGS}/ecology/viral_summary.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        r"""
        python workflow/scripts/py/aggregate_viral.py \
            --viral-dir {params.viral_dir} --samples {input.samples} --out {output.summary} > {log} 2>&1
        """


rule ecology_stats:
    """Phylogeny-controlled niche and viral-pressure associations (phyloglm per predictor,
    adjusted for completeness + genome size; BH). The primary run uses
    stats.confidence_tiers_primary; the same models are re-run on the A+B sensitivity tier
    set (stats.confidence_tiers_sensitivity) into ecology_associations_sensitivity.tsv."""
    input:
        master=rules.build_master_table.output.master,
        viral=rules.viral_summary.output.summary,
        tree=rules.species_tree.output.tree,
    output:
        assoc=f"{ECO}/ecology_associations.tsv",
        assoc_sens=f"{ECO}/ecology_associations_sensitivity.tsv",
    params:
        predictors=_eco_predictors(),
        tiers=ECO_TIERS,
        tiers_sens=ECO_TIERS_SENS,
        fdr=cfg("stats/fdr_method", "BH"),
        btol=cfg("stats/phyloglm_btol", 10),
        covariates=",".join(_eco_list(cfg("ecology/covariates", ["completeness", "genome_size"]))),
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/ecology/ecology_stats.log",
    conda:
        "../envs/rstats.yaml"
    threads: T_SMALL
    resources:
        mem_mb=16000,
        runtime=240,
    shell:
        r"""
        echo "== primary tiers {params.tiers}" > {log}
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --viral {input.viral} --tree {input.tree} \
            --predictors '{params.predictors}' --primary-tiers {params.tiers} \
            --covariates '{params.covariates}' \
            --fdr {params.fdr} --btol {params.btol} --seed {params.seed} \
            --out {output.assoc} >> {log} 2>&1
        echo "== sensitivity tiers {params.tiers_sens}" >> {log}
        Rscript workflow/scripts/R/ecology_stats.R \
            --master {input.master} --viral {input.viral} --tree {input.tree} \
            --predictors '{params.predictors}' --primary-tiers {params.tiers_sens} \
            --covariates '{params.covariates}' \
            --fdr {params.fdr} --btol {params.btol} --seed {params.seed} \
            --out {output.assoc_sens} >> {log} 2>&1
        """


ECO_MGEN_PREDICTORS = ["mgen_co2_reducing", "mgen_acetoclastic", "mgen_methylotrophic",
                       "mgen_methyl_dismutating", "mgen_methyl_reducing"]


def _eco_mgen(key, default):
    """config methanogenesis.stats.<key> (default when unset)."""
    return cfg(f"methanogenesis/stats/{key}", default)


def _eco_mgen_csv(key, default):
    """A comma list for methanogenesis_stats.R from a YAML list or a comma string."""
    return ",".join(_eco_list(_eco_mgen(key, default))) or ",".join(_eco_list(default))


def _eco_mgen_flag(key, default):
    v = _eco_mgen(key, default)
    s = str(v).strip().lower()
    if s in ("true", "1", "yes", "y", "t"):
        return "TRUE"
    if s in ("false", "0", "no", "n", "f"):
        return "FALSE"
    raise ValueError(f"config methanogenesis.stats.{key} must be true/false, got {v!r}")


rule methanogenesis_stats:
    """Design 6.2: within methanogens, does methanogenesis type (gene-based pathway
    capacity, the master's mgen_* columns) predict Prometheus presence, with type and
    habitat fitted jointly? Per predictor (methanogenesis_stats.R): marginal and
    habitat-adjusted (PRIMARY partial effect; habitat terms by |phi| under the
    events-per-variable rule, VIF guard) phyloglm, the stochastic-mapping replication
    guard, the corHMM hidden-rate correlation test and four sensitivity re-fits; BH over
    the primary rows (registry family `methanogenesis`, always exploratory). Each
    methanogenesis.stats.<key> overrides the script default; tier sets default to
    stats.confidence_tiers_primary / _sensitivity; FDR method, btol and seed are the
    pipeline-wide stats.fdr_method, stats.phyloglm_btol and seed. corHMM 2.8 is installed
    into the rstats env by rstats.post-deploy.sh (vendored, SHA-256 checked)."""
    input:
        master=rules.build_master_table.output.master,
        tree=rules.species_tree.output.tree,
    output:
        assoc=f"{ECO}/methanogenesis_associations.tsv",
    params:
        predictors=_eco_mgen_csv("predictors", ECO_MGEN_PREDICTORS),
        habitat=_eco_mgen_csv("habitat", ["auto"]),
        covariates=_eco_mgen_csv("covariates", ["completeness", "genome_size"]),
        tiers=_eco_mgen_csv("primary_tiers", ECO_TIERS),
        tiers_sens=_eco_mgen_csv("sensitivity_tiers", ECO_TIERS_SENS),
        min_n=int(_eco_mgen("min_n", 10)),
        min_class=int(_eco_mgen("min_class", 3)),
        min_changes=float(_eco_mgen("min_changes", 5)),
        nsim=int(_eco_mgen("nsim", 100)),
        corhmm=_eco_mgen_flag("corhmm", True),
        epv_min=float(_eco_mgen("epv_min", 10)),
        vif_max=float(_eco_mgen("vif_max", 5)),
        complete_threshold=float(_eco_mgen("complete_threshold", 95)),
        fdr=cfg("stats/fdr_method", "BH"),
        btol=cfg("stats/phyloglm_btol", 10),
        seed=cfg("seed", 12345),
    log:
        f"{LOGS}/ecology/methanogenesis_stats.log",
    conda:
        "../envs/rstats.yaml"
    threads: max(1, min(len(ECO_MGEN_PREDICTORS), int(T_MED)))
    resources:
        mem_mb=lambda wc, attempt: mem_mb(wc, attempt, base=16000, per_attempt=16000),
        runtime=1440,
    shell:
        r"""
        # nonfatal: methanogenesis_stats.R writes a header-only table (exit 0) only when the master table has no methanogen tip of the species tree (genuinely empty input); a missing column, package or crashed fit fails the rule
        Rscript workflow/scripts/R/methanogenesis_stats.R \
            --master {input.master} --tree {input.tree} \
            --predictors '{params.predictors}' --habitat '{params.habitat}' \
            --covariates '{params.covariates}' \
            --primary-tiers {params.tiers} --sensitivity-tiers {params.tiers_sens} \
            --min-n {params.min_n} --min-class {params.min_class} \
            --min-changes {params.min_changes} --nsim {params.nsim} --corhmm {params.corhmm} \
            --epv-min {params.epv_min} --vif-max {params.vif_max} \
            --complete-threshold {params.complete_threshold} \
            --fdr-method {params.fdr} --btol {params.btol} --seed {params.seed} \
            --threads {threads} --out {output.assoc} > {log} 2>&1
        """


rule plot_ecology:
    """Niche-association forest + viral / R-loop distributions by presence (phyloglm
    results in the titles)."""
    input:
        assoc=rules.ecology_stats.output.assoc,
        viral=rules.viral_summary.output.summary,
        master=rules.build_master_table.output.master,
    output:
        expand(f"{ECO}/plots/ecology_enrichment.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{ECO}/plots/viral_rloop_association.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{ECO}/plots",
        tiers=ECO_TIERS,
        formats=",".join(PLOT_FORMATS),
    log:
        f"{LOGS}/ecology/plot_ecology.log",
    conda:
        "../envs/plotting.yaml"
    threads: T_SMALL
    resources:
        mem_mb=4000,
        runtime=30,
    shell:
        r"""
        python workflow/scripts/py/plot_ecology.py \
            --associations {input.assoc} --viral {input.viral} --master {input.master} \
            --primary-tiers {params.tiers} --formats {params.formats} \
            --outdir {params.outdir} > {log} 2>&1
        """
