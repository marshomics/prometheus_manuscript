# =====================================================================
# Intrinsic-disorder module (OPT-IN: disorder.enabled: true).
# Consensus per-residue disorder (metapredict V3 + IUPred3, or the labelled, correlated
# metapredict + alphaPredict-pLDDT fallback) on intact, non-fragmented Prometheus loci,
# anchored to the UNtrimmed family alignment, compared across methanogens / niches with
# phylolm under the Contract-v2 tier rule, and (when the bacterial branch is also on)
# bacteria vs archaea. Built with the `disorder` target; reuses the archaeal validated
# proteins / catalog / alignment / species tree read-only.
# Every rule fails on a missing tool or a crashed predictor (Contract v2 item 4); only a
# genuinely empty input yields header-only tables / placeholder figures.
# =====================================================================

DIS = f"{RESULTS}/disorder"


def _dis_cfg(key, default=""):
    return cfg(f"disorder/{key}", default)


def _dis_host_background(wc):
    """Optional host-proteome background table, tracked as an input when configured."""
    p = _dis_cfg("host_background", "")
    return [p] if p else []


def _dis_tier_sets():
    """Contract v2 item 6 tier sets: the primary set, then the A+B sensitivity set."""
    prim = [str(t) for t in (cfg("stats/confidence_tiers_primary", ["A"]) or ["A"])]
    sens = [str(t) for t in (cfg("stats/confidence_tiers_sensitivity", ["A", "B"]) or [])]
    sets = [prim] + ([sens] if sens and set(sens) != set(prim) else [])
    return ";".join(",".join(s) for s in sets)


def _dis_predict_params():
    return dict(
        intact_only=str(_dis_cfg("intact_only", True)),
        meta_thr=_dis_cfg("metapredict_threshold", 0.5),
        meta_version=_dis_cfg("metapredict_version", "V3"),
        second=_dis_cfg("second_predictor", "auto"),
        iupred3_path=_dis_cfg("iupred3_path", ""),
        iupred3_mode=_dis_cfg("iupred3_mode", "long"),
        iupred3_thr=_dis_cfg("iupred3_threshold", 0.5),
        iupred3_smoothing=_dis_cfg("iupred3_smoothing", ""),
        plddt=_dis_cfg("alphapredict_plddt_cutoff", 50),
        min_idr=_dis_cfg("min_idr_len", 12),
        mask_lcr=str(_dis_cfg("mask_low_complexity", True)),
    )


rule disorder_predict:
    """Consensus per-residue disorder on the validated (intact) archaeal Prometheus loci."""
    input:
        proteins=rules.validate_prometheus.output.proteins,
        catalog=rules.validate_prometheus.output.catalog,
        host_bg=_dis_host_background,
    output:
        locus=f"{DIS}/disorder_per_locus.tsv",
        residue=f"{DIS}/disorder_per_residue.tsv",
    params:
        **_dis_predict_params(),
        host_bg=_dis_cfg("host_background", ""),
    log: f"{LOGS}/disorder/predict.log"
    conda: "../envs/disorder.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=240
    shell:
        r"""
        OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        python workflow/scripts/py/predict_disorder.py \
            --proteins {input.proteins} --catalog {input.catalog} \
            --intact-only {params.intact_only} --meta-threshold {params.meta_thr} \
            --metapredict-version {params.meta_version} --device cpu --threads {threads} \
            --second-predictor {params.second} \
            --iupred3-path "{params.iupred3_path}" --iupred3-mode {params.iupred3_mode} \
            --iupred3-threshold {params.iupred3_thr} --iupred3-smoothing "{params.iupred3_smoothing}" \
            --plddt-cutoff {params.plddt} \
            --min-idr-len {params.min_idr} --mask-lcr {params.mask_lcr} \
            --host-background "{params.host_bg}" \
            --out-locus {output.locus} --out-residue {output.residue} > {log} 2>&1
        """


rule disorder_map:
    """Map per-residue consensus disorder onto the UNtrimmed family alignment + domain track."""
    input:
        residue=rules.disorder_predict.output.residue,
        alignment=rules.family_align.output.aln,
        domtbl=rules.prometheus_domain_scan.output.domtbl,
    output:
        bycol=f"{DIS}/disorder_by_column.tsv",
    params:
        domains_json=lambda wc: json_domains(),
        min_support=_dis_cfg("min_domain_support", 0.5),
    log: f"{LOGS}/disorder/map.log"
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        python workflow/scripts/py/map_disorder_to_alignment.py \
            --residue {input.residue} --alignment {input.alignment} \
            --domtbl {input.domtbl} --domains '{params.domains_json}' \
            --min-domain-support {params.min_support} --out {output.bycol} > {log} 2>&1
        """


rule disorder_predict_bact:
    """Consensus disorder on the bacterial Prometheus loci (cross-domain comparison)."""
    input:
        proteins=rules.bact_validate.output.proteins,
        catalog=rules.bact_validate.output.catalog,
    output:
        locus=f"{DIS}/bacteria/disorder_per_locus.tsv",
        residue=f"{DIS}/bacteria/disorder_per_residue.tsv",
    params:
        **_dis_predict_params(),
    log: f"{LOGS}/disorder/predict_bact.log"
    conda: "../envs/disorder.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=240
    shell:
        r"""
        OMP_NUM_THREADS={threads} MKL_NUM_THREADS={threads} \
        python workflow/scripts/py/predict_disorder.py \
            --proteins {input.proteins} --catalog {input.catalog} \
            --intact-only {params.intact_only} --meta-threshold {params.meta_thr} \
            --metapredict-version {params.meta_version} --device cpu --threads {threads} \
            --second-predictor {params.second} \
            --iupred3-path "{params.iupred3_path}" --iupred3-mode {params.iupred3_mode} \
            --iupred3-threshold {params.iupred3_thr} --iupred3-smoothing "{params.iupred3_smoothing}" \
            --plddt-cutoff {params.plddt} \
            --min-idr-len {params.min_idr} --mask-lcr {params.mask_lcr} --host-background "" \
            --out-locus {output.locus} --out-residue {output.residue} > {log} 2>&1
        """


rule disorder_stats:
    """Disorder ~ methanogen / niche (+ lambda signal), phylolm, temperature-adjusted, tier rule."""
    input:
        master=rules.build_master_table.output.master,
        disorder=rules.disorder_predict.output.locus,
        tree=rules.species_tree.output.tree,
    output:
        assoc=f"{DIS}/disorder_associations.tsv",
    params:
        temp_cols=",".join(str(c) for c in (_dis_cfg("temperature_cols", []) or [])),
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
        tier_sets=_dis_tier_sets(),
        min_group=_dis_cfg("min_group_size", 3),
    log: f"{LOGS}/disorder/stats.log"
    conda: "../envs/rstats.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        Rscript workflow/scripts/R/disorder_stats.R \
            --master {input.master} --disorder {input.disorder} --tree {input.tree} \
            --temp-cols "{params.temp_cols}" --fdr {params.fdr} --seed {params.seed} \
            --tier-sets "{params.tier_sets}" --min-group {params.min_group} \
            --out {output.assoc} > {log} 2>&1
        """


def _dis_bact_locus(wc):
    return rules.disorder_predict_bact.output.locus if (DISORDER_ENABLED and BACT_ENABLED) else []


rule plot_disorder:
    """Disorder along the protein + by group + linker-vs-domain (+ cross-domain if on)."""
    input:
        bycol=rules.disorder_map.output.bycol,
        per_locus=rules.disorder_predict.output.locus,
        master=rules.build_master_table.output.master,
        bact=_dis_bact_locus,
    output:
        expand(f"{DIS}/plots/disorder_along_protein.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{DIS}/plots/disorder_by_group.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{DIS}/plots/disorder_linker_vs_domain.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{DIS}/plots",
        bact_path=(f"{DIS}/bacteria/disorder_per_locus.tsv" if (DISORDER_ENABLED and BACT_ENABLED) else ""),
        domains_json=lambda wc: json_domains(),
        formats=",".join(PLOT_FORMATS),
        tiers=",".join(str(t) for t in (cfg("stats/confidence_tiers_primary", ["A"]) or ["A"])),
        min_occ=_dis_cfg("min_column_occupancy", 0.5),
    log: f"{LOGS}/disorder/plot.log"
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_disorder.py \
            --by-column {input.bycol} --per-locus {input.per_locus} --master {input.master} \
            --bact-per-locus "{params.bact_path}" --outdir {params.outdir} \
            --domains '{params.domains_json}' --formats {params.formats} \
            --tiers "{params.tiers}" --min-occupancy {params.min_occ} > {log} 2>&1
        """
