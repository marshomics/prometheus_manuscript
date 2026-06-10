# =====================================================================
# Intrinsic-disorder module (OPT-IN: disorder.enabled: true).
# Consensus per-residue disorder (metapredict + IUPred3, or the alphaPredict
# fallback) on intact Prometheus loci, anchored to the family alignment, compared
# across methanogens/niches/lineages with phylolm, and (when the bacterial branch
# is also on) bacteria vs archaea. Built with the `disorder` target; reuses the
# archaeal validated proteins / catalog / alignment / species tree read-only.
# =====================================================================

DIS = f"{RESULTS}/disorder"


def _d(key, default=""):
    return cfg(f"disorder/{key}", default)


rule disorder_predict:
    """Consensus per-residue disorder on the validated (intact) archaeal Prometheus loci."""
    input:
        proteins=rules.validate_prometheus.output.proteins,
        catalog=rules.validate_prometheus.output.catalog,
    output:
        locus=f"{DIS}/disorder_per_locus.tsv",
        residue=f"{DIS}/disorder_per_residue.tsv",
    params:
        intact_only=str(_d("intact_only", True)),
        meta_thr=_d("metapredict_threshold", 0.5),
        iupred3_path=_d("iupred3_path", ""),
        iupred3_mode=_d("iupred3_mode", "long"),
        iupred3_thr=_d("iupred3_threshold", 0.5),
        plddt=_d("alphapredict_plddt_cutoff", 50),
        min_idr=_d("min_idr_len", 12),
        mask_lcr=str(_d("mask_low_complexity", True)),
        host_bg=_d("host_background", ""),
    log: f"{LOGS}/disorder/predict.log"
    conda: "../envs/disorder.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=240
    shell:
        r"""
        python workflow/scripts/py/predict_disorder.py \
            --proteins {input.proteins} --catalog {input.catalog} \
            --intact-only {params.intact_only} --meta-threshold {params.meta_thr} \
            --iupred3-path "{params.iupred3_path}" --iupred3-mode {params.iupred3_mode} \
            --iupred3-threshold {params.iupred3_thr} --plddt-cutoff {params.plddt} \
            --min-idr-len {params.min_idr} --mask-lcr {params.mask_lcr} \
            --host-background "{params.host_bg}" \
            --out-locus {output.locus} --out-residue {output.residue} > {log} 2>&1
        """


rule disorder_map:
    """Map per-residue consensus disorder onto the family alignment + domain track."""
    input:
        residue=rules.disorder_predict.output.residue,
        alignment=rules.family_align.output.aln,
        domtbl=rules.prometheus_domain_scan.output.domtbl,
    output:
        bycol=f"{DIS}/disorder_by_column.tsv",
    params:
        domains_json=lambda wc: json_domains(),
    conda: "../envs/annotation.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=60
    shell:
        r"""
        python workflow/scripts/py/map_disorder_to_alignment.py \
            --residue {input.residue} --alignment {input.alignment} \
            --domtbl {input.domtbl} --domains '{params.domains_json}' --out {output.bycol}
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
        intact_only=str(_d("intact_only", True)),
        meta_thr=_d("metapredict_threshold", 0.5),
        iupred3_path=_d("iupred3_path", ""),
        iupred3_mode=_d("iupred3_mode", "long"),
        iupred3_thr=_d("iupred3_threshold", 0.5),
        plddt=_d("alphapredict_plddt_cutoff", 50),
        min_idr=_d("min_idr_len", 12),
        mask_lcr=str(_d("mask_low_complexity", True)),
    log: f"{LOGS}/disorder/predict_bact.log"
    conda: "../envs/disorder.yaml"
    threads: T_MED
    resources: mem_mb=16000, runtime=240
    shell:
        r"""
        python workflow/scripts/py/predict_disorder.py \
            --proteins {input.proteins} --catalog {input.catalog} \
            --intact-only {params.intact_only} --meta-threshold {params.meta_thr} \
            --iupred3-path "{params.iupred3_path}" --iupred3-mode {params.iupred3_mode} \
            --iupred3-threshold {params.iupred3_thr} --plddt-cutoff {params.plddt} \
            --min-idr-len {params.min_idr} --mask-lcr {params.mask_lcr} --host-background "" \
            --out-locus {output.locus} --out-residue {output.residue} > {log} 2>&1
        """


rule disorder_stats:
    """Disorder ~ methanogen / niche, phylogeny-controlled (phylolm), temperature-aware."""
    input:
        master=rules.build_master_table.output.master,
        disorder=rules.disorder_predict.output.locus,
        tree=rules.species_tree.output.tree,
    output:
        assoc=f"{DIS}/disorder_associations.tsv",
    params:
        temp_cols=",".join(_d("temperature_cols", []) or []),
        fdr=cfg("stats/fdr_method", "BH"),
        seed=cfg("seed", 12345),
    log: f"{LOGS}/disorder/stats.log"
    conda: "../envs/rstats.yaml"
    threads: T_SMALL
    resources: mem_mb=8000, runtime=120
    shell:
        r"""
        Rscript workflow/scripts/R/disorder_stats.R \
            --master {input.master} --disorder {input.disorder} --tree {input.tree} \
            --temp-cols "{params.temp_cols}" --fdr {params.fdr} --seed {params.seed} \
            --out {output.assoc} > {log} 2>&1
        """


def _bact_disorder_locus(wc):
    return rules.disorder_predict_bact.output.locus if (DISORDER_ENABLED and BACT_ENABLED) else []


rule plot_disorder:
    """Disorder along the protein + by group + linker-vs-domain (+ cross-domain if on)."""
    input:
        bycol=rules.disorder_map.output.bycol,
        per_locus=rules.disorder_predict.output.locus,
        master=rules.build_master_table.output.master,
        bact=_bact_disorder_locus,
    output:
        expand(f"{DIS}/plots/disorder_along_protein.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{DIS}/plots/disorder_by_group.{{ext}}", ext=PLOT_FORMATS),
        expand(f"{DIS}/plots/disorder_linker_vs_domain.{{ext}}", ext=PLOT_FORMATS),
    params:
        outdir=f"{DIS}/plots",
        bact_path=(f"{DIS}/bacteria/disorder_per_locus.tsv" if (DISORDER_ENABLED and BACT_ENABLED) else ""),
    conda: "../envs/plotting.yaml"
    threads: T_SMALL
    resources: mem_mb=4000, runtime=30
    shell:
        r"""
        python workflow/scripts/py/plot_disorder.py \
            --by-column {input.bycol} --per-locus {input.per_locus} --master {input.master} \
            --bact-per-locus "{params.bact_path}" --outdir {params.outdir}
        """
