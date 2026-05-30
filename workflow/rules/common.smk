# =====================================================================
# common.smk — config, sample sheet, helpers, and the global target list.
# Kept deliberately static (no checkpoints) so `snakemake -n` validates the
# full DAG; data-dependent steps degrade gracefully on empty inputs.
# =====================================================================
import os
import json
import pandas as pd
from snakemake.utils import validate, min_version

min_version("7.25.0")

# ---- config ----------------------------------------------------------
configfile: "config/config.yaml"
try:
    validate(config, "../schemas/config.schema.yaml")
except Exception as e:
    print(f"[common] config schema warning: {e}")

RESULTS = config["results_dir"]
LOGS = config["logs_dir"]
GENOME_DIR = config["genomes_dir"]
EXT = config["genome_extension"]

def cfg(path, default=None):
    """Safe nested config getter: cfg('prometheus/min_len_aa', 800)."""
    node = config
    for key in path.split("/"):
        if isinstance(node, dict) and key in node:
            node = node[key]
        else:
            return default
    return node

# ---- samples ---------------------------------------------------------
samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index("genome_id", drop=False)
try:
    validate(samples, "../schemas/samples.schema.yaml")
except Exception as e:
    print(f"[common] samples schema warning: {e}")

GENOMES = samples["genome_id"].tolist()

def genome_fasta(wildcards):
    return os.path.join(GENOME_DIR, f"{wildcards.genome}.{EXT}")

def json_domains():
    """Serialize the Prometheus domain definitions for the validator CLI."""
    return json.dumps(cfg("prometheus/domains", {}))

# ---- resource helpers (Grid Engine via profiles/sge) -----------------
T_SMALL = cfg("threads/small", 1)
T_MED = cfg("threads/medium", 8)
T_LARGE = cfg("threads/large", 24)

def mem_mb(wildcards, attempt, base=4000, per_attempt=4000):
    return base + (attempt - 1) * per_attempt

def runtime_min(default=240):
    return default

# ---- which detection tools to run vs ingest --------------------------
RUN_DEFENSEFINDER = cfg("prometheus/run_defensefinder", True)
DEREPLICATE = cfg("dereplicate", True)
PER_DOMAIN_TREES = cfg("evolution/per_domain_trees", True)
DOMAINS = list(cfg("prometheus/domains", {"DUF4011": {}, "DUF3320": {}, "DNA2": {}, "nuclease": {}, "REase": {}}).keys())
RUN_GENERAX = cfg("evolution/reconciliation/run_generax", True)
RUN_ALE = cfg("evolution/reconciliation/run_ale", False)
RUN_SELECTION = cfg("evolution/selection/run", True)
SELECTION_METHODS = cfg("evolution/selection/methods", ["FEL", "MEME", "BUSTED", "aBSREL", "RELAX"])
PLOT_FORMATS = cfg("plotting/formats", ["png", "svg"])

# ---- final targets ---------------------------------------------------
def plot_targets():
    figs = [
        "qc/plots/qc_completeness_contamination",
        "qc/plots/qc_assembly_stats",
        "prometheus/plots/prometheus_detection_summary",
        "prometheus/plots/prometheus_architecture",
        "prometheus/plots/prometheus_validation_tree",
        "diversity/plots/subfamily_rarefaction",
        "diversity/plots/architecture_by_taxon",
        "distribution/plots/prevalence_by_taxon",
        "distribution/plots/association_forest",
        "distribution/plots/cooccurrence_heatmap",
        "evolution/plots/genetree_speciestree_tanglegram",
        "evolution/plots/ancestral_gain_loss",
        "evolution/plots/selection_by_domain",
        "ecology/plots/ecology_enrichment",
        "ecology/plots/viral_rloop_association",
        "synthesis/figures/master_tree",
        "synthesis/figures/summary_forest",
    ]
    out = []
    for f in figs:
        for ext in PLOT_FORMATS:
            out.append(f"{RESULTS}/{f}.{ext}")
    return out

def all_targets():
    tables = [
        f"{RESULTS}/qc/genome_metadata.tsv",
        f"{RESULTS}/prometheus/prometheus_catalog.tsv",
        f"{RESULTS}/prometheus/detection_sensitivity.tsv",
        f"{RESULTS}/diversity/subfamilies.tsv",
        f"{RESULTS}/diversity/architecture_typology.tsv",
        f"{RESULTS}/scaffolds/master_table.tsv",
        f"{RESULTS}/distribution/prevalence_by_taxon.tsv",
        f"{RESULTS}/distribution/association_tests.tsv",
        f"{RESULTS}/distribution/cooccurrence.tsv",
        f"{RESULTS}/evolution/ancestral/gains_losses.tsv",
        f"{RESULTS}/evolution/selection/selection_summary.tsv",
        f"{RESULTS}/ecology/ecology_associations.tsv",
        f"{RESULTS}/ecology/viral/viral_summary.tsv",
        f"{RESULTS}/synthesis/results_registry.tsv",
        f"{RESULTS}/provenance/versions.txt",
    ]
    return tables + plot_targets()
