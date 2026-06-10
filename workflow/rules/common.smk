# =====================================================================
# common.smk — config, sample sheet, helpers, and the global target list.
# Kept deliberately static (no checkpoints) so `snakemake -n` validates the
# full DAG; data-dependent steps degrade gracefully on empty inputs.
# =====================================================================
import os
import json
import math
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

# ---- raw input genome set + batching (Phase 1: preprocess) -----------
INPUT_METADATA = cfg("input/metadata", "input_genomes_metadata_with_niches.tsv")
INPUT_GENOMES_DIR = cfg("input/genomes_dir", "")
INPUT_GENOME_EXT = cfg("input/genome_extension", "fasta")
ID_COL = cfg("input/id_column", "genome_id")
FILE_COL = cfg("input/filename_column", "filename")
SRC_COL = cfg("input/source_column", "source")
BATCH_SIZE = int(cfg("preprocess/gtdbtk_batch_size", 20000))

def _count_rows(path):
    try:
        with open(path) as fh:
            return max(0, sum(1 for _ in fh) - 1)  # minus header
    except Exception:
        return 0

N_INPUT = _count_rows(INPUT_METADATA)
N_BATCHES = math.ceil(N_INPUT / BATCH_SIZE) if N_INPUT > 0 else 0
BATCHES = [f"{i:04d}" for i in range(N_BATCHES)]

# ---- analysis samples (species representatives; produced by Phase 1) --
def _load_samples(path):
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
        if df.empty or "genome_id" not in df.columns:
            return [], None
        return df["genome_id"].tolist(), df.set_index("genome_id", drop=False)
    except Exception:
        return [], None

GENOMES, samples = _load_samples(config["samples"])
if not GENOMES:
    print(f"[common] No analysis samples at {config['samples']} yet — run "
          f"`snakemake preprocess` first (Phase 1). Phase-2 targets stay empty until then.")
# Fail fast on the one ID assumption the whole pipeline rests on: combined protein
# ids are '{genome}__{orig}' and utils.genome_of splits on the first '__', so a
# genome_id must never itself contain '__' (would mis-recover the genome downstream).
_bad_ids = [g for g in GENOMES if "__" in str(g)]
if _bad_ids:
    raise ValueError(
        f"[common] {len(_bad_ids)} genome_id(s) contain '__' (e.g. {_bad_ids[:3]}); "
        "this breaks the '{genome}__{orig}' id scheme. Rename them in the input metadata.")

def genome_fasta(wildcards):
    """Prodigal & co. read the staged representative assembly. This is exactly the
    path recorded in the filtered table's `assembly_path` column (genomes_dir ==
    GENOME_DIR == results/preprocess/genomes), staged by filter_collapse."""
    return os.path.join(GENOME_DIR, f"{wildcards.genome}.{EXT}")

# ---- Bacterial cross-domain branch (opt-in; never alters the archaeal DAG) ------
BACT_ENABLED = bool(cfg("bacteria/enabled", False))

def _load_bacterial_genomes():
    """Bacterial genome_ids, from bacteria.qc (preferred) or bacteria.taxonomy."""
    for key in ("bacteria/qc", "bacteria/taxonomy"):
        p = cfg(key, "")
        if p and os.path.exists(p):
            try:
                df = pd.read_csv(p, sep="\t", dtype=str, comment="#")
                col = "genome_id" if "genome_id" in df.columns else (
                    "user_genome" if "user_genome" in df.columns else df.columns[0])
                ids = [str(g).strip() for g in df[col].dropna().tolist() if str(g).strip()]
                bad = [g for g in ids if "__" in g]
                if bad:
                    raise ValueError(f"[common] {len(bad)} bacterial genome_id(s) contain "
                                     f"'__' (e.g. {bad[:3]}); breaks the id scheme.")
                return sorted(dict.fromkeys(ids))
            except ValueError:
                raise
            except Exception as e:
                print(f"[common] could not read bacterial genomes from {p}: {e}")
    return []

BACT_GENOMES = _load_bacterial_genomes() if BACT_ENABLED else []
if BACT_ENABLED and not BACT_GENOMES:
    print("[common] bacteria.enabled but no bacterial genome list found yet "
          "(set bacteria.qc or bacteria.taxonomy). crossdomain targets stay empty.")

# ---- Bacterial provirus / R-loop module (activates when assemblies are provided) ----
BACT_VIRAL = BACT_ENABLED and bool(cfg("bacteria/genomes_dir", ""))

def bact_genome_fasta(wildcards):
    gd = cfg("bacteria/genomes_dir", "")
    ext = cfg("bacteria/genome_extension", "fasta")
    return os.path.join(gd, f"{wildcards.genome}.{ext}")

# ---- Intrinsic-disorder module (opt-in; consumes archaeal outputs read-only) ----
DISORDER_ENABLED = bool(cfg("disorder/enabled", False))

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
PER_DOMAIN_TREES = cfg("evolution/per_domain_trees", True)
DOMAINS = list(cfg("prometheus/domains", {"DUF4011": {}, "DUF3320": {}, "DNA2": {}, "nuclease": {}, "REase": {}}).keys())
RUN_GENERAX = cfg("evolution/reconciliation/run_generax", True)
RUN_ALE = cfg("evolution/reconciliation/run_ale", False)
RUN_SELECTION = cfg("evolution/selection/run", True)
SELECTION_METHODS = cfg("evolution/selection/methods", ["FEL", "MEME", "BUSTED", "aBSREL", "RELAX"])
PLOT_FORMATS = cfg("plotting/formats", ["png", "svg"])

# ---- Phase 1 (preprocess) targets ------------------------------------
def preprocess_targets():
    return [
        f"{RESULTS}/preprocess/annotated_metadata.tsv",
        f"{RESULTS}/preprocess/species_metadata.tsv",
        f"{RESULTS}/preprocess/representatives.txt",
        f"{RESULTS}/preprocess/genome_metadata.tsv",
        f"{RESULTS}/preprocess/taxonomy.tsv",
        f"{RESULTS}/preprocess/genomes/.staged",
        f"{RESULTS}/preprocess/plots/preprocess_quality.png",
    ]

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
        # --- extra figures (plots_extra.smk; added alongside the above) ---
        "qc/plots/prometheus_covariates",
        "distribution/plots/prevalence_caterpillar",
        "distribution/plots/cooccurrence_oddsratio",
        "diversity/plots/architecture_combinations",
        "diversity/plots/architecture_by_phylum",
        "diversity/plots/subfamily_identity_curve",
        "ecology/plots/viral_pressure_prometheus",
        "ecology/plots/niche_effect_forest",
        "ecology/plots/prevalence_niche_phylum",
        "evolution/plots/dnds_along_protein",
        "evolution/plots/reconciliation_events",
        "evolution/plots/gain_enrichment_null",
    ]
    out = []
    for f in figs:
        for ext in PLOT_FORMATS:
            out.append(f"{RESULTS}/{f}.{ext}")
    return out

def crossdomain_targets():
    """Opt-in bacterial / cross-domain outputs. Empty unless bacteria.enabled."""
    if not BACT_ENABLED:
        return []
    t = [
        f"{RESULTS}/bacteria/prometheus/bacterial_prometheus_catalog.tsv",
        f"{RESULTS}/bacteria/distribution/bacterial_prevalence_by_taxon.tsv",
        f"{RESULTS}/bacteria/diversity/bacterial_architecture_typology.tsv",
        f"{RESULTS}/crossdomain/crossdomain_summary.tsv",
        f"{RESULTS}/crossdomain/tree/crossdomain.treefile",
        f"{RESULTS}/bacteria/distribution/bacterial_cooccurrence.tsv",
        f"{RESULTS}/crossdomain/selection/crossdomain_selection_summary.tsv",
    ]
    for ext in PLOT_FORMATS:
        t.append(f"{RESULTS}/crossdomain/figures/crossdomain_tree.{ext}")
        t.append(f"{RESULTS}/crossdomain/figures/crossdomain_comparison.{ext}")
    # optional niche specificity control (activates only when niche labels are provided)
    if cfg("bacteria/niche_metadata", ""):
        t.append(f"{RESULTS}/bacteria/ecology/bacterial_niche_associations.tsv")
    # provirus / R-loop mechanism test (activates only when bacterial assemblies are provided)
    if BACT_VIRAL:
        t.append(f"{RESULTS}/bacteria/ecology/viral/bacterial_viral_summary.tsv")
        t.append(f"{RESULTS}/bacteria/ecology/bacterial_viral_associations.tsv")
        for ext in PLOT_FORMATS:
            t.append(f"{RESULTS}/crossdomain/figures/viral_rloop_crossdomain.{ext}")
    return t


def disorder_targets():
    """Opt-in disorder outputs. Empty unless disorder.enabled."""
    if not DISORDER_ENABLED:
        return []
    t = [
        f"{RESULTS}/disorder/disorder_per_locus.tsv",
        f"{RESULTS}/disorder/disorder_by_column.tsv",
        f"{RESULTS}/disorder/disorder_associations.tsv",
    ]
    for ext in PLOT_FORMATS:
        for fig in ("disorder_along_protein", "disorder_by_group", "disorder_linker_vs_domain"):
            t.append(f"{RESULTS}/disorder/plots/{fig}.{ext}")
    if BACT_ENABLED:
        t.append(f"{RESULTS}/disorder/bacteria/disorder_per_locus.tsv")
    return t


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
        f"{RESULTS}/evolution/reconciliation/reconciliation_events.tsv",
        f"{RESULTS}/evolution/ancestral/gains_losses.tsv",
        f"{RESULTS}/evolution/selection/selection_summary.tsv",
        f"{RESULTS}/ecology/ecology_associations.tsv",
        f"{RESULTS}/ecology/viral/viral_summary.tsv",
        f"{RESULTS}/synthesis/results_registry.tsv",
        f"{RESULTS}/provenance/versions.txt",
    ]
    # disorder_targets() is empty unless disorder.enabled, so the default run is
    # unaffected until you turn the module on.
    return tables + plot_targets() + disorder_targets()
