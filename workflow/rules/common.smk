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

# The shipped SGE workflow profile uses the Snakemake >= 8 executor-plugin interface.
min_version("8.0")

# ---- config ----------------------------------------------------------
configfile: "config/config.yaml"
# Fail fast on an invalid configuration (a schema violation is an error, not a warning).
validate(config, "../schemas/config.schema.yaml")

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

# The ONE genome-id rule (identical to batch_genomes.py ID_RE): a leading alphanumeric,
# then letters, digits, '_', '.', '-'. Such ids need no Newick quoting, are accepted by
# GTDB-Tk/IQ-TREE/MAFFT unchanged and equal utils.sanitize(id); '__' is reserved as the
# genome__protein separator.
import re as _re
_ID_RE = _re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")


def _check_genome_ids(ids, where):
    """Fail at DAG build on ids that break the pipeline (see _ID_RE), so a bad id never
    surfaces hours later as a run-time failure of batch_genomes.py or a tree tool."""
    bad_sep = [g for g in ids if "__" in str(g)]
    bad_chr = [g for g in ids if not _ID_RE.match(str(g))]
    if bad_sep or bad_chr:
        raise ValueError(
            f"[common] invalid genome_id(s) in {where}: "
            + (f"{len(bad_sep)} contain '__' (e.g. {bad_sep[:3]}); " if bad_sep else "")
            + (f"{len(bad_chr)} do not match [A-Za-z0-9][A-Za-z0-9_.-]* (e.g. {bad_chr[:3]}); "
               if bad_chr else "")
            + "rename them in the input metadata.")


def _check_no_inline_hash(path):
    """The input table is read with comment='#', which silently truncates a data row at
    an inline '#'. Only whole-line comments are allowed (same rule as batch_genomes.py)."""
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            s = line.strip()
            if s and not s.startswith("#") and "#" in line:
                raise ValueError(f"[common] {path}:{n}: '#' inside a data row would silently "
                                 f"truncate it; only whole-line '#' comments are allowed.")


def _count_input_genomes(path):
    """Number of genomes in the Phase-1 input table: data rows with a non-empty id
    (header, '#' comment lines and blank lines excluded). batch_genomes.py receives the
    resulting batch count from the rule, so the expected `.ready` files always match."""
    if not os.path.exists(path):
        return 0
    _check_no_inline_hash(path)
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#", skip_blank_lines=True)
    if ID_COL not in df.columns:
        raise ValueError(f"[common] input metadata {path} has no '{ID_COL}' column "
                         f"(config input.id_column).")
    ids = df[ID_COL].fillna("").astype(str).str.strip()
    ids = ids[ids != ""]
    dup = ids[ids.duplicated()].unique().tolist()
    if dup:
        raise ValueError(f"[common] {len(dup)} duplicated genome id(s) in {path} "
                         f"(e.g. {dup[:3]}); ids must be unique.")
    _check_genome_ids(ids.tolist(), path)
    return int(len(ids))

N_INPUT = _count_input_genomes(INPUT_METADATA)
N_BATCHES = math.ceil(N_INPUT / BATCH_SIZE) if N_INPUT > 0 else 0
BATCHES = [f"{i:04d}" for i in range(N_BATCHES)]

# ---- analysis samples (species representatives; produced by Phase 1) --
def _load_samples(path):
    """Species-representative ids from the Phase-1 table. A MISSING file means Phase 1
    has not run yet; a present-but-malformed table is an error (never silently empty)."""
    if not os.path.exists(path):
        return [], None
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if "genome_id" not in df.columns:
        raise ValueError(f"[common] samples table {path} has no 'genome_id' column.")
    ids = df["genome_id"].fillna("").astype(str).str.strip()
    if (ids == "").any():
        raise ValueError(f"[common] samples table {path} has {(ids == '').sum()} empty genome_id(s).")
    if ids.duplicated().any():
        raise ValueError(f"[common] samples table {path} has duplicated genome_id(s): "
                         f"{ids[ids.duplicated()].unique().tolist()[:5]}")
    df["genome_id"] = ids
    return ids.tolist(), df.set_index("genome_id", drop=False)

GENOMES, samples = _load_samples(config["samples"])
PHASE2_READY = bool(GENOMES)
if not PHASE2_READY:
    print(f"[common] No analysis samples at {config['samples']} yet: `snakemake` resolves to "
          f"the Phase-1 (preprocess) targets only. Re-run it after Phase 1 to build Phase 2.")
# Fail fast on the id assumptions the whole pipeline rests on (see _check_genome_ids):
# combined protein ids are '{genome}__{gene}' (utils.genome_of splits on the first '__'),
# GTDB-Tk rejects some characters, and tree tips must not need Newick quoting.
_check_genome_ids(GENOMES, config["samples"])

def genome_fasta(wildcards):
    """Prodigal & co. read the staged representative assembly. This is exactly the
    path recorded in the filtered table's `assembly_path` column (genomes_dir ==
    GENOME_DIR == results/preprocess/genomes), staged by filter_collapse."""
    return os.path.join(GENOME_DIR, f"{wildcards.genome}.{EXT}")

# ---- Methanogenesis pathway capacity (opt-in; design 3.3 / 6.2) -----------------
# The master table only carries the mgen_* columns, and methanogenesis_stats.R can only run,
# when this is on; every consumer must gate on it so the default DAG is unchanged.
MGEN_ENABLED = str(cfg("methanogenesis/enabled", False)).strip().lower() in ("true", "1", "yes")

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
                _check_genome_ids(ids, p)
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
DOMAINS = list(cfg("prometheus/domains", {"DUF4011": {}, "DNA2": {}, "REase": {}, "DUF3320": {}}).keys())
RUN_GENERAX = cfg("evolution/reconciliation/run_generax", True)
RUN_ALE = cfg("evolution/reconciliation/run_ale", False)
RUN_SELECTION = cfg("evolution/selection/run", True)
SELECTION_METHODS = cfg("evolution/selection/methods", ["FEL", "MEME", "BUSTED", "aBSREL", "RELAX"])
PLOT_FORMATS = cfg("plotting/formats", ["png", "svg"])
if isinstance(PLOT_FORMATS, str):
    PLOT_FORMATS = PLOT_FORMATS.split(",")
PLOT_FORMATS = [str(f).strip().lower().lstrip(".") for f in PLOT_FORMATS if str(f).strip()]
_PLOT_SUPPORTED = ("png", "svg", "pdf")  # what plotting_style.py and plotting_theme.R can write
if not PLOT_FORMATS or [f for f in PLOT_FORMATS if f not in _PLOT_SUPPORTED]:
    raise ValueError(f"config plotting.formats must be a non-empty subset of {list(_PLOT_SUPPORTED)}; "
                     f"got {cfg('plotting/formats', None)!r}")
# The figure helpers (plotting_style.py, plotting_theme.R) read plotting.* from these
# variables, so every figure rule writes exactly the formats it declares. Cluster jobs
# re-parse this file, so they receive the same values; the config is the source of truth.
os.environ["PROM_PLOT_FORMATS"] = ",".join(PLOT_FORMATS)
os.environ["PROM_DPI"] = str(int(cfg("plotting/dpi", 600)))
os.environ["PROM_FONT"] = str(cfg("plotting/font_family", "Arial"))
os.environ["PROM_FALLBACK_FONT"] = str(cfg("plotting/fallback_font", "DejaVu Sans"))
os.environ["PROM_FONT_SIZE"] = str(float(cfg("plotting/base_font_size", 7)))

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
    """Opt-in bacterial / cross-domain outputs. Empty unless bacteria.enabled.
    Needs the archaeal Phase-2 outputs, so before Phase 1 it resolves to Phase 1."""
    if not BACT_ENABLED:
        return []
    if not PHASE2_READY:
        return preprocess_targets()
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
    if not PHASE2_READY:
        return preprocess_targets()
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
    # Before Phase 1 has produced the species table there is no analysis genome set:
    # resolve `all` to the Phase-1 targets instead of building an EMPTY analysis
    # (which would later have to be invalidated and re-run).
    if not PHASE2_READY:
        return preprocess_targets()
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
    # Per-domain trees (domain shuffling test) were never requested by `all` before the
    # review (EVO-32): request them whenever evolution.per_domain_trees is on.
    if PER_DOMAIN_TREES:
        tables.append(f"{RESULTS}/evolution/domains/rf_distances.tsv")
    if MGEN_ENABLED:          # opt-in methanogenesis pathway capacity (design 3.3 / 6.2)
        tables += [f"{RESULTS}/scaffolds/metabolism/methanogenesis_pathways.tsv",
                   f"{RESULTS}/ecology/methanogenesis_associations.tsv"]
    # disorder_targets() is empty unless disorder.enabled, so the default run is
    # unaffected until you turn the module on.
    return tables + plot_targets() + disorder_targets()
