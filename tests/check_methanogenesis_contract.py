#!/usr/bin/env python3
"""Interface check for the methanogenesis-type integration (node-1.1:N2).

Verifies that the pieces agree on one contract, without running the pipeline:

1. the rules file loads, every KO it references has a KOfam-style profile name documented, and
   its component expressions still equal the KEGG module steps they cite;
2. the pathway table's column list (methanogenesis_pathways.out_columns) matches SCHEMA.md;
3. the master-table additive columns (build_master_table.MGEN_ADDITIVE) match SCHEMA.md and are
   additive (the existing backbone and additive columns are unchanged);
4. every config key the new rules read exists in config/config.yaml and validates against
   workflow/schemas/config.schema.yaml;
5. the statistics script's required master columns are the ones the master table provides;
6. README documents the new config keys.

Prints METHANOGENESIS CONTRACT OK on success.
"""
from __future__ import annotations

import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "workflow" / "scripts" / "py"
sys.path.insert(0, str(SCRIPTS))


def _python_with_jsonschema():
    """This interpreter if it has jsonschema, else the one running Snakemake (which requires
    it), else any `python` on PATH that has it. Returns None when there is none."""
    import importlib.util
    import os
    import shutil
    if importlib.util.find_spec("jsonschema"):
        return sys.executable
    cands = []
    smk = os.environ.get("PROM_SNAKEMAKE") or shutil.which("snakemake")
    if smk:
        cands.append(str(pathlib.Path(smk).resolve().parent / "python"))
    cands += [shutil.which("python3") or "", shutil.which("python") or ""]
    import subprocess
    for c in cands:
        if c and pathlib.Path(c).exists():
            r = subprocess.run([c, "-c", "import jsonschema"], capture_output=True)
            if r.returncode == 0:
                return c
    return None


def validate_config_against_schema(repo) -> str:
    """'' when config.yaml validates against the schema, else the problem.

    jsonschema is a hard requirement of Snakemake (which validates the config at run time),
    so a missing module here is a broken check, not a reason to skip the validation.
    """
    import subprocess
    py = _python_with_jsonschema()
    if py is None:
        return ("no interpreter with jsonschema found (set PROM_SNAKEMAKE or put snakemake on "
                "PATH): config.yaml was NOT validated against the schema")
    code = ("import sys, yaml, jsonschema\n"
            "cfg = yaml.safe_load(open(sys.argv[1])); sch = yaml.safe_load(open(sys.argv[2]))\n"
            "jsonschema.validate(cfg, sch); print('CONFIG VALID')\n")
    r = subprocess.run([py, "-c", code, str(repo / "config" / "config.yaml"),
                        str(repo / "workflow" / "schemas" / "config.schema.yaml")],
                       capture_output=True, text=True)
    if r.returncode != 0 or "CONFIG VALID" not in r.stdout:
        return f"config.yaml does not validate: {(r.stderr or r.stdout).strip()[-300:]}"
    return ""


def main() -> int:
    problems = []
    import kegg_module as km
    import methanogenesis_pathways as mpath
    import methanogenesis_search as msearch
    import build_master_table as bmt

    rules_path = REPO / "workflow" / "resources" / "methanogenesis" / "rules" / \
        "methanogenesis_rules.yaml"
    rules = msearch.load_yaml(rules_path)
    try:
        R = mpath.Rules(rules)
    except Exception as ex:                                   # noqa: BLE001
        print(f"PROBLEM rules file does not load: {ex}")
        print("1 problem(s)")
        return 1

    # 1. components still match the KEGG steps they cite
    for name, c in (rules.get("components") or {}).items():
        if "module" not in c or "step" not in c:
            continue
        step = km.steps(km.parse(rules["kegg_modules"][c["module"]]))[int(c["step"]) - 1]
        want = repr(km.parse(c["expr"]))
        alts = step.items if isinstance(step, km.Alt) else [step]
        ok = want in [repr(a) for a in alts] if c.get("part") else want == repr(step)
        if not ok:
            problems.append(f"component {name} no longer matches {c['module']} step {c['step']}")

    schema = (REPO / "workflow" / "scripts" / "SCHEMA.md").read_text()

    # 2. pathway table columns
    cols = mpath.out_columns()
    sec = schema.split("## scaffolds/metabolism/methanogenesis_pathways.tsv")
    if len(sec) != 2:
        problems.append("SCHEMA.md has no methanogenesis_pathways.tsv section")
    else:
        body = sec[1].split("\n## ")[0]
        for c in cols:
            if c not in body:
                problems.append(f"SCHEMA.md methanogenesis_pathways: column '{c}' undocumented")
    if "## scaffolds/metabolism/methanogenesis_ko_hits.tsv" not in schema:
        problems.append("SCHEMA.md has no methanogenesis_ko_hits.tsv section")

    # 3. master-table additive columns, and the interface stays additive
    master_sec = schema.split("## scaffolds/master_table.tsv")[1].split("\n## ")[0]
    for c in bmt.MGEN_ADDITIVE:
        if c not in master_sec:
            problems.append(f"SCHEMA.md master_table: column '{c}' undocumented")
    base = bmt.build_master_header(["group_X", "niche_Y"], with_pathways=False)
    full = bmt.build_master_header(["group_X", "niche_Y"], with_pathways=True)
    if full[:len(base)] != base:
        problems.append("master-table header is not additive: the pathway columns reorder it")
    if full[len(base):] != list(bmt.MGEN_ADDITIVE):
        problems.append("master-table pathway columns are not MGEN_ADDITIVE")

    # 4. config keys
    import yaml
    cfg = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    mg = cfg.get("methanogenesis")
    if not isinstance(mg, dict):
        problems.append("config.yaml has no methanogenesis block")
    else:
        for k in ("enabled", "rules", "readthrough", "absence_alpha", "stats"):
            if k not in mg:
                problems.append(f"config.yaml methanogenesis.{k} missing")
        for k in ("min_n", "min_class", "min_changes", "nsim", "corhmm", "epv_min", "vif_max",
                  "complete_threshold"):
            if k not in (mg.get("stats") or {}):
                problems.append(f"config.yaml methanogenesis.stats.{k} missing")
        if mg.get("enabled"):
            problems.append("config.yaml ships with methanogenesis.enabled true (must be opt-in)")
        if str(mg.get("rules", "")) != str(rules_path.relative_to(REPO)):
            problems.append("config.yaml methanogenesis.rules does not point at the shipped rules")
    if "kofam_profiles" not in (cfg.get("db") or {}):
        problems.append("config.yaml db.kofam_profiles missing")
    err = validate_config_against_schema(REPO)
    if err:
        problems.append(err)

    # 5. the statistics script reads the columns the master table writes
    stats_r = REPO / "workflow" / "scripts" / "R" / "methanogenesis_stats.R"
    if stats_r.exists():
        text = stats_r.read_text()
        for c in [f"mgen_{p}" for p in bmt.MGEN_PATHWAYS]:
            if c not in text:
                problems.append(f"methanogenesis_stats.R does not mention predictor '{c}'")
        for c in ("methanogen", "tree_tip"):
            if c not in text:
                problems.append(f"methanogenesis_stats.R does not read master column '{c}'")
        # the response must come from the Contract-item-6 tier rule, by either route
        if not any(c in text for c in ("prometheus_status_primary", "prometheus_tier")):
            problems.append("methanogenesis_stats.R does not use the Contract item 6 response "
                            "(prometheus_status_primary or prometheus_tier)")
    else:
        problems.append("methanogenesis_stats.R missing (leaf-1.1.3)")

    # 6. README documents the module and its keys
    readme = (REPO / "README.md").read_text()
    for k in ("methanogenesis.enabled", "db.kofam_profiles", "methanogenesis.absence_alpha",
              "methanogenesis.readthrough", "methanogenesis.rules"):
        if k.split(".")[-1] not in readme or "methanogenesis" not in readme:
            problems.append(f"README does not document {k}")

    if problems:
        for p in problems:
            print(f"PROBLEM {p}")
        print(f"{len(problems)} problem(s)")
        return 1
    print(f"rules={len(R.components)} components, {len(R.all_kos())} KOs; "
          f"pathway table {len(cols)} columns; master +{len(bmt.MGEN_ADDITIVE)} columns")
    print("METHANOGENESIS CONTRACT OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
