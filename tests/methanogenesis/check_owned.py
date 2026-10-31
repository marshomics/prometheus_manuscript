#!/usr/bin/env python3
"""Hygiene check for the files leaf-1.1.1 owns.

Every owned Python file must compile, carry a module docstring, and contain no placeholder
marker (TODO / FIXME / XXX / "not implemented" / a bare `pass  # ...` body), and the rules file
and its RULES.md must stay in step: every component, derived term, pathway and substrate flag of
the rules file needs a RULES.md section, every component needs a description, every KO the rules
reference must appear in RULES.md, and every DOI cited there must appear in citations.tsv with
CrossRef-resolved metadata. Prints OWNED FILES CLEAN on success.

    python tests/methanogenesis/check_owned.py
"""
from __future__ import annotations

import ast
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts" / "py"
RULES_DIR = REPO / "workflow" / "resources" / "methanogenesis" / "rules"
OWNED = ["kofam.py", "kegg_module.py", "methanogenesis_search.py",
         "methanogenesis_pathways.py", "fetch_kofam_subset.py", "mcr_search.py"]
PLACEHOLDER = re.compile(r"\b(TODO|FIXME|XXX|HACK)\b|not implemented|NotImplementedError")


def main() -> int:
    problems = []
    for name in OWNED:
        path = SCRIPTS / name
        if not path.exists():
            problems.append(f"{name}: missing")
            continue
        text = path.read_text()
        try:
            tree = ast.parse(text, filename=str(path))
        except SyntaxError as ex:
            problems.append(f"{name}: syntax error: {ex}")
            continue
        if not ast.get_docstring(tree):
            problems.append(f"{name}: no module docstring")
        for m in PLACEHOLDER.finditer(text):
            line = text[:m.start()].count("\n") + 1
            problems.append(f"{name}:{line}: placeholder marker {m.group(0)!r}")
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                body = [b for b in node.body if not isinstance(b, ast.Expr)
                        or not isinstance(b.value, ast.Constant)]
                if len(body) == 1 and isinstance(body[0], ast.Pass):
                    problems.append(f"{name}:{node.lineno}: empty function {node.name}")

    sys.path.insert(0, str(SCRIPTS))
    import methanogenesis_search as msearch                      # noqa: E402

    rules_path = RULES_DIR / "methanogenesis_rules.yaml"
    doc_path = RULES_DIR / "RULES.md"
    if not rules_path.exists():
        problems.append("methanogenesis_rules.yaml: missing")
    elif not doc_path.exists():
        problems.append("RULES.md: missing")
    else:
        rules = msearch.load_yaml(rules_path)
        doc = doc_path.read_text()
        names = (list(rules.get("components") or {}) + list(rules.get("derived") or {}) +
                 list(rules.get("pathways") or {}) + list(rules.get("substrates") or {}))
        for n in names:
            if not re.search(rf"(^|[^A-Za-z0-9_]){re.escape(n)}([^A-Za-z0-9_]|$)", doc, re.M):
                problems.append(f"RULES.md: no entry for rule '{n}'")
        for n, c in (rules.get("components") or {}).items():
            if not str(c.get("description", "")).strip():
                problems.append(f"methanogenesis_rules.yaml: component '{n}' has no description")
        for ko in msearch.rules_kos(rules):
            if ko not in doc:
                problems.append(f"RULES.md: KO {ko} of the rules is not documented")
        # Every DOI cited in RULES.md must appear in citations.tsv, whose rows were resolved
        # against CrossRef (first author / year / title recorded). A citation written from
        # memory and never resolved is the failure mode this check exists to prevent.
        cited = {d.rstrip(".") for d in re.findall(r"doi:([^\s\),;]+)", doc)}
        if not cited:
            problems.append("RULES.md: no citation (expected DOIs)")
        cit_path = RULES_DIR / "citations.tsv"
        if not cit_path.exists():
            problems.append("citations.tsv: missing (run the CrossRef verification)")
        else:
            import csv
            with open(cit_path) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            need = {"doi", "first_author", "year", "title", "verified", "retrieved"}
            if not rows or not need <= set(rows[0]):
                problems.append(f"citations.tsv: needs columns {sorted(need)}")
            else:
                known = {r["doi"].strip(): r for r in rows}
                for d in sorted(cited - set(known)):
                    problems.append(f"RULES.md: DOI {d} is not in citations.tsv (unverified)")
                for d, r in sorted(known.items()):
                    if r["verified"].strip() != "crossref":
                        problems.append(f"citations.tsv: {d} is not marked verified=crossref")
                    if not r["first_author"].strip() or not str(r["year"]).strip():
                        problems.append(f"citations.tsv: {d} has no resolved author/year")
                for d in sorted(set(known) - cited):
                    problems.append(f"citations.tsv: {d} is not cited in RULES.md")

    if problems:
        for p in problems:
            print(f"PROBLEM {p}")
        print(f"{len(problems)} problem(s)")
        return 1
    print(f"checked {len(OWNED)} owned script(s) and the rules file against RULES.md")
    print("OWNED FILES CLEAN")
    return 0


if __name__ == "__main__":
    sys.exit(main())
