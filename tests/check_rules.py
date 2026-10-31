"""Static hygiene check for Snakemake rule files (no Snakemake import needed).

Fails (exit 1) when, in any given .smk file:
  * a shell line swallows a failure with `|| true` / `|| :` and does not carry an
    explicit `# nonfatal: <reason>` justification on the same line;
  * the `$(grep -c ... || echo 0)` idiom is used (grep -c already prints 0 and exits 1,
    so `|| echo 0` yields "0\\n0" and breaks integer tests);
  * a rule (other than pure target rules with only `input:`) has no `log:` directive.
Prints RULES HYGIENE OK on success.
"""
import re
import sys

bad = []
for path in sys.argv[1:]:
    lines = open(path).read().splitlines()
    for i, ln in enumerate(lines, 1):
        code = ln.split("#", 1)[0] if not ln.lstrip().startswith("#") else ""
        if re.search(r"\|\|\s*(true|:)\b", code) and "# nonfatal:" not in ln:
            bad.append(f"{path}:{i}: swallowed failure without '# nonfatal:' justification: {ln.strip()}")
        if re.search(r"\|\|\s*echo\s+0", code):
            bad.append(f"{path}:{i}: '|| echo 0' count idiom (yields '0\\n0'): {ln.strip()}")
    # rule blocks
    text = "\n".join(lines)
    for m in re.finditer(r"^rule\s+(\w+)\s*:\s*\n(.*?)(?=^rule\s|\Z)", text, flags=re.S | re.M):
        name, body = m.group(1), m.group(2)
        keys = set(re.findall(r"^\s{4}(\w+)\s*:", body, flags=re.M))
        target_only = keys <= {"input"}
        if not target_only and "log" not in keys:
            bad.append(f"{path}: rule {name} has no log: directive")
if bad:
    print("\n".join(bad))
    print(f"RULES HYGIENE FAILED ({len(bad)} finding(s))")
    sys.exit(1)
print("RULES HYGIENE OK")
