"""Stage the KOfam profiles and thresholds the methanogenesis rules need.

The full KOfam profile set is ~1.5 GB compressed and holds ~28,500 profiles; the
methanogenesis rules reference 85 KOs (plus the three mcr KOs). This script extracts exactly
those profiles and writes a manifest recording the source URL, the archive's MD5 and the date,
so the staged subset is reproducible and auditable.

    python workflow/scripts/py/fetch_kofam_subset.py \\
        --rules workflow/resources/methanogenesis/rules/methanogenesis_rules.yaml \\
        --out-dir resources/db/kofam/methanogenesis_profiles \\
        --ko-list resources/db/kofam/ko_list

Sources (KEGG, https://www.genome.jp/ftp/db/kofam/): profiles.tar.gz and ko_list.gz. KOfam is
released under the KEGG license terms; the files are NOT redistributed with this repository.
With --archive / --ko-list-gz the script reads local copies instead of downloading (useful on
a cluster node without outbound access, and in the tests).

Every KO of the rules file must end up with a profile and a ko_list threshold, otherwise the
script fails: a silently incomplete profile set would turn "not searched" into "absent".
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils
import kofam

import argparse
import gzip
import hashlib
import json
import shutil
import tarfile
import tempfile
import time
import urllib.request

BASE = "https://www.genome.jp/ftp/db/kofam"
PROFILES_URL = f"{BASE}/profiles.tar.gz"
KO_LIST_URL = f"{BASE}/ko_list.gz"


def download(url: str, dest: str, retries: int = 5) -> str:
    """Download with retries; return the file's MD5."""
    h = hashlib.md5()
    for attempt in range(1, retries + 1):
        try:
            h = hashlib.md5()
            with urllib.request.urlopen(url, timeout=600) as r, open(dest, "wb") as out:
                while True:
                    chunk = r.read(1 << 22)
                    if not chunk:
                        break
                    h.update(chunk)
                    out.write(chunk)
            return h.hexdigest()
        except Exception as exc:                       # noqa: BLE001 - retry any transport error
            if attempt == retries:
                raise RuntimeError(f"download failed after {retries} attempts: {url}: {exc}")
            utils.log(f"download attempt {attempt} failed ({exc}); retrying")
            time.sleep(5 * attempt)
    raise AssertionError("unreachable")


def extract_profiles(archive: str, want: set, out_dir: str) -> list:
    os.makedirs(out_dir, exist_ok=True)
    got = []
    with tarfile.open(archive, "r:gz") as tf:
        for m in tf:
            stem = os.path.splitext(os.path.basename(m.name))[0]
            if m.isfile() and m.name.endswith((".hmm", ".HMM")) and stem in want:
                with open(os.path.join(out_dir, f"{stem}.hmm"), "wb") as out:
                    shutil.copyfileobj(tf.extractfile(m), out)
                got.append(stem)
    return sorted(got)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rules", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ko-list", default="", help="write the ko_list here (plain text)")
    ap.add_argument("--extra-kos", default="K00399,K00401,K00402",
                    help="comma list of further KOs to stage (default: the mcr subunits)")
    ap.add_argument("--archive", default="", help="local profiles.tar.gz (skip the download)")
    ap.add_argument("--ko-list-gz", default="", help="local ko_list.gz (skip the download)")
    ap.add_argument("--manifest", default="", help="default: <out-dir>/manifest.json")
    args = ap.parse_args(argv)

    import methanogenesis_search as msearch
    rules = msearch.load_yaml(args.rules)
    want = set(msearch.rules_kos(rules))
    want |= {k.strip() for k in args.extra_kos.split(",") if k.strip()}
    os.makedirs(args.out_dir, exist_ok=True)
    manifest = args.manifest or os.path.join(args.out_dir, "manifest.json")
    tmpd = tempfile.mkdtemp(prefix="kofam_")
    try:
        archive, md5 = args.archive, ""
        if not archive:
            archive = os.path.join(tmpd, "profiles.tar.gz")
            utils.log(f"downloading {PROFILES_URL} (~1.5 GB)")
            md5 = download(PROFILES_URL, archive)
        got = extract_profiles(archive, want, args.out_dir)
        missing = sorted(want - set(got))
        if missing:
            raise ValueError(f"{len(missing)} KO(s) of {args.rules} have no KOfam profile in "
                             f"{archive}: {', '.join(missing)}")
        thresholds = {}
        if args.ko_list:
            gz = args.ko_list_gz
            if not gz:
                gz = os.path.join(tmpd, "ko_list.gz")
                utils.log(f"downloading {KO_LIST_URL}")
                download(KO_LIST_URL, gz)
            os.makedirs(os.path.dirname(os.path.abspath(args.ko_list)) or ".", exist_ok=True)
            with gzip.open(gz, "rt") as fh, open(args.ko_list, "w") as out:
                for line in fh:
                    out.write(line)
            table = kofam.parse_ko_list(args.ko_list)
            no_thr = sorted(k for k in want
                            if table.get(k, {}).get("threshold") is None
                            and not kofam.parse_hmm_models(
                                os.path.join(args.out_dir, f"{k}.hmm"))[0]["ga"])
            if no_thr:
                raise ValueError(f"{len(no_thr)} KO(s) have no curated threshold (ko_list or GA): "
                                 f"{', '.join(no_thr)}")
            thresholds = {k: table[k] for k in want if k in table}
        with open(manifest, "w") as out:
            json.dump({"profiles_url": PROFILES_URL if not args.archive else args.archive,
                       "profiles_md5": md5, "ko_list_url": KO_LIST_URL if not args.ko_list_gz
                       else args.ko_list_gz, "rules": os.path.basename(args.rules),
                       "rules_sha256": hashlib.sha256(open(args.rules, "rb").read()).hexdigest(),
                       "n_profiles": len(got), "kos": got,
                       "thresholds": thresholds, "retrieved": time.strftime("%Y-%m-%d")},
                      out, indent=1)
        utils.log(f"staged {len(got)} KOfam profiles in {args.out_dir}; manifest {manifest}")
        return 0
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
