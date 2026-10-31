# Review — leaf 1.1.9, stage `disorder` (opt-in intrinsic-disorder module)

Scope (owned): `workflow/rules/disorder.smk`, `workflow/envs/disorder.yaml`,
`workflow/scripts/py/{predict_disorder,map_disorder_to_alignment,plot_disorder}.py`,
`workflow/scripts/R/disorder_stats.R`, `tests/disorder/**`.
Line numbers in "file:line" refer to the pristine copy (commit 56f5437).
Suite: `bash tests/run_leaf.sh disorder` → 28/28, ~38 s, ≤ 2 threads per tool. Real-tool
runs: metapredict 3.1.1 (V3 + pLDDT V1), segmasker (BLAST+ 2.17.0), MAFFT 7.526 on the real
ProA fixture + its real hmmscan domtbl, phylolm 2.6.5 / ape 5.8.1 (R 4.3.3), Snakemake 8.30
dry run of the `disorder` target. IUPred3 is licence-gated and was exercised through a
stub that mimics its library call (see U1).

Primary-source evidence collected for this review (all re-runnable):

* **PyPI JSON** (`pypi.org/pypi/<pkg>/<ver>/json`): metapredict releases 2.63, 2.65.0,
  2.65.1, 3.0.0–3.0.2, 3.1.1 (no 2.6.5). Requires: 2.63 → `alphaPredict==1.0`; 2.65.x →
  `alphaPredict>=1.3.0`; 3.0.x → `pytorch-lightning`; 3.1.1 → `torch>=2.3, numpy>=2.0,
  scipy>=1.13, protfasta>=0.1.25, getSequence>=2.2.1`. Linux ships as sdist only.
* **metapredict 3.1.1 sdist**: `setup.py` ("The extension is required, not optional: if it
  cannot be compiled the install fails"); `meta.py` `predict_disorder(inputs, version='V3',
  device=None, normalized=True, round_values=True, return_numpy=True, …)` — str → ndarray,
  dict → `{key: [seq, ndarray]}`; `predict_pLDDT(…, pLDDT_version='V2')`, V1 = alphaPredict
  network, 0–100 scale; `parameters.py` `METAPREDICT_V3_THRESHOLD=0.5`, `DEFAULT_NETWORK='V3'`;
  `backend/encode_sequence.py` raises `ValueError('Invalid amino acid detected')` for any
  non-canonical letter; `backend/network_parameters.py` — V2 trained on legacy metapredict +
  alphaPredict pLDDT, V3 "similar to v2 … real pLDDT scores based on AF2 V4 structures";
  changelog 3.0 "Removed alphaPredict as a dependency", 3.1.0 device-dependent float
  differences ("set device='cpu' if you need bit-for-bit reproducible scores").
* **Byte identity**: `metapredict-3.1.1/…/ppLDDT/networks/alpha_fold_networkV7_hs100_nL2_200e_all_prot.pt`
  and `alphapredict-1.2/…/networks/` same file, sha256 `354fb0005abfd70335e6db288b82e82ffaa8d9294650ab389197e5c08ad98ced`.
  Real run on p53: `alphaPredict.predict()` vs `metapredict.predict_pLDDT(pLDDT_version='V1')`
  max |Δ| = 1.03e-4 (alphaPredict rounds to 4 dp).
* **Real metapredict 3.1.1 (torch 2.15.0, numpy 2.1.3), p53 P04637**: fraction ≥ 0.5 =
  1.00 (TAD 1–60), 0.053 (DBD 101–290), 1.00 (CTD 361–393); batch-dict output equals the
  single-sequence call exactly (max |Δ| 0.0). pLDDT V1 < 50 marks only 31.7 % of the TAD and
  every pLDDT<50 residue is also metapredict-disordered (consensus = the pLDDT call).
* **protfasta 0.1.26** `_configs.STANDARD_CONVERSION`: B→N, U→C, X→G, Z→Q, `*`/`-`/space deleted.
* **segmasker (BLAST+ 2.17.0 osx-arm64; same interval format as 2.15)**: `-outfmt interval`
  prints `>` + the FASTA defline verbatim (no `lcl|` without `-parse_seqids`), then
  `start - end` **0-based inclusive** (checked against `-outfmt fasta` lower-casing: `30 - 45`
  = residues 31–46); a record without LCR still gets its header; a missing input exits 1;
  SEG defaults window 12 / locut 2.2 / hicut 2.5.
* **MAFFT 7.526**: drops a trailing `*`, upper-cases, keeps `X`.
* **micromamba linux-64 dry-run**: blast 2.15.0 (bioconda) requires `libzlib <1.3`; conda-forge
  pytorch 2.5.1 requires `libzlib >=1.3.1` → unsatisfiable; blast 2.16.0 / 2.17.0 solve.
* **phylolm 2.6.5**: `summary(fit)$coefficients` columns `Estimate, StdErr, t.value, p.value`;
  `fit$optpar` = λ; p = 2·pt(−|t|, n − d) (verified in the suite); no rootedness check (an
  unrooted tree is silently treated as rooted at its basal node); λ ≈ 0 fit logLik 39.3443 =
  closed-form λ = 0 GLS log-likelihood 39.34428.
* **IUPred** (Erdős & Dosztányi 2020 Curr Protoc Bioinform; Erdős et al. 2021 NAR 49:W297):
  library `iupred(seq, mode)` with `mode ∈ {long, short, glob}`; scores in [0, 1], > 0.5 =
  disordered; IUPred3 adds medium (Savitzky–Golay 19,5) / strong (MA 11) smoothing for the
  web server and the downloadable package; a re-packaging of the library returns
  `(scores, glob_text)` for `mode='long'`.

---

## Findings

### D-01 · critical · workflow/envs/disorder.yaml:15
**Issue.** `metapredict==2.6.5` does not exist on PyPI, so the env cannot be built and the
module can never run. **Evidence.** PyPI release list above. **Fix.** `metapredict==3.1.1`
(installed and API-tested here; V3 default network, deterministic batching, NumPy 2).
**Test.** `test_env_yaml_pins` (pin + installed version == 3.1.1), `test_real_metapredict_end_to_end`;
G3 pip check (`pip ok : metapredict==3.1.1`).

### D-02 · critical · workflow/envs/disorder.yaml:20
**Issue.** `alphaPredict==1.2` is unsatisfiable next to every metapredict 2.6x (2.63 needs
`==1.0`, 2.65.x needs `>=1.3.0`), so even a corrected 2.6x pin would not resolve.
**Evidence.** PyPI `requires_dist`. **Fix.** Separate package dropped; the fallback uses
metapredict's bundled pLDDT V1 network, byte-identical to alphaPredict's (sha256 above,
numerically equal to alphaPredict 1.2 within its rounding). The label `metapredict+alphapredict`
stays exact. **Test.** `test_env_yaml_pins`, `test_real_metapredict_end_to_end` (pLDDT V1 == direct call).

### D-03 · major · workflow/envs/disorder.yaml:12
**Issue.** `blast=2.15.0` (libzlib <1.3) makes the env unsolvable together with any current
conda-forge PyTorch. **Evidence.** micromamba conflict trace. **Fix.** `blast=2.17.0` (segmasker
interval format verified on 2.17.0). **Test.** micromamba linux-64 dry-run with a writable
cache: rc 0, 161 packages (python 3.11.18, pytorch 2.5.1 cpu_mkl, blast 2.17.0, c-compiler
2.0.0 → gcc 15.3); `test_real_segmasker_positions`.

### D-04 · major · workflow/envs/disorder.yaml (whole file)
**Issue.** metapredict has no linux wheel and its sdist must compile a required Cython
extension (no compiler in the env → install fails); pip would also pull the multi-GB CUDA
torch wheel; NumPy 1.26 pin is incompatible with metapredict 3.1. **Fix.** `c-compiler`,
`pytorch=2.5.1=cpu_*` (conda-forge CPU build satisfies `torch>=2.3`), `numpy=2.1.3`,
`scipy=1.14.1`, `matplotlib-base`, pinned `protfasta==0.1.26`, `getSequence==3.0.1`;
biopython (unused) removed. **Test.** env dry-run solve (U2), `test_env_yaml_pins`.

### D-05 · critical · predict_disorder.py:118-131, 280-286, 312-313
**Issue.** A missing metapredict produced header-only outputs (the rule "succeeded"), and any
per-sequence exception dropped that locus silently. metapredict raises on X/B/Z/U/`*`
(verified), so every protein carrying Prodigal's trailing `*` or an ambiguous `X` would have
vanished without a trace (see U4). **Fix.** `DisorderError` → exit 1 for a missing/crashing
predictor; one batch call; sequences sanitised (trailing `*` stripped, protfasta conversion,
any other character → error; `n_nonstandard` recorded). **Test.**
`test_predict_fails_on_missing_or_crashing_tools` (import error, crash, short/NaN tracks),
`test_sanitize_sequence`, `test_real_metapredict_end_to_end` (p53 + `*`).

### D-06 · major · predict_disorder.py:145-186
**Issue.** IUPred3: `res[0] if isinstance(res,(list,tuple))` took the first *score* of a bare
list (TypeError → silently switched to the CLI path); a failure for one locus returned None,
so that locus became metapredict-only while still labelled `metapredict+iupred3`; library
and CLI paths could use different smoothing defaults. **Fix.** Single library path, explicit
optional `--iupred3-smoothing`, tuple/bare-sequence handling with length/finite/[0,1]
validation, any failure → exit 1; mode restricted to long|short. **Test.**
`test_iupred_return_shapes`, `test_predict_iupred3_library_contract` (smoothing passed,
library default when unset, short/crash → fail).

### D-07 · major · predict_disorder.py:189-213
**Issue.** segmasker missing or failing → None → no LCR masking while `frac_disordered_nonlcr`
was still reported as the LCR-excluded "preferred" metric. **Fix.** Mask requested ⇒ segmasker
required (checked before metapredict runs); non-zero exit, missing records, or intervals
beyond the sequence → error; neutral temp ids (`q1..qN`) so FASTA-id parsing cannot break the
mapping. **Test.** `test_segmasker_interval_parser`, `test_real_segmasker_positions`
(positions equal segmasker's own lower-casing), stub crash/missing/beyond failures.

### D-08 · major · predict_disorder.py:217-234
**Issue.** Missing/unreadable catalog or missing columns returned None → the intact filter was
silently disabled; a missing `fragmented` value was kept. **Fix.** Strict columns
(`locus_id, catalytic_status, fragmented, confidence_tier`); intact ⇔ `catalytic_status ==
intact` AND `fragmented` explicitly false; protein absent from the catalog → error;
`confidence_tier` carried per locus. **Test.** `test_catalog_and_background_loading`,
failure cases "needs --catalog" / "absent from the catalog".

### D-09 · major · predict_disorder.py:61-71
**Issue.** Consensus silently truncated to the shorter track (`min(n, len)`) and treated a
missing second score as "ordered". **Fix.** Every track must have exactly L finite, in-range
values (L = sequence length) else exit 1. **Test.** `test_consensus_rules_and_thresholds`, stub
`short`/`nan`.

### D-10 · minor · predict_disorder.py:68
**Issue.** IUPred threshold applied as `>=`; IUPred's convention is score **> 0.5**.
**Fix.** Per-tool documented conventions: metapredict ≥ thr, IUPred > thr, pLDDT < cutoff.
**Test.** `test_consensus_rules_and_thresholds` (0.5 exactly → ordered for IUPred, disordered for metapredict).

### D-11 · major · predict_disorder.py:101
**Issue.** `frac_disordered_nonlcr` divided by all residues, so masked LCR residues still sat in
the denominator as "ordered" — not "disorder excluding LCR". **Fix.** Disordered non-LCR
residues / non-LCR residues (`''` when every residue is LCR). **Test.** `test_summarize_lcr_semantics`
(5/7), stub CLI (10/20).

### D-12 · major · predict_disorder.py:237-251; disorder.smk:34; disorder_stats.R (whole)
**Issue.** Background column auto-picked (first name containing "background" or "frac"), any
error → {}; the file was not a rule input; partial coverage mixed normalised and raw loci; and
the statistics never used the background at all. **Fix.** Strict `genome_id` +
`background_frac_disordered`|`background_frac` in [0,1], no duplicates, every scored genome
covered; `host_background` column carried and used as a covariate in every phylolm model
(no slope-1 assumption); tracked as a rule input when configured. **Test.**
`test_catalog_and_background_loading`, stub CLI (`disorder_minus_background`, partial table → fail).

### D-13 · minor · predict_disorder.py:126; disorder.smk:36-48
**Issue.** Network version implicit (metapredict's default moved V2→V3 in 3.0), device
auto-selected (CUDA/MPS scores differ slightly), torch threads unbounded on cluster nodes,
one call per sequence. **Fix.** `--metapredict-version V3`, `--device cpu`,
`torch.set_num_threads({threads})` + `OMP/MKL_NUM_THREADS`, batch dict call, provenance in the
new `predictor_detail` column. **Test.** stub asserts V3/cpu/batch and the torch thread cap;
real run checks the provenance line.

### D-14 · minor · predict_disorder.py:319-337
**Issue.** The per-residue table (~10⁷ rows for thousands of ~2,200-aa loci) was built in
memory; no residue letters for downstream identity checks. **Fix.** Streamed writing; additive
`aa` column. **Test.** stub/real known-answer tests read the residue table.

### D-15 · major · map_disorder_to_alignment.py:24-52, 100-101
**Issue.** Own substring name matching (Contract v2 item 3) and a domain track taken from a
single "reference" (longest) locus: when that locus lacks DUF3320 (≈54 % of Prometheus loci)
the DUF3320 columns are labelled linker for everyone. **Fix.** Per-locus categories via
`utils.load_domain_spec` / `utils.domain_intervals`; each column takes the plurality category of
the annotated aligned loci with support ≥ `--min-domain-support` (0.5); ties → NA.
**Test.** `test_map_known_answer` (vote, tie and support values), `test_map_real_fixture_with_mafft`
(real homologues: DUF4011 → DNA2 → REase → DUF3320 in canonical order).

### D-16 · major · map_disorder_to_alignment.py:98-122
**Issue.** No check that aligned rows are the scored sequences: a stale/re-ordered alignment, a
residue edited by the aligner, or a scored locus missing from the alignment shifted or dropped
data silently; ragged rows (family_align copies the raw FASTA for < 3 sequences) were treated
as an MSA. **Fix.** Per-locus residue-count and residue-identity checks; scored locus absent →
error; ragged rows → error unless < 3 sequences (header-only, documented). **Test.**
`test_map_validation_failures_and_empty`.

### D-17 · minor · map_disorder_to_alignment.py:92-96
**Issue.** Bad JSON → `{}` and an empty/unmatched domtbl silently produced an all-NA domain
track (and a "need both domain and linker" placeholder downstream). **Fix.** Empty spec with a
domtbl → error; no aligned locus with a categorised hit → error. **Test.** map failure test.

### D-18 · minor · map_disorder_to_alignment.py:135-136
**Issue.** `mean_disorder` and `frac_loci_disordered` are identical by construction (binary
consensus). **Fix.** Both kept (SCHEMA), documented; additive `occupancy`, `domain_support`,
`mean_meta_score` (continuous score). **Test.** `test_map_known_answer`.

### D-19 · major · plot_disorder.py:26-27
**Issue.** Hard-coded `DOMAIN_COLORS` including the non-existent `nuclease` category.
**Fix.** Categories and order from `--domains` (spec canonical order), Okabe–Ito palette; a
domain value outside the spec → error. **Test.** `test_plot_outputs_categories_and_tier_rule`
(spec categories in the legend, custom `XYZ` category, no `nuclease`, unknown → fail).

### D-20 · major · plot_disorder.py:96-118
**Issue.** Left join + `_bool(NaN)` coded genomes absent from the master table as
non-methanogens; ANME lumped with non-methanogens; no tier rule; `violinplot` on a 1-point or
constant group can fail. **Fix.** Inner join on pass_primary genomes; three groups
(non-methanogen / ANME / methanogen); Contract-v2 item 6 tier filter (`--tiers`) with the
exclusion count in the title; points always, violins only for ≥ 3 distinct values. **Test.**
`test_plot_outputs_categories_and_tier_rule`, `test_plot_tier_filter_regression`.

### D-21 · minor · plot_disorder.py:77, 124, 148; disorder.smk:146-149
**Issue.** Always wrote png+svg while the rule expects `plotting.formats`; one SVG rectangle per
alignment column; gappy low-occupancy columns plotted at full weight. **Fix.** `--formats`;
domain track drawn as contiguous runs; `--min-occupancy` (0.5) for both column figures.
**Test.** plot tests (`--formats png` / `svg`, file sets checked).

### D-22 · critical · disorder_stats.R:47, 86, 136-139
**Issue.** Top-level `tryCatch` wrote a header-only table and `quit(status = 0)` on any error;
a missing `phylolm` was "handled" the same way. **Fix.** No swallowing; required packages
checked with `stop()`; only genuinely empty inputs (no scored locus, `();`/< 3-tip tree) give a
header-only table; base R + ape + phylolm (readr/dplyr dropped). **Test.**
`test_stats_failures_are_errors` (missing phylolm via an empty `R_LIBS`, bad covariate, bad FDR,
missing column, non-binary predictor, non-numeric response, no branch lengths),
`test_stats_skips_low_prevalence_and_handles_empty`.

### D-23 · major · disorder_stats.R:118-123
**Issue.** Temperature covariates entered only the niche models, but methanogen status is itself
confounded with growth temperature (many methanogens are (hyper)thermophiles) and disordered
content tracks temperature. **Fix.** Configured temperature columns adjust every model
(constant/collinear covariates dropped with a note). **Test.** `test_stats_known_answer_and_tier_rule`
(`niche_hot` in the methanogen model; estimate equals an independent phylolm fit).

### D-24 · major · disorder_stats.R:70-75 (Contract v2 item 6)
**Issue.** Per-genome disorder averaged every scored locus regardless of tier; no exclusion
or count. **Fix.** `--tier-sets "A;A,B"` (primary + A+B sensitivity from config): a genome enters
iff its best validated tier ∈ T, only its tier-T loci are averaged, genomes outside T are
excluded and counted (`n_excluded_tier`, note); rows carry `tier_set`. **Test.**
`test_stats_known_answer_and_tier_rule` (B-only genomes excluded under A, present under A+B;
exact n), `test_plot_tier_filter_regression`.

### D-25 · major · disorder_stats.R:94-96
**Issue.** ANME (methanogen = False by construction) sat in the non-methanogen reference group,
so mcr-carrying methanotrophs diluted the methanogen contrast ("ANME separated" in the design).
**Fix.** ANME excluded from the methanogen contrast (counted in the note); shown as their own
group in the figure. **Test.** known-answer test (2 ANME excluded, n 18 → 16).

### D-26 · minor · disorder_stats.R:90, 99, 107
**Issue.** No minimum group size for methanogen; misspelt temperature columns silently ignored
(`intersect`); non-binary niche values coerced; the phylolm error text and λ-at-bound warnings
lost. **Fix.** `--min-group` (3) for every binary predictor; error for absent temperature
columns and non-binary values; notes keep fit errors and every phylolm warning. **Test.**
`test_stats_skips_low_prevalence_and_handles_empty`, failure test, known-answer test
(no warning leaks; boundary λ flagged).

### D-27 · minor · disorder_stats.R (missing test)
**Issue.** The README promises comparison across lineages, but no lineage / phylogenetic-signal
test existed. **Fix.** `phylogenetic_signal` row: Pagel's λ of the disorder trait, LRT against
λ = 0 using the closed-form λ = 0 GLS likelihood with phylolm's tip-variance-preserving
transform; boundary null 0.5·χ²₀ + 0.5·χ²₁. **Test.** `test_stats_lambda0_loglik_matches_phylolm`,
known-answer test (p recomputed from the recorded statistic).

### D-28 · minor · disorder_stats.R:107-112
**Issue.** No interval estimates. **Fix.** Additive `ci_low/ci_high`: 95 % t interval with
df = n − number of coefficients (the df phylolm's own p uses); `lambda`, `n_level1`, `response`.
**Test.** known-answer test (equal to an independent fit to 1e-8; df check).

### D-29 · minor · disorder_stats.R:53-58
**Issue.** Unrooted tree left unrooted when phangorn was absent (phylolm accepts it silently);
branch lengths never validated. **Fix.** Unrooted → phangorn midpoint (required then) and noted;
missing/negative/non-finite branch lengths → error. **Test.** failure test (no branch lengths).

### D-30 · major · disorder.smk:52-68, 135-157, 13, 131
**Issue.** `disorder_map` and `plot_disorder` had no `log:`; the host background was not an
input; tier sets, threads, network version and plot formats were not passed; helper names
`_d` / `_bact_disorder_locus` live in Snakemake's shared namespace. **Fix.** Logs for every rule
(stdout+stderr), `host_bg` input function, new params with `cfg()` defaults, namespaced
`_dis_*` helpers. **Test.** G2 `tests/check_rules.py` (RULES HYGIENE OK), `test_rule_hygiene`,
`test_snakemake_dry_run_disorder_target` (rendered commands).

---

## Verified correct, no change
* **V-1** `disorder_map` maps onto `rules.family_align.output.aln` — the UNtrimmed alignment
  (ClipKIT removes columns, so residue→column mapping must not use the trimmed one).
* **V-2** segmasker interval coordinates are 0-based inclusive; the pristine `+1` conversion was right.
* **V-3** `call_idrs` (runs of consensus 1s ≥ `min_idr_len`, 1-based inclusive) is correct;
  kept as plain runs (no gap closure; metapredict's own IDR caller smooths) — documented.
* **V-4** `utils.genome_of(locus)` on combined `genome__gene` ids; catalog values
  (`catalytic_status` ∈ intact|partial|degenerate|unknown, `fragmented` True/False from the
  Prodigal partial code) match leaf 1.1.2's current validator.
* **V-5** metapredict 0.5 threshold is the documented V2/V3 threshold; alphaPredict pLDDT is on
  the AF2 0–100 scale (pLDDT < 50 = very low confidence).
* **V-6** `results_registry.py` (synthesis leaf) reads `predictor, level, model, estimate,
  p_value` + `tier_set` (`canonical_tiers` → "A", "A,B"): compatible with the new table.

## Unresolved / for the driver
* **U-1** IUPred3 source is licence-gated and was not obtainable; the library contract
  (`iupred(seq, mode[, smoothing])` → `(scores, text)` or scores) is exercised with a stub and
  guarded at run time (length/range checks). The library's own smoothing default is unknown,
  hence the explicit optional `disorder.iupred3_smoothing`.
* **U-2** G3 cannot pass inside this sandbox: `.unlazy/toolchain.env` points micromamba at the
  driver's read-only cache ("Could not find any writable cache directory for repodata file").
  The identical command with a writable cache solves (rc 0); all three pip pins return HTTP 200.
* **U-3** The linux-64 build of the env (pip compiling metapredict's sdist with conda gcc) was
  not executed here (no linux host); locally the osx-arm64 wheel was used.
* **U-4** Whether the validated proteins still carry Prodigal's trailing `*` was not checked
  (the read was blocked by the safety check); the predictor handles both.
* **U-5** Repo-root `cache/` (≈23 MB of `mamba*` temp files) appears when check_env_solve.sh runs
  with a read-only `CONDA_PKGS_DIRS` (micromamba falls back to `./cache`); not deleted (outside
  owned paths).
* **U-6** The fallback consensus is not independent of metapredict (D-02 evidence): on p53 it
  collapses to the pLDDT < 50 call. Stage IUPred3, or evaluate AIUPred (IUPred3's successor,
  pip-installable from GitHub) as an installable independent second predictor.
