# Methanogenesis pathway rules — rationale and evidence

Every rule in `methanogenesis_rules.yaml` is listed here with what it asserts, why, and a
primary citation. A rule may be changed only together with its entry here and a re-run of the
real-genome known-answer test (`tests/methanogenesis`, fixture `real_genome_ko_calls.json`).

The calls are **genomic capacity**, not culture phenotype. A genome can carry a complete,
correct gene set for a substrate that the organism does not use in culture; §6 lists the known
cases in the reference set. Any statement derived from these columns must be phrased as
capacity ("encodes the methanol pathway"), never as growth.

## 1. What the caller decides

For each methanogen (credible Mcr placed with the methanogen references; `classify_methanogen.py`):

| Pathway | Meaning |
|---|---|
| `co2_reducing` | CO2 + H2 or formate -> CH4 through the full Wood-Ljungdahl methyl branch and Mtr |
| `acetoclastic` | acetate -> CH4 through CODH/ACS |
| `methylotrophic` | at least one methylated substrate is activated to methyl-CoM |
| `methyl_dismutating` | methylated substrate without an external electron donor (some methyl groups oxidised back to CO2 through Mtr run in reverse) |
| `methyl_reducing` | methylated substrate reduced with H2 or formate; the oxidative route is not available |

plus substrate flags `sub_h2`, `sub_formate`, `sub_methanol`, `sub_methylamine` (two-sided) and
`sub_methylthiol`, `sub_quatamine`, `sub_methoxyl` (one-sided; see §4).

Pathways are **independent indicators, not one category**: *Methanosarcina barkeri* is
CO2-reducing, acetoclastic and methylotrophic at once (Maeder et al. 2006,
doi:10.1128/JB.00810-06), so forcing one label per genome would make the trait a function of an
arbitrary precedence rule rather than of the biology.

## 2. KO detection

KO calls use the KofamScan rule exactly (`kofam.py`, shared with `mcr_search.py`): a hit is
accepted when its score on the KO's own `score_type` (full-sequence or best-domain) reaches the
curated KOfam threshold (Aramaki et al. 2019, doi:10.1093/bioinformatics/btz859). There is no
E-value fallback and no uncurated profile: a KO without a threshold is an error.

Three additions, each of which exists because a plain threshold rule produces a **false absence**:

* **Contig-edge fragments** (`partial_rescue`): a gene running off a contig end is accepted at
  `partial_score_fraction` x threshold, because truncation, not divergence, explains the lower
  score. Prodigal's partial flags come from `genes.tsv` (Contract v2 item 2).
* **Stop-codon read-through** (`readthrough_TAG`, `readthrough_TGA`), §3.
* **Near-threshold evidence** (`near_threshold`), §5.

## 3. Stop-codon read-through: pyrrolysine and selenocysteine

`readthrough.TAG` (residue `O`, pyrrolysine): `K16176` *mtmB*, `K16178` *mtbB*, `K14083` *mttB*.
The methylamine methyltransferases carry a single in-frame UAG read as pyrrolysine (Hao et al.
2002, doi:10.1126/science.1069556; Srinivasan et al. 2002, doi:10.1126/science.1069588).
Prodigal ends the gene at that codon, so the protein arrives as an N-terminal fragment that
fails a full-length threshold.

`readthrough.TGA` (residue `U`, selenocysteine): `K00201` *fwdB*, `K22516` *fdhA*, `K03388`
*hdrA*, `K00440` *frhA*, `K14126`/`K14127`/`K14128` *mvhADG*. Methanococcales and Methanopyrus
encode selenocysteine forms of these enzymes (Quitzke et al. 2018,
doi:10.1016/j.bbagen.2018.03.030).

A fragment is re-translated through the codon only when it is a homologous hit of one of those
KOs, its 3' end is complete (there is a stop codon to read through) and no internal stop
intervenes; the extended protein must then pass the **full** threshold. Prodigal genes lying
in-frame inside the extended span are recorded as `merged_proteins` and not counted again.

Evidence that this matters: in the 12-genome reference run, read-through recovered `K16176`
(MtmB) in *Ca.* Methanomethylophilus alvi, and `K00201` (FwdB) plus `K03388` (HdrA) in
*Methanocaldococcus jannaschii* and *Methanococcus maripaludis* — i.e. KOs whose absence would
otherwise have removed the Fwd and Hdr components of a textbook hydrogenotroph.

## 4. Components

`relax` makes a complex tolerant: the component is present when every `signature` entry is
satisfied (a string = that KO, a list = at least one of them) and at least `min_fraction` of
the mandatory subunit slots are detected. **A signature is always a set of interchangeable core
subunits, never a single subunit**: KEGG lists neither FwdB (`K00201`) nor HdrA (`K03388`) for
*M. jannaschii*, and a single-subunit signature therefore collapsed that genome to
`none_detected` during development.

`one_sided: true` marks a component whose absence is not evidence (`mts`, `quat`, `methoxyl`):
these functions have enzymes outside the KO set, so the flag is reported as `not_detected`, never
as `absent`.

| Component | KOs | Assertion and evidence |
|---|---|---|
| `fwd` | K00200 K00201 K00202 K00203 (K11261) (K00205 / K11260 / K00204) | formylmethanofuran dehydrogenase, M00567 step 1. Signature: FwdA (K00200) or FwdC (K00202) |
| `ftr` | K00672 | formyltransferase, M00567 step 2 |
| `mch` | K01499 | methenyl-H4MPT cyclohydrolase, M00567 step 3 |
| `mtd` | K00319 or K13942 | F420-dependent dehydrogenase Mtd, or the H2-forming Hmd, M00567 step 4 |
| `mer` | K00320 | methylene-H4MPT reductase, M00567 step 5 |
| `mtr` | K00577-K00584 (K00582, K00583 optional) | Na+-translocating methyl-H4MPT:CoM methyltransferase, M00567 step 6. Signature: MtrA or MtrH |
| `mcr` | K00399 K00401 K00402 | methyl-CoM reductase; informational only (methanogen status comes from the McrA placement) |
| `moco` | M00880 | molybdopterin / tungstopterin cofactor; required for a functional Fwd/Fmd, §5 |
| `hdr_abc`, `hdr_abc2` | K03388-K03390, K22480-K22482 | cytoplasmic heterodisulfide reductase (electron bifurcation; Kaster et al. 2011, doi:10.1073/pnas.1016761108) |
| `hdr_de` | K08264 K08265 | membrane-bound HdrDE with cytochrome b, the methanogens-with-cytochromes architecture (Thauer et al. 2008, doi:10.1038/nrmicro1931) |
| `frh`, `mvh`, `vht`, `ech` | K00440/K00441/K00443, K14126-K14128, K14068-K14070, K14086-K14091 | hydrogen uptake: F420-reducing, F420-non-reducing, methanophenazine-reducing, energy-converting |
| `fdh` | K22516 K00125 | F420-dependent formate dehydrogenase |
| `ackpta` | K00925 K00625 | acetate kinase then phosphotransacetylase, the low-affinity acetate activation of *Methanosarcina*. Relaxed on AckA, §5 |
| `acs` | K01895 | AMP-forming acetyl-CoA synthetase, the high-affinity route of *Methanothrix* (Smith & Ingram-Smith 2007, doi:10.1016/j.tim.2007.02.002) |
| `codh` | K00193 K00197 K00194 | CO dehydrogenase / acetyl-CoA synthase, which cleaves acetyl-CoA |
| `mta` | K14080 K04480 K14081 | methanol:CoM methyltransferase system MtaABC (Sauer et al. 1997, doi:10.1111/j.1432-1033.1997.t01-1-00670.x) |
| `methylamine_mt` | K14082 with (K16177/K16176), (K16179/K16178) or (K14084/K14083) | MtbA plus a mono-, di- or trimethylamine corrinoid protein; the methyltransferases are optional exactly as in KEGG M00563, because they are the pyrrolysine enzymes |
| `pyl` | K16180 K16181 K16182 (K11627) | pyrrolysine biosynthesis PylBCD, §5 |
| `mts` | K16954 K16955 | methylthiol:CoM methyltransferase (Tallant & Krzycki 1997, doi:10.1128/jb.179.22.6902-6911.1997). One-sided |
| `quat` | K25220, K25600 or K26847 | quaternary-amine methyltransferases: glycine betaine MtgB (Ticak et al. 2014, doi:10.1073/pnas.1409642111), carnitine MtcB, proline betaine MtpB. One-sided |
| `methoxyl` | K25916 K25917 | methoxylated-aromatic methyltransferases (*Methermicoccus*; Mayumi et al. 2016, doi:10.1126/science.aaf8821). One-sided |

Derived terms: `methyl_branch` = fwd+ftr+mch+mtd+mer; `functional_methyl_branch` = that plus
`moco`; `hdr` = any heterodisulfide reductase; `h2_uptake` = any uptake hydrogenase;
`electron_input` = `h2_uptake` or `fdh`; `oxidative_methyl_route` =
`functional_methyl_branch` + `mtr` + `hdr_de`; `methylamine` = `methylamine_mt` + `pyl`;
`methyl_substrate` = any methylated-substrate system; `acetate_core` = `codh` + `mtr` + `hdr_de`.

## 5. The five decisions that keep the calls honest

1. **Fwd without its cofactor is not CO2 reduction.** `co2_reducing` requires
   `functional_methyl_branch`, i.e. the methyl branch **and** molybdopterin biosynthesis
   (M00880). *Methanosphaera stadtmanae* carries a complete `fwdABCDG` operon yet cannot reduce
   CO2: it lacks molybdopterin biosynthesis, so the enzyme cannot be matured, and it is
   restricted to methanol + H2 (Fricke et al. 2006, doi:10.1128/JB.188.2.642-658.2006). A
   KO-presence rule calls that genome hydrogenotrophic; this one does not.
2. **A methylamine call needs the pyrrolysine machinery.** `methylamine` = the KEGG M00563
   methyltransferase set **and** `pyl`. The MttB superfamily (COG5598) is large and most members
   are *not* methylamine enzymes: they demethylate quaternary amines such as glycine betaine,
   carnitine and proline betaine and lack pyrrolysine (Ticak et al. 2014,
   doi:10.1073/pnas.1409642111). Those hits set `sub_quatamine`, not `sub_methylamine`.
3. **Acetyl-CoA synthetase is also anabolic.** CODH/ACS and Acs serve carbon assimilation in
   hydrogenotrophs, so `acetoclastic` is called only with the whole acetate-utilising
   architecture (`acetate_core` = CODH/ACS + Mtr + HdrDE) plus acetate activation: AckA/Pta, or
   Acs in a genome that is neither CO2-reducing nor methylotrophic. An Acs-only methylotroph
   therefore comes out **unresolved** (`maybe(...)`) rather than falsely present or absent —
   the honest answer, because its Acs may be anabolic (*Methanococcoides burtonii*,
   *Methanomethylovorans hollandica* in the reference run).
4. **Absence is gated on completeness.** A missing component contributes
   p = (1 - completeness)^k for its k missing loci, combined over a rule with product (`all`),
   maximum (`any`) and complement (`not`); a pathway is called absent only when the resulting
   upper bound is <= `absence_alpha` (default 0.05) and present only when the lower bound is
   >= 1 - alpha. One missing locus in a 60%-complete MAG gives p = 0.4 and stays `unresolved`;
   the same locus in a 99%-complete genome gives 0.01 and is an absence. Completeness is
   CheckM2's estimate, and the bound assumes loci are missed independently.
5. **A threshold just above a divergent orthologue is not an absence.** A KO whose best hit
   reaches >= `near_threshold_fraction` (0.9) of the threshold is recorded as `near_threshold`:
   never presence, but enough to withdraw a confident absence (the component becomes
   `equivocal`, interval [0, 1]). Real case: KOfam `K00625` (phosphotransacetylase) scores 379.0
   in *Methanosarcina acetivorans* C2A against a curated threshold of 387.57, a 2% shortfall in
   a textbook acetoclastic methanogen. `ackpta` is additionally relaxed on AckA (`K00925`) for
   the same reason; AckA together with `acetate_core` is the diagnostic combination.

## 6. Reference set and known capacity-vs-phenotype discordances

The rules reproduce the published physiology of 12 cultured methanogens
(`tests/methanogenesis/fixtures/real_genome_ko_calls.json`; assemblies -> Prodigal 2.6.3 ->
`methanogenesis_search.py` with the curated KOfam subset): 34 of 36 asserted calls resolved and
correct, 0 mismatches, 2 unresolved (the Acs-only methylotrophs of decision 3).

Two reference genomes carry a pathway they do not use in culture. Both are reported as
gene-level positives, which is the correct behaviour for a capacity call:

* *Methanobrevibacter smithii* ATCC 35061 encodes `mtaABC` yet grows only on H2/formate + CO2;
  it uses neither methanol nor methylamines (Samuel et al. 2007, doi:10.1073/pnas.0704189104).
* *Methanosarcina acetivorans* C2A encodes the methyl branch, Mtr and hydrogenases yet does not
  grow on H2 + CO2, because it does not use H2 as an electron donor (Guss et al. 2009,
  doi:10.1128/JB.00563-08).

## 7. What these rules do not decide

No rule here distinguishes methanogenesis from anaerobic methane oxidation: ANME genomes carry a
near-complete reverse pathway, so they are excluded upstream by the McrA placement
(`classify_methanogen.py`) and reported as `not_applicable`. Methoxydotrophic methanogenesis has
no KEGG module; `sub_methoxyl` uses the *Methermicoccus* KOs and is one-sided.

## 8. Every KO used, with its component

KO names are KEGG's own (`https://rest.kegg.jp/list/ko`, retrieved 2026-10-31).
`RT` marks a KO whose gene may be read through a stop codon (§3).

| KO | gene / name | component(s) | RT |
|---|---|---|---|
| `K00125` | fdhB — formate dehydrogenase (coenzyme F420) beta subunit | fdh |  |
| `K00193` | cdhC — acetyl-CoA decarbonylase/synthase, CODH/ACS complex subunit beta | codh |  |
| `K00194` | cdhD, acsD — acetyl-CoA decarbonylase/synthase, CODH/ACS complex subunit delta | codh |  |
| `K00197` | cdhE, acsC — acetyl-CoA decarbonylase/synthase, CODH/ACS complex subunit gamma | codh |  |
| `K00200` | fwdA, fmdA — formylmethanofuran dehydrogenase subunit A | fwd |  |
| `K00201` | fwdB, fmdB — formylmethanofuran dehydrogenase subunit B | fwd | TGA |
| `K00202` | fwdC, fmdC — formylmethanofuran dehydrogenase subunit C | fwd |  |
| `K00203` | fwdD, fmdD — formylmethanofuran dehydrogenase subunit D | fwd |  |
| `K00204` | fwdH — 4Fe-4S ferredoxin | fwd |  |
| `K00205` | fwdF, fmdF — 4Fe-4S ferredoxin | fwd |  |
| `K00319` | mtd — methylenetetrahydromethanopterin dehydrogenase | mtd |  |
| `K00320` | mer — 5,10-methylenetetrahydromethanopterin reductase | mer |  |
| `K00399` | mcrA — methyl-coenzyme M reductase alpha subunit | mcr |  |
| `K00401` | mcrB — methyl-coenzyme M reductase beta subunit | mcr |  |
| `K00402` | mcrG — methyl-coenzyme M reductase gamma subunit | mcr |  |
| `K00440` | frhA — coenzyme F420 hydrogenase subunit alpha | frh | TGA |
| `K00441` | frhB — coenzyme F420 hydrogenase subunit beta | frh |  |
| `K00443` | frhG — coenzyme F420 hydrogenase subunit gamma | frh |  |
| `K00577` | mtrA — tetrahydromethanopterin S-methyltransferase subunit A | mtr |  |
| `K00578` | mtrB — tetrahydromethanopterin S-methyltransferase subunit B | mtr |  |
| `K00579` | mtrC — tetrahydromethanopterin S-methyltransferase subunit C | mtr |  |
| `K00580` | mtrD — tetrahydromethanopterin S-methyltransferase subunit D | mtr |  |
| `K00581` | mtrE — tetrahydromethanopterin S-methyltransferase subunit E | mtr |  |
| `K00582` | mtrF — tetrahydromethanopterin S-methyltransferase subunit F | mtr |  |
| `K00583` | mtrG — tetrahydromethanopterin S-methyltransferase subunit G | mtr |  |
| `K00584` | mtrH — tetrahydromethanopterin S-methyltransferase subunit H | mtr |  |
| `K00625` | pta — phosphate acetyltransferase | ackpta |  |
| `K00672` | ftr — formylmethanofuran--tetrahydromethanopterin N-formyltransferase | ftr |  |
| `K00925` | ackA — acetate kinase | ackpta |  |
| `K01499` | mch — methenyltetrahydromethanopterin cyclohydrolase | mch |  |
| `K01895` | ACSS1_2, acs — acetyl-CoA synthetase | acs |  |
| `K03388` | hdrA2 — heterodisulfide reductase subunit A2 | hdr_abc | TGA |
| `K03389` | hdrB2 — heterodisulfide reductase subunit B2 | hdr_abc |  |
| `K03390` | hdrC2 — heterodisulfide reductase subunit C2 | hdr_abc |  |
| `K03635` | MOCS2B, moaE — molybdopterin synthase catalytic subunit | moco |  |
| `K03637` | moaC, CNX3 — cyclic pyranopterin monophosphate synthase | moco |  |
| `K03638` | moaB — molybdopterin adenylyltransferase | moco |  |
| `K03639` | moaA, CNX2 — GTP 3',8-cyclase | moco |  |
| `K03750` | moeA — molybdopterin molybdotransferase | moco |  |
| `K03831` | mogA — molybdopterin adenylyltransferase | moco |  |
| `K04480` | mtaB — methanol---5-hydroxybenzimidazolylcobamide Co-methyltransferase | mta |  |
| `K08264` | hdrD — heterodisulfide reductase subunit D | hdr_de |  |
| `K08265` | hdrE — heterodisulfide reductase subunit E | hdr_de |  |
| `K11260` | fwdG — 4Fe-4S ferredoxin | fwd |  |
| `K11261` | fwdE, fmdE — formylmethanofuran dehydrogenase subunit E | fwd |  |
| `K11627` | pylS — pyrrolysyl-tRNA synthetase | pyl |  |
| `K13942` | hmd — 5,10-methenyltetrahydromethanopterin hydrogenase | mtd |  |
| `K14068` | vhoA, vhtA — methanophenazine hydrogenase, large subunit | vht |  |
| `K14069` | vhoC, vhtC — methanophenazine hydrogenase, cytochrome b subunit | vht |  |
| `K14070` | vhoG, vhtG — methanophenazine hydrogenase | vht |  |
| `K14080` | mtaA | mta |  |
| `K14081` | mtaC — methanol corrinoid protein | mta |  |
| `K14082` | mtbA | methylamine_mt |  |
| `K14083` | mttB — trimethylamine---corrinoid protein Co-methyltransferase | methylamine_mt | TAG |
| `K14084` | mttC — trimethylamine corrinoid protein | methylamine_mt |  |
| `K14086` | echA — ech hydrogenase subunit A | ech |  |
| `K14087` | echB — ech hydrogenase subunit B | ech |  |
| `K14088` | echC — ech hydrogenase subunit C | ech |  |
| `K14089` | echD — ech hydrogenase subunit D | ech |  |
| `K14090` | echE — ech hydrogenase subunit E | ech |  |
| `K14091` | echF — ech hydrogenase subunit F | ech |  |
| `K14126` | mvhA, vhuA, vhcA — F420-non-reducing hydrogenase large subunit | mvh | TGA |
| `K14127` | mvhD, vhuD, vhcD — F420-non-reducing hydrogenase iron-sulfur subunit | mvh | TGA |
| `K14128` | mvhG, vhuG, vhcG — F420-non-reducing hydrogenase small subunit | mvh | TGA |
| `K15376` | GPHN — gephyrin | moco |  |
| `K16176` | mtmB — methylamine---corrinoid protein Co-methyltransferase | methylamine_mt | TAG |
| `K16177` | mtmC — monomethylamine corrinoid protein | methylamine_mt |  |
| `K16178` | mtbB — dimethylamine---corrinoid protein Co-methyltransferase | methylamine_mt | TAG |
| `K16179` | mtbC — dimethylamine corrinoid protein | methylamine_mt |  |
| `K16180` | pylB — methylornithine synthase | pyl |  |
| `K16181` | pylC — 3-methylornithine--L-lysine ligase | pyl |  |
| `K16182` | pylD — 3-methylornithyl-N6-L-lysine dehydrogenase | pyl |  |
| `K16954` | mtsA — methylthiol:coenzyme M methyltransferase | mts |  |
| `K16955` | mtsB — methylated-thiol--corrinoid protein | mts |  |
| `K20967` | MOCS1 — GTP 3',8-cyclase / cyclic pyranopterin monophosphate synthase | moco |  |
| `K21142` | moaX — MoaE-MoaD fusion protein | moco |  |
| `K22480` | hdrA1 — heterodisulfide reductase subunit A1 | hdr_abc2 |  |
| `K22481` | hdrB1 — heterodisulfide reductase subunit B1 | hdr_abc2 |  |
| `K22482` | hdrC1 — heterodisulfide reductase subunit C1 | hdr_abc2 |  |
| `K22516` | fdhA — formate dehydrogenase (coenzyme F420) alpha subunit | fdh | TGA |
| `K25220` | mtgB — glycine betaine---corrinoid protein Co-methyltransferase | quat |  |
| `K25600` | mtcB — L-carnitine---corrinoid protein Co-methyltransferase | quat |  |
| `K25916` | mtvB, mtoB — methoxylated aromatic compound---corrinoid protein Co-methyltransferase | methoxyl |  |
| `K25917` | mtoA | methoxyl |  |
| `K26847` | mtpB — proline betaine --- corrinoid protein Co-methyltransferase | quat |  |
