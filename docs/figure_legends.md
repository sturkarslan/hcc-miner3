# Figure explanations and legends

Figures: `results/09_figures/figure1.{pdf,png}`, `figure2.{pdf,png}` (183 mm, Nature double column; Liberation Sans
embedded as TrueType). Source data: `figure1_source_data.xlsx`, `figure2_source_data.xlsx`. Built by
`scripts/09_publication_figures.py` from saved results only; program names from `config/program_labels.tsv`
(each with its evidence).

> **Pending re-run (2026-10-01).** Published classes and signatures now come from the reference panel
> (`config/reference_panel.yaml`): Montironi 2023 immune classes, Sia 2017 immune class, Haber 2023 IFNAP,
> Zhu 2022 ABRS, Gao 2019 proteogenomic axes, Désert 2017 zonation classes and immune-cell sets first, with
> Hoshida 2009 / Boyault 2007 / Chiang 2008 kept as references. The NTP null distribution was also corrected
> (step 07b). Panel descriptions below are updated; **numbers and the narrative claims about classes and the
> head-to-head comparison are from the previous run** and must be refreshed after
> `bash scripts/slurm/rerun_reference_panel.sh` and curation of `config/program_labels.tsv`.

---

## Figure 1

### What the figure shows

We built a regulatory network of hepatocellular carcinoma from 929 tumours in three cohorts (US, Chinese and French)
and asked what it captures. The network organises about 8,400 genes into 76 transcriptional programs. Those programs
fall into three familiar blocks of liver-cancer biology: proliferation and progenitor features; hepatocyte
differentiation and WNT/β-catenin signalling; and immune and stromal content. Tumour states defined by the network
line up with published HCC classes: states high in proliferation programs are enriched for proliferative classes and
TP53 mutations, while states high in differentiation programs are enriched for hepatocyte-like and CTNNB1-mutant
tumours. Causal inference connects driver mutations to specific regulators and programs. For example, CTNNB1 mutation
acts through the Wnt effectors LEF1 and TCF7 to switch on the β-catenin and low-recurrence programs. The network also separates
two drivers that are usually grouped together as WNT-pathway mutations. CTNNB1-mutant tumours switch on the
WNT/β-catenin programs and the classic liver targets GLUL, AXIN2 and NKD1. AXIN1-mutant tumours do not: they show
only a partial rise in a few targets (LGR5, TBX3, RNF43, ZNRF3). Both groups are low in immune and stromal programs,
but the detail differs: AXIN1 mutants have lower MHC class I antigen-presentation genes (B2M, HLA-A/B, TAP1, PSMB9),
which CTNNB1 mutants do not, whereas T-cell genes are lower mainly in CTNNB1 mutants. Ten programs move in opposite
directions in the two groups, and the causal layer places LEF1 and TCF7 up under CTNNB1 and down under AXIN1. Low WNT
output in AXIN1-mutant HCC and reduced immune infiltration in both groups have been reported before; the split
between the two immune-low patterns and the opposite programs is the part that needs independent validation.

### Legend

**Fig. 1 | A causal and mechanistic regulatory network of hepatocellular carcinoma.**
**a**, Study design. Expression from 929 primary HCCs (TCGA-LIHC, n = 366; CLCA, n = 239; LICA-FR, n = 324) was
harmonized to 13,866 genes and batch-corrected (ComBat). MINER inferred regulons, transcriptional programs and
states; regulons derived from four co-expression modules that tracked a suspected technical signal (186 of 4,294)
were excluded. Causal inference linked 113 genomic features to regulators and
regulons. Ridge risk models on program activity were trained in TCGA or CLCA and tested in the other, and in three
external cohorts (GSE14520, n = 221; LIRI-JP, n = 203; GSE76427, n = 115) with fixed weights.
**b**, Network size at each level (log scale).
**c**, Mean regulon dysregulation (over- minus under-expressed membership) of each program (rows) in each state
(columns). States are hierarchically clustered (average linkage, correlation distance); programs are grouped into
biology blocks (anchored on the reference panel: proliferation/progenitor, differentiated/WNT, immune/stromal).
Tracks: mean risk score per state; fraction of tumours in published HCC classes, labelled by biology with the class
code or source in grey: Montironi 2023 immune classes (Inflamed 20-gene signature by nearest-template prediction, then
Sia 2017 immune class within inflamed tumours and CTNNB1 mutation within non-inflamed tumours), and the reference
Hoshida S1–S3 and Boyault G1–G6 classes by nearest-template prediction; mutation frequencies.
**d**, Pearson correlation across tumours between program activity and reference-panel signature scores (immune and
immunotherapy-response signatures from Montironi 2023, Sia 2017, Haber 2023 and Zhu 2022; WNT/β-catenin activation;
Gao 2019 proteogenomic axes; Désert 2017 zonation classes; selected hallmarks; Hoshida classes as references), with
each tumour's mean expression regressed out of both. Program names are coloured by the sign of their risk-model
weight (red adverse, blue protective). Tags give the source; H, MSigDB hallmark.
**e**, Causal flows from CTNNB1 mutation (n = 222 tumours) to regulators and programs: the regulators with the largest
effects (one per regulon family) plus LEF1 and TCF7. Edges: red, up in mutant tumours; blue, down; dashed, the
regulator represses its regulon; width proportional to |Cohen's d|. Flows shown are high confidence: Benjamini–Hochberg
q ≤ 0.1 for the regulon and driver–regulator tests, same direction in every cohort tested, |d| ≥ 0.5 and ≥ 30
altered tumours.
**f**, AXIN1 and CTNNB1 mutations have different program effects. Left, effect of AXIN1-only (n = 77) and
CTNNB1-only (n = 218) mutation on each of the 76 programs relative to tumours wild type for both (n = 620), in s.d.
of program activity (linear model adjusted for TP53 mutation, proliferation program activity and cohort; four tumours
with both mutations excluded). Colours: significant for CTNNB1 only, AXIN1 only, both in the same direction (shared),
both in opposite directions, or neither (Benjamini–Hochberg q ≤ 0.05 over programs). Right, the same model for the
expression of WNT target, MHC class I and T-cell genes.

---

## Figure 2

### What the figure shows

A risk score built from the network's programs separates patients by recurrence and survival across six independent
test settings. Its largest weights fall on proliferation, MYC and progenitor programs (adverse) and on
hepatocyte-differentiation, low-recurrence and immune programs (protective). The score reproduces known HCC biology:
it is highest in proliferative and TP53-related classes and lowest in hepatocyte-like and WNT/β-catenin classes. That
holds even against the author-assigned classes of a cohort it was never trained on. Network states ordered by this
score show a matching gradient of observed recurrence. The model was validated across cohorts in discovery and in
three external cohorts spanning Affymetrix, RNA-seq and Illumina platforms and hepatitis B, hepatitis C and mixed
aetiologies. Validation is strongest in LIRI-JP, where the top-risk fifth had a 5.3-fold higher death rate. The score
performs comparably to, but not better than, published signatures such as the proliferation class, so its value lies
in explaining risk mechanistically rather than in predicting it more accurately. When the network is rebuilt without
each discovery cohort, regulators and key programs are recovered and the risk model performs just as well in the
unseen cohort.

### Legend

**Fig. 2 | Program-based risk, its biology and its validation.**
**a**, Largest adverse (red) and protective (blue) weights of the ridge model trained on TCGA recurrence.
**b**, Risk score (within-cohort z) by published HCC class, labelled by biology: Montironi 2023 immune classes and
Hoshida classes (calls for all 929 tumours, as in Fig. 1c), and the LICA-FR authors' Boyault labels (n = 324).
Bars, median and interquartile range; Kruskal–Wallis P.
**c**, States ordered by mean risk score. Top, class and mutation fractions; middle, risk score; bottom, observed recurrence as GuanRank computed within each cohort (TCGA and CLCA; 1,
earliest event). Boxes, median and interquartile range; whiskers, 1.5× interquartile range. State mean risk versus
median GuanRank: Spearman ρ = 0.64, P = 2.2 × 10⁻⁴ (28 states with ≥ 5 tumours).
**d**, **e**, Kaplan–Meier curves for the within-cohort top 20% of risk scores versus the rest. **d**, Cross-cohort
validation in discovery (TCGA-trained model in CLCA; CLCA-trained model in TCGA). **e**, External validation with the
pre-specified TCGA-trained models, never refit. C, Harrell's C-index; HR, Cox hazard ratio of top 20% versus rest.
**f**, Cox hazard ratio per s.d. of risk score (95% confidence interval) in every test cohort, with random-effects
(DerSimonian–Laird) pooled estimates over external cohorts.
**g**, C-index of the network model compared with Cox models on published signature scores, trained in the same cohort:
Hoshida classes (2009); Chiang proliferation signature (2008); immune and WNT classes (Montironi 2023 Inflamed,
Sia 2017 immune class, Haber 2023 IFNAP, WNT/β-catenin activation); Gao 2019 proteogenomic prognosis axes (PYCR2,
ADH1A); and all published scores together (also TP53-mutant proteome, Désert stem-like, Hoshida S1–S3, Chiang
proliferation and CTNNB1, KRT19, plus tumour mean expression). Adding the network score to all known scores: likelihood-ratio P ≥ 0.10 in every
cohort.
**h**, Leave-one-cohort-out stability: the network, causal inference and risk model were rebuilt without each
discovery cohort. Bars show regulators recovered; median correlation, in the held-out cohort, between each program's
activity and its best leave-out match (all programs; the quarter with the largest risk weights); and high-confidence
driver–regulator edges recovered (CTNNB1; all drivers). Text, C-index of models built entirely
without the test cohort.

---

## Panel → source

| Panel | Result files (under `results/`) |
|---|---|
| 1b | `04_miner/combat/{coexpr,mechinf}`, `05_causal/combat/{highConfidenceCausalResults.csv,regulon_families.tsv}` |
| 1c | `04_miner/combat/subtypes_filtered/`, `07_post/subtype_mapping/subtypes_filtered/ntp_calls_*.tsv`, `03_genomics_clinical/genomic_features.csv` |
| 1d | `07_post/subtype_mapping/subtypes_filtered/program_signature_correlation.tsv` |
| 1e | `05_causal/combat/highConfidenceCausalResults.csv` |
| 1f | `07_post/axin1_ctnnb1/{program_effects.tsv,gene_effects.tsv,groups.tsv}` (step 07f) |
| 2a | `06_risk/combat/predictor_ridge_programs_TCGA_RFS_h36m/weights.tsv` |
| 2b, 2c | `06_risk/.../predictions.tsv`, NTP calls, `01_harmonized/samples.tsv`, `03_genomics_clinical/survival_*_h36m_miner.csv` |
| 2d | `06_risk/combat/predictor_ridge_programs_{TCGA,CLCA}_RFS_h36m/` |
| 2e–g | `08_validation/<cohort>/{scores,evaluation,head_to_head}.tsv`, `07_post/figures/risk_head_to_head.tsv` |
| 2h | `08_validation/loco/{loco_summary.tsv,loco_programs.tsv}` |
