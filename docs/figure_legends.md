# Figure legends (draft)

Figures: `results/09_figures/figure1.{pdf,png}`, `figure2.{pdf,png}` (183 mm, Nature double column; Liberation Sans,
fonts embedded as TrueType). Source data: `figure1_source_data.xlsx`, `figure2_source_data.xlsx`.
Built by `scripts/09_publication_figures.py` from saved results only; program names from `config/program_labels.tsv`
(each with its evidence).

## Figure 1 | A causal and mechanistic regulatory network of hepatocellular carcinoma

**a**, Study design. The network was inferred from 929 primary HCCs in three cohorts (TCGA-LIHC, n = 366; CLCA,
n = 239; LICA-FR, n = 324): expression was harmonized to 13,866 genes and batch-corrected (ComBat), and a technical
RNA-quality axis was removed. MINER then inferred mechanistic regulons, transcriptional programs and states. Causal
inference linked 113 genomic features (driver mutations, pathways, TERT promoter, focal and arm-level copy number)
to regulators and regulons. Risk models (ridge regression on program activity; 36-month horizon) were trained in
TCGA or CLCA and tested in the other. Three external cohorts (GSE14520, n = 221; LIRI-JP, n = 203; GSE76427, n = 115)
were scored with fixed weights.

**b**, Network size at each level (log scale).

**c**, Mean regulon dysregulation (over- minus under-expressed membership) of each program in each transcriptional
state. States are clustered (average linkage, correlation distance); programs are grouped by the biology block of
their best-matching signatures. Tracks show each state's mean risk score, the fraction of tumours assigned to
published HCC classes by nearest-template prediction (labelled by the biology each class represents; class codes
in grey: Hoshida S1–S3, Boyault G1–G6, Chiang (C)) and driver alteration frequencies.

**d**, Correlation across tumours between program activity and published signature scores (per-sample global mean
regressed out of both). Rows are curated program names, coloured by the direction of the program's weight in the
risk model (red adverse, blue protective). Columns are grouped by biology: proliferation/progenitor,
differentiated/WNT, immune/stromal. H, MSigDB hallmark.

**e**, Causal flows from CTNNB1 mutation (n = 222): the regulators with the largest effects (one per regulon family,
plus the canonical Wnt effectors LEF1 and TCF7) and the programs they reach. Edge colour: up (red) or down (blue) in
mutant tumours; dashed: the regulator represses its regulon; width ∝ |Cohen's d|. Program box border: adverse (red)
or protective (blue) in the risk model. High-confidence flows: BH q ≤ 0.1 for the regulon and driver–regulator
tests, consistent direction in every cohort tested, |d| ≥ 0.5, and ≥ 30 altered tumours.

**f**, Each driver's "net risk push" (sum over its causal regulon families of sign(d) × the program's risk weight;
no survival data used) against its observed association with recurrence (Cox z per cohort, Stouffer meta-z over
TCGA and CLCA; copy-number features TCGA only). Spearman ρ = 0.64, P = 0.001, 23 drivers.

## Figure 2 | Program-based risk, its biology and its validation

**a**, The eight largest adverse (red) and protective (blue) program weights in the TCGA-trained recurrence model.

**b**, Risk score (within-cohort z) by published HCC class, labelled by the biology each class represents
(Kruskal–Wallis P). Hoshida and Chiang classes are nearest-template-prediction calls in all 929 tumours; Boyault
classes are the LICA-FR authors' labels (n = 324), so the model is tested against labels it was never trained on.

**c**, States ordered by mean risk score. Top: fraction of each state in selected classes and with TP53 or CTNNB1
mutations. Middle: risk score. Bottom: observed recurrence as GuanRank computed within each cohort (TCGA + CLCA; 1 =
earliest event). State mean risk vs median GuanRank: Spearman ρ = 0.64, P = 2.2 × 10⁻⁴ (28 states with ≥ 5
tumours).

**d**, Cross-cohort validation in discovery: Kaplan–Meier curves of recurrence for the within-cohort top 20%
predicted risk vs the rest, for the TCGA-trained model in CLCA and the CLCA-trained model in TCGA.

**e**, External validation with the pre-specified TCGA-trained models (never refit): GSE14520 OS, LIRI-JP OS,
GSE76427 recurrence.

**f**, Hazard ratio per s.d. of risk score (95% CI) in every test cohort, and random-effects pooled estimates over the
external cohorts (OS: 1.47, 0.92–2.36, heterogeneity P = 0.005; recurrence: 1.25, 1.05–1.49, heterogeneity P = 0.65).

**g**, C-index of the MINER program model vs Cox models on known subtype scores, all trained in the same cohort
(TCGA, or CLCA for the TCGA test): Hoshida classes, the Chiang proliferation signature, and all known scores
(Hoshida S1–S3, proliferation, CTNNB1, KRT19, stem cell, global mean). Adding the MINER score to all known scores:
likelihood-ratio P ≥ 0.10 in every cohort. The model performs comparably to published signatures but doesn't
improve on them.

**h**, Leave-one-cohort-out stability. The network was rebuilt without each discovery cohort (ComBat, MINER and causal
inference on the remaining two cohorts only; technical regulons removed by gene content). Bars show the fraction of
full-network regulators recovered, the median correlation between each full-network program's activity and its best
leave-out match in the held-out cohort (all programs, and the 25% with the largest risk weights), and the fraction of
high-confidence driver→regulator edges (same direction) recovered among the leave-out MINER-filtered flows (CTNNB1
and all drivers). C-index of a ridge program model built entirely without the test cohort (network and training):
TCGA RFS 0.62; CLCA RFS 0.60, OS 0.67.

## Panel → source

| Panel | Result files |
|---|---|
| 1b | `04_miner/combat/{coexpr,mechinf}`, `05_causal/combat/{highConfidenceCausalResults.csv,regulon_families.tsv}` |
| 1c | `04_miner/combat/subtypes_filtered/`, `07_post/subtype_mapping/subtypes_filtered/ntp_calls_*.tsv`, `03_genomics_clinical/genomic_features.csv` |
| 1d | `07_post/subtype_mapping/subtypes_filtered/program_signature_correlation.tsv` |
| 1e | `05_causal/combat/highConfidenceCausalResults.csv` |
| 1f | `07_post/figures/causal_driver_summary.tsv` |
| 2a | `06_risk/combat/predictor_ridge_programs_TCGA_RFS_h36m/weights.tsv` |
| 2b, 2c | `06_risk/.../predictions.tsv`, NTP calls, `01_harmonized/samples.tsv`, `03_genomics_clinical/survival_*_h36m_miner.csv` |
| 2d | `06_risk/combat/predictor_ridge_programs_{TCGA,CLCA}_RFS_h36m/` |
| 2e–g | `08_validation/<cohort>/{scores,evaluation,head_to_head}.tsv`, `07_post/figures/risk_head_to_head.tsv` |
