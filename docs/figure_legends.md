# Figure explanations and legends

Figures: `results/09_figures/figure1.{pdf,png}`, `figure2.{pdf,png}` (183 mm, Nature double column; Liberation Sans
embedded as TrueType). Source data: `figure1_source_data.xlsx`, `figure2_source_data.xlsx`. Built by
`scripts/09_publication_figures.py` from saved results only; program names from `config/program_labels.tsv`
(each with its evidence).

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
acts through the Wnt effectors LEF1 and TCF7 to switch on the β-catenin and low-recurrence programs. We then asked
whether a driver's causal effects on the programs anticipate its association with recurrence, testing this strictly
out of sample: program risk weights learned in one cohort were compared with driver outcome in the other. With weights
from CLCA and outcome in TCGA the association is positive but modest (ρ = 0.34), and it weakens when correlated
copy-number arms are counted once (ρ = 0.16). The reverse direction has only ten testable drivers and shows no
association. The signal that exists sits in the programs each driver causally reaches (ρ = 0.32) rather than in the
programs it does not reach (ρ = −0.01). We therefore treat this as suggestive, not as evidence that the causal layer
predicts driver prognosis. An earlier version of this panel reported ρ = 0.64; that estimate reused the training
cohort's outcome and was inflated.

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
biology blocks. Tracks: mean risk score per state; fraction of tumours in published HCC classes by nearest-template
prediction, labelled by biology (class codes in grey: Hoshida S1–S3, Boyault G1–G6, Chiang (C)); mutation frequencies.
**d**, Pearson correlation across tumours between program activity and published signature scores, with each tumour's
mean expression regressed out of both. Program names are coloured by the sign of their risk-model weight (red adverse,
blue protective). H, MSigDB hallmark.
**e**, Causal flows from CTNNB1 mutation (n = 222 tumours) to regulators and programs: the regulators with the largest
effects (one per regulon family) plus LEF1 and TCF7. Edges: red, up in mutant tumours; blue, down; dashed, the
regulator represses its regulon; width proportional to |Cohen's d|. Flows shown are high confidence: Benjamini–Hochberg
q ≤ 0.1 for the regulon and driver–regulator tests, same direction in every cohort tested, |d| ≥ 0.5 and ≥ 30
altered tumours.
**f**, Out-of-sample test of whether a driver's causal effects on programs track its association with recurrence.
Risk push: mean over the driver's high-confidence regulon families of Cohen's d times the risk weight of the family's
program, with weights from the ridge model trained in one cohort. y axis: Cox z for recurrence (36 months) of the
driver in the other cohort (≥ 10 altered and ≥ 10 wild-type tumours). Left, CLCA weights and TCGA outcome (71 drivers
in 33 clusters of overlapping alterations): Spearman ρ = 0.34 (cluster-bootstrap 95% interval 0.08–0.53), P = 0.041
by permutation of program weights, P = 0.070 against random regulon families; ρ = 0.16 with one driver per cluster.
Right, TCGA weights and CLCA outcome (10 drivers, 8 clusters; CLCA has no copy-number data): ρ = 0.14 (−0.75 to 0.78),
permutation P = 0.74. Point size, number of altered tumours in the test cohort. Predictor and endpoint were fixed
before the analysis.

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
**b**, Risk score (within-cohort z) by published HCC class, labelled by biology. Hoshida and Chiang classes are
nearest-template-prediction calls (all 929 tumours); Boyault classes are the LICA-FR authors' labels (n = 324).
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
Hoshida classes; Chiang proliferation signature; all known scores (Hoshida S1–S3, proliferation, CTNNB1, KRT19, stem
cell and tumour mean expression). Adding the network score to all known scores: likelihood-ratio P ≥ 0.10 in every
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
| 1f | `07_post/driver_push/{driver_table.tsv,summary.tsv}` (step 07e) |
| 2a | `06_risk/combat/predictor_ridge_programs_TCGA_RFS_h36m/weights.tsv` |
| 2b, 2c | `06_risk/.../predictions.tsv`, NTP calls, `01_harmonized/samples.tsv`, `03_genomics_clinical/survival_*_h36m_miner.csv` |
| 2d | `06_risk/combat/predictor_ridge_programs_{TCGA,CLCA}_RFS_h36m/` |
| 2e–g | `08_validation/<cohort>/{scores,evaluation,head_to_head}.tsv`, `07_post/figures/risk_head_to_head.tsv` |
| 2h | `08_validation/loco/{loco_summary.tsv,loco_programs.tsv}` |
