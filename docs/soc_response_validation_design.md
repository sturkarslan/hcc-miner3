# Do network activity patterns recapitulate response to standard-of-care HCC drugs? Validation design

Status: design only (2026-10-02). Nothing below has been run. Numbers marked *(verify)* were not re-checked against
the primary paper in this session; all others were checked against the trial report or a review of it on 2026-10-02.

## 1. Standard of care and the response rates to reproduce

Current guidelines (BCLC 2022, AASLD 2023, ESMO, ASCO 2024, NCCN) agree on this ordering. Response rates are
confirmed objective response rate (ORR); RECIST 1.1 by independent review unless stated. ORR depends strongly on the
criteria (mRECIST gives higher rates) and on the reviewer (investigator vs independent), so every comparison below
must be like for like.

### First line, advanced HCC (BCLC C, or B not suitable for locoregional therapy)

| Regimen | Trial | ORR, experimental | ORR, control (sorafenib or len/sor) | Median OS (months) |
|---|---|---|---|---|
| Atezolizumab + bevacizumab | IMbrave150 (updated) | 30% (mRECIST 35%) | 11% (mRECIST 14%) | 19.2 vs 13.4 *(verify)* |
| Nivolumab + ipilimumab | CheckMate 9DW | 36% (CR 7%) | 13% (lenvatinib or sorafenib) | 23.7 vs 20.6 |
| Camrelizumab + rivoceranib | CARES-310 (mostly Asian) | 25.4% | 5.9% | 22.1 vs 15.2 *(verify)* |
| Durvalumab + tremelimumab (STRIDE) | HIMALAYA | 20.1% | 5.1% | 16.4 vs 13.8 *(verify)* |
| Durvalumab alone | HIMALAYA | 17.0% *(verify)* | 5.1% | 16.6 *(verify)* |
| Tislelizumab | RATIONALE-301 | 14.3% | 5.4% | 15.9 vs 14.1 *(verify)* |
| Lenvatinib | REFLECT | 18.8% (independent RECIST); 24.1% investigator mRECIST; 40.6% independent mRECIST | sorafenib 9.2% (investigator mRECIST) | 13.6 vs 12.3 *(verify)* |
| Sorafenib | SHARP / Asia-Pacific | 2% / 3.3% *(verify)* | placebo | 10.7 vs 7.9 *(verify)* |

Atezolizumab + bevacizumab is the preferred first line. STRIDE is preferred when the bleeding risk is high, and a TKI
(lenvatinib or sorafenib) when immunotherapy is contraindicated.

### Second line (mostly after sorafenib)

| Regimen | Trial | ORR | Notes |
|---|---|---|---|
| Regorafenib | RESORCE | 11% (mRECIST) | sorafenib-tolerant patients |
| Cabozantinib | CELESTIAL | 4% | |
| Ramucirumab | REACH-2 | 5% | AFP ≥ 400 ng/mL only |
| Pembrolizumab | KEYNOTE-240 / -394 | 18.3% / 12.7% *(verify)* | |
| Nivolumab + ipilimumab | CheckMate 040 | 32% *(verify)* | |

### Earlier stages

- BCLC 0–A: resection, ablation or transplant. No adjuvant therapy is established; IMbrave050 (adjuvant
  atezolizumab + bevacizumab) lost its early recurrence-free survival benefit at the updated analysis *(verify)*.
- BCLC B: TACE; objective response by mRECIST about 50% in pooled series *(verify)*. TACE + immunotherapy
  combinations (EMERALD-1, LEAP-012) improve progression-free survival *(verify)*.

### Real-world response rates

- Atezolizumab + bevacizumab: pooled over 47 studies / 5,400 patients, ORR 26.7% (95% CI 24.6–29.1) RECIST and
  34.0% (30.3–37.8) mRECIST, close to the trial.
- Lenvatinib: 19–42% by mRECIST across Japanese and Korean multicentre series.
- STRIDE: 16–52% in small real-world series. Too variable to use as a target.

### Subgroup patterns the model should reproduce

These are the strongest tests, because they do not depend on matching the trial population as a whole.

1. **Etiology.** In a meta-analysis of three PD-(L)1 trials (> 1,600 patients), immunotherapy did not improve survival in
   non-viral HCC. NASH-HCC did worse on anti-PD-(L)1 in two more cohorts (Pfister, Nature 2021).
2. **WNT/β-catenin.** CTNNB1-mutant tumours were refractory to checkpoint inhibitors (Harding, Clin Cancer Res
   2019). The "immune exclusion" (CTNNB1) class was predicted to resist immunotherapy (Montironi 2023).
3. **Pre-existing immunity and oncofetal genes, under atezolizumab + bevacizumab** (Zhu, Nat Med 2022; 358 patients
   from GO30140 and IMbrave150):
   - Better outcomes with high CD274, T-effector signature and CD8 density, and with the 10-gene ABRS.
   - Worse outcomes with a high Treg : T-effector ratio and with GPC3/AFP expression.
   - The bevacizumab contribution is linked to VEGFR2, Tregs and myeloid inflammation.
4. **Anti-PD-1 and interferon.** The IFNAP signature (interferon / antigen presentation) and the Sia immune class
   predict anti-PD-1 response (Haber, Gastroenterology 2023; 83 tumours with expression).
5. **Ramucirumab.** Benefit only when AFP ≥ 400 ng/mL (REACH-2).
6. **Adjuvant sorafenib (STORM).** Tumours that responded to sorafenib were enriched in CD4+ T, B and cytolytic NK
   cells (Pinyol, Gut 2019; GSE109211).

## 2. What "recapitulate" means: three levels, tested in this order

| Level | Question | Data | Strength |
|---|---|---|---|
| L3, patient | In treated patients, does a pre-specified network score separate responders from non-responders, and does it add to the published biomarkers? | treated cohorts with pre-treatment expression and response | strongest; needs treated data |
| L2, subgroup | Does the network predict response differences between clinical or genomic subgroups in the direction and size reported by the trials? | our 929 discovery tumours + external cohorts; trial subgroup ORR/HR | good; no treated data needed |
| L1, population | Weighted to a trial's population, does the predicted share of responders match the observed ORR, and does the ranking of regimens match? | L3 calibration transported to our cohorts | weakest; only valid after L3 calibration |

**Important distinction.** The step-06 risk score is **prognostic**: a high-risk tumour does worse on any therapy.
"Recapitulating drug response" needs a **predictive** signal. That means either:
- a drug × score interaction in randomized data (IMbrave150: atezolizumab + bevacizumab vs sorafenib), or
- scores that rank regimens differently from each other.

The risk score alone is therefore a negative control, not the test.

## 3. Pre-specified drug → network mapping (freeze before any outcome data is loaded)

One drug-class score = mean of the standardized activities of the listed programs and panel signatures, with the
given signs and equal weights. No fitting on response. Program IDs are from `config/program_labels.tsv`.

| Drug class | Up (predicts response) | Down (predicts resistance) | Network-specific hypothesis |
|---|---|---|---|
| PD-(L)1 ± CTLA-4 (atezo-bev, STRIDE, nivo-ipi, tislelizumab, camrelizumab, pembrolizumab) | P11 immune class (T/NK), P60 IFN-γ / IFNAP, P16 IFN-α; ABRS; Sia immune class | WNT/β-catenin P2, P7, P9, P39; oncofetal / progenitor P22, P23 (GPC3/AFP); Treg : Teff ratio | AXIN1-mutant tumours: MHC class I low (B2M, HLA-B, TAP1, PSMB9 down; step 07f) by a different route from CTNNB1 → predicted ICI-resistant, untested |
| Anti-VEGF component (bevacizumab, ramucirumab; TKIs) | P33 endothelial / angiogenesis; VEGFA; myeloid inflammation P50 | — | Ramucirumab: AFP / progenitor-high (P22, P23) |
| Lenvatinib (VEGFR1–3, FGFR1–4) | P33; FGF19 / FGFR4 axis; 11q13 (CCND1/FGF19) amplification flows (step 05: E2F up, 10 adverse programs) | — | 11q13-amplified tumours predicted lenvatinib-favoured |
| Sorafenib (RAF, VEGFR, PDGFR) | MAPK activity; VEGFA (6p21) amplification; immune P11 (STORM responders) | — | |
| Cabozantinib / regorafenib (MET, AXL, VEGFR2) | MET / HGF activity, P33 | — | |
| TACE (BCLC B) | — | hypoxia / glycolysis P75, proliferation P55 | |

**Gaps to fill before freezing:**
- Zhu 2022 T-effector, Treg, myeloid-inflammation and angiogenesis signatures are in the paper's Methods, not the
  supplement. They are still missing from `subtype_signatures_custom.tsv` (open item in PROJECT_LOG).
- A MAPK activity set (e.g. PROGENy-style), FGF19/FGFR4 targets and MET targets need adding to the panel.

The mapping goes into `config/drug_network_map.yaml`. Its commit hash is recorded in PROJECT_LOG **before** any
treated cohort is scored. This is the protection against tuning the mapping to the outcome.

## 4. Treated cohorts (patient level, L3)

| Cohort | Treatment, setting | n | Response label | Access | Use |
|---|---|---|---|---|---|
| GO30140 + IMbrave150 (Zhu 2022) | atezo-bev, atezo alone, sorafenib; 1L advanced | 358 | RECIST response, PFS, OS | controlled (EGA via Roche data access committee); application needed | primary ICI test; **predictive interaction** vs sorafenib arm |
| Haber 2023 | anti-PD-1, advanced | 83 expression / 72 mutations | response, OS | check (13 centres; may be on request) | ICI replication; mutation subgroups |
| Samsung pembrolizumab phase 2 (Genome Med 2021) | pembrolizumab after sorafenib | 60 | response | check | ICI replication, 2L |
| GSE109211 (STORM, Pinyol 2019) | adjuvant sorafenib vs placebo | 140 (67 sorafenib: 21 responders / 46 non-responders) | RFS-based responder | public (GEO) | sorafenib; **placebo arm = interaction control** |
| GSE104580 | TACE | 147 (81 responders / 66 non-responders) | mRECIST within 3 months | public (GEO) | TACE |
| CheckMate 040 (Sangro 2020) | nivolumab | — | response | not public | — |

Scoring uses the step-08 portable method unchanged:
- regulon score = mean z of its genes within the cohort;
- program activity = mean regulon score;
- no refitting.

Report network gene coverage per cohort; it was 75–99% in the step-08 cohorts. Array cohorts (GSE109211,
GSE104580) are loaded with the generic GEO loader; one probe per symbol, highest mean.

**Patient-level analyses:**
1. Each drug-class score vs response:
   - Discrimination: AUC with bootstrap CI, and odds ratio per SD.
   - Calibration, after a logistic link fitted in a calibration cohort: calibration slope and intercept, and observed
     vs expected responders.
2. Head-to-head against published biomarkers: ABRS, IFNAP, Sia immune class, CD274, Montironi classes, plus the
   step-06 risk score. Likelihood-ratio test of the network score added to them, as in step 08. Prognosis was
   "comparable, not better", so expect the same unless the causal layer adds something.
3. Predictive interaction:
   - IMbrave150 atezo-bev vs sorafenib: Cox for OS / PFS and logistic for response, with treatment × score.
   - GSE109211: sorafenib vs placebo, recurrence-free survival × score.
4. Train / test across cohorts if any refitting is wanted: fit on GO30140, test on IMbrave150 and Haber; never in
   the same cohort.

## 5. Subgroup structure (L2): no treated data needed

Compute each drug-class score in the 929 discovery tumours, and in external cohorts where the covariate exists.
Compare subgroup means with the direction of the trial evidence:

| Contrast | Expected (from trials / literature) | Our data |
|---|---|---|
| Viral (HBV / HCV) vs non-viral / NASH, ICI score | viral higher | etiology in TCGA, LICA-FR (CLCA mostly HBV), LIRI-JP |
| CTNNB1-mutant vs wild type, ICI score | mutant lower | 03b genomic features (all three discovery cohorts) |
| AXIN1-mutant vs wild type, ICI score | **prediction: lower, via MHC-I** (new) | 03b; test in Haber / Zhu mutation data |
| AFP-high vs low, ramucirumab / anti-VEGF score | higher | AFP in TCGA, CLCA *(check fields)* |
| Asian vs non-Asian, ICI score | small or no difference | cohort / ancestry |
| 11q13 amplification, lenvatinib score | higher (hypothesis) | TCGA copy number |

Outcomes:
- Direction concordance across contrasts (sign test).
- Effect sizes, with within-cohort standardization.
- The causal layer (step 05) supplies the mechanism for the genomic contrasts: a driver → regulator → program path
  for each predicted difference. This is what the network adds over a signature lookup.

## 6. Population transport (L1): predicted vs observed ORR per regimen

(2026-10-05: superseded in mechanics by the trial-emulation workflow, `docs/trial_emulation_design.md`. There, each
trial arm is emulated by 1,000 synthetic cohorts matched to its size, eligibility and baseline, and responders are
called from drug DCNA with a threshold calibrated on other trials. Arms are listed in `config/trials.yaml`, which
replaces the `soc_regimens.yaml` planned in section 8. The weighting idea below is kept as the population-matching
step.)

Only for drug classes with an L3 calibration (realistically ICI from Zhu / Haber, sorafenib from STORM, TACE from
GSE104580). Steps:

1. **Calibrate.** Fit a logistic model, P(response) = f(score), in the calibration cohort (e.g. the atezo-bev arm of
   GO30140 + IMbrave150).
2. **Target population.** Weight our tumours (or, better, the advanced-stage external cohorts) to each trial's
   baseline by inverse-odds weighting:
   - etiology (HBV / HCV / non-viral);
   - Asia;
   - BCLC C or macrovascular invasion;
   - AFP ≥ 400;
   - ECOG and Child–Pugh A, which are not available and are assumed.
3. **Predict.** Predicted ORR = weighted mean of P(response), with a bootstrap CI that propagates calibration
   uncertainty. Compare with the observed trial ORR (exact binomial CI). Do the same for the real-world pooled
   estimate.
4. **Across regimens.** Rank correlation of predicted vs observed ORR; only about 6–8 points, so report it but don't
   lean on it. The sorafenib control arms (5–13% depending on criteria) are a built-in low anchor.

**Population caveat.**
- Most discovery tumours are resected and earlier-stage: TCGA stage I–II is 74%; CLCA is the exception, with 57%
  BCLC C and 57% microvascular invasion.
- Advanced and metastatic tumours, and biopsies, differ in immune / stromal content.
- So: primary L1 target population = the treated cohorts' own baselines. Sensitivity analyses: CLCA BCLC C /
  microvascular invasion, and TCGA stage III–IV.

## 7. Pre-specified success criteria (write into PROJECT_LOG before running)

1. **L3:**
   - ICI score AUC ≥ 0.65, with the 95% CI excluding 0.5, in at least one independent ICI cohort.
   - Treatment × score interaction p < 0.05 in IMbrave150 (atezo-bev vs sorafenib).
   - No such interaction for the step-06 risk score (negative control).
2. **L2:** direction concordant for ≥ 4 of the 5 literature-backed contrasts. The AXIN1 contrast is reported
   separately as a new prediction.
3. **L1:** predicted ORR within the trial's 95% CI for the calibrated regimens. Report the absolute error.
4. **Added value:** likelihood-ratio p < 0.05 for the network score over ABRS + IFNAP + Sia immune class, in a cohort
   not used for any fitting.

**Negative controls:**
- the step-06 risk score;
- random program sets of the same size, with the same sign mix;
- the drug → program mapping permuted across drug classes;
- the sorafenib or placebo arms for the ICI score.

## 8. Implementation (new step 10)

- `config/soc_regimens.yaml`: one entry per regimen and trial arm with:
  - setting, n;
  - ORR (criteria, reviewer), DCR, median OS / PFS;
  - population descriptors (etiology %, Asia %, BCLC C / macrovascular invasion / extrahepatic spread %, AFP ≥ 400 %);
  - source.

  Plus the real-world pooled estimates. The table in section 1 is the seed; fill the *(verify)* cells from the
  papers.
- `config/drug_network_map.yaml`: section 3, frozen.
- `scripts/10a_score_treated_cohorts.py`: step-08 loaders and portable scoring for GSE109211, GSE104580, then the
  controlled cohorts when approved. Writes program activities and panel scores.
- `scripts/10b_response_validation.py`: L3 analyses (AUC, OR, calibration, head-to-head, interaction).
- `scripts/10c_subgroups_transport.py`: L2 contrasts with causal-path annotation; L1 weighting and predicted vs
  observed ORR.
- **Figure:** predicted vs observed ORR per regimen (forest style, trial CI and prediction CI side by side); subgroup
  concordance panel.

**Order of work:**
1. Fill the signature gaps (Zhu Teff / Treg / myeloid / angiogenesis; MAPK, FGF19 / FGFR4, MET sets) and freeze the
   map.
2. Run the public cohorts (GSE109211, GSE104580) and the L2 contrasts in the discovery data.
3. Apply for Zhu 2022 (EGA) and ask Haber / Samsung for data.
4. L3 ICI analyses and L1 transport once those arrive.

## Sources

- Guidelines: [ASCO 2024 guideline summary](https://www.cancer.fr/professionnels-de-sante/veille/nota-bene-cancer/bulletin-n-469/systemic-therapy-for-advanced-hepatocellular-carcinoma-asco-guideline), [guideline comparison](https://www.guidelinecentral.com/insights/hcc-side-by-side/)
- Trials:
  - IMbrave150 response: [Kudo et al., PMC8366100](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8366100/)
  - CheckMate 9DW: [CancerNetwork 2025](https://www.cancernetwork.com/view/efficacy-data-support-nivolumab-ipilimumab-as-first-line-hcc-treatment), [ASCO Post](https://ascopost.com/news/april-2025/fda-approves-nivolumab-plus-ipilimumab-for-unresectable-or-metastatic-hcc/)
  - CARES-310: [ESMO](https://esmo.org/oncology-news/camrelizumab-plus-rivoceranib-improves-survival-compared-with-sorafenib-in-first-line-treatment-for-patients-with-unresectable-hepatocellular-carcinoma)
  - REFLECT: [OncLive](https://onclive.com/view/independent-review-backs-key-reflect-findings-with-lenvatinib-in-hcc)
  - RATIONALE-301: [Novartis](https://www.novartis.com/news/media-releases/novartis-announces-tislelizumab-demonstrated-efficacy-and-tolerability-first-line-advanced-liver-cancer-phase-iii-trial)
  - Second line: [PMC9747269 table](https://pmc.ncbi.nlm.nih.gov/articles/PMC9747269/table/T9)
- Real world:
  - atezo-bev systematic review and Japanese cohort: [PubMed 37680945](https://pubmed.ncbi.nlm.nih.gov/37680945/), [Anticancer Res 2022](https://ar.iiarjournals.org/content/42/11/5465)
  - lenvatinib: [Liver Cancer 2021](https://doaj.org/article/b281b7e657a440fba18cf97545a1e8bc)
  - STRIDE: [Karger 2024](https://karger.com/article/doi/10.1159/000542517)
- Biomarkers:
  - Zhu 2022 ([Nat Med 28:1599](https://scholars.lib.ntu.edu.tw/handle/123456789/617262))
  - Haber 2023 ([Gastroenterology 164:72](https://repositori.upf.edu/items/7070966c-9a60-43e2-bd7b-1f163834a33d))
  - Pfister 2021 ([Nature](https://www.doi.org/10.1038/s41586-021-03362-0))
  - Harding / Wnt exclusion ([PMC6445700](https://pmc.ncbi.nlm.nih.gov/articles/PMC6445700))
  - Samsung pembrolizumab ([PMC8734300](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8734300/))
- Treated cohorts: [GSE109211](https://omicsdi.org/dataset/geo/GSE109211), [GSE104580](https://omicsdi.org/dataset/geo/GSE104580)
