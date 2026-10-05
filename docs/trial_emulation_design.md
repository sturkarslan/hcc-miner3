# In-silico trial emulation: predicted vs published response rates from network DCNA

Status: design (2026-10-05); nothing run. The same workflow serves the HCC (`sturkarslan/hcc-miner3`) and ICC
(`sturkarslan/icc-miner3`) networks. This file is kept identical in both repos; disease-specific inputs are in each
repo's `config/trials.yaml` and `docs/soc_response_validation_design.md`.

## Idea

For each published trial arm:
1. Read its eligibility, size, baseline characteristics and ORR.
2. Draw 1,000 synthetic cohorts of the same size and make-up from our tumours.
3. Call each tumour responder / non-responder from the drug's **DCNA** (drug-constrained network activity).
4. Compare the distribution of predicted ORR with the published ORR.

## 1. DCNA

DCNA definition (Baliga lab, multiple myeloma; [PMC12081869](https://pmc.ncbi.nlm.nih.gov/articles/PMC12081869/)):
the mean activity of all regulons that contain, or are regulated by, the drug's target(s), in one patient. For an
antagonist / inhibitor, high DCNA means likely sensitive; for an agonist, low DCNA.

Choices for our networks:
- **Targets** come from ChEMBL mechanisms (target gene, action type), fetched by a script and never typed by hand.
  DrugBank is the fallback.
- **"Regulated by"** has two sources:
  - MINER regulator → regulon edges, when the target is a regulator;
  - step-05 causal flows, when the target is a mutated or fused driver. Example: FGFR2 fusion → the regulators of
    ICC programs 49, 51, 52, 54.
- **Activity** = regulon activity used everywhere else (mean z of regulon genes within cohort), so external cohorts
  are scored the same way.
- **Combinations: independent drug action.** A tumour responds if it responds to at least one component. This is
  the null model under which most combination benefit in trials is explained (Palmer & Sorger, Cell 2017
  *(verify citation)*). It needs no extra parameters, and deviations from it are themselves a result.
- **Drugs without a protein target in tumour cells** need a stated proxy, flagged as such:
  - cisplatin and oxaliplatin (DNA);
  - FOLFOX components;
  - checkpoint antibodies, whose targets sit mostly on immune cells: CD274 / PDCD1 / CTLA4 regulons pick up
    immune-infiltrate programs. That is the intended meaning, but it should be stated.

## 2. Trial table (`config/trials.yaml`)

One entry per **arm**:
- NCT ID; setting; drugs;
- n;
- ORR, with criteria (RECIST 1.1 / mRECIST) and reviewer (independent / investigator);
- DCR; median PFS / OS;
- eligibility: disease site, stage, biomarker selection (FGFR2 fusion, IDH1, HER2, AFP ≥ 400, …), prior lines,
  Child–Pugh;
- baseline characteristics: sex, age, region / Asia, etiology, site (intrahepatic / extrahepatic / gallbladder),
  macrovascular invasion, extrahepatic spread.

Each value carries a source and `verified: true / false`.

Retrieval (`scripts/tools/fetch_trials.py`, run on the login node):
1. ClinicalTrials.gov API v2 (`/api/v2/studies/{NCT}`): eligibility text, enrolment, arms, and, for trials with
   posted results, baseline-characteristic and outcome-measure tables including ORR.
2. Fallback: the publication (Table 1 and response table) via the PMC open-data bucket. The same route was used for
   the classifier supplements.

Eligibility text is mapped to filters by hand once per trial and recorded in the YAML. Free text is too ambiguous
to parse automatically.

## 3. Synthetic cohorts

For each arm, b = 1 … 1,000:
1. **Eligible pool.** Tumours that pass the trial's filters available in our data: biomarker status from step 03b
   (FGFR2 fusion, IDH1, ERBB2 amplification, …), site, stage. A filter we cannot evaluate (ECOG, Child–Pugh, prior
   lines) is listed as unmatched; it is never silently ignored.
2. **Matching weights.** Raking (iterative proportional fitting) of pool weights so that the weighted marginals match
   the trial's baseline: sex, age group, region, etiology, site, stage / vascular invasion. Report effective sample
   size and the largest weight; flag arms where the effective sample size is below ~3 × the number of distinct
   tumours drawn.
3. **Draw** n tumours with replacement, probability ∝ weight.
4. **Predicted ORR_b** = fraction of the drawn tumours called responders (section 4).

The 1,000 ORR_b values give the prediction interval. This combines trial-size sampling noise (≈ binomial) and
composition uncertainty. It does not include calibration uncertainty, which is added by bootstrapping the threshold
fit (section 4).

**Pool size limits what can be emulated.** Biomarker-selected trials draw from tiny pools: e.g. FGFR2 fusion 28
(FU-iCCA) + 7 (OEP002768), ERBB2 amplification ~11–13. A synthetic cohort of 108 drawn from 28 tumours repeats each
tumour about four times. The interval then reflects those 28 tumours, not the population. Report the pool size
next to every result.

## 4. Responder threshold: the one decision that makes or breaks the test

**Problem.** If the DCNA threshold for a trial is set so that the predicted responder fraction equals **that trial's**
ORR, then predicted ORR = observed ORR by construction, and the comparison tests nothing. The threshold must come
from data that are not the trial being predicted.

Options, all expressed as a DCNA quantile or z in a fixed reference population (FU-iCCA, or the HCC discovery set),
so that a threshold carries over between trials:

| | Threshold from | Tests |
|---|---|---|
| A. Leave-one-trial-out (recommended) | the other trials of the same drug class (e.g. pemigatinib → futibatinib; TOPAZ-1 → KEYNOTE-966; IMbrave150 → STRIDE / nivo-ipi), refit with each held out | whether population and eligibility differences between trials explain their different ORRs |
| B. Control-arm anchor | the control arm of the same trial (e.g. GemCis in TOPAZ-1, sorafenib in IMbrave150), then predict the experimental arm with the drug-specific DCNA | whether the network predicts the **added** effect of the new drug |
| C. One global rule | no trial data: e.g. DCNA z > 1 in the reference population, or a cut fitted on a treated cohort with patient-level response (GSE255058 for ICC; GSE109211 / GSE104580 / Zhu 2022 for HCC) | absolute calibration; the strictest test |
| D. Per-trial fit (the original proposal) | the trial itself | nothing about ORR; still useful as a **descriptive** step (which tumours and programs make up the "responders" at the trial's ORR), and for the comparisons in section 5 that do not use ORR |

Recommendation:
- **A** as the primary test, **B** as secondary, **C** where a treated cohort exists.
- **D** only for description, never as evidence.

## 5. Outputs and how to judge them

Per arm:
- predicted ORR median and 95% interval;
- observed ORR with its exact binomial CI;
- coverage (is the observed ORR inside the interval?);
- absolute error.

Across arms:
- Spearman ρ and mean absolute error between predicted and observed ORR;
- for B, predicted vs observed **difference** between arms.

Pre-specified success criteria (fix before running):
1. Under A, observed ORR inside the 95% prediction interval for ≥ 70% of arms, and ρ > 0 with permutation p < 0.05.
2. **Biomarker contrast reproduced.** Without being told the ORR, the network ranks biomarker-selected trials above
   unselected ones for the matching drug. For ICC: FGFR inhibitors in FGFR2-fusion tumours (37–42%) vs the same DCNA
   in all-comers; for HCC: ramucirumab in AFP ≥ 400.

**Negative controls:**
- random target sets of the same size: random-regulon DCNA;
- drugs permuted between arms;
- the step-06 risk score in place of DCNA.

The risk score is prognostic, so it should not reproduce drug-specific differences.

**Known biases** (state in every result):
- Our tumours are mostly resected, earlier-stage and treatment-naive, while trial patients are advanced and often
  pretreated.
- Biliary trials mix intrahepatic, extrahepatic and gallbladder tumours; use the intrahepatic subgroup ORR where
  published, otherwise flag.
- RECIST vs mRECIST, and investigator vs independent review, change ORR by up to 2-fold (lenvatinib: 19% → 41%).
  Compare like with like, and never pool criteria.

## 6. Implementation (step 10, shared code shape)

- `scripts/tools/fetch_drug_targets.py`: ChEMBL mechanisms → `config/drug_targets.tsv` (drug, target gene, action,
  source). Proxies for target-less drugs are listed explicitly.
- `scripts/tools/fetch_trials.py`: ClinicalTrials.gov v2 + manual curation → `config/trials.yaml`.
- `scripts/10a_dcna.py`:
  - DCNA per drug × tumour for the discovery set and every external cohort (portable scoring);
  - regulon sets per drug written out, so each DCNA is auditable.
- `scripts/10b_trial_emulation.py`: eligibility pools, raking, 1,000 draws per arm, thresholds by A / B / C (and D
  for description), the outputs of section 5.
- `scripts/10c_trial_figures.py`:
  - predicted vs observed ORR per arm (interval vs CI);
  - calibration plot;
  - per-arm composition of predicted responders (programs, states, published classes).

Seed 12 for all draws. Results go to `results/10_trials/`.
