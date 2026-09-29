# HCC regulatory network model (MINER3)

A causal and mechanistic regulatory network model of hepatocellular carcinoma built with MINER3 from three cohorts: TCGA-LIHC, CLCA 2024 and LICA-FR.

Data is not in this repository. It stays on the server under `data/`.

## Pipeline

| Step | What | Output |
|---|---|---|
| 01 | Harmonize expression: TPM for all cohorts (TCGA from GDC STAR), gene IDs to Ensembl, primary HCC tumors only, log2(TPM+1), gene filtering | `results/01_harmonized/` |
| 02 | Cross-cohort batch correction (ComBat, with per-cohort z-score as comparison) and QC | `results/02_batch_corrected/` |
| 03 | Build mutation and CNA matrices; harmonize clinical and survival tables | `results/03_genomics_clinical/` |
| 04 | MINER3: `coexpr` → `mechinf` → `subtypes` (regulons, programs, states, member matrices) | `results/04_miner/` |
| 05 | MINER3 causal inference (mutation/CNA → regulator → regulon) | `results/05_causal/` |
| 06 | MINER3 risk prediction and survival (`miner3-riskpredict` / `miner3-survival`, plus the library's predictor) | `results/06_risk/` |
| 07 | Post-analysis: states and programs, hallmark and functional enrichment, subtype concordance | `results/07_post/` |

## Layout

```
config/     parameters and sample/cohort tables (small, versioned)
scripts/    numbered pipeline scripts + SLURM submit wrappers
docs/       method notes
data/       raw cohort data           (not versioned)
results/    all generated outputs     (not versioned)
envs/       project-local conda envs  (not versioned)
```

## Environment

MINER3 runs from the `miner3` conda env (`/users/sturkars/mambaforge/envs/miner3`).

## Project log

Decisions, issues and open questions are tracked in [PROJECT_LOG.md](PROJECT_LOG.md).
