# HCC regulatory network model (MINER3)

A causal and mechanistic regulatory network model of hepatocellular carcinoma built with MINER3 from three cohorts: TCGA-LIHC, CLCA 2024 and LICA-FR.

Data is not in this repository. It stays on the server under `data/`.

## Pipeline

| Step | What | Output |
|---|---|---|
| 01 | Harmonize expression: gene IDs to Ensembl, primary HCC tumors only, log2(x+1), gene filtering | `results/01_harmonized/` |
| 02 | Cross-cohort batch correction (ComBat, with per-cohort z-score as comparison) and QC | `results/02_batch_corrected/` |
| 03 | Build mutation and CNA matrices; harmonize clinical and survival tables | `results/03_genomics_clinical/` |
| 04 | MINER3: coexpression → mechanistic inference → bicluster members → subtypes | `results/04_miner/` |
| 05 | MINER3 causal inference (mutation/CNA → regulator → regulon) | `results/05_causal/` |
| 06 | Risk modeling and survival (train TCGA OS; test CLCA RFS and external cohorts) | `results/06_risk/` |
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
