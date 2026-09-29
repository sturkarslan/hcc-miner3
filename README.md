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
| 08 | Validation: leave-one-cohort-out stability of regulons and states; external cohorts ICGC LIRI-JP (RNA-seq) and GSE14520/LCI (Affymetrix) | `results/08_validation/` |

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
Steps 01–02 run from a project-local env with ComBat (`inmoose`):

```
mamba env create -p envs/hcc-prep -f config/environment_prep.yml
```

## Running

All parameters are in `config/params.yaml`. Entries marked `VERIFY` are column or file names that have not
been checked against the real files yet; scripts stop and list the available columns when one is wrong.

```
mkdir -p logs
sbatch scripts/slurm/01_harmonize_expression.sbatch
sbatch scripts/slurm/02_batch_correct.sbatch
sbatch scripts/slurm/03_clinical_survival.sbatch
```

Step 02 writes `results/02_batch_corrected/expression_combat_z.csv` (and `expression_cohort_z_z.csv` for
comparison), genes × samples with Ensembl IDs, already z-scored for `miner3-coexpr --skip_tpm`.
Check the QC tables and figures before step 04.

### QC figures

Each step writes PNGs to `results/<step>/qc/`:

| Figure | Shows | Look for |
|---|---|---|
| `h1_sample_gene_flow` | samples and genes kept at each harmonization step | unexpected losses |
| `h2_id_mapping` | how identifiers reached Ensembl (tier), ambiguous and unmapped | large unmapped share |
| `h3_per_sample` | per-sample median/IQR of log2(TPM+1), genes detected, % TPM in shared genes | outlier samples, cohort-wide shifts |
| `h4_gene_detection` | per-gene detection fraction by cohort, with the filter threshold | cohort-specific detection |
| `h5_gene_means_between_cohorts` | per-gene mean expression, cohort vs cohort | gene-specific cohort effects (why batch correction is needed) |
| `h6_sample_correlation` | sample–sample correlation, uncorrected | cohort blocks |
| `b1_pca_by_cohort` | PC1/PC2 by cohort for each matrix, with silhouette and kNN mixing | cohorts overlapping after correction |
| `b2_pca_by_<label>` | LICA-FR labels, one level highlighted per panel | label structure kept after correction |
| `b3_pc_association` | variance of PC1–10 explained by cohort and by labels | cohort η² → 0, label η² kept |
| `b4_gene_cohort_variance` | per-gene variance explained by cohort (ECDF) | shift to ~0 after correction |
| `b5_qc_metrics` | all summary metrics side by side | see figure subtitle |
| `b6_programs_by_cohort`, `b7_programs_by_<label>` | CTNNB1, proliferation, hepatocyte program scores | label differences kept; cohort differences removed (see composition caveat in the log) |
| `b8_rle` | relative log expression per sample, before vs after ComBat | medians centred on 0 |
| `b9_sample_correlation_corrected` | sample–sample correlation after correction | cohort blocks gone |
| `s1_kaplan_meier` | OS and RFS by cohort with numbers at risk | plausible event rates and follow-up |

## Project log

Decisions, issues and open questions are tracked in [PROJECT_LOG.md](PROJECT_LOG.md).
