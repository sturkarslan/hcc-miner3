# HCC MINER3 Project Log

Running record of decisions, issues and open questions. Newest entries at the top of each section.
Keep entries short: date, what, why.

---

## Open questions

- **[2026-09-29] LICA-FR survival data.** `clinical_molecular_annotations.xlsx` has no OS/RFS columns. Is follow-up data available from the LICA-FR authors (Zucman-Rossi lab)? Without it, LICA-FR contributes to network inference and subtype/state analysis only, not to survival modeling.
- **[2026-09-29] CLCA overall-survival time.** CLCA has `OS_STATUS` but no OS time; only `RFS` (days) + `RFS_STATUS`. Is OS time available from the CLCA portal or supplementary tables?
- **[2026-09-29] GitHub remote.** `gh` CLI is not installed on the server. Need a repo URL (user creates it on github.com) or a token to create one.

---

## Decisions

- **[2026-09-29] All project files stay under `/proj/omics4tb2/sturkarslan/HCC`.** Includes this log (used instead of the Claude memory dir), temp files, and any new conda env (`envs/`).
- **[2026-09-29] Data never goes to git.** `data/`, `results/`, `envs/`, and binary/large formats are in `.gitignore`. Only code, configs and docs are pushed.
- **[2026-09-29] Explicit cross-cohort batch correction is required.** MINER3's `correct_batch_effects()` only z-scores each gene across all samples (plus optional TPM quantile normalization when sample means vary). It does not model cohort, so it won't remove TCGA/CLCA/LICA-FR platform effects. Plan: log2(x+1) → gene filtering → cohort correction (ComBat, with per-cohort z-scoring as a comparison) → evaluate with PCA/silhouette by cohort and preservation of known biology (Boyault G1–G6, CTNNB1 program, proliferation class) → z-scored matrix to MINER with `--skip_tpm`.
- **[2026-09-29] Proposed validation (pending user sign-off).** Build the network on all three cohorts (maximizes power for regulon discovery). Validate by (a) leave-one-cohort-out stability of regulons and states, and (b) external cohorts not used in training: ICGC LIRI-JP (RNA-seq, OS) and GSE14520/LCI (Affymetrix, OS + RFS). Risk model trained on TCGA OS; tested on CLCA RFS and external sets.

---

## Data inventory (2026-09-29)

| Cohort | Expression | Units / IDs | n expr | Mutations | CNA | Survival |
|---|---|---|---|---|---|---|
| TCGA-LIHC PanCan 2018 | `data_mrna_seq_v2_rsem.txt` | RSEM normalized counts; Hugo + Entrez | 366 tumors (+ 50 normals in `normals/`) | MAF (`data_mutations.txt`) | GISTIC `data_cna.txt`, log2, seg, arm-level | OS, DSS, DFS, PFS (months) |
| CLCA 2024 (Chinese, HBV-dominant) | `data_mrna_seq_tpm.txt` | TPM; Hugo | 239 | MAF (WGS, 494 cases) | not in folder | OS status only; RFS days + status |
| LICA-FR | `RNAseq_tpm_457s.RData` (`tpm_`) | TPM; gene symbols | 457 (includes relapse/metastasis; NTL rows in annotation) | `mutation.dat_revised` (3.7M rows, WES/WGS) + driver summary table (529) | `wxs_seg` segments (529) | none found |

Other notes:
- LICA-FR annotation IDs are `#1010T`; expression/mutation IDs are `CHC1010T`. Map `#` → `CHC`.
- LICA-FR annotation includes 99 NTL (non-tumor liver), 50 relapse, 9 metastases, 11 HCC-on-HCA. Restrict to primary HCC tumors for the network; keep others for later analyses.
- TCGA includes 3 fibrolamellar carcinomas. Exclude them.
- LICA-FR has rich labels (molecular_group, Boyault G1–G6, immune class, etiology) useful as ground truth for checking batch correction and states.
- Etiology is confounded with cohort (CLCA ≈ HBV; LICA-FR mixed alcohol/NASH/HBV/HCV; TCGA mixed). Batch correction can partly remove real etiologic signal. Record this as a caveat and check etiology-specific genes after correction.

---

## Environment / tooling

- **[2026-09-29]** `miner3` conda env: `/users/sturkars/mambaforge/envs/miner3` (isb_miner3 1.2.4). CLIs available: `miner3-coexpr`, `-mechinf`, `-subtypes`, `-causalinference`, `-survival`, `-riskpredict`. **`miner3-bcmembers` has no CLI entry point**, although `miner/bcmembers.py` exists. Need a wrapper (`python -m` or a small script) to produce `coherentMembers.csv` for causal inference.
- **[2026-09-29]** No R env has `sva` or `readxl`. Base mambaforge Python has `openpyxl`; the `miner3` env does not. A project-local env is needed for ComBat (R `sva` or Python `inmoose`), under `envs/`.
- **[2026-09-29]** Compute: login node has 16 cores and 15 GB RAM. SLURM partitions `active`/`slow`/`urgent` on nodes baliga[1-5]. Run MINER via sbatch.
- Reference pipeline: `/proj/omics4tb2/sturkarslan/GBM-15370004/analysis/scripts/network_analysis/02_run_miner3.sh` (coexpr `-mg 6` → mechinf `-mc 0.1` → bcmembers → subtypes, all `--skip_tpm`).

---

## Issues

_(none yet)_
