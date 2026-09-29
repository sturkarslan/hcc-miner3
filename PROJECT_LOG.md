# HCC MINER3 Project Log

Running record of decisions, issues and open questions. Newest entries at the top of each section.
Keep entries short: date, what, why.

---

## Open questions

- **[2026-09-29] Steps 01–02 written, not yet run on real data.** `scripts/01_harmonize_expression.py`, `scripts/02_batch_correct.py`, `config/params.yaml`. Tested end to end only on synthetic data (cloud session; no server access). Before the first run:
  - Check the `VERIFY` entries in `config/params.yaml`: CLCA and LICA-FR folder names, and the LICA-FR annotation columns (ID, sample type and its HCC value, Boyault/molecular/immune/etiology labels).
  - Fill `config/tcga_exclude.tsv` with the 3 fibrolamellar TCGA cases.
  - Download the HGNC complete set to `data/reference/hgnc_complete_set.txt`. Without it, symbols map through MINER's `identifier_mappings.txt` only.
  - If pyreadr loses the row names of `tpm_`, the script stops. Export the matrix to TSV from R and set `format: tsv_tpm`; the error message gives the command.
  - Create `envs/hcc-prep` from `config/environment_prep.yml`.
- **[2026-09-29] LICA-FR survival data.** Still not found.
  - `clinical_molecular_annotations.xlsx` has no OS/RFS columns.
  - Jia & Tang 2021 (J Clin Transl Hepatol, PMC9039713) Table S2 (`JCTH-10-273-s002.csv`, provided by user) pools ICGC LIRI-JP (232), LICA-FR (161) and LIHC-US/TCGA (294) by ICGC donor ID (`DO…`). Grouping by ID block and censoring pattern, the block `DO228883…, DO231795–DO231932, DO44634–DO44868, DO50743–DO50748, DO50810–DO50974` looks like LICA-FR. Only 6 of those donors have vital status + time. **Inference, not confirmed**: no ICGC donor-ID → CHC-ID mapping is available (ICGC DCC portal retired).
  - The HCC-subtypes Shiny app (http://51.159.169.17:3838/HCC-subtypes/) is interactive only; its data can't be retrieved with a plain HTTP fetch. Need its underlying data or to ask the authors.
  - Next option: request follow-up from the LICA-FR authors (Zucman-Rossi lab).
- **[2026-09-29] CLCA overall-survival time.** CLCA has `OS_STATUS` but no OS time; only `RFS` (days) + `RFS_STATUS`. Is OS time available from the CLCA portal or supplementary tables?
- **[2026-09-29] GitHub push blocked.** Remote set to `git@github.com:sturkarslan/hcc-miner3.git`. Neither `~/.ssh/id_ed25519` nor `~/.ssh/id_rsa` is accepted by GitHub (`Permission denied (publickey)`), and `gh` is not installed. Fix: add `~/.ssh/id_ed25519.pub` (fingerprint `SHA256:+ft1jzbaIpoPKP0I62/teOsZzXZPrfcYcT8dTSNLtkU`) at github.com → Settings → SSH keys, then `git push -u origin main`.

---

## Decisions

- **[2026-09-29] Step 01 harmonization rules.** Symbol → Ensembl tiers, first hit wins: HGNC approved → MINER Gene Name → HGNC previous → aliases (MINER Synonym + HGNC alias). A symbol matching more than one Ensembl ID in its tier is dropped and logged in `id_mapping_<cohort>.tsv`. Rows mapping to the same Ensembl ID are summed, since TPM is additive. Only genes shared by all cohorts are kept, and TPM is renormalized to 1e6 over them. A gene is kept if TPM ≥ 1 in ≥ 20% of samples in **every** cohort, so cohort-specific detection can't drive the batch effect. Output is log2(TPM+1). For TCGA, only `-01A` samples, one per patient.
- **[2026-09-29] Step 02: ComBat primary, per-cohort z as comparison.** ComBat (inmoose, parametric, cohort as batch, no covariates, because the LICA-FR labels exist in only one cohort). Genes with zero variance in any cohort are dropped. The final matrix is the gene z-score across all samples, clipped below −4 as in MINER's `preProcessTPM`. QC per matrix: silhouette and kNN mixing by cohort on the top 20 PCs, within-cohort structure preservation (Spearman of gene-centered sample-sample correlations), LICA-FR label silhouettes, and marker-program scores (CTNNB1, proliferation, hepatocyte) with per-cohort mean/SD and Kruskal–Wallis H by label.
  - **Caveat, seen in the synthetic test.** Correction removes real differences in subtype composition between cohorts, along with the batch effect. For example, a CTNNB1 program enriched in LICA-FR is pulled to the cross-cohort mean. The within-LICA-FR Boyault association is preserved. The same applies to etiology (CLCA ≈ HBV). Keep this in mind when comparing state frequencies across cohorts.
- **[2026-09-29] Use GDC STAR TPM for TCGA, not the cBioPortal file.** The cBioPortal `data_mrna_seq_v2_rsem.txt` is RSEM upper-quartile-normalized counts, not TPM (column sums ≈ 2.8×10⁷). Downloaded `TCGA-LIHC.star_tpm.tsv.gz` from UCSC Xena GDC hub into `data/LIHC-TCGA-GDC-Xena/`. **Values are log2(TPM + 0.001)**; back-transform before combining with CLCA/LICA-FR TPM. Ensembl IDs are versioned (GENCODE v36). 371 primary tumors (`-01A`), 3 recurrent (`-02`), 50 normals (`-11`).
- **[2026-09-29] Risk modeling uses MINER3's own framework**, not a custom plan (user direction). Details under Environment/tooling.
- **[2026-09-29] Follow the working GBM MINER3 run** in `/proj/omics4tb2/sturkarslan/GBM-15370004/analysis_stringent_sct/`: `miner3-coexpr -mg 6` → `miner3-mechinf -mc 0.1` → `miner3-subtypes`, all with `--skip_tpm` on pre-z-scored input. `miner3-subtypes` writes `coherentMembers.csv`, `overExpressedMembers.csv`, `transcriptional_programs.json`, `transcriptional_states.json`, so no separate bcmembers step is needed.
- **[2026-09-29] All project files stay under `/proj/omics4tb2/sturkarslan/HCC`.** Includes this log (used instead of the Claude memory dir), temp files, cloned references (`external/`), and any new conda env (`envs/`).
- **[2026-09-29] Data never goes to git.** `data/`, `results/`, `envs/`, `external/`, and binary/large formats are in `.gitignore`. Only code, configs and docs are pushed.
- **[2026-09-29] Explicit cross-cohort batch correction is required.** Verified in `miner.py` (see Environment/tooling): MINER3 has no TMM step and no cohort-aware batch correction. Plan: TPM → log2(TPM+1) → gene filtering → cohort correction (ComBat; per-cohort z-scoring as comparison) → QC (PCA/silhouette by cohort; preservation of Boyault G1–G6, CTNNB1 program, proliferation class) → z-scored matrix to MINER with `--skip_tpm`.
- **[2026-09-29] Proposed validation (pending user sign-off).** Build the network on all three cohorts. Validate by (a) leave-one-cohort-out stability of regulons and states, and (b) external cohorts: ICGC LIRI-JP (RNA-seq; survival available in the user-provided Table S2) and GSE14520/LCI (Affymetrix, OS + RFS).

---

## Data inventory (2026-09-29)

| Cohort | Expression | Units / IDs | n expr | Mutations | CNA | Survival |
|---|---|---|---|---|---|---|
| TCGA-LIHC (GDC via Xena) | `LIHC-TCGA-GDC-Xena/TCGA-LIHC.star_tpm.tsv.gz` | log2(TPM+0.001); Ensembl v36 | 371 primary + 50 normal | cBioPortal MAF | cBioPortal GISTIC, log2, seg, arm-level | Xena survival + cBioPortal OS/DSS/DFS/PFS |
| TCGA-LIHC PanCan 2018 (cBioPortal) | `data_mrna_seq_v2_rsem.txt` | RSEM UQ-normalized counts; Hugo + Entrez | 366 | (same as above) | (same) | (same) |
| CLCA 2024 (Chinese, HBV-dominant) | `data_mrna_seq_tpm.txt` | TPM; Hugo | 239 | MAF (WGS, 494 cases) | not in folder | OS status only; RFS days + status |
| LICA-FR | `RNAseq_tpm_457s.RData` (`tpm_`) | TPM; gene symbols | 457 (includes relapse/metastasis; NTL rows in annotation) | `mutation.dat_revised` (3.7M rows, WES/WGS) + driver summary table (529) | `wxs_seg` segments (529) | none (see open question) |

Other notes:
- LICA-FR annotation IDs are `#1010T`; expression/mutation IDs are `CHC1010T`. Map `#` → `CHC`.
- LICA-FR annotation includes 99 NTL (non-tumor liver), 50 relapse, 9 metastases, 11 HCC-on-HCA. Restrict to primary HCC tumors for the network; keep others for later analyses.
- TCGA includes 3 fibrolamellar carcinomas. Exclude them.
- LICA-FR has rich labels (molecular_group, Boyault G1–G6, immune class, etiology) useful as ground truth for checking batch correction and states.
- Etiology is confounded with cohort (CLCA ≈ HBV; LICA-FR mixed alcohol/NASH/HBV/HCV; TCGA mixed). Batch correction can partly remove real etiologic signal. Record this as a caveat and check etiology-specific genes after correction.
- Gene IDs: TCGA is Ensembl; CLCA and LICA-FR are symbols. Map symbols → Ensembl (MINER's `identifier_mappings.txt` + HGNC aliases) and log unmapped/ambiguous genes.

---

## Environment / tooling

- **[2026-09-29] MINER3 preprocessing, verified in `miner.py`.** `preprocess()` = `remove_null_rows` → `correct_batch_effects()` → ID conversion.
  - `correct_batch_effects()` z-scores each gene across all samples. Only when `--skip_tpm` is absent **and** the SD of per-sample means is ≥ 0.15 does it call `preProcessTPM()`.
  - `preProcessTPM()`: keep genes non-zero in ≥ 50% of samples → scale each sample so its positive-value median = 1023 → log2(x+1) → quantile-normalize samples, then genes → gene z-score (clipped below −4).
  - **No TMM, no ComBat, no cohort term.** This is a median-scaling + quantile normalization, so it cannot remove cohort effects when cohorts differ in gene-specific ways (platform, pipeline, gene annotation).
  - `zscore()` returns the input unchanged if the mean of positive values is < 0.1, so pre-z-scored input passes through untouched.
- **[2026-09-29] MINER3 risk prediction.**
  - `miner3-riskpredict` CLI takes a JSON spec (`exp`, `idmap`, `coexpression_dictionary`, `coexpression_modules`, `regulon_modules`, `mechanistic_output`, `regulon_df`, `overexpressed_members`, `underexpressed_members`, `eigengenes`, `filtered_causal_results`, `transcriptional_programs`, `transcriptional_states`, `primary_survival_data`, optional `translocations`). It runs GuanRank on KM estimates + Cox per regulon and per program, KM plots, and activity heatmaps. It needs causal-inference output first.
  - `miner3-survival` is similar and takes a network dir + survival CSV (first two columns: duration, observed) + translocation CSV.
  - The full risk classifier (`generatePredictor`, `riskStratification`, `iAUC`, `predictRisk`; xgboost or decision tree, 20% high-risk cutoff) is in the library, not the CLI. It's used in the paper notebooks (`external/miner3_mattwall/miner/src/MINER_Figure_5.ipynb`). **`xgboost` is not installed in the `miner3` env.**
- **[2026-09-29] Known MINER3 patch.** `miner.py` `decompose()` is patched (empty-array guard) in the shared `miner3` env; backup at `miner.py.backup`. Patch script: `GBM-15370004/analysis_stringent_sct/scripts/network_analysis/patch_miner3.py`.
- **[2026-09-29]** `miner3` conda env: `/users/sturkars/mambaforge/envs/miner3` (isb_miner3 1.2.4). CLIs: `miner3-coexpr`, `-mechinf`, `-subtypes`, `-causalinference`, `-survival`, `-riskpredict`. `miner3-bcmembers` has no entry point, but it isn't needed (see Decisions).
- **[2026-09-29]** No R env has `sva` or `readxl`. Base mambaforge Python has `openpyxl`; the `miner3` env does not. A project-local env is needed for ComBat (R `sva` or Python `inmoose`), under `envs/`.
- **[2026-09-29]** Compute: login node has 16 cores and 15 GB RAM. SLURM partitions `active`/`slow`/`urgent` on nodes baliga[1-5]. Run MINER via sbatch.
- **[2026-09-29]** Reference clone of MattWallScientist/miner3 at `external/miner3_mattwall` (gitignored).

---

## Issues

- **[2026-09-29] Resolved.** `git push` was blocked because GitHub didn't recognize the server's SSH key. User added `~/.ssh/id_ed25519.pub` to GitHub; `main` now pushes to `git@github.com:sturkarslan/hcc-miner3.git`.
