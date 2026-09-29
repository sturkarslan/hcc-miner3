# HCC MINER3 Project Log

Running record of decisions, issues and open questions. Newest entries at the top of each section.
Keep entries short: date, what, why.

---

## Open questions

- **[2026-09-29] Steps 01–02 run on the server (final: SLURM 14944, 14945); step 03 waits for the CLCA supplementary table.** Pre-run checklist from the cloud session, done on the server:
  - VERIFY entries resolved against the real files: CLCA folder is `data/HCC-CLCA-2024`; LICA-FR annotation columns are `Sample`, `Sample_type`, `primary_diagnosis`, `Per_patient_analysis`, `Boyault_transcriptomic_group`, `molecular_group`, `immune_class_RNAseq_Fluidigm`, `Etiology`. MINER idmap points to the `miner3` env copy.
  - `config/tcga_exclude.tsv` has the 3 fibrolamellar cases (TCGA-DD-A4NB, TCGA-MR-A8JO, TCGA-RC-A6M5).
  - HGNC complete set downloaded to `data/reference/hgnc_complete_set.txt`.
  - `envs/hcc-prep` created. The conda and pip caches point into `envs/`, so nothing is written outside the project.
- **[2026-09-29] LICA-FR survival data.** Still not found.
  - `clinical_molecular_annotations.xlsx` has no OS/RFS columns.
  - Jia & Tang 2021 (J Clin Transl Hepatol, PMC9039713) Table S2 (`JCTH-10-273-s002.csv`, provided by user) pools ICGC LIRI-JP (232), LICA-FR (161) and LIHC-US/TCGA (294) by ICGC donor ID (`DO…`). Grouping by ID block and censoring pattern, the block `DO228883…, DO231795–DO231932, DO44634–DO44868, DO50743–DO50748, DO50810–DO50974` looks like LICA-FR. Only 6 of those donors have vital status + time. **Inference, not confirmed**: no ICGC donor-ID → CHC-ID mapping is available (ICGC DCC portal retired).
  - The HCC-subtypes Shiny app (http://51.159.169.17:3838/HCC-subtypes/) is interactive only; its data can't be retrieved with a plain HTTP fetch. Need its underlying data or to ask the authors.
  - Next option: request follow-up from the LICA-FR authors (Zucman-Rossi lab).
- **[2026-09-29] CLCA survival: partly resolved, censoring times assumed.** The user provided the Nature 2024 Supplementary Table 1 (`41586_2024_7054_MOESM3_ESM.csv`; **not yet on the server**: copy to `data/HCC-CLCA-2024/`). Table 1a, 494 patients: operation date (day), recurrence and death status, recurrence and death dates (**month only**), and RFS days for recurrences.
  - **There is no last-follow-up date.** Event-free patients are censored at an assumed cutoff of `2021-09-30`, because the latest recorded recurrence and death are both in 2021-09. The planned cross-check against cBioPortal can't work: the cBioPortal CLCA `RFS` column has times only for the 169 recurred patients, so it gives no censoring times. The cutoff stays an assumption unless the CLCA authors give a last-follow-up date.
  - OS for deaths = death month (mid-month) − operation date, so ±15 days. RFS for recurrences uses the table's day counts.
  - Usable: OS 384 patients / 74 deaths. 4 deaths have no date and are dropped; 106 have unknown status. RFS 386 / 169 recurrences; 108 unknown. 1 patient has a death month before the recurrence month (rounding), so OS is set to RFS.
  - Still to check on the server: that CLCA expression sample IDs match `CLCA_####` (step 03 logs the overlap).

---

## Decisions

- **[2026-09-29] Step 02 result: ComBat is the MINER input** (`results/02_batch_corrected/expression_combat_z.csv`, 13,866 × 929). Before correction, CLCA separates on PC2 (TCGA–CLCA gene means r = 0.81; 38% of genes differ > 2-fold). After ComBat, cohort silhouette 0.18 → 0.00 and kNN mixing 0.04 → 0.73, with within-cohort sample structure preserved (ρ 0.999 TCGA, 0.993 CLCA, 1.000 LICA-FR; per-cohort z gives 0.95–0.98). LICA-FR label signal is unchanged: Kruskal–Wallis H for CTNNB1/proliferation/hepatocyte programs by Boyault class 149/143/124 (uncorrected 149/143/129). The CTNNB1 program is highest in G5/G6 and proliferation in G1–G3, as expected.
  - Composition caveat confirmed on real data: before correction, the CTNNB1 program mean is higher in LICA-FR (+0.15) than TCGA (−0.15), consistent with more CTNNB1-mutant tumors in European cohorts. ComBat sets all cohort means to ~0.
- **[2026-09-29] Bug fixed: TCGA Xena values are log2(TPM + 1), not log2(TPM + 0.001).** My earlier note was wrong. Evidence: zero-TPM entries are stored as 0.0, and back-transforming with 0.001 gave column sums of 1.06×10⁶ (10⁶ + ~1 per gene × 60,616 genes). The first step 01 run (14942) had every shared gene "detected" in all TCGA samples. `pseudocount: 1` now, and TCGA column sums are 1.00×10⁶.
- **[2026-09-29] Symbol mapping uses GENCODE v36 names first** (Xena probemap, the release TCGA uses), then HGNC. HGNC's current Ensembl IDs differ from v36 for some genes (e.g. SOD2 is ENSG00000291237 in HGNC, ENSG00000112096 in v36), and those genes were silently dropped from the shared set. Names duplicated in v36 fall through to HGNC. LICA-FR rows mapped went from 71% to 96%; shared genes went from 22,224 to 22,374.
- **[2026-09-29] CLCA TPM contains non-polyA small RNAs.** RPPH1 and RMRP take ~34% of CLCA TPM on average (range ~5–70% per sample), consistent with a total-RNA library. They aren't in the shared gene set, and renormalizing TPM over shared genes removes the effect.
- **[2026-09-29] To check: 5 LICA-FR samples with only 10–20% of TPM in shared genes** (`results/01_harmonized/sample_stats.tsv`, `tpm_frac_shared_genes`). They may be outliers; check them in PCA before MINER.
- **[2026-09-29] Step 01 result (real data): 14,657 genes × 929 samples.** TCGA 366 (371 primary minus 3 fibrolamellar and duplicate/non-`01A` aliquots), CLCA 239, LICA-FR 324. 22,224 Ensembl genes shared by all cohorts; 14,657 pass TPM ≥ 1 in ≥ 20% of every cohort. Symbol mapping: CLCA 92% of rows, LICA-FR 71% (warning), but unmapped rows hold only ~0.5% (CLCA) and ~1% (LICA-FR) of each sample's TPM. They are mostly clone-name lncRNAs, so the loss is negligible. 3 LICA-FR expression columns (`BCB307T`, `BCM257T`, `BCM269T`) are not in the annotation and are dropped.
  - LICA-FR: 440 annotated samples pass the filter, but only 324 have RNA-seq (the other 116 are WES-only).
- **[2026-09-29] LICA-FR sample selection: primary, HCC, one tumor per patient.** The annotation has one row per tumor (628). Step 01 now takes a multi-column filter (`keep_filters`): `Sample_type == primary` & `primary_diagnosis == HCC` & `Per_patient_analysis == yes` (the authors' one-tumor-per-patient flag). This drops NTL, relapse, metastasis, HCC-on-HCA and second tumors of the same patient.
- **[2026-09-29] TCGA survival in step 03.** Built from the cBioPortal PanCan patient file (TCGA-CDR times, months × 30.44 → days). OS = `OS_MONTHS`/`OS_STATUS`. The RFS-like endpoint is `PFS_MONTHS`/`PFS_STATUS` (CDR PFI, which TCGA-CDR recommends for LIHC); it's not strictly RFS, so keep this in mind when comparing with CLCA RFS. Step 03 skips a cohort whose input file is missing, instead of stopping.
- **[2026-09-29] `envs/hcc-prep` pins numpy ≥ 2 and matplotlib ≥ 3.9 through conda.** Otherwise pip's `inmoose` → `fastcluster` pulls numpy 2 over a numpy-1 matplotlib build, and imports fail (`_ARRAY_API not found`).
- **[2026-09-29] SLURM wrappers** source conda from `${CONDA_BASE:-/users/sturkars/mambaforge}` (conda may not be on PATH in batch jobs) and set `PYTHONUTF8=1`, because the server locale is ISO-8859-1 and Python would otherwise read UTF-8 text files with the wrong encoding.
- **[2026-09-29] QC figures for steps 01–03** (`scripts/qc_plots.py`, written to `results/<step>/qc/`; list in README).
  - Colour follows the entity. Cohorts use categorical slots 1–3, the only slots validated for all-pairs use in scatters. Correction methods use gray/violet/green, so they never reuse a cohort colour.
  - Labels with more than 3 levels (Boyault G1–G6 etc.) are shown as small multiples, one level highlighted per panel, not as extra colours.
  - Both palettes were run through the colour validator: they pass colourblind separation. The uncorrected gray fails the chroma check on purpose, as a neutral baseline.
- **[2026-09-29] Step 01 harmonization rules.** Symbol → Ensembl tiers, first hit wins: HGNC approved → MINER Gene Name → HGNC previous → aliases (MINER Synonym + HGNC alias). A symbol matching more than one Ensembl ID in its tier is dropped and logged in `id_mapping_<cohort>.tsv`. Rows mapping to the same Ensembl ID are summed, since TPM is additive. Only genes shared by all cohorts are kept, and TPM is renormalized to 1e6 over them. A gene is kept if TPM ≥ 1 in ≥ 20% of samples in **every** cohort, so cohort-specific detection can't drive the batch effect. Output is log2(TPM+1). For TCGA, only `-01A` samples, one per patient.
- **[2026-09-29] Step 02: ComBat primary, per-cohort z as comparison.** ComBat (inmoose, parametric, cohort as batch, no covariates, because the LICA-FR labels exist in only one cohort). Genes with zero variance in any cohort are dropped. The final matrix is the gene z-score across all samples, clipped below −4 as in MINER's `preProcessTPM`. QC per matrix: silhouette and kNN mixing by cohort on the top 20 PCs, within-cohort structure preservation (Spearman of gene-centered sample-sample correlations), LICA-FR label silhouettes, and marker-program scores (CTNNB1, proliferation, hepatocyte) with per-cohort mean/SD and Kruskal–Wallis H by label.
  - **Caveat, seen in the synthetic test.** Correction removes real differences in subtype composition between cohorts, along with the batch effect. For example, a CTNNB1 program enriched in LICA-FR is pulled to the cross-cohort mean. The within-LICA-FR Boyault association is preserved. The same applies to etiology (CLCA ≈ HBV). Keep this in mind when comparing state frequencies across cohorts.
- **[2026-09-29] Use GDC STAR TPM for TCGA, not the cBioPortal file.** The cBioPortal `data_mrna_seq_v2_rsem.txt` is RSEM upper-quartile-normalized counts, not TPM (column sums ≈ 2.8×10⁷). Downloaded `TCGA-LIHC.star_tpm.tsv.gz` from UCSC Xena GDC hub into `data/LIHC-TCGA-GDC-Xena/`. **Values are log2(TPM + 0.001)**; back-transform before combining with CLCA/LICA-FR TPM. Ensembl IDs are versioned (GENCODE v36). 371 primary tumors (`-01A`), 3 recurrent (`-02`), 50 normals (`-11`).
- **[2026-09-29] Risk modeling uses MINER3's own framework**, not a custom plan (user direction). Details under Environment/tooling.
- **[2026-09-29] Follow the working GBM MINER3 run** in `/proj/omics4tb2/sturkarslan/GBM-15370004/analysis_stringent_sct/`: `miner3-coexpr -mg 6` → `miner3-mechinf -mc 0.1` → `miner3-subtypes`, all with `--skip_tpm` on pre-z-scored input. `miner3-subtypes` writes `coherentMembers.csv`, `overExpressedMembers.csv`, `transcriptional_programs.json`, `transcriptional_states.json`, so no separate bcmembers step is needed.
- **[2026-09-29] All project files stay under `/proj/omics4tb2/sturkarslan/HCC`.** Includes this log (used instead of the Claude memory dir), temp files, cloned references (`external/`), and any new conda env (`envs/`).
- **[2026-09-29] Data never goes to git.** `data/`, `results/`, `envs/`, `external/`, and binary/large formats are in `.gitignore`. Only code, configs and docs are pushed.
- **[2026-09-29] Explicit cross-cohort batch correction is required.** Verified in `miner.py` (see Environment/tooling): MINER3 has no TMM step and no cohort-aware batch correction. Plan: TPM → log2(TPM+1) → gene filtering → cohort correction (ComBat; per-cohort z-scoring as comparison) → QC (PCA/silhouette by cohort; preservation of Boyault G1–G6, CTNNB1 program, proliferation class) → z-scored matrix to MINER with `--skip_tpm`.
- **[2026-09-29] Validation approved by user.** Build the network on all three cohorts. Validate by (a) leave-one-cohort-out stability of regulons and states, and (b) external cohorts: ICGC LIRI-JP (RNA-seq; survival available in the user-provided Table S2) and GSE14520/LCI (Affymetrix, OS + RFS). Added as step 08.

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
