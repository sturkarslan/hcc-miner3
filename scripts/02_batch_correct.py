#!/usr/bin/env python
"""Step 02: cross-cohort batch correction and QC.

Input: results/01_harmonized/expression_log2tpm1.csv and samples.tsv.

Methods (config batch.methods):
  combat    ComBat on log2(TPM+1) with cohort as batch (inmoose.pycombat_norm). Primary.
  cohort_z  gene z-score within each cohort. Comparison.

Each corrected matrix is then gene z-scored across all samples and clipped below
batch.z_clip_low, so it can go to MINER3 with --skip_tpm (MINER's zscore()
passes pre-z-scored input through unchanged).

Outputs (results/02_batch_corrected/):
  expression_<method>_log2.csv   corrected, before the final z-score (combat only)
  expression_<method>_z.csv      MINER3 input, genes x samples, Ensembl IDs
  qc_metrics.tsv                 one row per matrix (uncorrected + each method)
  qc_programs.tsv                marker-program scores: per-cohort mean/SD, LICA-FR label association
  qc/b1..b9_*.png                PCA by cohort and by label, PC-covariate association, per-gene
                                 cohort variance, metric summary, program scores, RLE, correlation
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from hcc_common import load_params, p, setup_logging

IN, OUT = "01_harmonized", "02_batch_corrected"


def zscore_rows(df):
    sd = df.std(axis=1, ddof=1).replace(0, np.nan)
    return df.sub(df.mean(axis=1), axis=0).div(sd, axis=0).fillna(0.0)


def combat(expr, batch, parametric, log):
    from inmoose.pycombat import pycombat_norm
    out = pycombat_norm(expr, batch.tolist(), par_prior=parametric)
    out = pd.DataFrame(np.asarray(out), index=expr.index, columns=expr.columns)
    log.info("ComBat done (parametric=%s)", parametric)
    return out


def cohort_z(expr, batch):
    return pd.concat([zscore_rows(expr.loc[:, (batch == b).values]) for b in batch.unique()],
                     axis=1)[expr.columns]


# ------------------------------------------------------------------ QC

def pcs(expr, n):
    x = expr.sub(expr.mean(axis=1), axis=0).T.values
    n = min(n, min(x.shape) - 1)
    return PCA(n_components=n, random_state=0).fit_transform(x)


def knn_mixing(X, labels, k=15):
    """Mean fraction of each sample's k neighbours from other cohorts, divided by the
    fraction expected under perfect mixing. 1 = well mixed, 0 = fully separated."""
    labels = np.asarray(labels)
    k = min(k, len(labels) - 1)
    idx = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X, return_distance=False)[:, 1:]
    other = (labels[idx] != labels[:, None]).mean(axis=1)
    freq = pd.Series(labels).map(pd.Series(labels).value_counts(normalize=True)).values
    return float(np.mean(other / (1 - freq)))


def within_cohort_preservation(before, after, batch):
    """Spearman correlation between sample-sample correlation structures before and
    after correction, per cohort (1 = within-cohort structure unchanged). Genes are
    centered within the cohort first, so removing a gene's cohort mean does not count."""
    res = {}
    for b in batch.unique():
        cols = batch.index[batch == b]
        x0, x1 = before[cols], after[cols]
        c0 = np.corrcoef(x0.sub(x0.mean(axis=1), axis=0).T.values)
        c1 = np.corrcoef(x1.sub(x1.mean(axis=1), axis=0).T.values)
        iu = np.triu_indices_from(c0, k=1)
        res[b] = stats.spearmanr(c0[iu], c1[iu])[0]
    return res


def program_scores(z, genes, programs, log):
    sym2ens = genes.reset_index().dropna(subset=["symbol"]).groupby("symbol")["ensembl"].first()
    out = {}
    for prog, syms in programs.items():
        ens = [sym2ens[s] for s in syms if s in sym2ens.index and sym2ens[s] in z.index]
        if len(ens) < 2:
            log.warning("program %s: only %d of %d genes present; skipped", prog, len(ens), len(syms))
            continue
        out[prog] = z.loc[ens].mean()
    return pd.DataFrame(out)


def label_association(scores, labels):
    """Kruskal-Wallis H of each program score across a categorical label."""
    ok = labels.dropna()
    groups = [g.index for _, g in ok.groupby(ok) if len(g) >= 3]
    if len(groups) < 2:
        return {}
    return {prog: stats.kruskal(*[scores.loc[g, prog] for g in groups])[0] for prog in scores.columns}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    args = ap.parse_args()

    P = load_params(args.params)
    B = P["batch"]
    indir = p(os.path.join(P["paths"]["results"], IN))
    outdir = p(os.path.join(P["paths"]["results"], OUT))
    log = setup_logging(outdir, "02_batch_correct")

    expr = pd.read_csv(os.path.join(indir, "expression_log2tpm1.csv"), index_col=0)
    samples = pd.read_csv(os.path.join(indir, "samples.tsv"), sep="\t", index_col="sample").loc[expr.columns]
    genes = pd.read_csv(os.path.join(indir, "genes.tsv"), sep="\t", index_col="ensembl")
    batch = samples["cohort"]
    log.info("Input: %d genes x %d samples; cohorts %s", *expr.shape, batch.value_counts().to_dict())

    # ComBat cannot handle genes with zero variance inside a batch
    zero_var = pd.concat([expr.loc[:, batch == b].var(axis=1) == 0 for b in batch.unique()], axis=1).any(axis=1)
    if zero_var.any():
        log.info("Dropping %d genes with zero variance in at least one cohort", zero_var.sum())
        expr = expr.loc[~zero_var]

    corrected = {"uncorrected": expr}
    for m in B["methods"]:
        if m == "combat":
            corrected[m] = combat(expr, batch, B.get("combat_parametric", True), log)
            corrected[m].to_csv(os.path.join(outdir, f"expression_{m}_log2.csv"))
        elif m == "cohort_z":
            corrected[m] = cohort_z(expr, batch)
        else:
            raise ValueError(f"unknown batch method {m!r}")

    finals = {}
    for m, mat in corrected.items():
        z = zscore_rows(mat)
        if B.get("z_clip_low") is not None:
            z = z.clip(lower=B["z_clip_low"])
        finals[m] = z
        if m != "uncorrected":
            z.to_csv(os.path.join(outdir, f"expression_{m}_z.csv"))
            log.info("Wrote MINER3 input expression_%s_z.csv (%d x %d)", m, *z.shape)

    # ---- QC
    Q = B.get("qc", {})
    label_cols = [c for c in samples.columns if c not in ("cohort", "patient")]
    rows, prog_rows, prog_scores = [], [], {}
    for m, z in finals.items():
        X = pcs(z, Q.get("n_pcs", 20))
        row = {"matrix": m,
               "silhouette_cohort": silhouette_score(X, batch.values),
               "knn_mixing_cohort": knn_mixing(X, batch.values)}
        if m != "uncorrected":
            for b, r in within_cohort_preservation(corrected["uncorrected"], corrected[m], batch).items():
                row[f"within_structure_rho_{b}"] = r
        for lab in label_cols:
            has = samples[lab].notna()
            if samples.loc[has, lab].nunique() >= 2:
                # label separation inside the samples that carry the label (LICA-FR)
                Xl = X[has.values]
                row[f"silhouette_{lab}"] = silhouette_score(Xl, samples.loc[has, lab].astype(str).values)
        rows.append(row)

        sc = prog_scores[m] = program_scores(z, genes, Q.get("programs", {}), log)
        for prog in sc.columns:
            pr = {"matrix": m, "program": prog}
            for b in batch.unique():
                pr[f"mean_{b}"] = sc.loc[batch == b, prog].mean()
                pr[f"sd_{b}"] = sc.loc[batch == b, prog].std()
            for lab in label_cols:
                h = label_association(sc, samples[lab])
                if prog in h:
                    pr[f"kruskalH_{lab}"] = h[prog]
            prog_rows.append(pr)

    qc = pd.DataFrame(rows)
    qc.to_csv(os.path.join(outdir, "qc_metrics.tsv"), sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(prog_rows).to_csv(os.path.join(outdir, "qc_programs.tsv"), sep="\t", index=False,
                                   float_format="%.4f")
    log.info("QC metrics (silhouette_cohort: lower is better; knn_mixing_cohort: 1 = mixed; "
             "label silhouettes and within_structure_rho: should stay close to uncorrected)\n%s",
             qc.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    import qc_plots
    written = qc_plots.batch_report(outdir, finals, corrected, samples, genes, qc, prog_scores, label_cols)
    log.info("QC figures in %s/qc: %s", outdir, ", ".join(written))

if __name__ == "__main__":
    main()
