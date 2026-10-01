#!/usr/bin/env python
"""Step 07f: AXIN1 versus CTNNB1 mutation, two WNT-pathway drivers with different programs.

Both genes are counted as WNT/beta-catenin pathway drivers in HCC. This step asks, for every program,
whether AXIN1-only and CTNNB1-only mutant tumours differ from double wild-type tumours and from each
other. All programs are tested (no selection); BH q-values are over programs.

Model, per program (activity standardized over all profiled tumours):
    activity ~ AXIN1_only + CTNNB1_only + covariates + cohort
Tumours with both mutations are excluded (too few to estimate). Covariates are in
post.axin1_ctnnb1.adjust (features from the genomic matrix; programs by id). A program used as a
covariate is not adjusted for itself. The same model is fitted to marker genes (z expression) and,
within each cohort without the cohort term, to check that signs agree.

Causal layer: high-confidence driver -> regulator edges of the two drivers (step 05), one row per
regulator, with the direction under each driver.

Inputs: results/04_miner/<matrix>/<post.subtypes_dir>/{eigengenes.csv, transcriptional_programs.json},
        results/03_genomics_clinical/genomic_features.csv, results/01_harmonized/{samples,genes}.tsv,
        results/02_batch_corrected/expression_<matrix>_z.csv,
        results/05_causal/<matrix>/highConfidenceCausalResults.csv
Outputs (results/07_post/axin1_ctnnb1/):
  groups.tsv             sample -> group, cohort
  program_effects.tsv    per program: beta, t, p, q for each driver and for AXIN1 - CTNNB1; per-cohort betas
  gene_effects.tsv       same for marker genes
  regulator_contrast.tsv regulators reached by either driver, direction and largest |d| under each
  qc/g1_program_effects.png
"""

import argparse
import json
import os

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import qc_plots as Q  # noqa: E402
from hcc_common import load_params, p, setup_logging  # noqa: E402

OUT = os.path.join("07_post", "axin1_ctnnb1")


def bh(pv):
    pv = np.asarray(pv, float)
    o = np.argsort(pv)
    q = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def ols(y, X):
    """OLS with intercept already in X. Returns beta, covariance, residual df."""
    XtX_inv = np.linalg.pinv(X.T @ X)
    b = XtX_inv @ X.T @ y
    df = len(y) - np.linalg.matrix_rank(X)
    s2 = ((y - X @ b) ** 2).sum() / df
    return b, s2 * XtX_inv, df


def fit(y, a, c, covars, cohort):
    """Effects of AXIN1-only (a) and CTNNB1-only (c) and their difference. cohort=None: no cohort term."""
    cols = [np.ones(len(y)), a, c] + [v for v in covars]
    if cohort is not None:
        cols += [(cohort == k).astype(float).values for k in sorted(cohort.unique())[1:]]
    X = np.column_stack(cols)
    b, V, df = ols(y, X)
    out = {}
    for name, L in (("AXIN1", [0, 1, 0]), ("CTNNB1", [0, 0, 1]), ("diff", [0, 1, -1])):
        L = np.array(L + [0] * (X.shape[1] - 3), float)
        est, se = L @ b, np.sqrt(L @ V @ L)
        out[f"beta_{name}"], out[f"t_{name}"] = est, est / se
        out[f"p_{name}"] = 2 * stats.t.sf(abs(est / se), df)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    C = P["post"]["axin1_ctnnb1"]
    matrix = args.matrix or P["miner"]["matrix"]
    res = p(P["paths"]["results"])
    outdir = os.path.join(res, OUT)
    os.makedirs(os.path.join(outdir, "qc"), exist_ok=True)
    log = setup_logging(outdir, "07f_axin1_ctnnb1")

    # ---- program activity as steps 06 / 07e, then standardized per program
    sdir = os.path.join(res, "04_miner", matrix, P["post"]["subtypes_dir"])
    E = pd.read_csv(os.path.join(sdir, "eigengenes.csv"), index_col=0)
    E.index = E.index.astype(str)
    progs = {str(k): [str(r) for r in v] for k, v in json.load(open(os.path.join(sdir, "transcriptional_programs.json"))).items()}
    Ez = E.sub(E.mean(1), axis=0)
    act = pd.DataFrame({k: Ez.loc[[r for r in v if r in Ez.index]].mean() for k, v in progs.items()}).T
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    F = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)

    fa, fc = C["features"]["AXIN1"], C["features"]["CTNNB1"]
    feat_cov = [x for x in C["adjust"] if x in F.index]
    prog_cov = [str(x) for x in C["adjust"] if str(x) in act.index and x not in F.index]
    keep = [s for s in act.columns if s in F.columns and F.loc[[fa, fc] + feat_cov, s].notna().all()]
    Fk = F.loc[:, keep]
    both = (Fk.loc[fa] == 1) & (Fk.loc[fc] == 1)
    keep = [s for s in keep if not both[s]]
    Fk, cohort = F.loc[:, keep], samples["cohort"].reindex(keep)
    a, c = Fk.loc[fa].values.astype(float), Fk.loc[fc].values.astype(float)
    grp = pd.Series(np.where(a == 1, "AXIN1", np.where(c == 1, "CTNNB1", "WT")), index=keep, name="group")
    pd.concat([grp, cohort], axis=1).to_csv(os.path.join(outdir, "groups.tsv"), sep="\t")
    log.info("Profiled tumours %d (both mutated, excluded: %d). %s", len(keep), int(both.sum()),
             grp.groupby(cohort).value_counts().unstack(fill_value=0).to_dict("index"))
    log.info("Adjusted for features %s, programs %s, cohort", feat_cov, prog_cov)

    A = act.loc[:, keep]
    A = A.sub(A.mean(1), axis=0).div(A.std(1), axis=0)

    def run(M, drop_self):
        rows = {}
        for k in M.index:
            y = M.loc[k].values.astype(float)
            cov = [Fk.loc[f].values.astype(float) for f in feat_cov] + \
                  [A.loc[q].values for q in prog_cov if not (drop_self and q == k)]
            r = fit(y, a, c, cov, cohort)
            for coh in sorted(cohort.unique()):
                m = (cohort == coh).values
                if a[m].sum() >= C["min_group_cohort"] and c[m].sum() >= C["min_group_cohort"]:
                    rc = fit(y[m], a[m], c[m], [v[m] for v in cov], None)
                    r[f"beta_AXIN1_{coh}"], r[f"beta_CTNNB1_{coh}"] = rc["beta_AXIN1"], rc["beta_CTNNB1"]
            rows[k] = r
        T = pd.DataFrame(rows).T
        for n in ("AXIN1", "CTNNB1", "diff"):
            T[f"q_{n}"] = bh(T[f"p_{n}"].astype(float))
        for n in ("AXIN1", "CTNNB1"):
            cc = [x for x in T.columns if x.startswith(f"beta_{n}_")]
            T[f"sign_consistent_{n}"] = [bool((np.sign(T.loc[i, cc].astype(float)) == np.sign(T.loc[i, f"beta_{n}"])).all())
                                         for i in T.index]
        return T

    # ---- programs
    T = run(A, drop_self=True)
    T.insert(0, "n_regulons", [len(progs[k]) for k in T.index])
    alpha = C["q_max"]
    sa, sc, sd = T["q_AXIN1"] <= alpha, T["q_CTNNB1"] <= alpha, T["q_diff"] <= alpha
    same = np.sign(T["beta_AXIN1"].astype(float)) == np.sign(T["beta_CTNNB1"].astype(float))
    T["pattern"] = np.select([sa & sc & ~same, sa & sc & same, sc & ~sa, sa & ~sc], ["opposite", "shared", "CTNNB1 only", "AXIN1 only"], "neither")
    T.index.name = "program"
    T.sort_values("p_diff").to_csv(os.path.join(outdir, "program_effects.tsv"), sep="\t")
    log.info("Programs (q <= %.2f): %s; differ between drivers: %d of %d", alpha, T["pattern"].value_counts().to_dict(),
             int(sd.sum()), len(T))
    log.info("Largest differences:\n%s", T.sort_values("p_diff")[["beta_AXIN1", "t_AXIN1", "beta_CTNNB1", "t_CTNNB1", "q_diff", "pattern"]]
             .head(15).astype({"beta_AXIN1": float, "t_AXIN1": float, "beta_CTNNB1": float, "t_CTNNB1": float, "q_diff": float}).round(3).to_string())

    # ---- marker genes
    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t")
    sym = genes[genes["kept"]].drop_duplicates("symbol").set_index("symbol")["ensembl"]
    want = [(g, s) for g, ss in C["genes"].items() for s in ss]
    found = [(g, s) for g, s in want if s in sym.index]
    X = pd.read_csv(os.path.join(res, "02_batch_corrected", f"expression_{matrix}_z.csv"), index_col=0)
    Xg = X.loc[[sym[s] for _, s in found], keep]
    Xg.index = [s for _, s in found]
    G = run(Xg, drop_self=False)
    G.insert(0, "set", [g for g, _ in found])
    G.index.name = "gene"
    G.to_csv(os.path.join(outdir, "gene_effects.tsv"), sep="\t")
    log.info("Marker genes (missing: %s):\n%s", [s for _, s in want if s not in sym.index],
             G[["set", "beta_AXIN1", "t_AXIN1", "beta_CTNNB1", "t_CTNNB1", "q_diff"]].to_string(float_format="%.2f"))

    # ---- causal layer: regulators under each driver
    hc = pd.read_csv(os.path.join(res, "05_causal", matrix, "highConfidenceCausalResults.csv"), index_col=0)
    hc = hc[hc["Mutation"].isin([fa, fc])].assign(ad=lambda d: d["cohen_d"].abs())
    R = (hc.groupby(["regulator_symbol", "Mutation"]).agg(edge=("MutationRegulatorEdge", "first"), n_regulons=("Regulon_ID", "nunique"),
                                                          max_abs_d=("ad", "max")).unstack("Mutation"))
    R.columns = [f"{x}_{'AXIN1' if m == fa else 'CTNNB1'}" for x, m in R.columns]
    R["relation"] = np.where(R["edge_AXIN1"].isna(), "CTNNB1 only", np.where(R["edge_CTNNB1"].isna(), "AXIN1 only",
                             np.where(R["edge_AXIN1"] == R["edge_CTNNB1"], "same direction", "opposite direction")))
    R.to_csv(os.path.join(outdir, "regulator_contrast.tsv"), sep="\t")
    log.info("Regulators: %s. Opposite: %s", R["relation"].value_counts().to_dict(),
             sorted(R.index[R["relation"] == "opposite direction"]))

    # ---- QC figure: every program, CTNNB1 effect vs AXIN1 effect
    Q.setup() if hasattr(Q, "setup") else None
    fig, ax = plt.subplots(figsize=(4.2, 4))
    col = {"opposite": Q.SLOTS[1], "shared": Q.SLOTS[0], "CTNNB1 only": Q.SLOTS[2], "AXIN1 only": Q.SLOTS[3], "neither": "#bbbbbb"}
    for k, d in T.groupby("pattern"):
        ax.scatter(d["beta_CTNNB1"].astype(float), d["beta_AXIN1"].astype(float), s=14, color=col[k], label=f"{k} ({len(d)})")
    for k in T.sort_values("p_diff").index[:12]:
        ax.annotate(f"P{k}", (float(T.loc[k, "beta_CTNNB1"]), float(T.loc[k, "beta_AXIN1"])), fontsize=6, xytext=(2, 2), textcoords="offset points")
    ax.axhline(0, color="k", lw=0.4)
    ax.axvline(0, color="k", lw=0.4)
    ax.set_xlabel("CTNNB1-mutant effect on program (s.d.)")
    ax.set_ylabel("AXIN1-mutant effect on program (s.d.)")
    ax.legend(fontsize=6, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "qc", "g1_program_effects.png"), dpi=200)
    log.info("Done")


if __name__ == "__main__":
    main()
