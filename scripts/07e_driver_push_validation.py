#!/usr/bin/env python
"""Step 07e: unbiased test of "causal net risk push explains driver prognosis" (Figure 1f).

Figure 1f (step 07d) correlates each driver's net risk push, computed with program weights of the
ridge model trained in TCGA, with the driver's Cox z meta-analysed over TCGA + CLCA. Problems this
step addresses:
  1. Circularity: the weights were learned from TCGA recurrence, and TCGA recurrence is half of the
     y-axis. Here weights from the model trained in cohort A are only ever compared with driver Cox z
     computed in the other cohort B (both directions). The old in-sample pairing is reported as a
     reference row so its inflation is visible.
  2. Breadth: the original push is a sum over causal families, so drivers touching many families
     (arm events) get large values regardless of per-family effect. Mean and effect-size-weighted
     versions are added; a strict version takes Cohen's d from the training cohort only.
  3. Does the causal layer add anything? The model's total predicted risk shift of a driver,
     sum_k w_k * Delta_k (Delta_k = altered - wild-type difference of standardized program activity in
     the test cohort), is split into the programs the driver causally reaches and the rest.
  4. Significance: permutation nulls that keep the dependence between drivers:
       weights  - program weights shuffled across programs (do the specific weights matter?)
       families - each driver's causal families (or reached programs) replaced by the same number of
                  random ones, with the driver's own effect on them (does causal selection matter?)
  5. Non-independent drivers (TP53 / p53 pathway / 17p loss; CTNNB1 / WNT pathway ...): drivers are
     clustered by overlap of altered tumours; effective n, rho on cluster representatives and a
     cluster-bootstrap CI are reported.
  6. Pre-specified analysis: primary = config post.driver_push.primary; everything else is labelled
     sensitivity. All eligible drivers are used (no hand-picked arm events); the Figure 1f set is a
     sensitivity subset.

No survival information enters the predictors except through the training cohort's weights; the
test cohort's survival is used only for the y-axis.

Inputs: results/04_miner/<matrix>/<risk.subtypes_dir>/{eigengenes.csv, transcriptional_programs.json},
        results/05_causal/<matrix>/{highConfidenceCausalResults.csv, regulon_families.tsv},
        results/03_genomics_clinical/{genomic_features.csv, genomic_features_info.tsv, survival_*},
        results/06_risk/<matrix>/predictor_ridge_programs_<A>_<EP>_h36m/{weights.tsv, predictions.tsv}
Outputs (results/07_post/driver_push/):
  driver_table.tsv     one row per driver x pairing x endpoint: test-cohort Cox z (raw, stage-adjusted),
                       all predictors, cluster
  summary.tsv          one row per pairing x endpoint x adjustment x driver set x predictor: n, n_clusters,
                       Spearman rho, naive p, permutation p (weights, families), cluster-representative rho,
                       cluster-bootstrap 95% CI; `primary` flags the pre-specified test
  driver_clusters.tsv  driver -> cluster; jaccard.tsv
  qc/e1-e4 PNGs
"""

import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from hcc_common import load_params, p, setup_logging

warnings.filterwarnings("ignore")
OUT = os.path.join("07_post", "driver_push")
PREDICTORS = ["sign_sum", "sign_mean", "d_mean", "d_mean_strict", "delta_total", "delta_causal", "delta_noncausal"]
FAMILY_NULL = {"sign_sum", "sign_mean", "d_mean", "delta_causal"}


# ------------------------------------------------------------------ helpers

def cohen_d(M, a, b):
    """Cohen's d of every row of M (DataFrame) between sample lists a and b."""
    x, y = M[a].values, M[b].values
    sp = np.sqrt(((len(a) - 1) * x.var(1, ddof=1) + (len(b) - 1) * y.var(1, ddof=1)) / (len(a) + len(b) - 2))
    return pd.Series((x.mean(1) - y.mean(1)) / np.where(sp > 0, sp, np.nan), index=M.index)


def cox_z(x, srv, stage=None):
    from lifelines import CoxPHFitter
    d = srv.assign(x=x)
    if stage is not None:
        d = d.assign(stage=stage).dropna()
        if d["stage"].nunique() < 2:
            return np.nan
    if d["x"].sum() < 3 or (d["x"] == 0).sum() < 3 or d["observed"].sum() < 3:
        return np.nan
    try:
        return float(CoxPHFitter().fit(d, "duration", "observed").summary.loc["x", "z"])
    except Exception:  # noqa: BLE001 - separation / convergence
        return np.nan


def stage_ordinal(df, cohort):
    if cohort == "CLCA" and "BCLC" in df:
        return df["BCLC"].astype(str).str.strip().map({"0": 0, "A": 1, "B": 2, "C": 3, "D": 4}).astype(float)
    if cohort == "TCGA" and "stage" in df:
        return df["stage"].astype(str).str.extract(r"STAGE (I{1,3}V?|IV)")[0].map(
            {"I": 1, "II": 2, "III": 3, "IV": 4}).astype(float)
    return pd.Series(np.nan, index=df.index)


def spearman(x, y):
    ok = ~(np.isnan(x) | np.isnan(y))
    if ok.sum() < 5:
        return np.nan, np.nan, int(ok.sum())
    r, pv = stats.spearmanr(x[ok], y[ok])
    return float(r), float(pv), int(ok.sum())


def driver_clusters(F, drivers, thr):
    """Average-linkage clusters on Jaccard overlap of altered samples (samples profiled for both)."""
    J = pd.DataFrame(np.eye(len(drivers)), index=drivers, columns=drivers)
    for i, a in enumerate(drivers):
        for b in drivers[i + 1:]:
            both = F.loc[[a, b]].dropna(axis=1)
            A, B = both.loc[a] == 1, both.loc[b] == 1
            u = (A | B).sum()
            J.loc[a, b] = J.loc[b, a] = (A & B).sum() / u if u else 0.0
    if len(drivers) < 2:
        return J, pd.Series(1, index=drivers)
    lab = fcluster(linkage(squareform(1 - J.values, checks=False), "average"), 1 - thr, "distance")
    return J, pd.Series(lab, index=drivers)


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    Rk, Q = P["risk"], P["post"]["driver_push"]
    matrix = args.matrix or P["miner"]["matrix"]
    res = p(P["paths"]["results"])
    outdir = os.path.join(res, OUT)
    log = setup_logging(outdir, "07e_driver_push_validation")
    rng = np.random.default_rng(int(Q.get("seed", 7)))
    sfx = Rk["survival_suffix"]

    # ---- network: program activity exactly as step 06 (mean of regulon eigengenes centred over samples)
    sdir = os.path.join(res, "04_miner", matrix, Rk["subtypes_dir"])
    E = pd.read_csv(os.path.join(sdir, "eigengenes.csv"), index_col=0)
    E.index = E.index.astype(str)
    progs = {str(k): [str(r) for r in v] for k, v in json.load(open(os.path.join(sdir, "transcriptional_programs.json"))).items()}
    Ez = E.sub(E.mean(1), axis=0)
    act = pd.DataFrame({k: Ez.loc[[r for r in v if r in Ez.index]].mean() for k, v in progs.items()}).T
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    cohort = samples["cohort"].reindex(E.columns)

    # ---- causal flows: one row per driver x family (largest |d|), as step 07d
    cdir = os.path.join(res, "05_causal", matrix)
    hc = pd.read_csv(os.path.join(cdir, "highConfidenceCausalResults.csv"), index_col=0)
    hc["program"] = pd.to_numeric(hc["program"], errors="coerce").astype("Int64").astype(str)
    hc["Regulon_ID"] = hc["Regulon_ID"].astype(str)
    fams = (hc.assign(ad=hc["cohen_d"].abs()).sort_values("ad", ascending=False)
            .drop_duplicates(["Mutation", "family"]).drop(columns="ad"))
    fam_df = pd.read_csv(os.path.join(cdir, "regulon_families.tsv"), sep="\t", index_col=0)
    fam_df.index = fam_df.index.astype(str)
    fam_df["program"] = pd.to_numeric(fam_df["program"], errors="coerce").astype("Int64").astype(str)
    fam_df = fam_df.loc[fam_df.index.intersection(E.index)]

    F = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
    F = F.reindex(columns=E.columns)
    info = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features_info.tsv"), sep="\t", index_col=0)
    drivers_all = [f for f in fams["Mutation"].unique() if f in F.index]
    fig_arms = set(P["post"].get("causal_figures", {}).get("arm_events", []))
    fig_set = {f for f in drivers_all if info.loc[f, "type"] != "arm_cna"} | (fig_arms & set(drivers_all))

    def surv_file(c, ep):
        s = pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_{c}_{ep}{sfx}_miner.csv"), index_col=0)
        s.columns = ["duration", "observed"]
        return s.loc[s.index.intersection(E.columns)]

    clin = {}
    for c in Rk["survival_cohorts"]:
        f = os.path.join(res, "03_genomics_clinical", f"survival_{c}.tsv")
        if os.path.exists(f):
            t = pd.read_csv(f, sep="\t", index_col=0)
            pt = samples.loc[samples["cohort"] == c, "patient"]
            clin[c] = stage_ordinal(t.reindex(pt.values).set_axis(pt.index), c)

    # pooled per-driver d of every family (largest-|d| regulon), for the random-family null
    fam_ids = fam_df["family"].values
    d_all_fam = {}
    for f in drivers_all:
        row = F.loc[f].dropna()
        mut, wt = list(row.index[row == 1]), list(row.index[row == 0])
        d = cohen_d(E, mut, wt)
        t = pd.DataFrame({"family": fam_ids, "program": fam_df["program"].values, "d": d.reindex(fam_df.index).values})
        t = t.dropna(subset=["d"]).assign(ad=lambda x: x["d"].abs()).sort_values("ad", ascending=False)
        d_all_fam[f] = t.drop_duplicates("family").set_index("family")[["program", "d"]]

    pairings = [dict(x, kind="out_of_sample") for x in Q["pairings"]]
    if Q.get("in_sample_reference", True):
        pairings += [dict(x, test=x["train"], kind="in_sample") for x in Q["pairings"]]

    rows, summ, nulls = [], [], {}
    clusters_out = {}
    for pr in pairings:
        A, B, ep_w = pr["train"], pr["test"], pr["train_endpoint"]
        tag = f"ridge_programs_{A}_{ep_w}{sfx}"
        wdir = os.path.join(res, "06_risk", matrix, f"predictor_{tag}")
        w = pd.read_csv(os.path.join(wdir, "weights.tsv"), sep="\t", index_col=0).iloc[:, 0]
        w.index = w.index.astype(str)
        prog_ids = [k for k in act.index if k in w.index]
        w = w.reindex(prog_ids)
        trainS = surv_file(A, ep_w).index
        mu, sd = act.loc[prog_ids, trainS].mean(1), act.loc[prog_ids, trainS].std(1).replace(0, 1)
        Xs = act.loc[prog_ids].sub(mu, axis=0).div(sd, axis=0)
        pf = os.path.join(wdir, "predictions.tsv")
        if os.path.exists(pf):
            pred = pd.read_csv(pf, sep="\t", index_col=0)["risk_score"]
            mine = (Xs.T @ w).reindex(pred.index)
            log.info("%s: reconstructed score vs step-06 predictions r = %.4f", tag, np.corrcoef(mine, pred)[0, 1])
        pname = f"{A}->{B}" + (" (in-sample)" if pr["kind"] == "in_sample" else "")
        k_index = {k: i for i, k in enumerate(prog_ids)}

        for ep in Q["endpoints"]:
            srvB = surv_file(B, ep)
            stageB = clin.get(B, pd.Series(dtype=float)).reindex(srvB.index)
            drivers, recs = [], []
            S_sign, S_d, S_ds, Delta, Reach, nfam, nfam_s = [], [], [], [], [], [], []
            for f in drivers_all:
                rowB = F.loc[f, srvB.index].dropna()
                if (rowB == 1).sum() < Q["min_altered_test"] or (rowB == 0).sum() < Q["min_wt_test"]:
                    continue
                ff = fams[(fams["Mutation"] == f) & fams["program"].isin(prog_ids)]
                if ff.empty:
                    continue
                x = rowB.astype(float)
                z = cox_z(x.values, srvB.loc[x.index])
                zs = cox_z(x.values, srvB.loc[x.index], stageB.loc[x.index].values) if Q.get("stage_adjust") else np.nan
                # Delta_k in the test cohort (expression only; all profiled test samples)
                exB = F.loc[f, cohort.index[cohort == B]].dropna()
                mB, wB = list(exB.index[exB == 1]), list(exB.index[exB == 0])
                dl = (Xs[mB].mean(1) - Xs[wB].mean(1)).values
                # strict d: family's regulon d within the training cohort only
                exA = F.loc[f, cohort.index[cohort == A]].dropna()
                mA, wA = list(exA.index[exA == 1]), list(exA.index[exA == 0])
                dA = cohen_d(E.loc[ff["Regulon_ID"].unique()], mA, wA) if min(len(mA), len(wA)) >= 5 else None
                ss, sdd, sds, rc = (np.zeros(len(prog_ids)) for _ in range(4))
                ns = 0
                for _, r in ff.iterrows():
                    k = k_index[r["program"]]
                    ss[k] += np.sign(r["cohen_d"])
                    sdd[k] += r["cohen_d"]
                    rc[k] = 1
                    if dA is not None and not np.isnan(dA.get(r["Regulon_ID"], np.nan)):
                        sds[k] += dA[r["Regulon_ID"]]
                        ns += 1
                drivers.append(f)
                S_sign.append(ss)
                S_d.append(sdd)
                S_ds.append(sds)
                Delta.append(dl)
                Reach.append(rc)
                nfam.append(len(ff))
                nfam_s.append(ns)
                recs.append({"pairing": pname, "kind": pr["kind"], "train": A, "test": B, "endpoint": ep,
                             "driver": f, "type": info.loc[f, "type"], "n_altered_test": int((rowB == 1).sum()),
                             "n_families": len(ff), "n_programs_reached": int(rc.sum()), "z": z, "z_stage": zs,
                             "in_figure_1f": f in fig_set})
            if len(drivers) < 5:
                log.warning("%s %s: only %d eligible drivers; skipped", pname, ep, len(drivers))
                continue
            S_sign, S_d, S_ds, Delta, Reach = map(np.array, (S_sign, S_d, S_ds, Delta, Reach))
            nfam, nfam_s = np.array(nfam, float), np.array(nfam_s, float)

            def predictors(wv, Ss=S_sign, Sd=S_d, Rc=Reach, nf=nfam):
                tot = Delta @ wv
                cau = (Delta * Rc) @ wv
                with np.errstate(invalid="ignore", divide="ignore"):
                    strict = np.where(nfam_s > 0, (S_ds @ wv) / nfam_s, np.nan)
                return {"sign_sum": Ss @ wv, "sign_mean": (Ss @ wv) / nf, "d_mean": (Sd @ wv) / nf,
                        "d_mean_strict": strict, "delta_total": tot, "delta_causal": cau, "delta_noncausal": tot - cau}

            wv = w.values
            obs = predictors(wv)
            tab = pd.DataFrame(recs)
            for k, v in obs.items():
                tab[k] = v
            J, cl = driver_clusters(F, drivers, Q["cluster_jaccard"])
            tab["cluster"] = cl.reindex(drivers).values
            clusters_out[(pname, ep)] = (J, cl)
            rows.append(tab)

            # ---- nulls (shared across driver sets / adjustments; recomputed on each subset)
            nperm = int(Q["n_perm"])
            null_w = [predictors(rng.permutation(wv)) for _ in range(nperm)]
            null_f = []
            for _ in range(int(Q["n_perm_families"])):
                Ss, Sd, Rc = np.zeros_like(S_sign), np.zeros_like(S_d), np.zeros_like(Reach)
                for i, f in enumerate(drivers):
                    pool = d_all_fam[f][d_all_fam[f]["program"].isin(prog_ids)]
                    pick = pool.iloc[rng.choice(len(pool), size=min(int(nfam[i]), len(pool)), replace=False)]
                    for prg, d in zip(pick["program"], pick["d"]):
                        Ss[i, k_index[prg]] += np.sign(d)
                        Sd[i, k_index[prg]] += d
                    Rc[i, rng.choice(len(prog_ids), size=int(Reach[i].sum()), replace=False)] = 1
                null_f.append(predictors(wv, Ss, Sd, Rc))
            nulls[(pname, ep)] = (null_w, null_f)

            sets = {"all": np.ones(len(drivers), bool),
                    "non_arm": (tab["type"] != "arm_cna").values,
                    "arm": (tab["type"] == "arm_cna").values,
                    "figure_1f_set": tab["in_figure_1f"].values}
            for adj, ycol in (("none", "z"), ("stage", "z_stage")):
                y = tab[ycol].values.astype(float)
                if np.isnan(y).all():
                    continue
                for sname, m in sets.items():
                    if m.sum() < 5:
                        continue
                    clm = tab["cluster"].values[m]
                    reps = (tab[m].assign(i=np.arange(len(tab))[m]).sort_values("n_altered_test", ascending=False)
                            .drop_duplicates("cluster")["i"].values)
                    for pk in PREDICTORS:
                        x = obs[pk]
                        rho, pv, n = spearman(x[m], y[m])
                        if np.isnan(rho):
                            continue
                        pw = np.mean([abs(spearman(nw[pk][m], y[m])[0]) >= abs(rho) for nw in null_w])
                        pfam = (np.mean([abs(spearman(nf_[pk][m], y[m])[0]) >= abs(rho) for nf_ in null_f])
                                if pk in FAMILY_NULL else np.nan)
                        rr, _, nr = spearman(x[reps], y[reps])
                        # cluster bootstrap
                        ucl = np.unique(clm)
                        idx_by = {c: np.where(m & (tab["cluster"].values == c))[0] for c in ucl}
                        boots = []
                        for _ in range(int(Q["n_boot"])):
                            ii = np.concatenate([idx_by[c] for c in rng.choice(ucl, size=len(ucl), replace=True)])
                            boots.append(spearman(x[ii], y[ii])[0])
                        lo, hi = np.nanpercentile(boots, [2.5, 97.5]) if np.isfinite(boots).sum() > 10 else (np.nan, np.nan)
                        summ.append({"pairing": pname, "kind": pr["kind"], "endpoint": ep, "adjustment": adj,
                                     "driver_set": sname, "predictor": pk, "n": n, "n_clusters": len(ucl),
                                     "rho": rho, "p_naive": pv, "p_perm_weights": pw, "p_perm_families": pfam,
                                     "rho_cluster_reps": rr, "n_cluster_reps": nr, "boot_lo": lo, "boot_hi": hi})
            log.info("%s %s: %d drivers in %d clusters", pname, ep, len(drivers), tab["cluster"].nunique())

    if not rows:
        raise RuntimeError("no pairing produced eligible drivers")
    table = pd.concat(rows, ignore_index=True)
    table.to_csv(os.path.join(outdir, "driver_table.tsv"), sep="\t", index=False, float_format="%.4g")
    S = pd.DataFrame(summ)
    pr_ = Q["primary"]
    S["primary"] = ((S["kind"] == "out_of_sample") & (S["endpoint"] == pr_["endpoint"]) & (S["adjustment"] == pr_["adjustment"])
                    & (S["driver_set"] == pr_["driver_set"]) & (S["predictor"] == pr_["predictor"]))
    S.to_csv(os.path.join(outdir, "summary.tsv"), sep="\t", index=False, float_format="%.4g")
    cl_rows = []
    for (pn, ep), (J, cl) in clusters_out.items():
        cl_rows += [{"pairing": pn, "endpoint": ep, "driver": d, "cluster": c} for d, c in cl.items()]
    pd.DataFrame(cl_rows).to_csv(os.path.join(outdir, "driver_clusters.tsv"), sep="\t", index=False)
    key0 = next(iter(clusters_out))
    clusters_out[key0][0].to_csv(os.path.join(outdir, "jaccard.tsv"), sep="\t", float_format="%.3f")

    show = S[(S["endpoint"] == pr_["endpoint"]) & (S["adjustment"] == "none") & (S["driver_set"] == "all")]
    log.info("Spearman rho, all drivers, %s, unadjusted (p_perm_weights / p_perm_families; cluster-bootstrap CI):\n%s",
             pr_["endpoint"], show[["pairing", "predictor", "n", "n_clusters", "rho", "p_naive", "p_perm_weights",
                                    "p_perm_families", "rho_cluster_reps", "boot_lo", "boot_hi"]]
             .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    log.info("PRIMARY (pre-specified):\n%s", S[S["primary"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))

    import qc_plots
    written = qc_plots.driver_push_report(outdir, table, S, nulls, clusters_out, pr_)
    log.info("QC figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
