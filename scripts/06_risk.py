#!/usr/bin/env python
"""Step 06: risk analysis with MINER3's framework (miner3 env + project xgboost 1.7.6).

Survival is never pooled across cohorts (PROJECT_LOG): every model is fit within a cohort, on the
36-month-horizon files from step 03, and evaluated on another cohort by within-cohort discrimination.

Part A - prognostic regulons, programs and states (MINER's survival analysis, per cohort):
  Cox per regulon eigengene / program activity (mean regulon eigengene) / state membership, fit
  separately in each cohort with survival; Stouffer meta-z over cohorts (weighted by sqrt(events)),
  BH over units, and a consistency flag (same sign in every cohort).
Part B - MINER risk predictor (generatePredictor):
  features = regulon dysregulation (overExpressedMembers - underExpressedMembers, -1/0/1), label =
  top risk.class1_proportion by GuanRank within the training cohort. generatePredictor picks the
  best of risk.iterations random splits by test AUC, so its internal test AUC is optimistic; the
  honest estimate is the external cohort, which the model never sees. External evaluation:
  MINER riskStratification (integrated AUC over early-event cutoffs, as in MINER), Harrell's C of the
  predicted probability, Cox HR of predicted high- vs low-risk, and the same HR adjusted for stage
  (BCLC for CLCA, AJCC for TCGA). The trained model also scores LICA-FR (no survival) for step 07.

Part C - MINER ridge risk model (optimize_ridge_model): features = program activity (primary) or
  regulon eigengenes, standardized with training-cohort statistics; label = top class1_proportion vs
  bottom (1 - train_cut_low) by GuanRank in the training cohort; alpha chosen by MINER's bootstrap
  out-of-bag AUC within the training cohort. External evaluation: C-index, Cox HR per SD, HR of the
  within-cohort top class1_proportion, stage-adjusted HR per SD, MINER integrated AUC.

Outputs (results/06_risk/<matrix>/):
  prognostic_{regulons,programs,states}_<EP>.tsv
  predictor_<train>_<EP>/  xgboost: model.json, predictions.tsv (all cohorts), evaluation.tsv, MINER PDFs
  predictor_ridge_<features>_<train>_<EP>/  weights.tsv, predictions.tsv, evaluation.tsv, alpha_optimization.pdf
  risk_summary.tsv, qc/r1_*.png, qc/r2_*.png
"""

import argparse
import json
import os
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

from hcc_common import load_params, p, setup_logging

warnings.filterwarnings("ignore")
OUT = "06_risk"
_G = {}


def bh(pv):
    pv = np.asarray(pv, float)
    q = np.full_like(pv, np.nan)
    ok = ~np.isnan(pv)
    n = ok.sum()
    if n == 0:
        return q
    o = np.argsort(pv[ok])
    r = pv[ok][o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(n)
    qq[o] = np.minimum(r, 1)
    q[ok] = qq
    return q


def cox_z(x, dur, obs):
    """Cox z and p for one covariate (lifelines; MINER's survivalMedianAnalysisDirect uses z)."""
    from lifelines import CoxPHFitter
    d = pd.DataFrame({"x": x, "duration": dur, "observed": obs})
    try:
        s = CoxPHFitter().fit(d, "duration", "observed").summary
        return s.loc["x", "z"], s.loc["x", "p"], s.loc["x", "exp(coef)"]
    except Exception:  # noqa: BLE001 - separation / no variance
        return np.nan, np.nan, np.nan


def _cox_task(key):
    X, srv = _G["X"], _G["srv"]
    x = X.loc[key, srv.index].values.astype(float)
    if x.std() == 0:
        return key, np.nan, np.nan, np.nan
    x = (x - x.mean()) / x.std() if _G["standardize"] else x
    return (key, *cox_z(x, srv["duration"].values, srv["observed"].values))


def per_cohort_cox(X, srv_by_cohort, standardize, cores, log, what):
    """Cox per row of X within each cohort; Stouffer meta over cohorts."""
    out = pd.DataFrame(index=X.index)
    zs, ws = [], []
    for c, srv in srv_by_cohort.items():
        srv = srv.loc[srv.index.intersection(X.columns)]
        _G.update(X=X, srv=srv, standardize=standardize)
        with Pool(cores) as pool:
            res = pool.map(_cox_task, list(X.index))
        r = pd.DataFrame(res, columns=["key", "z", "p", "hr"]).set_index("key")
        out[f"z_{c}"], out[f"p_{c}"], out[f"hr_{c}"] = r["z"], r["p"], r["hr"]
        out[f"q_{c}"] = bh(r["p"].values)
        zs.append(r["z"])
        ws.append(np.sqrt(srv["observed"].sum()))
        log.info("%s Cox in %s: n=%d, events=%d; q<=0.1: %d of %d", what, c, len(srv), int(srv["observed"].sum()),
                 int((out[f"q_{c}"] <= 0.1).sum()), len(X))
    Z = pd.concat(zs, axis=1).values
    w = np.array(ws)
    ok = ~np.isnan(Z).any(axis=1)
    meta = np.full(len(X), np.nan)
    meta[ok] = (Z[ok] * w).sum(1) / np.sqrt((w ** 2).sum())
    out["meta_z"] = meta
    out["meta_p"] = 2 * stats.norm.sf(np.abs(meta))
    out["meta_q"] = bh(out["meta_p"].values)
    out["consistent"] = ok & (np.sign(Z).min(axis=1) == np.sign(Z).max(axis=1))
    out["direction"] = np.where(out["meta_z"] > 0, "adverse", "protective")
    return out.sort_values("meta_p")


def stage_ordinal(df, cohort):
    if cohort == "CLCA":
        return df["BCLC"].astype(str).str.strip().map({"0": 0, "A": 1, "B": 2, "C": 3, "D": 4}).astype(float)
    if cohort == "TCGA":
        return df["stage"].astype(str).str.extract(r"STAGE (I{1,3}V?|IV)")[0].map(
            {"I": 1, "II": 2, "III": 3, "IV": 4}).astype(float)
    return pd.Series(np.nan, index=df.index)


def evaluate_external(clf, mtrx, guan, clin, cohort, tag, outdir, log):
    """MINER riskStratification + C-index + (stage-adjusted) Cox HR of predicted class."""
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    from miner import miner
    lbls = clf.predict(np.array(mtrx.T))
    prob = pd.Series(clf.predict_proba(np.array(mtrx.T))[:, 1], index=mtrx.columns)
    aucs, cutoffs, _, _, miner_hr, _, _ = miner.riskStratification(
        lbls, mtrx, guan, tag, clf, guan_rank=False, resultsDirectory=None, plot_all=False, plot_any=True)
    import matplotlib.pyplot as plt
    plt.savefig(os.path.join(outdir, f"{tag}_external_{cohort}.pdf"), bbox_inches="tight")
    plt.close("all")
    srv = guan[["duration", "observed"]].loc[mtrx.columns.intersection(guan.index)]
    r = {"cohort": cohort, "n": len(srv), "events": int(srv["observed"].sum()),
         "pct_predicted_high": 100 * lbls.mean(),
         "miner_integrated_auc": float(np.mean(aucs)), "miner_hr_statistic": float(miner_hr),
         "c_index": concordance_index(srv["duration"], -prob[srv.index], srv["observed"])}
    d = srv.assign(high=pd.Series(lbls, index=mtrx.columns)[srv.index].astype(int))
    if d["high"].nunique() == 2:
        s = CoxPHFitter().fit(d, "duration", "observed").summary.loc["high"]
        r.update(hr=s["exp(coef)"], hr_lo=s["exp(coef) lower 95%"], hr_hi=s["exp(coef) upper 95%"], hr_p=s["p"])
        st = stage_ordinal(clin.reindex(srv.index), cohort)
        da = d.assign(stage=st).dropna()
        if da["stage"].nunique() > 1 and da["high"].nunique() == 2:
            sa = CoxPHFitter().fit(da, "duration", "observed").summary
            r.update(hr_stage_adj=sa.loc["high", "exp(coef)"], hr_stage_adj_lo=sa.loc["high", "exp(coef) lower 95%"],
                     hr_stage_adj_hi=sa.loc["high", "exp(coef) upper 95%"], hr_stage_adj_p=sa.loc["high", "p"],
                     n_stage_adj=len(da))
    log.info("%s -> %s: n=%d events=%d, C=%.3f, MINER iAUC=%.3f, HR high vs low %.2f (%.2f-%.2f) p=%.2g; "
             "stage-adjusted HR %.2f p=%.2g", tag, cohort, r["n"], r["events"], r["c_index"],
             r["miner_integrated_auc"], r.get("hr", np.nan), r.get("hr_lo", np.nan), r.get("hr_hi", np.nan),
             r.get("hr_p", np.nan), r.get("hr_stage_adj", np.nan), r.get("hr_stage_adj_p", np.nan))
    return r, prob, pd.Series(lbls, index=mtrx.columns)


class _ScoreClassifier:
    """Adapter so MINER's riskStratification can evaluate a continuous risk score."""
    def __init__(self, score):
        self.s = (score - score.min()) / (score.max() - score.min())

    def predict_proba(self, X):
        v = self.s.values
        return np.column_stack([1 - v, v])


def evaluate_score(score, srv, guan, clin, cohort, tag, outdir, top, log):
    """External evaluation of a continuous risk score within one cohort."""
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    from miner import miner
    s = score[srv.index]
    z = (s - s.mean()) / s.std()
    high = (s >= s.quantile(1 - top)).astype(int)
    r = {"cohort": cohort, "n": len(srv), "events": int(srv["observed"].sum()),
         "c_index": concordance_index(srv["duration"], -s, srv["observed"])}
    d = srv.assign(x=z)
    c = CoxPHFitter().fit(d, "duration", "observed").summary.loc["x"]
    r.update(hr_per_sd=c["exp(coef)"], hr_per_sd_lo=c["exp(coef) lower 95%"], hr_per_sd_hi=c["exp(coef) upper 95%"],
             hr_per_sd_p=c["p"])
    c = CoxPHFitter().fit(srv.assign(high=high), "duration", "observed").summary.loc["high"]
    r.update(hr=c["exp(coef)"], hr_lo=c["exp(coef) lower 95%"], hr_hi=c["exp(coef) upper 95%"], hr_p=c["p"])
    st = stage_ordinal(clin.reindex(srv.index), cohort)
    da = d.assign(stage=st).dropna()
    if da["stage"].nunique() > 1:
        c = CoxPHFitter().fit(da, "duration", "observed").summary.loc["x"]
        r.update(hr_per_sd_stage_adj=c["exp(coef)"], hr_per_sd_stage_adj_lo=c["exp(coef) lower 95%"],
                 hr_per_sd_stage_adj_hi=c["exp(coef) upper 95%"], hr_per_sd_stage_adj_p=c["p"], n_stage_adj=len(da))
    mt = pd.DataFrame([s.values], columns=s.index)  # riskStratification only needs the sample columns
    aucs, *_ = miner.riskStratification(high.values, mt, guan.loc[s.index], tag, _ScoreClassifier(s),
                                        guan_rank=False, resultsDirectory=None, plot_all=False, plot_any=False)
    r["miner_integrated_auc"] = float(np.mean(aucs))
    log.info("%s -> %s: n=%d events=%d, C=%.3f, HR/SD %.2f (%.2f-%.2f) p=%.1e, top %.0f%% HR %.2f p=%.1e, "
             "stage-adj HR/SD %.2f p=%.1e, MINER iAUC %.3f", tag, cohort, r["n"], r["events"], r["c_index"],
             r["hr_per_sd"], r["hr_per_sd_lo"], r["hr_per_sd_hi"], r["hr_per_sd_p"], 100 * top, r["hr"], r["hr_p"],
             r.get("hr_per_sd_stage_adj", np.nan), r.get("hr_per_sd_stage_adj_p", np.nan), r["miner_integrated_auc"])
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    ap.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--skip", default="", help="comma list: prognostic,predictor,ridge")
    args = ap.parse_args()
    P = load_params(args.params)
    M, R = P["miner"], P["risk"]
    matrix = args.matrix or M["matrix"]
    res = p(P["paths"]["results"])
    sdir = os.path.join(res, "04_miner", matrix, R["subtypes_dir"])
    outdir = os.path.join(res, OUT, matrix)
    log = setup_logging(outdir, "06_risk")
    from miner import miner

    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    E = pd.read_csv(os.path.join(sdir, "eigengenes.csv"), index_col=0)
    E.index = E.index.astype(str)
    over = pd.read_csv(os.path.join(sdir, "overExpressedMembers.csv"), index_col=0)
    under = pd.read_csv(os.path.join(sdir, "underExpressedMembers.csv"), index_col=0)
    over.index, under.index = over.index.astype(str), under.index.astype(str)
    diff = (over - under).astype(int)
    progs = json.load(open(os.path.join(sdir, "transcriptional_programs.json")))
    states = json.load(open(os.path.join(sdir, "transcriptional_states.json")))
    log.info("Regulons %d, programs %d, states %d, samples %d", E.shape[0], len(progs), len(states), E.shape[1])

    # survival per cohort and endpoint (sample-indexed MINER files, common horizon)
    sfx = R["survival_suffix"]
    surv = {}
    for ep in R["endpoints"]:
        for c in R["survival_cohorts"]:
            f = os.path.join(res, "03_genomics_clinical", f"survival_{c}_{ep}{sfx}_miner.csv")
            s = pd.read_csv(f, index_col=0)
            s.columns = ["duration", "observed"]
            surv[(c, ep)] = s.loc[s.index.intersection(E.columns)]
    clin = {c: pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_{c}.tsv"), sep="\t", index_col=0)
            for c in R["survival_cohorts"]}
    # clinical tables are patient-indexed; map to expression samples
    for c in clin:
        pt = samples.loc[samples["cohort"] == c, "patient"]
        clin[c] = clin[c].reindex(pt.values).set_axis(pt.index)

    summary = []
    # ---------------- Part A: prognostic regulons, programs, states
    if "prognostic" not in args.skip:
        fam = pd.read_csv(os.path.join(res, "05_causal", matrix, "regulon_families.tsv"), sep="\t", index_col=0)
        fam.index = fam.index.astype(str)
        Ez = E.sub(E.mean(1), axis=0)
        prog_act = pd.DataFrame({k: Ez.loc[[r for r in v if r in Ez.index]].mean() for k, v in progs.items()}).T
        state_mem = pd.DataFrame({k: E.columns.isin(v).astype(int) for k, v in states.items()}, index=E.columns).T
        for ep in R["endpoints"]:
            sbc = {c: surv[(c, ep)] for c in R["survival_cohorts"]}
            for name, X, std in (("regulons", E, True), ("programs", prog_act, True), ("states", state_mem, False)):
                t = per_cohort_cox(X, sbc, std, args.cores, log, f"{name} {ep}")
                if name == "regulons":
                    t = t.join(fam[["family", "program", "regulator_symbol", "n_genes"]])
                if name == "states":
                    t["n_samples"] = state_mem.sum(1)
                t.index.name = name[:-1]
                t.to_csv(os.path.join(outdir, f"prognostic_{name}{sfx}_{ep}.tsv"), sep="\t")
                sig = t[(t["meta_q"] <= R["q_max"]) & t["consistent"]]
                log.info("%s %s: %d of %d prognostic (meta q <= %.2f, same sign in every cohort): %d adverse, %d protective",
                         name, ep, len(sig), len(t), R["q_max"], (sig["direction"] == "adverse").sum(),
                         (sig["direction"] == "protective").sum())
                summary.append({"part": "prognostic", "unit": name, "endpoint": ep, "n_tested": len(t),
                                "n_prognostic": len(sig), "adverse": int((sig["direction"] == "adverse").sum()),
                                "protective": int((sig["direction"] == "protective").sum())})

    # ---------------- Part B: MINER predictor, train in one cohort, test in the others
    if "predictor" not in args.skip:
        for tr in R["train_test"]:
            train, ep = tr["train"], tr["endpoint"]
            tag = f"{train}_{ep}{sfx}"
            pdir = os.path.join(outdir, f"predictor_{tag}")
            os.makedirs(pdir, exist_ok=True)
            srv_tr = surv[(train, ep)]
            guan_tr = miner.guanRank(miner.kmAnalysis(srv_tr.copy(), "duration", "observed"))
            mem_tr = diff[srv_tr.index]
            log.info("Training %s: %d samples, %d events, top %.0f%% by GuanRank = high risk",
                     tag, len(srv_tr), int(srv_tr["observed"].sum()), 100 * R["class1_proportion"])
            clf, _, _, mean_aucs, mean_hrs, _, _, _ = miner.generatePredictor(
                [mem_tr], [guan_tr], [tag], iterations=R["iterations"], method=R["method"],
                n_estimators=R["n_estimators"], output_directory=pdir, test_only=True, metric="roc_auc",
                class1_proportion=R["class1_proportion"], test_proportion=R["test_proportion"])
            import matplotlib.pyplot as plt
            plt.close("all")
            if R["method"] == "xgboost":
                clf.save_model(os.path.join(pdir, "model.json"))
            rows = [{"cohort": f"{train} (internal test, best of {R['iterations']} splits; optimistic)",
                     "miner_integrated_auc": float(np.max(mean_aucs)) if len(mean_aucs) else np.nan}]
            preds = []
            for c in tr["test"]:
                srv_te = surv[(c, ep)]
                guan_te = miner.guanRank(miner.kmAnalysis(srv_te.copy(), "duration", "observed"))
                r, prob, lbl = evaluate_external(clf, diff[srv_te.index], guan_te, clin.get(c, pd.DataFrame()), c,
                                                 tag, pdir, log)
                rows.append(r)
            # score every network sample (LICA-FR has no survival; used in step 07)
            prob_all = pd.Series(clf.predict_proba(np.array(diff.T))[:, 1], index=diff.columns)
            pred = pd.DataFrame({"cohort": samples["cohort"].reindex(diff.columns), "prob_high_risk": prob_all,
                                 "predicted_high_risk": clf.predict(np.array(diff.T)).astype(int)})
            pred["used_for_training"] = pred.index.isin(srv_tr.index)
            pred.to_csv(os.path.join(pdir, "predictions.tsv"), sep="\t")
            ev = pd.DataFrame(rows)
            ev.to_csv(os.path.join(pdir, "evaluation.tsv"), sep="\t", index=False)
            # features the model uses
            imp = pd.Series(clf.feature_importances_, index=diff.index).sort_values(ascending=False)
            imp = imp[imp > 0].to_frame("importance")
            imp.to_csv(os.path.join(pdir, "feature_importance.tsv"), sep="\t")
            log.info("%s: model uses %d regulons", tag, len(imp))
            for r in rows[1:]:
                summary.append({"part": "predictor", "unit": tag, "endpoint": ep, **{k: r[k] for k in r}})

    # ---------------- Part C: MINER ridge risk model (optimize_ridge_model) on program / regulon activity
    if "ridge" not in args.skip:
        from sklearn.linear_model import Ridge
        Ez = E.sub(E.mean(1), axis=0)
        feats = {"programs": pd.DataFrame({k: Ez.loc[[r for r in v if r in Ez.index]].mean()
                                           for k, v in progs.items()}).T,
                 "regulons": E}
        for fname in R["ridge"]["features"]:
            X = feats[fname]
            for tr in R["train_test"]:
                train, ep = tr["train"], tr["endpoint"]
                tag = f"ridge_{fname}_{train}_{ep}{sfx}"
                pdir = os.path.join(outdir, f"predictor_{tag}")
                os.makedirs(pdir, exist_ok=True)
                srv_tr = surv[(train, ep)]
                guan_tr = miner.guanRank(miner.kmAnalysis(srv_tr.copy(), "duration", "observed"))
                mu, sd = X[srv_tr.index].mean(1), X[srv_tr.index].std(1).replace(0, 1)
                Xs = X.sub(mu, axis=0).div(sd, axis=0)          # standardized with training-cohort statistics
                np.random.seed(R["ridge"]["seed"])
                alpha, mean_auc, _, a_range = miner.optimize_ridge_model(
                    Xs[guan_tr.index], guan_tr, n_iter=R["ridge"]["n_iter"],
                    train_cut_high=R["class1_proportion"], train_cut_low=R["ridge"]["train_cut_low"],
                    max_range=R["ridge"]["max_alpha"], range_step=R["ridge"]["alpha_step"],
                    savefile=os.path.join(pdir, "alpha_optimization.pdf"))
                import matplotlib.pyplot as plt
                plt.close("all")
                n = len(guan_tr)
                hr_ids = guan_tr.index[:round(n * R["class1_proportion"])]
                lr_ids = guan_tr.index[round(n * R["ridge"]["train_cut_low"]):]
                y = np.r_[np.ones(len(hr_ids)), np.zeros(len(lr_ids))]
                model = Ridge(alpha=alpha, random_state=0).fit(Xs[list(hr_ids) + list(lr_ids)].T.values, y)
                score = pd.Series(model.predict(Xs.T.values), index=Xs.columns)
                log.info("%s: alpha %d (bootstrap OOB AUC %.3f in %s)", tag, alpha, mean_auc.max(), train)
                rows = [{"cohort": f"{train} (training, apparent)",
                         **evaluate_score(score, srv_tr, guan_tr, clin.get(train, pd.DataFrame()), train, tag, pdir,
                                          R["class1_proportion"], log)}]
                for c in tr["test"]:
                    srv_te = surv[(c, ep)]
                    guan_te = miner.guanRank(miner.kmAnalysis(srv_te.copy(), "duration", "observed"))
                    r = evaluate_score(score, srv_te, guan_te, clin.get(c, pd.DataFrame()), c, tag, pdir,
                                       R["class1_proportion"], log)
                    rows.append(r)
                    summary.append({"part": "ridge", "unit": tag, "endpoint": ep, "alpha": alpha, **r})
                pd.DataFrame(rows).to_csv(os.path.join(pdir, "evaluation.tsv"), sep="\t", index=False)
                pd.Series(model.coef_, index=X.index, name="weight").sort_values().to_csv(
                    os.path.join(pdir, "weights.tsv"), sep="\t")
                pred = pd.DataFrame({"cohort": samples["cohort"].reindex(score.index), "risk_score": score})
                pred["used_for_training"] = pred.index.isin(srv_tr.index)
                for c in pred["cohort"].dropna().unique():   # within-cohort top fraction = predicted high risk
                    m = pred["cohort"] == c
                    pred.loc[m, "predicted_high_risk"] = (pred.loc[m, "risk_score"] >= pred.loc[m, "risk_score"].quantile(
                        1 - R["class1_proportion"])).astype(int)
                pred.to_csv(os.path.join(pdir, "predictions.tsv"), sep="\t")

    if summary:
        pd.DataFrame(summary).to_csv(os.path.join(outdir, "risk_summary.tsv"), sep="\t", index=False)
    try:
        import qc_plots
        written = qc_plots.risk_report(outdir, R, surv)
        log.info("QC figures: %s", ", ".join(written))
    except Exception as e:  # noqa: BLE001
        log.warning("figures failed: %s", e)


if __name__ == "__main__":
    main()
