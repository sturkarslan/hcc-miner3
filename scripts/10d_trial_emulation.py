#!/usr/bin/env python
"""Step 10d: emulate published HCC trials with synthetic cohorts drawn from the discovery tumours, predict each arm's
objective response rate (ORR) from drug-constrained regulon activity (DCRA, step 10b), and compare with the observed ORR.

Trials: config/clinical_trials.yaml; per-arm ORR, denominator, mean / median age, % female, % Asian from ClinicalTrials.gov
posted results (scripts/tools/fetch_trials.py). Primary ORR measure per arm: RECIST 1.1 by independent review in the global
population when posted, else the trial's single posted ORR (its criteria are recorded). Placebo arms, regional / PD-L1 /
expansion sub-cohorts and safety run-ins are excluded. Count units are converted to percentages.

Synthetic cohorts: entropy-balancing weights over the discovery patients so that the weighted % female, % Asian and mean age
equal the arm's baseline (only the characteristics the trial posted); then n_cohorts cohorts of the arm's size are drawn
with replacement using those weights. Pool: all discovery tumours with DCRA, or (pool: advanced) only stage III-IV /
BCLC B-C tumours, closer to the trials' eligibility (unresectable / advanced).

Patient response: DCRA of the arm's drugs. Monotherapy: responder if DCRA > tau. Combination: independent drug action
(responder if any drug's DCRA > tau). Drugs with no mapped regulon are skipped; an arm with none is not predicted.

Thresholds (tau), from the published ORR:
  A. per-arm tau that reproduces each arm's own ORR (the calibration itself; circular, reported only to show how tau
     varies across drugs and across trials of the same drug)
  B. leave-one-trial-out global tau: one tau shared by all drugs, fitted on the other trials' arms (least squares on ORR),
     then the held-out trial's arms are predicted -> the real test
  C. same-drug transfer: tau for a drug fitted on its arms in the other trials, applied to this trial's arm
Comparison: predicted ORR (mean and 95% interval over the synthetic cohorts) vs observed ORR (exact binomial 95% CI);
Pearson / Spearman across arms, mean absolute error, coverage. Null for B: drug labels permuted across arms (same taus).
Outputs: results/10_response/trials/{emulation_arms.tsv, emulation_predictions.tsv, emulation_summary.tsv} + figure
"""

import os
import re

import numpy as np
import pandas as pd
import yaml
from scipy import optimize, stats

from hcc_common import load_params, p, setup_logging

EXCLUDE = re.compile(r"placebo|safety run-in|mee cohort|open label|china|pd-l1", re.I)


def pick_orr(o):
    """One ORR per arm: RECIST 1.1 independent review, global population; else the only posted value."""
    o = o[~o["measure"].str.contains("China|PD-L1|Duration", case=False)]
    o = o[~o["group"].str.contains("China", case=False)]
    out = []
    for (trial, grp), g in o.groupby(["trial", "group"]):
        def score(r):
            t = (r["measure"] + " " + str(r["class"])).lower()
            return (("mrecist" in t or "modified" in t) * -2 + ("investigator" in t or "orr-inv" in t) * -1 +
                    ("irf" in t or "birc" in t or "bicr" in t or "independent" in t) * 1 + ("recist 1.1" in t or "recist v1.1" in t) * 1)
        g = g.assign(s=g.apply(score, axis=1)).sort_values("s", ascending=False)
        r = g.iloc[0]
        v = r["value"] / r["n"] * 100 if str(r["unit"]).lower() == "participants" else r["value"]
        out.append({"trial": trial, "group": grp, "orr": v / 100, "orr_n": r["n"], "orr_measure": (r["measure"] + " " + str(r["class"])).strip()})
    return pd.DataFrame(out)


def entropy_weights(X, target):
    """Weights w (sum 1) minimizing KL to uniform subject to X' w = target (columns already centred on target)."""
    Z = (X - target).values
    def dual(lam):
        e = np.exp(Z @ lam)
        return np.log(e.mean())
    def grad(lam):
        e = np.exp(Z @ lam)
        return (Z * e[:, None]).sum(0) / e.sum()
    lam = optimize.minimize(dual, np.zeros(Z.shape[1]), jac=grad, method="BFGS").x
    w = np.exp(Z @ lam)
    return w / w.sum(), np.abs(grad(lam)).max()


def main():
    P = load_params(None)
    res = p(P["paths"]["results"])
    outdir = os.path.join(res, "10_response", "trials")
    log = setup_logging(outdir, "10d_trial_emulation")
    cfg = yaml.safe_load(open(p("config/clinical_trials.yaml")))
    E = cfg.get("emulation", {})
    ncoh = int(E.get("n_cohorts", 1000))
    rng = np.random.default_rng(int(E.get("seed", 7)))

    # ---------------- trial arms
    A = pd.read_csv(os.path.join(outdir, "trial_arms.tsv"), sep="\t")
    O = pick_orr(pd.read_csv(os.path.join(outdir, "orr_measures.tsv"), sep="\t"))
    A = A.merge(O, on=["trial", "group"], how="inner")
    A = A[~A["group"].str.contains(EXCLUDE)]
    A.loc[A["female_frac"] > 0.6, "female_frac"] = np.nan           # implausible posted value (Asia-Pacific): not used
    drugs = []
    for _, r in A.iterrows():
        t = next(t for t in cfg["trials"] if t["name"] == r["trial"])
        m = next((v for k, v in t["arms"].items() if re.search(k, r["group"], re.I)), None)
        drugs.append(m)
    A["drugs"] = drugs
    A = A.dropna(subset=["drugs", "orr"]).reset_index(drop=True)

    # ---------------- patient pool: DCRA + covariates
    dc = pd.read_csv(os.path.join(res, "10_response", "dcna", "dcra_discovery.tsv"), sep="\t", index_col=0)
    sm = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    cov = []
    t = pd.read_csv(os.path.join(res, "03_genomics_clinical", "survival_TCGA.tsv"), sep="\t", index_col=0)
    tc = pd.read_csv(p("data/LIHC-TCGA-PanCancerAtlas_2018/data_clinical_patient.txt"), sep="\t", comment="#", index_col=0)
    for s in sm.index[sm["cohort"] == "TCGA"]:
        pt = sm.loc[s, "patient"]
        st = str(t["stage"].get(pt, ""))
        cov.append({"sample": s, "age": t["age"].get(pt), "female": float(str(t["sex"].get(pt, "")).lower().startswith("f")),
                    "asian": float(str(tc["RACE"].get(pt, "")).upper() == "ASIAN") if pt in tc.index else np.nan,
                    "advanced": float("III" in st or "IV" in st) if st not in ("", "nan") else np.nan})
    c = pd.read_csv(os.path.join(res, "03_genomics_clinical", "survival_CLCA.tsv"), sep="\t", index_col=0)
    for s in sm.index[sm["cohort"] == "CLCA"]:
        pt = sm.loc[s, "patient"]
        cov.append({"sample": s, "age": c["age"].get(pt), "female": 1.0 - c["sex_male"].get(pt, np.nan), "asian": 1.0,
                    "advanced": float(str(c["BCLC"].get(pt)) in ("B", "C"))})
    la = pd.read_excel(p("data/LICA-FR/clinical_molecular_annotations.xlsx"))
    la["sample"] = la["Sample"].astype(str).str.replace("^#", "CHC", regex=True)
    la = la.set_index("sample")
    for s in sm.index[sm["cohort"] == "LICA_FR"]:
        if s in la.index:
            r = la.loc[s] if la.loc[[s]].shape[0] == 1 else la.loc[[s]].iloc[0]
            cov.append({"sample": s, "age": pd.to_numeric(r["Age_at_Sampling"], errors="coerce"),
                        "female": float(str(r["Gender"]).upper().startswith("F")), "asian": 0.0,
                        "advanced": float(str(r["BCLC"]).strip() in ("B", "C")) if pd.notna(r["BCLC"]) else np.nan})
    C = pd.DataFrame(cov).set_index("sample")
    C = C.loc[C.index.intersection(dc.columns)]
    C["age"] = pd.to_numeric(C["age"], errors="coerce")
    log.info("Pool: %d tumours with DCRA and covariates; female %.2f, Asian %.2f, mean age %.1f, advanced %.2f (n %d)", len(C),
             C["female"].mean(), C["asian"].mean(), C["age"].mean(), C["advanced"].mean(), int(C["advanced"].sum()))

    def weights(pool, r):
        cols, target = [], []
        for col, key in (("female", "female_frac"), ("asian", "asian_frac"), ("age", "age")):
            if pd.notna(r[key]):
                cols.append(col)
                target.append(r[key])
        X = pool[cols].astype(float)
        ok = X.notna().all(1)
        w = np.zeros(len(pool))
        if not cols:
            w[:] = 1 / len(pool)
            return w, 0.0, cols
        ww, err = entropy_weights(X[ok], np.array(target))
        w[ok.values] = ww
        return w, err, cols

    # ---------------- simulate
    pools = {"all": C.dropna(subset=["female", "age"]), "advanced": C[(C["advanced"] == 1)].dropna(subset=["female", "age"])}
    taus = np.round(np.linspace(-1.0, 1.0, 401), 3)          # DCRA (mean trinary activity) lies in [-1, 1]
    sims = {}
    rows = []
    for pname, pool in pools.items():
        for i, r in A.iterrows():
            dr = [d for d in r["drugs"] if d in dc.index]
            if not dr:
                continue
            w, err, cols = weights(pool, r)
            n = int(r["orr_n"])
            idx = rng.choice(len(pool), size=(ncoh, n), p=w)               # synthetic cohorts
            V = dc.loc[dr, pool.index].values.max(axis=0)                     # independent action: best drug
            resp = (V[idx][:, :, None] > taus[None, None, :]).mean(axis=1)   # cohorts x taus
            sims[(pname, i)] = resp
            ess = 1 / np.sum(w ** 2)
            rows.append({"pool": pname, "arm": i, "trial": r["trial"], "group": r["group"], "drugs": "+".join(dr), "n": n, "orr": r["orr"],
                         "balance_vars": ",".join(cols), "balance_max_err": err, "effective_n": ess})
    AR = pd.DataFrame(rows)
    A.assign(drugs=A["drugs"].map("+".join)).to_csv(os.path.join(outdir, "emulation_arms.tsv"), sep="\t", float_format="%.4g")

    def mean_curve(pname, i):
        return sims[(pname, i)].mean(0)

    preds = []
    for pname in pools:
        ar = AR[AR["pool"] == pname]
        # A. per-arm tau
        for _, r in ar.iterrows():
            mc = mean_curve(pname, r["arm"])
            k = int(np.argmin(np.abs(mc - r["orr"])))
            preds.append({"pool": pname, "mode": "A per-arm tau (calibration)", "trial": r["trial"], "group": r["group"], "drugs": r["drugs"],
                          "orr": r["orr"], "n": r["n"], "tau": taus[k], "pred": mc[k], "lo": np.percentile(sims[(pname, r["arm"])][:, k], 2.5),
                          "hi": np.percentile(sims[(pname, r["arm"])][:, k], 97.5)})
        # B. leave-one-trial-out global tau; C. same-drug transfer
        for trial in ar["trial"].unique():
            tr, te = ar[ar["trial"] != trial], ar[ar["trial"] == trial]
            sse = sum((mean_curve(pname, a) - o) ** 2 for a, o in zip(tr["arm"], tr["orr"]))
            k = int(np.argmin(sse))
            naive_all = tr["orr"].mean()
            for _, r in te.iterrows():
                s = sims[(pname, r["arm"])][:, k]
                preds.append({"pool": pname, "mode": "B leave-one-trial-out global tau", "trial": trial, "group": r["group"], "drugs": r["drugs"],
                              "orr": r["orr"], "n": r["n"], "tau": taus[k], "pred": s.mean(), "lo": np.percentile(s, 2.5), "hi": np.percentile(s, 97.5),
                              "naive": naive_all})
                same = tr[tr["drugs"] == r["drugs"]]
                if len(same):
                    sse2 = sum((mean_curve(pname, a) - o) ** 2 for a, o in zip(same["arm"], same["orr"]))
                    k2 = int(np.argmin(sse2))
                    s2 = sims[(pname, r["arm"])][:, k2]
                    preds.append({"pool": pname, "mode": "C same-drug transfer", "trial": trial, "group": r["group"], "drugs": r["drugs"],
                                  "orr": r["orr"], "n": r["n"], "tau": taus[k2], "pred": s2.mean(), "lo": np.percentile(s2, 2.5),
                                  "hi": np.percentile(s2, 97.5), "n_train_arms": len(same), "naive": same["orr"].mean()})
        # null for B: permute drug assignment across arms (each arm simulated with another arm's drugs' DCRA), same procedure
        null_r = []
        for _ in range(200):
            perm = rng.permutation(ar["arm"].values)
            mp = dict(zip(ar["arm"], perm))
            pr, ob = [], []
            for trial in ar["trial"].unique():
                tr, te = ar[ar["trial"] != trial], ar[ar["trial"] == trial]
                sse = sum((mean_curve(pname, mp[a]) - o) ** 2 for a, o in zip(tr["arm"], tr["orr"]))
                k = int(np.argmin(sse))
                for _, r in te.iterrows():
                    pr.append(mean_curve(pname, mp[r["arm"]])[k])
                    ob.append(r["orr"])
            null_r.append(stats.pearsonr(pr, ob).statistic if np.std(pr) > 0 else 0)
        AR.loc[AR["pool"] == pname, "null_r_mean"] = np.mean(null_r)
        AR.loc[AR["pool"] == pname, "null_r_95"] = np.percentile(null_r, 95)
    PR = pd.DataFrame(preds)
    ci = PR.apply(lambda r: stats.binomtest(int(round(r["orr"] * r["n"])), int(r["n"])).proportion_ci(), axis=1)
    PR["orr_lo"], PR["orr_hi"] = [c.low for c in ci], [c.high for c in ci]
    PR["covered"] = (PR["orr"] >= PR["lo"]) & (PR["orr"] <= PR["hi"])
    PR.to_csv(os.path.join(outdir, "emulation_predictions.tsv"), sep="\t", index=False, float_format="%.4g")
    summ = []
    for (pname, mode), g in PR.groupby(["pool", "mode"]):
        summ.append({"pool": pname, "mode": mode, "arms": len(g), "pearson": stats.pearsonr(g["pred"], g["orr"]).statistic if g["pred"].std() > 0 else np.nan,
                     "spearman": stats.spearmanr(g["pred"], g["orr"]).statistic, "mae": (g["pred"] - g["orr"]).abs().mean(),
                     "coverage": g["covered"].mean(), "tau_range": f"{g['tau'].min():.2f} to {g['tau'].max():.2f}",
                     "naive_mae": (g["naive"] - g["orr"]).abs().mean() if "naive" in g and g["naive"].notna().any() else np.nan,
                     "null_pearson_mean": AR.loc[AR["pool"] == pname, "null_r_mean"].iloc[0] if mode.startswith("B") else np.nan,
                     "null_pearson_95": AR.loc[AR["pool"] == pname, "null_r_95"].iloc[0] if mode.startswith("B") else np.nan})
    S = pd.DataFrame(summ)
    S.to_csv(os.path.join(outdir, "emulation_summary.tsv"), sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 250, "display.max_rows", 200, "display.max_colwidth", 45):
        log.info("Arms: %d (balance on the posted characteristics; effective n of the weighted pool):\n%s", AR["arm"].nunique(),
                 AR[AR["pool"] == "all"][["trial", "group", "drugs", "n", "orr", "balance_vars", "effective_n"]].round(3).to_string(index=False))
        log.info("Per-arm tau (mode A, pool all):\n%s", PR[(PR["pool"] == "all") & PR["mode"].str.startswith("A")][["trial", "drugs", "orr", "tau"]].round(3).to_string(index=False))
        log.info("Predictions (pool all, modes B / C):\n%s", PR[(PR["pool"] == "all") & ~PR["mode"].str.startswith("A")][
            ["mode", "trial", "drugs", "orr", "pred", "lo", "hi", "tau", "covered"]].round(3).to_string(index=False))
        log.info("Summary:\n%s", S.round(3).to_string(index=False))

    # ---------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 6, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), layout="constrained")
    g = PR[(PR["pool"] == "all") & PR["mode"].str.startswith("B")]
    ax = axes[0]
    ax.errorbar(100 * g["orr"], 100 * g["pred"], yerr=[100 * (g["pred"] - g["lo"]), 100 * (g["hi"] - g["pred"])],
                xerr=[100 * (g["orr"] - g["orr_lo"]), 100 * (g["orr_hi"] - g["orr"])], fmt="o", ms=3, lw=0.5, color="#2a78d6")
    for _, r in g.iterrows():
        ax.annotate(f"{r['trial']}: {r['drugs'].split('+')[0]}", (100 * r["orr"], 100 * r["pred"]), fontsize=3.8, xytext=(2, 2), textcoords="offset points")
    m = 100 * max(g["orr_hi"].max(), g["hi"].max()) + 2
    ax.plot([0, m], [0, m], color="grey", lw=0.5, ls="--")
    s = S[(S["pool"] == "all") & S["mode"].str.startswith("B")].iloc[0]
    ax.set_title(f"leave-one-trial-out global threshold\nr = {s['pearson']:.2f} (null 95% {s['null_pearson_95']:.2f}), MAE {100 * s['mae']:.1f} pts", fontsize=6)
    ax.set_xlabel("observed ORR (%)")
    ax.set_ylabel("predicted ORR (%)")
    ax = axes[1]
    g = PR[(PR["pool"] == "all") & PR["mode"].str.startswith("C")]
    if len(g):
        ax.errorbar(100 * g["orr"], 100 * g["pred"], yerr=[100 * (g["pred"] - g["lo"]), 100 * (g["hi"] - g["pred"])], fmt="o", ms=3, lw=0.5, color="#e34948")
        for _, r in g.iterrows():
            ax.annotate(f"{r['trial']}: {r['drugs']}", (100 * r["orr"], 100 * r["pred"]), fontsize=3.8, xytext=(2, 2), textcoords="offset points")
        m = 100 * max(g["orr"].max(), g["hi"].max()) + 2
        ax.plot([0, m], [0, m], color="grey", lw=0.5, ls="--")
        s = S[(S["pool"] == "all") & S["mode"].str.startswith("C")].iloc[0]
        ax.set_title(f"same-drug threshold transfer\nr = {s['pearson']:.2f}, MAE {100 * s['mae']:.1f} pts, coverage {100 * s['coverage']:.0f}%", fontsize=6)
    ax.set_xlabel("observed ORR (%)")
    ax = axes[2]
    g = PR[(PR["pool"] == "all") & PR["mode"].str.startswith("A")].sort_values("tau")
    ax.scatter(g["tau"], range(len(g)), s=8, color=["#e34948" if "sorafenib" == d else "#2a78d6" for d in g["drugs"]])
    ax.set_yticks(range(len(g)))
    ax.set_yticklabels([f"{t}: {d}" for t, d in zip(g["trial"], g["drugs"])], fontsize=4)
    ax.set_xlabel("DCRA threshold reproducing the arm's ORR")
    ax.set_title("per-arm threshold (calibration only)\nred = sorafenib arms", fontsize=6)
    fig.savefig(os.path.join(res, "10_response", "figures", "trial_emulation.png"), dpi=300)
    fig.savefig(os.path.join(res, "10_response", "figures", "trial_emulation.pdf"))


if __name__ == "__main__":
    main()
