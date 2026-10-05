#!/usr/bin/env python
"""Step 10b: drug-constrained network activity (DCNA), as in the GBM (gbmMINER) drug-response work.

  1. Regulon activity per tumour: MINER generateRegulonActivity (trinary -1 / 0 / +1: regulon genes coherently
     under / not / over-expressed relative to the tumour's own expression distribution, binomial p = 0.05).
     Discovery tumours: MINER's own over- minus under-expressed membership (subtypes_filtered). Other cohorts:
     computed from the within-cohort gene z-scores.
  2. Drug -> regulon set R(D): regulons whose regulator is a drug target, or that contain a drug target as a member
     (Open Targets targets, data/reference/ot_drug_targets.tsv; no expansion).
  3. DCRA(D, patient) = mean trinary activity over R(D) (regulon level).
     DCPA(D, patient) = mean activity over the programs holding any regulon of R(D), program activity = mean trinary
     activity of its regulons (program level).
  4. Prediction: inhibitors / antagonists -> responder if DCNA > 0 (agonists: <= 0).
Tests:
  - GSE109211 (STORM): sorafenib DCNA in the sorafenib arm (placebo arm = control): predicted vs observed responder
    (Fisher), AUC of continuous DCNA, and treatment x DCNA interaction.
  - GSE104580 (TACE): doxorubicin (TOP2A) as the TACE chemotherapy proxy (the agent is not given in GEO).
  - Discovery (929): predicted-responder fraction per drug vs reported trial ORR (monotherapy rows of
    config/soc_regimens.yaml), Spearman; exploratory (resected tumours, not trial populations).
Outputs: results/10_response/dcna/{drug_regulons.tsv, dcna_<cohort>.tsv, dcna_tests.tsv, dcna_vs_orr.tsv}
"""

import importlib.util
import json
import os

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from hcc_common import load_params, miner_id_backmap, p, setup_logging


def main():
    from miner import miner
    P = load_params(None)
    res = p(P["paths"]["results"])
    mx = P["miner"]["matrix"]
    outdir = os.path.join(res, "10_response", "dcna")
    os.makedirs(outdir, exist_ok=True)
    log = setup_logging(outdir, "10b_dcna")
    spec = importlib.util.spec_from_file_location("s10", os.path.join(os.path.dirname(os.path.abspath(__file__)), "10_response_validation.py"))
    s10 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s10)

    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col=0)
    back = miner_id_backmap(p(P["miner"]["idmap"]), genes.index)
    mdir = os.path.join(res, "04_miner", mx)
    regs_miner = json.load(open(os.path.join(mdir, "mechinf", "regulons_filtered.json")))
    regs = {k: [back.get(g, g) for g in v] for k, v in regs_miner.items()}
    rdf = pd.read_csv(os.path.join(mdir, "mechinf", "regulonDf.csv"))
    regulator = rdf.groupby(rdf["Regulon_ID"].astype(str))["Regulator"].first().map(lambda g: back.get(g, g))
    progs = {str(k): [str(r) for r in v] for k, v in json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_programs.json"))).items()}
    reg2prog = {r: k for k, v in progs.items() for r in v}

    T = pd.read_csv(p("data/reference/ot_drug_targets.tsv"), sep="\t")
    rows, R = [], {}
    for d, g in T.groupby("drug"):
        tg = set(g["target_ensembl"])
        rs = sorted({r for r in regs if regulator.get(r) in tg} | {r for r, v in regs.items() if tg & set(v)})
        act = g["action_type"].str.upper().iloc[0]
        R[d] = (rs, -1 if act in ("AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR") else 1)
        rows.append({"drug": d, "action": act, "targets": ",".join(sorted(g["target_symbol"].unique())), "n_regulons": len(rs),
                     "as_regulator": sum(regulator.get(r) in tg for r in rs), "programs": len({reg2prog.get(r) for r in rs} - {None})})
    DR = pd.DataFrame(rows).set_index("drug")
    DR.to_csv(os.path.join(outdir, "drug_regulons.tsv"), sep="\t")
    log.info("Drug -> regulons:\n%s", DR.to_string())

    def dcna(A):
        """A: regulons x samples trinary. Returns DCRA and DCPA (drugs x samples)."""
        A.index = A.index.astype(str)
        PA = pd.DataFrame({k: A.loc[[r for r in v if r in A.index]].mean() for k, v in progs.items() if any(r in A.index for r in v)}).T
        dr, dp = {}, {}
        for d, (rs, sign) in R.items():
            rr = [r for r in rs if r in A.index]
            if not rr:
                continue
            dr[d] = A.loc[rr].mean()
            pp = sorted({reg2prog[r] for r in rr if r in reg2prog and reg2prog[r] in PA.index})
            if pp:
                dp[d] = PA.loc[pp].mean()
        return pd.DataFrame(dr).T, pd.DataFrame(dp).T

    def trinary(z):
        rm = {k: [g for g in v if g in z.index] for k, v in regs.items()}
        rm = {k: v for k, v in rm.items() if len(v) > 1}
        return miner.generateRegulonActivity(rm, z, p=0.05)

    out, tests = {}, []
    # ---- discovery: MINER's own trinary membership
    over = pd.read_csv(os.path.join(mdir, "subtypes_filtered", "overExpressedMembers.csv"), index_col=0)
    under = pd.read_csv(os.path.join(mdir, "subtypes_filtered", "underExpressedMembers.csv"), index_col=0)
    out["discovery"] = dcna(over - under)
    # ---- treated cohorts
    cohorts = {"GSE109211": ("data/treated/GSE109211/matrix/GSE109211_series_matrix.txt.gz", "data/treated/GPL13938.txt", "Symbol", 0,
                             ("ch:outcome", "responder"), "ch:treatment", "sorafenib"),
               "GSE104580": ("data/treated/GSE104580/matrix/GSE104580_series_matrix.txt.gz", "data/treated/GPL570.annot.gz", "Gene symbol", 27,
                             ("ch:subject subgroup", "TACE responders"), None, "doxorubicin")}

    class _S:   # symbol -> Ensembl only
        sym2ens = genes.reset_index().drop_duplicates("symbol").set_index("symbol")["ensembl"]
        to_ens = s10.Scorer.to_ens
    for name, (path, annot, sym, skip, (lc, lv), arm, drug) in cohorts.items():
        x = s10.load_series(p(path), p(annot), sym, log, skip)
        z = s10.zrows(_S.to_ens(_S, x)).dropna(how="all")
        out[name] = dcna(trinary(z))
        smp = pd.read_csv(p(f"data/treated/{name}/samples_geo.tsv"), sep="\t", index_col=0)
        for level, M in (("regulon (DCRA)", out[name][0]), ("program (DCPA)", out[name][1])):
            if drug not in M.index:
                continue
            s = M.loc[drug]
            sign = R[drug][1]
            pred = ((sign * s) > 0).astype(int)
            y = (smp.loc[s.index, lc] == lv).astype(int)
            arms = [("all", s.index)] + ([(a, s.index[smp.loc[s.index, arm] == a]) for a in ("Sor", "Plac")] if arm else [])
            for an, idx in arms:
                ct = pd.crosstab(pred[idx], y[idx]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
                orr, pf = stats.fisher_exact(ct.values)
                au = stats.mannwhitneyu(s[idx][y[idx] == 1], s[idx][y[idx] == 0]).statistic / (y[idx].sum() * (1 - y[idx]).sum())
                tests.append({"cohort": name, "drug": drug, "level": level, "arm": an, "n": len(idx), "pred_responders": int(pred[idx].sum()),
                              "response_rate_pred_resp": ct.loc[1, 1] / max(ct.loc[1].sum(), 1), "response_rate_pred_nonresp": ct.loc[0, 1] / max(ct.loc[0].sum(), 1),
                              "fisher_or": orr, "fisher_p": pf, "auc_continuous": au})
            if arm:
                import statsmodels.formula.api as smf
                d = pd.DataFrame({"y": y, "s": s, "sor": (smp.loc[s.index, arm] == "Sor").astype(int)})
                m = smf.logit("y ~ s * sor", d).fit(disp=0)
                tests.append({"cohort": name, "drug": drug, "level": level, "arm": "interaction", "n": len(d), "interaction_p": m.pvalues["s:sor"],
                              "or_placebo": np.exp(m.params["s"]), "or_sorafenib": np.exp(m.params["s"] + m.params["s:sor"])})
    for name, (a, b) in out.items():
        a.to_csv(os.path.join(outdir, f"dcra_{name}.tsv"), sep="\t", float_format="%.4g")
        b.to_csv(os.path.join(outdir, f"dcpa_{name}.tsv"), sep="\t", float_format="%.4g")
    Tt = pd.DataFrame(tests)
    Tt.to_csv(os.path.join(outdir, "dcna_tests.tsv"), sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 220):
        log.info("Treated cohorts:\n%s", Tt.round(3).to_string(index=False))

    # ---- discovery: predicted-responder fraction vs trial ORR (monotherapy rows)
    soc = yaml.safe_load(open(p("config/soc_regimens.yaml")))["regimens"]
    mono = {r["regimen"].split(" (")[0]: r["orr"] for r in soc if " + " not in r["regimen"]}
    rows = []
    for level, M in (("regulon (DCRA)", out["discovery"][0]), ("program (DCPA)", out["discovery"][1])):
        for d, orr in mono.items():
            if d in M.index:
                rows.append({"level": level, "drug": d, "trial_orr": orr, "pred_responder_frac": float(((R[d][1] * M.loc[d]) > 0).mean()),
                             "mean_dcna": float(M.loc[d].mean()), "n_regulons": DR.loc[d, "n_regulons"]})
    O = pd.DataFrame(rows)
    O.to_csv(os.path.join(outdir, "dcna_vs_orr.tsv"), sep="\t", index=False, float_format="%.4g")
    for level, g in O.groupby("level"):
        rho, pv = stats.spearmanr(g["trial_orr"], g["pred_responder_frac"])
        log.info("%s: predicted-responder fraction vs trial ORR, %d drugs: Spearman %.2f (p %.2f)\n%s", level, len(g), rho, pv,
                 g.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
