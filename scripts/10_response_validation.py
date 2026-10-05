#!/usr/bin/env python
"""Step 10: do network activity patterns recapitulate response to standard-of-care drugs?
(docs/soc_response_validation_design.md; drug -> network map frozen in config/drug_network_map.yaml before scoring.)

Scoring (same portable method as step 08, no refitting): regulon = mean z of its genes within the cohort; program =
mean of its regulons; signatures = mean z (signed sets: up minus down); each component standardized within the cohort;
drug-class score = mean of components x sign. The step-06 risk score (TCGA RFS ridge) is a negative control.

L3, patient level (public treated cohorts):
  GSE109211  STORM adjuvant sorafenib vs placebo (FFPE, Illumina DASL); responder label in both arms -> AUC per arm and
             the treatment x score interaction (logistic). A predictive sorafenib score responds in the sorafenib arm
             only.
  GSE104580  TACE (Affymetrix U133 Plus 2, biopsies); responder vs non-responder (mRECIST).
  For each drug-class score: AUC (responder higher) with bootstrap 95% CI; null = the same class with its program
  components replaced by random programs of the same signs (n random draws), signatures kept.
L2, subgroups in the 929 discovery tumours (no treated data): ICI score by etiology (LICA-FR: viral vs non-viral),
  CTNNB1 and AXIN1 mutation; lenvatinib score by CCND1/FGF19 amplification. OLS on the within-cohort z with cohort
  terms; expected directions from the literature (design section 5); sign concordance.
L1 (population transport) needs an ICI-treated calibration cohort and is not run.

Outputs: results/10_response/{scores_<cohort>.tsv, l3_auc.tsv, l3_interaction.tsv, l2_subgroups.tsv}
"""

import json
import os

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from hcc_common import load_params, miner_id_backmap, p, setup_logging

OUT = "10_response"


def zrows(x):
    return x.sub(x.mean(1), axis=0).div(x.std(1).replace(0, np.nan), axis=0)


def zcol(v):
    return (v - v.mean()) / v.std()


def load_series(path, annot, sym_col, log, skip_annot=0, log2=None):
    x = pd.read_csv(path, sep="\t", comment="!", index_col=0)
    x = x.apply(pd.to_numeric, errors="coerce")
    if log2 is None:
        log2 = x.max().max() > 100
    if log2:
        x = np.log2(x.clip(lower=1))
    if annot.endswith(".txt"):
        lines = open(annot).read().split("\n")
        i = lines.index("!platform_table_begin")
        a = pd.read_csv(annot, sep="\t", skiprows=i + 1, usecols=["ID", sym_col], low_memory=False, comment=None)
        a = a[a["ID"] != "!platform_table_end"]
    else:
        a = pd.read_csv(annot, sep="\t", skiprows=skip_annot, usecols=["ID", sym_col], low_memory=False)
    a = a.dropna()
    a = a[~a[sym_col].astype(str).str.contains("///")]
    x = x.loc[x.index.intersection(a["ID"])]
    x["symbol"] = a.set_index("ID").loc[x.index, sym_col].values
    x = x.assign(m=x.drop(columns="symbol").mean(1)).sort_values("m", ascending=False).drop_duplicates("symbol")
    x = x.drop(columns="m").set_index("symbol")
    log.info("%s: %d genes x %d samples (log2 applied: %s)", os.path.basename(path), *x.shape, log2)
    return x


class Scorer:
    def __init__(self, P, res, log):
        mx = P["miner"]["matrix"]
        genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col=0)
        back = miner_id_backmap(p(P["miner"]["idmap"]), genes.index)
        mdir = os.path.join(res, "04_miner", mx)
        self.reg = {k: [back.get(g, g) for g in v] for k, v in json.load(open(os.path.join(mdir, "mechinf", "regulons_filtered.json"))).items()}
        self.progs = {str(k): [str(r) for r in v] for k, v in json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_programs.json"))).items()}
        self.sig = pd.read_csv(os.path.join(res, "07_post", "signatures", "signatures.tsv"), sep="\t").dropna(subset=["ensembl"])
        self.sym2ens = genes.reset_index().drop_duplicates("symbol").set_index("symbol")["ensembl"]
        self.genes = genes
        self.w = pd.read_csv(os.path.join(res, "06_risk", mx, f"predictor_ridge_programs_{'TCGA_RFS_h36m'}", "weights.tsv"),
                             sep="\t", index_col=0)["weight"].rename(index=str)
        self.log = log

    def to_ens(self, x):
        x = x.loc[x.index.intersection(self.sym2ens.index)]
        x.index = self.sym2ens[x.index].values
        return x.groupby(level=0).mean()

    def programs(self, z):
        reg = {k: z.loc[[g for g in v if g in z.index]].mean() for k, v in self.reg.items()
               if sum(g in z.index for g in v) >= max(3, 0.5 * len(v))}
        reg = pd.DataFrame(reg).T
        prog = pd.DataFrame({k: reg.loc[[r for r in v if r in reg.index]].mean() for k, v in self.progs.items()
                             if any(r in reg.index for r in v)}).T
        cov = len(set(z.index) & set().union(*map(set, self.reg.values()))) / len(set().union(*map(set, self.reg.values())))
        return zrows(prog), cov, len(reg)

    def signature(self, z, entry):
        def m(name):
            g = self.sig.loc[self.sig["set"] == name, "ensembl"]
            g = [x for x in g.unique() if x in z.index]
            return z.loc[g].mean() if len(g) >= 3 else None
        if " - " in entry:
            a, b = [e.strip() for e in entry.split(" - ")]
            va, vb = m(a), m(b)
            return None if va is None else (va - vb if vb is not None else va)
        return m(entry)

    def drug_scores(self, z, prog, dmap):
        out, comps = {}, {}
        for d, spec in dmap.items():
            parts = []
            for k, sgn in (spec.get("programs") or {}).items():
                if str(k) in prog.index:
                    parts.append(sgn * prog.loc[str(k)])
            for e, sgn in (spec.get("signatures") or {}).items():
                v = self.signature(z, e)
                if v is not None:
                    parts.append(sgn * zcol(v))
                    comps[(d, e)] = zcol(v)
            out[d] = zcol(pd.concat(parts, axis=1).mean(1)) if parts else np.nan
        risk = prog.fillna(0).T @ self.w.reindex(prog.index).fillna(0)
        out["risk_score (negative control)"] = zcol(risk)
        return pd.DataFrame(out), comps


def auc(score, y):
    a, b = score[y == 1], score[y == 0]
    return stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))


def boot_ci(score, y, rng, n=1000):
    idx = np.arange(len(y))
    v = []
    for _ in range(n):
        i = rng.choice(idx, len(idx))
        if y.iloc[i].nunique() == 2:
            v.append(auc(score.iloc[i], y.iloc[i]))
    return np.percentile(v, [2.5, 97.5])


def main():
    P = load_params(None)
    res = p(P["paths"]["results"])
    outdir = os.path.join(res, OUT)
    os.makedirs(outdir, exist_ok=True)
    log = setup_logging(outdir, "10_response_validation")
    dm = yaml.safe_load(open(p("config/drug_network_map.yaml")))
    dmap = dm["drug_classes"]
    log.info("Drug map version %s; classes %s", dm["version"], list(dmap))
    S = Scorer(P, res, log)
    rng = np.random.default_rng(10)
    nrand = int(dm["negative_controls"]["random_program_sets"])

    cohorts = {
        "GSE109211": dict(path="data/treated/GSE109211/matrix/GSE109211_series_matrix.txt.gz", annot="data/treated/GPL13938.txt",
                          sym="Symbol", label=("ch:outcome", "responder"), arm="ch:treatment"),
        "GSE104580": dict(path="data/treated/GSE104580/matrix/GSE104580_series_matrix.txt.gz", annot="data/treated/GPL570.annot.gz",
                          sym="Gene symbol", skip=27, label=("ch:subject subgroup", "TACE responders"), arm=None),
    }
    l3, inter = [], []
    for name, c in cohorts.items():
        x = load_series(p(c["path"]), p(c["annot"]), c["sym"], log, c.get("skip", 0))
        z = zrows(S.to_ens(x)).dropna(how="all")
        prog, cov, nreg = S.programs(z)
        D, comps = S.drug_scores(z, prog, dmap)
        smp = pd.read_csv(p(f"data/treated/{name}/samples_geo.tsv"), sep="\t", index_col=0).reindex(D.index)
        y = (smp[c["label"][0]] == c["label"][1]).astype(int)
        D["responder"] = y
        if c["arm"]:
            D["arm"] = smp[c["arm"]]
        # published biomarkers as references
        for e in ("custom:ZHU2022_ABRS", "custom:HABER2023_IFNAP", "custom:SIA2017_IMMUNE_CLASS_UP"):
            v = S.signature(z, e)
            if v is not None:
                D[e.replace("custom:", "ref:")] = zcol(v)
        D.to_csv(os.path.join(outdir, f"scores_{name}.tsv"), sep="\t", float_format="%.4g")
        log.info("%s: %d samples (%d responders); network genes present %.0f%%; regulons scored %d; programs %d", name, len(D), y.sum(),
                 100 * cov, nreg, len(prog))
        arms = [("all", D.index)] + ([(a, D.index[D["arm"] == a]) for a in sorted(D["arm"].dropna().unique())] if c["arm"] else [])
        for arm, idx in arms:
            yy = y[idx]
            for col in [k for k in D.columns if k not in ("responder", "arm")]:
                s = D.loc[idx, col]
                if s.isna().all():
                    continue
                a = auc(s, yy)
                lo, hi = boot_ci(s, yy, rng)
                r = {"cohort": name, "arm": arm, "score": col, "n": len(idx), "responders": int(yy.sum()), "auc": a, "ci_lo": lo, "ci_hi": hi,
                     "p_mwu": stats.mannwhitneyu(s[yy == 1], s[yy == 0]).pvalue}
                if col in dmap and dmap[col].get("programs"):
                    sp = dmap[col]["programs"]
                    null = []
                    pool = [k for k in prog.index if k not in {str(k) for k in sp}]
                    for _ in range(nrand):
                        ks = rng.choice(pool, len(sp), replace=False)
                        parts = [sg * prog.loc[k] for k, sg in zip(ks, sp.values())] + [sg * comps[(col, e)] for e, sg in
                                                                                         (dmap[col].get("signatures") or {}).items() if (col, e) in comps]
                        null.append(auc(zcol(pd.concat(parts, axis=1).mean(1))[idx], yy))
                    r["null_auc_mean"] = float(np.mean(null))
                    r["p_vs_random_programs"] = (1 + np.sum(np.abs(np.array(null) - 0.5) >= abs(a - 0.5))) / (1 + nrand)
                l3.append(r)
        if c["arm"]:
            import statsmodels.formula.api as smf
            for col in [k for k in D.columns if k not in ("responder", "arm")]:
                d = D[["responder", "arm", col]].dropna().rename(columns={col: "s"})
                d["sor"] = (d["arm"] == "Sor").astype(int)
                try:
                    m = smf.logit("responder ~ s * sor", d).fit(disp=0)
                    inter.append({"cohort": name, "score": col, "or_per_sd_placebo": np.exp(m.params["s"]),
                                  "or_per_sd_sorafenib": np.exp(m.params["s"] + m.params["s:sor"]), "interaction_p": m.pvalues["s:sor"]})
                except Exception as e:  # noqa: BLE001
                    log.warning("%s %s: interaction model failed (%s)", name, col, e)
    L3 = pd.DataFrame(l3)
    L3.to_csv(os.path.join(outdir, "l3_auc.tsv"), sep="\t", index=False, float_format="%.4g")
    I = pd.DataFrame(inter)
    I.to_csv(os.path.join(outdir, "l3_interaction.tsv"), sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 220, "display.max_rows", 200):
        log.info("L3 AUC (responder higher):\n%s", L3.round(3).to_string(index=False))
        log.info("GSE109211 treatment x score interaction (logistic):\n%s", I.round(3).to_string(index=False))

    # ---------------- L2: subgroup contrasts in the discovery tumours
    z = pd.read_csv(os.path.join(res, "02_batch_corrected", f"expression_{P['miner']['matrix']}_z.csv"), index_col=0)
    sm = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    G = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
    Dd = []
    for c in sm["cohort"].unique():
        zc = zrows(z[sm.index[sm["cohort"] == c]])
        prog, _, _ = S.programs(zc)
        Dd.append(S.drug_scores(zc, prog, dmap)[0])
    Dd = pd.concat(Dd)
    Dd.to_csv(os.path.join(outdir, "scores_discovery.tsv"), sep="\t", float_format="%.4g")
    et = sm["etiology"].astype(str)
    viral = pd.Series(np.where(et.str.contains("HBV|HCV"), 1, np.where(et.isin(["nan", "", "None"]), np.nan, 0)), index=sm.index)
    contrasts = [("ici", "viral vs non-viral (LICA-FR)", viral, +1),
                 ("ici", "CTNNB1-mutant vs wild type", G.loc["MUT_CTNNB1"].reindex(sm.index), -1),
                 ("ici", "AXIN1-mutant vs wild type (new prediction)", G.loc["MUT_AXIN1"].reindex(sm.index), -1),
                 ("lenvatinib", "CCND1/FGF19 amplified vs not (hypothesis)",
                  G.loc["AMP_CCND1_FGF19"].reindex(sm.index) if "AMP_CCND1_FGF19" in G.index else None, +1),
                 ("risk_score (negative control)", "CTNNB1-mutant vs wild type", G.loc["MUT_CTNNB1"].reindex(sm.index), 0)]
    import statsmodels.formula.api as smf
    l2 = []
    for d, lab, g, exp in contrasts:
        if g is None:
            continue
        dd = pd.DataFrame({"s": Dd[d], "g": g, "cohort": sm["cohort"]}).dropna()
        m = smf.ols("s ~ g + C(cohort)" if dd["cohort"].nunique() > 1 else "s ~ g", dd).fit()
        l2.append({"score": d, "contrast": lab, "n": len(dd), "n_group": int(dd["g"].sum()), "effect_sd": m.params["g"],
                   "p": m.pvalues["g"], "expected_sign": exp,
                   "concordant": (np.sign(m.params["g"]) == exp) if exp else np.nan})
    L2 = pd.DataFrame(l2)
    L2.to_csv(os.path.join(outdir, "l2_subgroups.tsv"), sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 220):
        log.info("L2 subgroup contrasts (discovery; effect in within-cohort SD):\n%s", L2.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
