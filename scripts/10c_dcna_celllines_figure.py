#!/usr/bin/env python
"""Step 10c: DCNA in HCC cell lines against measured drug sensitivity (GDSC), and the step-10 drug-response figure.

Cell lines: HCC models of the Sanger Cell Model Passports (model_list) with RNA-seq (rnaseq_all, TPM) and GDSC IC50
(GDSC2 preferred; GDSC1 for drugs not in GDSC2). Expression -> log2(TPM + 1) -> Ensembl -> gene z across the HCC lines ->
MINER trinary regulon activity -> DCRA / DCPA per drug (Open Targets targets, as step 10b). Inhibitors / antagonists are
predicted sensitive when DCNA > 0.
Tests (all drugs with mapped regulons and >= 8 lines):
  - pooled: drug-standardized ln IC50 (z across lines within drug) of predicted responders vs non-responders (lower =
    more sensitive), Mann-Whitney; per drug, the difference in mean z and Spearman(DCNA, ln IC50);
  - null: line labels permuted within drug (n_perm) for the pooled difference.
Figure (results/10_response/figures/drug_response.{png,pdf}):
  a GSE109211 response rate by DCNA-predicted class, sorafenib and placebo arms   b GSE104580 TACE (doxorubicin proxy)
  c cell lines, all drugs pooled: z ln IC50 by predicted class (regulon / program level)
  d cell lines, standard-of-care drugs   e DCRA of standard-of-care drugs across the 929 discovery tumours
Outputs: results/10_response/dcna/{celllines_dcra.tsv, celllines_dcpa.tsv, celllines_tests.tsv, celllines_by_drug.tsv}
"""

import json
import os
import zipfile

import numpy as np
import pandas as pd
from scipy import stats

from hcc_common import load_params, miner_id_backmap, p, setup_logging

SOC = ["sorafenib", "regorafenib", "cabozantinib", "doxorubicin", "epirubicin", "5-fluorouracil", "fluorouracil", "cisplatin",
       "oxaliplatin", "gemcitabine", "axitinib", "pazopanib", "mitomycin-c"]


def load_cmp_expression(zpath, models, log):
    z = zipfile.ZipFile(zpath)
    name = [n for n in z.namelist() if n.endswith(".csv")][0]
    parts = []
    with z.open(name) as fh:
        for ch in pd.read_csv(fh, chunksize=2_000_000, low_memory=False):
            ch = ch[ch["model_id"].isin(models)]
            if len(ch):
                parts.append(ch)
    d = pd.concat(parts)
    log.info("CMP RNA-seq: %d rows for %d models; sources %s; columns %s", len(d), d["model_id"].nunique(),
             d["data_source"].value_counts().to_dict() if "data_source" in d else "-", list(d.columns))
    if "data_source" in d:   # one source per model: Sanger first
        d["rank"] = d["data_source"].map({"Sanger": 0}).fillna(1)
        best = d.groupby("model_id")["rank"].min()
        d = d[d["rank"] == d["model_id"].map(best)]
    x = d.pivot_table(index="gene_symbol", columns="model_id", values="tpm", aggfunc="mean")
    return np.log2(x.fillna(0) + 1)


def main():
    from miner import miner
    P = load_params(None)
    res = p(P["paths"]["results"])
    mx = P["miner"]["matrix"]
    outdir = os.path.join(res, "10_response", "dcna")
    figdir = os.path.join(res, "10_response", "figures")
    os.makedirs(figdir, exist_ok=True)
    log = setup_logging(outdir, "10c_dcna_celllines")
    rng = np.random.default_rng(3)

    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col=0)
    back = miner_id_backmap(p(P["miner"]["idmap"]), genes.index)
    mdir = os.path.join(res, "04_miner", mx)
    regs = {k: [back.get(g, g) for g in v] for k, v in json.load(open(os.path.join(mdir, "mechinf", "regulons_filtered.json"))).items()}
    rdf = pd.read_csv(os.path.join(mdir, "mechinf", "regulonDf.csv"))
    regulator = rdf.groupby(rdf["Regulon_ID"].astype(str))["Regulator"].first().map(lambda g: back.get(g, g))
    progs = {str(k): [str(r) for r in v] for k, v in json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_programs.json"))).items()}
    reg2prog = {r: k for k, v in progs.items() for r in v}
    T = pd.read_csv(p("data/reference/ot_drug_targets.tsv"), sep="\t")
    R = {}
    for d, g in T.groupby("drug"):
        tg = set(g["target_ensembl"])
        rs = sorted({r for r in regs if regulator.get(r) in tg} | {r for r, v in regs.items() if tg & set(v)})
        if rs:
            act = str(g["action_type"].iloc[0]).upper()
            R[d] = (rs, -1 if act in ("AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR") else 1)
    log.info("Drugs with Open Targets targets: %d; with >= 1 mapped regulon: %d", T["drug"].nunique(), len(R))

    # ---- cell lines
    ml = pd.read_csv(p("data/celllines/model_list_20240110.csv"), low_memory=False)
    hcc = ml[ml["cancer_type"].astype(str).str.contains("Hepatocellular", case=False)]
    G = pd.read_csv(p("data/celllines/gdsc_hcc_ic50.tsv"), sep="\t")
    G["drug"] = G["DRUG_NAME"].astype(str).str.lower()
    G = G.sort_values("ds", ascending=False).drop_duplicates(["drug", "SANGER_MODEL_ID"])     # GDSC2 first
    x = load_cmp_expression(p("data/celllines/rnaseq_all_20220624.zip"), set(hcc["model_id"]), log)
    sym2ens = genes.reset_index().drop_duplicates("symbol").set_index("symbol")["ensembl"]
    x = x.loc[x.index.intersection(sym2ens.index)]
    x.index = sym2ens[x.index].values
    x = x.groupby(level=0).mean()
    x = x.loc[x.index.intersection(genes.index[genes["kept"]])]
    x = x.loc[x.var(1) > 0]
    zc = x.sub(x.mean(1), axis=0).div(x.std(1), axis=0)
    log.info("HCC lines with RNA-seq: %d (with GDSC data: %d); genes %d", zc.shape[1], len(set(zc.columns) & set(G["SANGER_MODEL_ID"])), zc.shape[0])
    rm = {k: [g for g in v if g in zc.index] for k, v in regs.items()}
    A = miner.generateRegulonActivity({k: v for k, v in rm.items() if len(v) > 1}, zc, p=0.05)
    A.index = A.index.astype(str)
    PA = pd.DataFrame({k: A.loc[[r for r in v if r in A.index]].mean() for k, v in progs.items() if any(r in A.index for r in v)}).T
    dra, dpa = {}, {}
    for d, (rs, sign) in R.items():
        rr = [r for r in rs if r in A.index]
        if rr:
            dra[d] = A.loc[rr].mean()
            pp = sorted({reg2prog[r] for r in rr if r in reg2prog and reg2prog[r] in PA.index})
            if pp:
                dpa[d] = PA.loc[pp].mean()
    DRA, DPA = pd.DataFrame(dra).T, pd.DataFrame(dpa).T
    DRA.to_csv(os.path.join(outdir, "celllines_dcra.tsv"), sep="\t", float_format="%.4g")
    DPA.to_csv(os.path.join(outdir, "celllines_dcpa.tsv"), sep="\t", float_format="%.4g")

    tests, bydrug, pooled = [], [], {}
    for level, M in (("regulon (DCRA)", DRA), ("program (DCPA)", DPA)):
        rows = []
        for d in M.index:
            g = G[G["drug"] == d].set_index("SANGER_MODEL_ID")["LN_IC50"]
            lines = g.index.intersection(M.columns)
            if len(lines) < 8:
                continue
            s, ic = M.loc[d, lines], g[lines]
            zi = (ic - ic.mean()) / ic.std()
            pred = (R[d][1] * s > 0).astype(int)
            for l in lines:
                rows.append({"drug": d, "line": l, "dcna": s[l], "pred": pred[l], "z_ln_ic50": zi[l], "ln_ic50": ic[l]})
            rho = stats.spearmanr(s, ic).statistic if s.nunique() > 1 else np.nan
            bydrug.append({"level": level, "drug": d, "n_lines": len(lines), "n_pred_resp": int(pred.sum()), "spearman_dcna_lnic50": rho,
                           "delta_z": zi[pred == 1].mean() - zi[pred == 0].mean() if 0 < pred.sum() < len(pred) else np.nan})
        X = pd.DataFrame(rows)
        pooled[level] = X
        a, b = X.loc[X["pred"] == 1, "z_ln_ic50"], X.loc[X["pred"] == 0, "z_ln_ic50"]
        obs = a.mean() - b.mean()
        null = []
        for _ in range(1000):
            Xp = X.copy()
            Xp["pred"] = Xp.groupby("drug")["pred"].transform(lambda v: rng.permutation(v.values))
            null.append(Xp.loc[Xp["pred"] == 1, "z_ln_ic50"].mean() - Xp.loc[Xp["pred"] == 0, "z_ln_ic50"].mean())
        bd = pd.DataFrame([r for r in bydrug if r["level"] == level])
        tests.append({"level": level, "drugs": X["drug"].nunique(), "pairs": len(X), "pred_resp_pairs": int(X["pred"].sum()),
                      "mean_z_pred_resp": a.mean(), "mean_z_pred_nonresp": b.mean(), "delta": obs,
                      "mwu_p": stats.mannwhitneyu(a, b).pvalue, "perm_p_one_sided": (1 + np.sum(np.array(null) <= obs)) / 1001,
                      "drugs_with_negative_rho": int((bd["spearman_dcna_lnic50"] < 0).sum()), "drugs_with_rho": int(bd["spearman_dcna_lnic50"].notna().sum()),
                      "sign_test_p": stats.binomtest(int((bd["spearman_dcna_lnic50"] < 0).sum()), int(bd["spearman_dcna_lnic50"].notna().sum())).pvalue})
    Tt, Bd = pd.DataFrame(tests), pd.DataFrame(bydrug)
    Tt.to_csv(os.path.join(outdir, "celllines_tests.tsv"), sep="\t", index=False, float_format="%.4g")
    Bd.to_csv(os.path.join(outdir, "celllines_by_drug.tsv"), sep="\t", index=False, float_format="%.4g")
    with pd.option_context("display.width", 220):
        log.info("Cell lines, pooled (lower z ln IC50 = more sensitive):\n%s", Tt.round(3).to_string(index=False))
        log.info("Standard-of-care drugs:\n%s", Bd[Bd["drug"].isin(SOC)].round(3).to_string(index=False))

    # ---------------- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 6, "axes.titlesize": 6.5, "axes.spines.top": False, "axes.spines.right": False,
                         "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"]})
    RED, BLUE = "#e34948", "#2a78d6"
    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 1.1], hspace=0.65, wspace=1.1)
    tt = pd.read_csv(os.path.join(outdir, "dcna_tests.tsv"), sep="\t")
    ax = fig.add_subplot(gs[0, 0:2])
    t = tt[(tt["cohort"] == "GSE109211") & (tt["level"] == "regulon (DCRA)") & tt["arm"].isin(["Sor", "Plac"])].set_index("arm")
    xs = np.arange(2)
    for j, (k, col, lab) in enumerate((("response_rate_pred_nonresp", BLUE, "predicted non-responder"), ("response_rate_pred_resp", RED, "predicted responder"))):
        ax.bar(xs + (j - 0.5) * 0.36, 100 * t.loc[["Sor", "Plac"], k], 0.34, color=col, label=lab)
    for i, a_ in enumerate(["Sor", "Plac"]):
        ax.text(i, 100 * t.loc[a_, ["response_rate_pred_resp", "response_rate_pred_nonresp"]].max() + 3, f"OR {t.loc[a_, 'fisher_or']:.1f}\nP {t.loc[a_, 'fisher_p']:.2f}",
                ha="center", fontsize=5)
    ax.set_xticks(xs)
    ax.set_xticklabels(["sorafenib arm", "placebo arm"])
    ax.set_ylabel("observed responders (%)")
    ax.set_ylim(0, 60)
    ax.legend(fontsize=4.8, frameon=False, loc="upper left")
    ax.set_title("a  STORM (GSE109211): sorafenib DCRA", loc="left")
    ax = fig.add_subplot(gs[0, 2:4])
    t = tt[(tt["cohort"] == "GSE104580")].set_index("level")
    for j, (k, col) in enumerate((("response_rate_pred_nonresp", BLUE), ("response_rate_pred_resp", RED))):
        ax.bar(xs + (j - 0.5) * 0.36, 100 * t.loc[["regulon (DCRA)", "program (DCPA)"], k], 0.34, color=col)
    for i, l_ in enumerate(["regulon (DCRA)", "program (DCPA)"]):
        ax.text(i, 100 * t.loc[l_, ["response_rate_pred_resp", "response_rate_pred_nonresp"]].max() + 3, f"OR {t.loc[l_, 'fisher_or']:.2f}\nP {t.loc[l_, 'fisher_p']:.0e}",
                ha="center", fontsize=5)
    ax.set_xticks(xs)
    ax.set_xticklabels(["regulon level", "program level"])
    ax.set_ylim(0, 95)
    ax.set_title("b  TACE (GSE104580): doxorubicin DCNA", loc="left")
    for j, level in enumerate(("regulon (DCRA)", "program (DCPA)")):
        ax = fig.add_subplot(gs[0, 4 + j])
        X = pooled[level]
        dd = [X.loc[X["pred"] == 0, "z_ln_ic50"], X.loc[X["pred"] == 1, "z_ln_ic50"]]
        bp = ax.boxplot(dd, widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color="k"))
        for b_, c_ in zip(bp["boxes"], (BLUE, RED)):
            b_.set_facecolor(c_)
            b_.set_alpha(0.6)
        for i, v in enumerate(dd):
            ax.scatter(i + 1 + rng.uniform(-0.18, 0.18, len(v)), v, s=1, color="k", alpha=0.25, lw=0)
        r = Tt.set_index("level").loc[level]
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["non-resp.", "resp."], fontsize=5)
        ax.set_title(f"{'c' if j == 0 else ''}  {level.split(' ')[0]}\n{int(r['drugs'])} drugs, {len(set(X['line']))} lines\nperm P {r['perm_p_one_sided']:.2f}",
                     loc="left", fontsize=5.5)
        if j == 0:
            ax.set_ylabel("ln IC50 (z within drug)")
    ax = fig.add_subplot(gs[1, 0:3])
    X = pooled["regulon (DCRA)"]
    soc = [d for d in SOC if d in set(X["drug"])]
    for i, d in enumerate(soc):
        v = X[X["drug"] == d]
        for k, col in ((0, BLUE), (1, RED)):
            w = v[v["pred"] == k]["ln_ic50"]
            ax.scatter(np.full(len(w), i + (k - 0.5) * 0.3) + rng.uniform(-0.06, 0.06, len(w)), w, s=5, color=col, lw=0)
            if len(w):
                ax.plot([i + (k - 0.5) * 0.3 - 0.12, i + (k - 0.5) * 0.3 + 0.12], [w.median()] * 2, color="k", lw=0.8)
    ax.set_xticks(range(len(soc)))
    ax.set_xticklabels(soc, rotation=40, ha="right")
    ax.set_ylabel("ln IC50 (µM)")
    ax.set_title("d  HCC cell lines, standard-of-care drugs (regulon DCNA; red = predicted responder)", loc="left")
    ax = fig.add_subplot(gs[1, 3:6])
    D = pd.read_csv(os.path.join(outdir, "dcra_discovery.tsv"), sep="\t", index_col=0)
    show = [d for d in ["sorafenib", "lenvatinib", "regorafenib", "cabozantinib", "atezolizumab", "durvalumab", "pembrolizumab",
                        "ipilimumab", "bevacizumab", "doxorubicin", "fluorouracil"] if d in D.index]
    M = D.loc[show]
    o = np.argsort(M.loc["sorafenib"].values) if "sorafenib" in M.index else np.arange(M.shape[1])
    im = ax.imshow(M.values[:, o], aspect="auto", cmap="RdBu_r", vmin=-0.6, vmax=0.6, interpolation="none")
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(show)
    ax.set_xticks([])
    ax.set_xlabel(f"{M.shape[1]} discovery tumours (ordered by sorafenib DCRA)")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("DCRA", fontsize=5)
    ax.set_title("e  DCRA of standard-of-care drugs", loc="left")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(figdir, f"drug_response.{ext}"), dpi=300, bbox_inches="tight")
    log.info("Figure written")


if __name__ == "__main__":
    main()
