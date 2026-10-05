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
    """Cell Model Passports rnaseq_tpm (wide: 4 header rows model_id / model_name / dataset_name / data_source, then
    gene_id, symbol, one column per model; one source per model already chosen by the Passports)."""
    z = zipfile.ZipFile(zpath)
    with z.open("rnaseq_tpm_20220624.csv") as fh:
        hdr = pd.read_csv(fh, nrows=4, header=None, low_memory=False)
    ids = hdr.iloc[0, 2:].tolist()
    keep = [i for i, m in enumerate(ids) if m in models]
    with z.open("rnaseq_tpm_20220624.csv") as fh:
        x = pd.read_csv(fh, skiprows=5, header=None, usecols=[1] + [k + 2 for k in keep], low_memory=False)
    x.columns = ["symbol"] + [ids[k] for k in keep]
    x = x.dropna(subset=["symbol"]).set_index("symbol").apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    log.info("CMP RNA-seq TPM: %d genes x %d HCC models (sources %s)", *x.shape,
             pd.Series(hdr.iloc[3, [k + 2 for k in keep]].values).value_counts().to_dict())
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
            rho = stats.spearmanr(s, ic)[0] if s.nunique() > 1 else np.nan
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
    # null 2: each drug's regulon set replaced by random regulons of the same size (keeps line-level sensitivity axes)
    allr = list(A.index)
    X0 = pooled["regulon (DCRA)"]
    ic_z = X0.set_index(["drug", "line"])["z_ln_ic50"]
    obs = X0.loc[X0["pred"] == 1, "z_ln_ic50"].mean() - X0.loc[X0["pred"] == 0, "z_ln_ic50"].mean()
    nullr = []
    for _ in range(200):
        dl = []
        for d in X0["drug"].unique():
            k = len([r for r in R[d][0] if r in A.index])
            rr = rng.choice(allr, k, replace=False)
            sc = A.loc[rr].mean()
            lines = X0.loc[X0["drug"] == d, "line"]
            pr = (R[d][1] * sc[lines] > 0).astype(int).values
            dl.append(pd.DataFrame({"pred": pr, "z": ic_z.loc[d].loc[lines].values}))
        dl = pd.concat(dl)
        nullr.append(dl.loc[dl["pred"] == 1, "z"].mean() - dl.loc[dl["pred"] == 0, "z"].mean())
    tests[0]["random_regulon_null_mean"] = float(np.nanmean(nullr))
    tests[0]["random_regulon_p_one_sided"] = (1 + np.sum(np.array(nullr) <= obs)) / (1 + len(nullr))
    log.info("Regulon-level pooled delta %.3f; random same-size regulon sets: mean %.3f, P %.3f", obs, np.nanmean(nullr),
             tests[0]["random_regulon_p_one_sided"])
    for level, X in pooled.items():
        X.to_csv(os.path.join(outdir, f"celllines_pairs_{level.split(' ')[0]}.tsv"), sep="\t", index=False, float_format="%.4g")
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
    from statsmodels.stats.proportion import proportion_confint

    def rate_line(ax, score, y, groups, col, lab, dx=0.0):
        """Observed response rate (95% Wilson CI) per DCNA group; groups = Series of ordered group labels."""
        xs, rs, lo, hi, ns = [], [], [], [], []
        for k, g in enumerate(groups.cat.categories):
            m = groups == g
            n, e = int(m.sum()), int(y[m].sum())
            if n == 0:
                continue
            l_, h_ = proportion_confint(e, n, method="wilson")
            xs.append(k + dx); rs.append(100 * e / n); lo.append(100 * (e / n - l_)); hi.append(100 * (h_ - e / n)); ns.append(n)
        ax.errorbar(xs, rs, yerr=[lo, hi], color=col, marker="o", ms=3, lw=1, capsize=1.5, label=lab)
        for x_, r_, n_ in zip(xs, rs, ns):
            ax.text(x_, -6, f"{n_}", ha="center", va="top", fontsize=4.2, color=col)
        ax.text(-0.6, -6, "n", ha="right", va="top", fontsize=4.2, color="k")

    # a: STORM, sorafenib DCRA quartiles (within arm), observed response rate per quartile, sorafenib vs placebo arm
    ax = fig.add_subplot(gs[0, 0:2])
    smp = pd.read_csv(p("data/treated/GSE109211/samples_geo.tsv"), sep="\t", index_col=0)
    sc = pd.read_csv(os.path.join(outdir, "dcra_GSE109211.tsv"), sep="\t", index_col=0).loc["sorafenib"]
    y = (smp.loc[sc.index, "ch:outcome"] == "responder").astype(int)
    labels = ["Q1\n(lowest)", "Q2", "Q3", "Q4\n(highest)"]
    for arm, col, lab, dx in (("Sor", RED, "sorafenib arm", -0.08), ("Plac", "#6b6b6b", "placebo arm", 0.08)):
        idx = smp.index[(smp["ch:treatment"] == arm)].intersection(sc.index)
        q = pd.qcut(sc[idx].rank(method="first"), 4, labels=labels)
        rate_line(ax, sc[idx], y[idx], q, col, lab, dx)
    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=5)
    ax.set_xlabel("sorafenib DCRA quartile (within arm)", labelpad=8)
    ax.set_ylabel("observed responders (%)")
    ax.set_ylim(-12, 80)
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.legend(fontsize=4.8, frameon=False, loc="upper left")
    ax.set_title("a  STORM (GSE109211)\n    sorafenib vs placebo", loc="left")
    # b: TACE, doxorubicin DCRA (3 TOP2A regulons: discrete values), observed response rate per level
    ax = fig.add_subplot(gs[0, 2:4])
    smp = pd.read_csv(p("data/treated/GSE104580/samples_geo.tsv"), sep="\t", index_col=0)
    sc = pd.read_csv(os.path.join(outdir, "dcra_GSE104580.tsv"), sep="\t", index_col=0).loc["doxorubicin"]
    y = (smp.loc[sc.index, "ch:subject subgroup"] == "TACE responders").astype(int)
    cats = ["−1", "−0.67 to −0.33", "0", "0.33 to 0.67", "+1"]
    g = pd.cut(sc, [-1.01, -0.99, -0.01, 0.01, 0.99, 1.01], labels=cats)
    rate_line(ax, sc, y, g, "#2a78d6", "TACE")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(cats, fontsize=5)
    ax.set_xlabel("doxorubicin (TOP2A) DCRA", labelpad=8)
    ax.set_ylim(-12, 100)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_title("b  TACE (GSE104580)", loc="left")
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
        pp = r["perm_p_one_sided"]
        extra = f"\nrandom regulons P {r['random_regulon_p_one_sided']:.3f}" if pd.notna(r.get("random_regulon_p_one_sided", np.nan)) else ""
        ax.set_title(f"{'c' if j == 0 else ''}  {level.split(' ')[0]} level\n{int(r['drugs'])} drugs, {len(set(X['line']))} lines\n"
                     f"label perm P {'< 0.001' if pp <= 0.001 else f'{pp:.3f}'}{extra}", loc="left", fontsize=5.2)
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
    ax.set_title("d  HCC lines, SOC drugs (red = predicted responder)", loc="left")
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
