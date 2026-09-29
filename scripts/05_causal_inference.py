#!/usr/bin/env python
"""Step 05: MINER3 causal inference (genomic feature -> regulator -> regulon), run in the miner3 env.

Calls the MINER library directly instead of miner3-causalinference, because the CLI re-derives
regulons from coregulationModules.json and would ignore the step 04c filter. Additions to
MINER's causalNetworkAnalysis():

  1. Profiled samples only. MINER treats every non-1 column as wild type, so each feature is run
     on the samples profiled for it (genomic_features.csv: empty = not profiled).
  2. Multiple testing. MINER filters on nominal p < 0.05 and writes only hits. Here the same two
     tests are recomputed for ALL regulons and regulators of each feature (regulon eigengene Welch
     t-test; mutation-regulator edge = mean t over the regulons that contain the regulator gene, else
     the regulator's own expression t-test, as in MINER), and Benjamini-Hochberg q-values are
     attached per feature.
  3. Cohort consistency. MINER's test pools cohorts. For every flow the same statistics are computed
     within each cohort with >= causal.min_cohort_altered mutants and wild types, and the flow is
     marked consistent when the regulon and edge directions agree with the pooled ones in every
     tested cohort.
  4. Redundancy. Regulons are grouped into families (average-linkage on eigengene correlation,
     r >= causal.family_r) and mapped to transcriptional programs; results are summarized per
     feature x family, so ~100 near-identical regulons count as one signal.

Inputs:  results/02_batch_corrected/expression_<matrix>_z.csv, results/04_miner/<matrix>/mechinf/
         {regulons_filtered.json, regulonDf.csv}, subtypes_filtered/{transcriptional_programs.json,
         coherentMembers.csv}, results/03_genomics_clinical/genomic_features.csv
Outputs (results/05_causal/<matrix>/):
  causal_results/pooled/*.csv   MINER per-feature output
  completeCausalResults.csv     all MINER flows, annotated (symbols, q-values, cohort stats, family, program)
  filteredCausalResults.csv     MINER's CLI filter (same four criteria as miner3-causalinference)
  highConfidenceCausalResults.csv  filtered + q <= causal.q_max for both tests + cohort-consistent + |d| >= min_abs_d
                                   + feature altered in >= min_altered_hc samples; sorted by |d|
  causal_by_family.tsv          one row per feature x regulon family (best flow = largest |d|, ranked by |d|)
  causal_by_feature.tsv         per feature: flows, regulators, families, programs at each level
  regulon_families.tsv          regulon -> family, program, regulator symbol
  wiring_diagram.csv            MINER wiringDiagram on the filtered results (if it runs)
  qc/c1_causal_by_feature.png, qc/c2_feature_program_heatmap.png
"""

import argparse
import json
import os
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from hcc_common import load_params, miner_id_backmap, p, setup_logging

OUT = "05_causal"
_G = {}  # shared state for worker processes (set before the pool forks)


def bh(pv):
    pv = np.asarray(pv, float)
    n = np.sum(~np.isnan(pv))
    q = np.full_like(pv, np.nan)
    ok = ~np.isnan(pv)
    order = np.argsort(pv[ok])
    ranked = pv[ok][order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q_ok = np.empty(n)
    q_ok[order] = np.minimum(ranked, 1)
    q[ok] = q_ok
    return q


def welch(E, mut, wt):
    t, pv = stats.ttest_ind(E[:, mut], E[:, wt], equal_var=False, axis=1)
    return t, pv


def feature_stats(feature, cols_mut, cols_wt):
    """Regulon t/p for all regulons and edge t/p for all regulators (MINER's definitions).
    For a regulator that is itself a regulon gene, MINER uses the mean t and mean -log10 p over those
    regulons; the edge p here is the matching geometric-mean p (10^-mean(-log10 p)), used for BH."""
    E, X, reg_ids, reg_members, regulators = _G["E"], _G["X"], _G["reg_ids"], _G["reg_members"], _G["regulators"]
    si = _G["sample_index"]
    mut = [si[s] for s in cols_mut]
    wt = [si[s] for s in cols_wt]
    rt, rp = welch(E, mut, wt)
    rt_s = pd.Series(rt, index=reg_ids)
    edge_t, edge_p = {}, {}
    for r in regulators:
        member_of = reg_members.get(r)
        if member_of:   # regulator is a gene in some regulons: mean t / mean -log10 p over them
            ts = rt_s[member_of].values
            ps = rp[[_G["reg_pos"][i] for i in member_of]]
            edge_t[r] = np.mean(ts)
            edge_p[r] = 10 ** (-np.mean(-np.log10(np.clip(ps, 1e-300, 1))))
        elif r in _G["gene_index"]:
            x = X[_G["gene_index"][r]]
            t, pv = stats.ttest_ind(x[mut], x[wt], equal_var=False)
            edge_t[r], edge_p[r] = t, pv
    return (pd.DataFrame({"t": rt, "p": rp}, index=reg_ids),
            pd.DataFrame({"t": pd.Series(edge_t), "p": pd.Series(edge_p)}))


def gene_arms(probemap, cytoband, back):
    """MINER gene ID -> chromosome arm (hg38; GENCODE v36 positions, UCSC cytoband centromeres)."""
    pm = pd.read_csv(probemap, sep="\t")
    pm["ens"] = pm["id"].str.split(".").str[0]
    pm = pm.drop_duplicates("ens").set_index("ens")
    cb = pd.read_csv(cytoband, sep="\t", header=None, names=["c", "s", "e", "b", "st"])
    cen = cb[cb["st"] == "acen"].groupby("c")["s"].min()
    arm = {}
    for g, r in pm.iterrows():
        if r["chrom"] in cen.index:
            arm[g] = r["chrom"][3:] + ("p" if r["chromStart"] < cen[r["chrom"]] else "q")
    inv = {v: k for k, v in back.items()}  # our ID -> MINER ID
    out = dict(arm)
    for ours, a_ in arm.items():
        if ours in inv:
            out[inv[ours]] = a_
    return out


def run_miner_feature(args):
    """MINER causalNetworkAnalysis for one feature on its profiled samples (worker)."""
    from miner import causal_inference as causalinf
    feature, row = args
    mm = pd.DataFrame([row.astype(int).values], index=[feature], columns=row.index)
    causalinf.causalNetworkAnalysis(regulon_matrix=_G["regulon_df"], expression_matrix=_G["exp_data"],
                                    reference_matrix=_G["eigengenes"], mutation_matrix=mm,
                                    resultsDirectory=_G["result_dir"], minRegulons=1,
                                    significance_threshold=_G["alpha"], causalFolder="pooled")
    return feature


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    ap.add_argument("--cores", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    ap.add_argument("--features", default=None, help="comma-separated subset (testing)")
    args = ap.parse_args()
    P = load_params(args.params)
    M, C = P["miner"], P["causal"]
    matrix = args.matrix or M["matrix"]
    res = p(P["paths"]["results"])
    mdir = os.path.join(res, "04_miner", matrix)
    outdir = os.path.join(res, OUT, matrix)
    log = setup_logging(outdir, "05_causal_inference")
    from miner import miner, causal_inference as causalinf

    # ---- expression and regulons, in MINER IDs (as miner3-causalinference does)
    expfile = os.path.join(res, "02_batch_corrected", f"expression_{matrix}_z.csv")
    exp_data, _ = miner.preprocess(expfile, p(M["idmap"]), do_preprocess_tpm=False)
    regulons = json.load(open(os.path.join(mdir, "mechinf", "regulons_filtered.json")))
    rdf = pd.read_csv(os.path.join(mdir, "mechinf", "regulonDf.csv"), index_col=0)
    rdf["Regulon_ID"] = rdf["Regulon_ID"].astype(str)
    regulon_df = rdf[rdf["Regulon_ID"].isin(regulons)].reset_index(drop=True)
    eigengenes = miner.getEigengenes(regulons, exp_data, regulon_dict=None, saveFolder=None)
    eigengenes = eigengenes * (np.percentile(exp_data, 95) / np.percentile(eigengenes, 95))
    eigengenes.index = np.array(eigengenes.index).astype(str)
    log.info("Expression %d x %d; %d regulons, %d regulators", *exp_data.shape, len(regulons),
             regulon_df["Regulator"].nunique())

    # ---- genomic features (columns = expression sample IDs)
    F = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
    F = F.loc[:, F.columns.isin(exp_data.columns)]
    if args.features:
        F = F.loc[args.features.split(",")]
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    cohort = samples["cohort"]
    log.info("%d genomic features", len(F))

    # ---- symbols and IDs
    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    back = miner_id_backmap(p(M["idmap"]), genes.index)
    def sym(g):
        s_ = genes["symbol"].get(back.get(g, g))
        return s_ if isinstance(s_, str) and s_ else g

    # ---- shared state for workers
    reg_ids = list(eigengenes.index)
    reg_members = regulon_df.groupby("Gene")["Regulon_ID"].apply(lambda s: sorted(set(s) & set(reg_ids))).to_dict()
    _G.update(dict(E=eigengenes.values, X=exp_data.values, reg_ids=reg_ids,
                   reg_pos={r: i for i, r in enumerate(reg_ids)}, reg_members=reg_members,
                   regulators=sorted(regulon_df["Regulator"].unique()),
                   gene_index={g: i for i, g in enumerate(exp_data.index)},
                   sample_index={s: i for i, s in enumerate(exp_data.columns)},
                   regulon_df=regulon_df, exp_data=exp_data, eigengenes=eigengenes,
                   result_dir=os.path.join(outdir, "causal_results"), alpha=C["alpha"]))
    os.makedirs(os.path.join(outdir, "causal_results", "pooled"), exist_ok=True)

    # ---- 1. MINER causal analysis, one feature at a time on profiled samples
    t0 = time.time()
    tasks = [(f, F.loc[f].dropna()) for f in F.index]
    with Pool(args.cores) as pool:
        for i, f in enumerate(pool.imap_unordered(run_miner_feature, tasks), 1):
            if i % 10 == 0 or i == len(tasks):
                log.info("MINER causal analysis: %d / %d features (%.1f min)", i, len(tasks), (time.time() - t0) / 60)
    files = os.listdir(os.path.join(outdir, "causal_results", "pooled"))
    if not files:
        log.warning("MINER found no causal flows for any feature")
        return
    cr = causalinf.readCausalFiles(os.path.join(outdir, "causal_results"))
    cr = cr.reset_index(drop=True)  # MINER indexes rows by regulon ID, which repeats across features
    num = [c for c in cr.columns if c not in ("Mutation", "Regulator", "Regulon")]
    cr[num] = cr[num].apply(pd.to_numeric, errors="coerce")
    cr["Regulon_ID"] = cr["Regulon"].str.replace("R-", "", regex=False)
    log.info("MINER flows: %d rows over %d features", len(cr), cr["Mutation"].nunique())

    # ---- 2 + 3. pooled q-values, per-cohort statistics, effect size, adjusted test
    feat_type_all = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features_info.tsv"),
                                sep="\t", index_col=0)["type"]
    Fa = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
    burden = (Fa.loc[[x for x in Fa.index if x.startswith("ARM_")]] == 1).sum(axis=0)
    Emat = eigengenes

    def cohen_d(mut, wt):
        x, y = Emat[mut].values, Emat[wt].values
        sp = np.sqrt(((len(mut) - 1) * x.var(1, ddof=1) + (len(wt) - 1) * y.var(1, ddof=1)) / (len(mut) + len(wt) - 2))
        return pd.Series((x.mean(1) - y.mean(1)) / sp, index=Emat.index)

    def adjusted_test(f, row, ftype, burden, cohort):
        """OLS per regulon: eigengene ~ feature + cohort dummies [+ other arm events for CNA features]."""
        S = list(row.index)
        x = row.values.astype(float)
        cols = [np.ones(len(S)), x]
        label = "cohort"
        if ftype in ("arm_cna", "focal_cna"):
            other = burden.reindex(S).values - (x if ftype == "arm_cna" else 0)
            cols.append(other)
            label = "cohort+arm_burden"
        coh = pd.get_dummies(cohort.reindex(S)).values.astype(float)
        if coh.shape[1] > 1:
            cols += list(coh[:, 1:].T)
        X = np.column_stack(cols)
        Y = Emat[S].values
        XtX = np.linalg.pinv(X.T @ X)
        B = Y @ X @ XtX
        R = Y - B @ X.T
        df = X.shape[0] - np.linalg.matrix_rank(X)
        t = B[:, 1] / np.sqrt((R ** 2).sum(1) / df * XtX[1, 1])
        return t, 2 * stats.t.sf(np.abs(t), df), label

    ann = []
    for f in cr["Mutation"].unique():
        row = F.loc[f].dropna()
        mut, wt = list(row.index[row == 1]), list(row.index[row == 0])
        rs, es = feature_stats(f, mut, wt)
        rs["q"] = bh(rs["p"].values)
        es["q"] = bh(es["p"].values)
        sub = cr[cr["Mutation"] == f]
        a = pd.DataFrame(index=sub.index)
        a["q_regulon"] = rs["q"].reindex(sub["Regulon_ID"]).values
        a["q_edge"] = es["q"].reindex(sub["Regulator"]).values
        pooled_rt = rs["t"].reindex(sub["Regulon_ID"]).values
        pooled_et = es["t"].reindex(sub["Regulator"]).values
        tested, agree = np.zeros(len(sub), int), np.zeros(len(sub), int)
        for c in ["TCGA", "CLCA", "LICA_FR"]:
            cm = [s for s in mut if cohort.get(s) == c]
            cw = [s for s in wt if cohort.get(s) == c]
            if len(cm) < C["min_cohort_altered"] or len(cw) < C["min_cohort_altered"]:
                a[f"t_regulon_{c}"] = np.nan
                a[f"t_edge_{c}"] = np.nan
                continue
            crs, ces = feature_stats(f, cm, cw)
            trc = crs["t"].reindex(sub["Regulon_ID"]).values
            tec = ces["t"].reindex(sub["Regulator"]).values
            a[f"t_regulon_{c}"], a[f"t_edge_{c}"] = trc, tec
            ok = (np.sign(trc) == np.sign(pooled_rt)) & (np.sign(tec) == np.sign(pooled_et))
            tested += 1
            agree += ok.astype(int)
        a["n_cohorts_tested"], a["n_cohorts_consistent"] = tested, agree
        a["cohort_consistent"] = (tested >= C["min_cohorts"]) & (agree == tested)
        # effect size (pooled Cohen's d of the regulon eigengene, altered vs wild type)
        a["cohen_d"] = cohen_d(mut, wt).reindex(sub["Regulon_ID"]).values
        # adjusted regulon test: cohort (+ arm-event burden for CNA features, excluding the feature itself)
        adj_t, adj_p, adj_label = adjusted_test(f, row, feat_type_all.get(f, ""), burden, cohort)
        adj_q = pd.Series(bh(adj_p), index=reg_ids)
        a["q_regulon_adjusted"] = adj_q.reindex(sub["Regulon_ID"]).values
        a["t_regulon_adjusted"] = pd.Series(adj_t, index=reg_ids).reindex(sub["Regulon_ID"]).values
        a["adjusted_for"] = adj_label
        ann.append(a)
    cr = cr.join(pd.concat(ann))

    # ---- 4. families and programs
    E = eigengenes.values
    D = np.clip(1 - np.corrcoef(E), 0, 2)
    np.fill_diagonal(D, 0)
    fam = fcluster(linkage(squareform(D, checks=False), "average"), 1 - C["family_r"], "distance")
    progs = json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_programs.json")))
    r2p = {r: k for k, v in progs.items() for r in v}
    reg_of = regulon_df.groupby("Regulon_ID")["Regulator"].first()
    fam_df = pd.DataFrame({"family": [f"F{x}" for x in fam], "program": [r2p.get(r, "") for r in reg_ids],
                           "regulator": [reg_of.get(r, "") for r in reg_ids]}, index=pd.Index(reg_ids, name="regulon"))
    fam_df["regulator_symbol"] = fam_df["regulator"].map(sym)
    fam_df["n_genes"] = [len(regulons[r]) for r in reg_ids]
    fam_df.to_csv(os.path.join(outdir, "regulon_families.tsv"), sep="\t")
    log.info("Regulon families at r >= %.2f: %d (from %d regulons)", C["family_r"], fam_df["family"].nunique(), len(fam_df))

    cr["regulator_symbol"] = cr["Regulator"].map(sym)
    cr["family"] = cr["Regulon_ID"].map(fam_df["family"])
    cr["program"] = cr["Regulon_ID"].map(fam_df["program"])
    feat_type = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features_info.tsv"),
                            sep="\t", index_col=0)["type"]
    cr["feature_type"] = cr["Mutation"].map(feat_type)
    cr["n_altered"] = cr["Mutation"].map((F == 1).sum(axis=1))
    # gene-dosage flag for arm features: >= 50% of the regulon's genes lie on the altered arm
    garm = gene_arms(p(C["gene_positions"]), p(C["cytoband"]), back)
    is_arm = cr["feature_type"] == "arm_cna"
    arm_name = cr["Mutation"].str.split("_").str[1]
    cr["frac_genes_on_arm"] = np.nan
    cr.loc[is_arm, "frac_genes_on_arm"] = [np.mean([garm.get(g) == a_ for g in regulons[r]])
                                           for r, a_ in zip(cr.loc[is_arm, "Regulon_ID"], arm_name[is_arm])]
    cr["cis_dosage"] = cr["frac_genes_on_arm"] >= 0.5
    cr.to_csv(os.path.join(outdir, "completeCausalResults.csv"))

    # MINER's CLI filter (same four criteria as miner3-causalinference)
    a = C["alpha"]
    filt = cr[(cr["-log10(p)_Regulon_stratification"] >= -np.log10(a))
              & (cr["Fraction_of_edges_correctly_aligned"] >= 0.5)
              & (cr["RegulatorRegulon_Spearman_p-value"] <= a)
              & (cr["-log10(p)_MutationRegulatorEdge"] >= -np.log10(a))]
    filt.to_csv(os.path.join(outdir, "filteredCausalResults.csv"))
    hc = filt[(filt["q_regulon"] <= C["q_max"]) & (filt["q_edge"] <= C["q_max"]) & filt["cohort_consistent"]
              & (filt["cohen_d"].abs() >= C["min_abs_d"]) & (filt["n_altered"] >= C["min_altered_hc"])]
    # rank by effect size: strongest shift of the regulon first
    hc = hc.assign(abs_d=hc["cohen_d"].abs()).sort_values("abs_d", ascending=False).drop(columns="abs_d")
    hc.to_csv(os.path.join(outdir, "highConfidenceCausalResults.csv"))
    log.info("Features below min_altered_hc (%d altered), excluded from high-confidence: %s", C["min_altered_hc"],
             sorted(set(filt.loc[filt["n_altered"] < C["min_altered_hc"], "Mutation"])))
    log.info("Flows: MINER %d, MINER-filtered %d, high-confidence %d (q <= %.2f both tests, consistent in all "
             "tested cohorts, >= %d cohorts, |d| >= %.2f, >= %d altered); of these, adjusted q <= %.2f: %d; cis-dosage: %d",
             len(cr), len(filt), len(hc), C["q_max"], C["min_cohorts"], C["min_abs_d"], C["min_altered_hc"], C["q_max"],
             int((hc["q_regulon_adjusted"] <= C["q_max"]).sum()), int(hc["cis_dosage"].sum()))

    # ---- summaries
    def by_family(df, level):
        g = df.assign(abs_d=df["cohen_d"].abs()).sort_values("abs_d", ascending=False).groupby(["Mutation", "family"])
        out = g.agg(program=("program", "first"), n_flows=("Regulon", "size"), n_regulons=("Regulon", "nunique"),
                    regulators=("regulator_symbol", lambda s: ",".join(pd.unique(s)[:10])),
                    n_regulators=("Regulator", "nunique"),
                    best_regulon=("Regulon", "first"), best_regulator=("regulator_symbol", "first"),
                    direction=("Regulon_stratification_t-statistic", lambda s: "up" if s.iloc[0] > 0 else "down"),
                    best_neglog10p=("-log10(p)_Regulon_stratification", "first"),
                    best_q_regulon=("q_regulon", "min"),
                    best_cohen_d=("cohen_d", "first"),
                    frac_adjusted_q_ok=("q_regulon_adjusted", lambda s: float((s <= C["q_max"]).mean())),
                    frac_cohort_consistent=("cohort_consistent", "mean"),
                    frac_cis_dosage=("cis_dosage", "mean"))
        out["level"] = level
        out = out.reset_index()
        out["rank_in_feature"] = out.groupby("Mutation")["best_cohen_d"].transform(
            lambda v: v.abs().rank(ascending=False, method="min")).astype(int)
        return out.sort_values(["Mutation", "rank_in_feature"])

    fam_tab = pd.concat([by_family(filt, "miner_filtered"), by_family(hc, "high_confidence")] if len(hc)
                        else [by_family(filt, "miner_filtered")])
    fam_tab.to_csv(os.path.join(outdir, "causal_by_family.tsv"), sep="\t", index=False)
    rows = []
    for f in F.index:
        r = {"feature": f, "type": feat_type.get(f), "n_altered": int((F.loc[f] == 1).sum()),
             "n_profiled": int(F.loc[f].notna().sum())}
        for name, df in (("miner", cr), ("filtered", filt), ("high_conf", hc)):
            d = df[df["Mutation"] == f]
            r[f"{name}_flows"] = len(d)
            r[f"{name}_regulators"] = d["Regulator"].nunique()
            r[f"{name}_families"] = d["family"].nunique()
            r[f"{name}_programs"] = d["program"].nunique()
        rows.append(r)
    feat_tab = pd.DataFrame(rows).sort_values("high_conf_families", ascending=False)
    feat_tab.to_csv(os.path.join(outdir, "causal_by_feature.tsv"), sep="\t", index=False)
    log.info("Top features by high-confidence families:\n%s", feat_tab.head(25).to_string(index=False))

    # ---- MINER wiring diagram on the filtered flows
    try:
        coher = pd.read_csv(os.path.join(mdir, "subtypes_filtered", "coherentMembers.csv"), index_col=0, header=0)
        # wiringDiagram reads the row index as the regulon ID and columns by position (MINER's order)
        causalinf.wiringDiagram(filt.set_index("Regulon_ID"), regulons, coher, include_genes=False,
                                savefile=os.path.join(outdir, "wiring_diagram.csv"))
    except Exception as e:  # noqa: BLE001 - MINER's wiringDiagram is fragile; results above stand on their own
        log.warning("wiringDiagram failed: %s", e)

    try:
        import qc_plots
        written = qc_plots.causal_report(outdir, feat_tab, fam_tab)
        log.info("QC figures: %s", ", ".join(written))
    except Exception as e:  # noqa: BLE001
        log.warning("figures failed: %s", e)


if __name__ == "__main__":
    main()
