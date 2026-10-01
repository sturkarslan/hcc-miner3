#!/usr/bin/env python
"""Step 07b: map MINER transcriptional states and programs onto published liver-cancer subtypes.

1. Sample calls: Nearest Template Prediction (Hoshida 2010) for each classifier in
   results/07_post/signatures/classifiers.json. Cosine similarity of each sample (template genes,
   z-scored) to each class template (+1 up, -1 dn); p-value from permuting gene labels; BH FDR over
   samples; FDR >= post.ntp.fdr -> "unassigned". Calls are checked against the LICA-FR labels.
2. States (sample groups): contingency with every classifier call, LICA-FR label and cohort;
   one-sided Fisher enrichment per state x level (BH within annotation); adjusted Rand index.
3. Programs and regulons (gene sets): hypergeometric overlap with every library signature
   (background = genes in MINER coexpression modules); BH over all tests.
4. Program activity vs signature score: Pearson r across samples between program activity (mean
   regulon eigengene) and signature score (mean z of its genes), after regressing each sample's mean z
   over all genes out of both (a per-sample level that otherwise correlates everything with everything);
   the unadjusted matrix is kept as program_signature_correlation_raw.tsv.

MINER IDs are mapped back to project Ensembl IDs (hcc_common.miner_id_backmap).
Inputs: results/04_miner/<matrix>/{<post.subtypes_dir>,mechinf}/, results/02_batch_corrected/,
        results/01_harmonized/, results/07_post/signatures/
Outputs (results/07_post/subtype_mapping/<subtypes_dir>/): ntp_calls_<classifier>.tsv,
        ntp_vs_labels.tsv, state_enrichment.tsv, state_annotation.tsv, state_ari.tsv,
        program_signature_overlap.tsv, regulon_signature_overlap.tsv,
        program_signature_correlation.tsv, qc/p1-p4 PNGs
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import adjusted_rand_score

from hcc_common import load_params, miner_id_backmap, p, setup_logging


def bh(pv):
    pv = np.asarray(pv, float)
    out = np.full_like(pv, np.nan)
    ok = ~np.isnan(pv)
    n = ok.sum()
    if n:
        order = np.argsort(pv[ok])
        ranked = pv[ok][order] * n / np.arange(1, n + 1)
        q = np.minimum.accumulate(ranked[::-1])[::-1]
        tmp = np.empty(n)
        tmp[order] = np.minimum(q, 1)
        out[ok] = tmp
    return out


def zscore_rows(df):
    sd = df.std(axis=1).replace(0, np.nan)
    return df.sub(df.mean(axis=1), axis=0).div(sd, axis=0).fillna(0)


def ntp(z, classes, nperm, seed, fdr_cut):
    """z: genes x samples. classes: {cls: {up: [...], dn: [...]}}."""
    genes = sorted({g for c in classes.values() for g in c["up"] + c["dn"]} & set(z.index))
    T = pd.DataFrame(0.0, index=list(classes), columns=genes)
    for cls, c in classes.items():
        T.loc[cls, [g for g in c["up"] if g in genes]] += 1
        T.loc[cls, [g for g in c["dn"] if g in genes]] -= 1
    X = z.loc[genes].values
    Xn = X / (np.linalg.norm(X, axis=0, keepdims=True) + 1e-12)
    Tn = T.values / (np.linalg.norm(T.values, axis=1, keepdims=True) + 1e-12)
    S = Tn @ Xn                                     # classes x samples
    best = S.max(axis=0)
    # Null as in Hoshida's NTP: the template scored against random gene sets of the same size drawn from
    # all genes (permuting only the template genes gives no null for one-direction signatures, whose
    # similarity is invariant to that permutation)
    rng = np.random.default_rng(seed)
    Z = z.values
    exceed = np.zeros(X.shape[1])
    for _ in range(nperm):
        Xr = Z[rng.choice(Z.shape[0], size=len(genes), replace=False)]
        Xr = Xr / (np.linalg.norm(Xr, axis=0, keepdims=True) + 1e-12)
        exceed += (Tn @ Xr).max(axis=0) >= best
    pval = (exceed + 1) / (nperm + 1)
    fdr = bh(pval)
    best_cls = np.array(T.index)[S.argmax(axis=0)]
    out = pd.DataFrame(S.T, index=z.columns, columns=[f"sim_{c}" for c in T.index])
    out.insert(0, "call", np.where(fdr < fdr_cut, best_cls, "unassigned"))
    out.insert(1, "nearest", best_cls)
    out.insert(2, "p", pval)
    out.insert(3, "fdr", fdr)
    return out, len(genes)


def enrichment_table(groups, labels, name):
    """groups, labels: Series indexed by sample. One-sided Fisher per group x level."""
    d = pd.DataFrame({"g": groups, "l": labels}).dropna()
    d = d[d["l"] != "unassigned"]
    rows = []
    for g in sorted(d["g"].unique(), key=lambda v: (len(str(v)), str(v))):
        for lv in sorted(d["l"].unique()):
            a = int(((d.g == g) & (d.l == lv)).sum())
            b = int(((d.g == g) & (d.l != lv)).sum())
            c = int(((d.g != g) & (d.l == lv)).sum())
            e = int(((d.g != g) & (d.l != lv)).sum())
            orr, pv = stats.fisher_exact([[a, b], [c, e]], alternative="greater")
            rows.append((name, g, lv, a + b, a + c, a, a / (a + b) if a + b else np.nan, orr, pv))
    t = pd.DataFrame(rows, columns=["annotation", "state", "level", "n_state", "n_level", "n_both",
                                    "frac_of_state", "odds_ratio", "p"])
    t["fdr"] = bh(t["p"])
    ari = adjusted_rand_score(d["g"].astype(str), d["l"].astype(str)) if len(d) else np.nan
    return t, ari, len(d)


def overlap_tests(gene_sets, sigs, background, min_genes, what):
    """gene_sets: {id: set(genes)}; sigs: {set: set(genes)} (already within background)."""
    N = len(background)
    rows = []
    for k, gs in gene_sets.items():
        gs = gs & background
        if not gs:
            continue
        for s, sg in sigs.items():
            if len(sg) < min_genes:
                continue
            ov = gs & sg
            if not ov:
                continue
            pv = stats.hypergeom.sf(len(ov) - 1, N, len(sg), len(gs))
            expct = len(gs) * len(sg) / N
            rows.append((k, s, len(gs), len(sg), len(ov), len(ov) / expct, pv, ov))
    t = pd.DataFrame(rows, columns=[what, "signature", f"{what}_genes", "signature_genes", "overlap",
                                    "fold_enrichment", "p", "overlap_ids"])
    t["fdr"] = bh(t["p"])
    return t.sort_values("p")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    ap.add_argument("--subtypes-dir", default=None, help="override post.subtypes_dir")
    args = ap.parse_args()
    P = load_params(args.params)
    Q = P["post"]
    matrix = args.matrix or P["miner"]["matrix"]
    sub = args.subtypes_dir or Q.get("subtypes_dir", "subtypes")
    res = p(P["paths"]["results"])
    mdir = os.path.join(res, "04_miner", matrix)
    outdir = os.path.join(res, "07_post", "subtype_mapping", sub)
    log = setup_logging(outdir, "07b_subtype_mapping")
    sigdir = os.path.join(res, "07_post", "signatures")

    z = pd.read_csv(os.path.join(res, "02_batch_corrected", f"expression_{matrix}_z.csv"), index_col=0)
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample").loc[z.columns]
    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    sym = genes["symbol"].fillna("")
    label_cols = [c for c in samples.columns if c not in ("cohort", "patient")]
    back = miner_id_backmap(p(P["miner"]["idmap"]), z.index)
    fix = lambda gs: {back.get(g, g) for g in gs}  # noqa: E731

    # ---- 1. NTP calls
    N = Q.get("ntp", {})
    zz = z
    if N.get("center", "all") == "cohort":
        zz = pd.concat([zscore_rows(z.loc[:, samples["cohort"] == c]) for c in samples["cohort"].unique()], axis=1)[z.columns]
    classifiers = json.load(open(os.path.join(sigdir, "classifiers.json")))
    calls = {}
    val_rows = []
    for cl, classes in classifiers.items():
        t, ng = ntp(zz, classes, int(N.get("permutations", 1000)), int(N.get("seed", 1)), float(N.get("fdr", 0.05)))
        t.to_csv(os.path.join(outdir, f"ntp_calls_{cl}.tsv"), sep="\t", float_format="%.4g")
        calls[cl] = t["call"]
        log.info("NTP %s (%d template genes): %s", cl, ng,
                 pd.crosstab(samples["cohort"], t["call"]).to_string().replace("\n", "\n    "))
        for lab in label_cols:
            d = pd.DataFrame({"call": t["call"], "label": samples[lab]}).dropna()
            d = d[d["call"] != "unassigned"]
            if d["label"].nunique() < 2 or len(d) < 10:
                continue
            same = d["call"].astype(str) == d["label"].astype(str)
            val_rows.append({"classifier": cl, "label": lab, "n": len(d),
                             "ari": adjusted_rand_score(d["label"].astype(str), d["call"]),
                             "exact_match": same.mean() if same.any() else np.nan})
    # Montironi 2023 immune classes (Suppl Fig 19): Inflamed signature, then Sia immune class within
    # inflamed tumours and CTNNB1 mutation within non-inflamed tumours
    MC = Q.get("montironi")
    if MC and MC["inflamed_classifier"] in calls and MC["immune_classifier"] in calls:
        F = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
        ctn = F.loc[MC["ctnnb1_feature"]].reindex(z.columns) if MC["ctnnb1_feature"] in F.index else \
            pd.Series(np.nan, index=z.columns)
        infl = calls[MC["inflamed_classifier"]] == MC["inflamed_class"]
        imm = calls[MC["immune_classifier"]] == MC["immune_class"]
        mc = pd.Series(np.where(infl, np.where(imm, "Immune", "Immune-like"),
                                np.where(ctn == 1, "Excluded", np.where(ctn == 0, "Intermediate", None))),
                       index=z.columns)
        calls["montironi"] = mc
        mc.rename("call").to_frame().assign(cohort=samples["cohort"], ctnnb1=ctn).to_csv(
            os.path.join(outdir, "montironi_classes.tsv"), sep="\t")
        # same file layout as the NTP classifiers, so downstream steps load every class the same way
        mc.fillna("unassigned").rename("call").to_frame().to_csv(os.path.join(outdir, "ntp_calls_montironi.tsv"), sep="\t")
        log.info("Montironi classes (CTNNB1-unprofiled non-inflamed tumours left missing):\n%s",
                 pd.crosstab(samples["cohort"], mc.fillna("missing")).to_string())
        for lab in label_cols:
            d = pd.DataFrame({"call": mc, "label": samples[lab]}).dropna()
            if d["label"].nunique() >= 2 and len(d) >= 10:
                val_rows.append({"classifier": "montironi", "label": lab, "n": len(d),
                                 "ari": adjusted_rand_score(d["label"].astype(str), d["call"]), "exact_match": np.nan})
    val = pd.DataFrame(val_rows)
    val.to_csv(os.path.join(outdir, "ntp_vs_labels.tsv"), sep="\t", index=False, float_format="%.3f")
    if len(val):
        log.info("NTP calls vs LICA-FR labels:\n%s", val.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

    # ---- 2. states
    states = json.load(open(os.path.join(mdir, sub, "transcriptional_states.json")))
    state_of = pd.Series({s: k for k, v in states.items() for s in v}).reindex(z.columns)
    log.info("%d states; %d of %d samples assigned", len(states), state_of.notna().sum(), len(state_of))
    annots = dict(**{f"ntp_{k}": v for k, v in calls.items()},
                  **{lab: samples[lab] for lab in label_cols}, cohort=samples["cohort"])
    enr, ari = [], []
    for name, lab in annots.items():
        t, a, n = enrichment_table(state_of, lab, name)
        enr.append(t)
        ari.append({"annotation": name, "n_samples": n, "ari": a})
    enr = pd.concat(enr, ignore_index=True)
    enr.to_csv(os.path.join(outdir, "state_enrichment.tsv"), sep="\t", index=False, float_format="%.4g")
    ari = pd.DataFrame(ari)
    ari.to_csv(os.path.join(outdir, "state_ari.tsv"), sep="\t", index=False, float_format="%.3f")
    log.info("States vs annotations (ARI):\n%s", ari.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    sig_enr = enr[(enr["fdr"] < 0.05) & (enr["odds_ratio"] > 1)]
    sig_enr = sig_enr.sort_values("fdr").assign(
        txt=lambda d: [f"{lv} ({f:.0%})" for lv, f in zip(d["level"], d["frac_of_state"])])
    ann = sig_enr.groupby(["state", "annotation"])["txt"].agg(";".join).unstack("annotation")
    ann = ann.reindex(sorted(states, key=int))
    ann.insert(0, "n_samples", [len(states[s]) for s in ann.index])
    ann.index.name = "state"
    ann.to_csv(os.path.join(outdir, "state_annotation.tsv"), sep="\t")

    # ---- 3. programs / regulons vs signature gene sets
    reg_file = "regulons_filtered.json" if sub.endswith("filtered") else "regulons.json"
    regulons = {k: fix(v) for k, v in json.load(open(os.path.join(mdir, "mechinf", reg_file))).items()}
    programs = json.load(open(os.path.join(mdir, sub, "transcriptional_programs.json")))
    prog_genes = {k: set().union(*(regulons[str(r)] for r in v if str(r) in regulons)) for k, v in programs.items()}
    background = set().union(*(fix(v) for v in json.load(open(os.path.join(mdir, "mechinf", "coexpressionDictionary.json"))).values()))
    sig = pd.read_csv(os.path.join(sigdir, "signatures.tsv"), sep="\t").dropna(subset=["ensembl"])
    sigs = {s: set(d["ensembl"]) & background for s, d in sig.groupby("set")}
    mg = int(Q.get("min_set_genes", 5))
    log.info("Background %d genes; %d programs; %d regulons (%s); %d signatures with >= %d background genes",
             len(background), len(programs), len(regulons), reg_file, sum(len(v) >= mg for v in sigs.values()), mg)
    to_sym = lambda ids: ",".join(sorted(sym.get(g, g) or g for g in ids))  # noqa: E731
    pt = overlap_tests(prog_genes, sigs, background, mg, "program")
    pt["overlap_genes"] = pt.pop("overlap_ids").map(to_sym)
    pt.to_csv(os.path.join(outdir, "program_signature_overlap.tsv"), sep="\t", index=False, float_format="%.4g")
    rt = overlap_tests(regulons, sigs, background, mg, "regulon")
    rt["overlap_genes"] = rt.pop("overlap_ids").map(to_sym)
    rt.to_csv(os.path.join(outdir, "regulon_signature_overlap.tsv"), sep="\t", index=False, float_format="%.4g")
    top = pt[pt["fdr"] < 0.05].sort_values("fdr").groupby("program").head(3)
    log.info("Programs: %d of %d have >= 1 signature at FDR < 0.05. Top per program:\n%s",
             pt.loc[pt.fdr < 0.05, "program"].nunique(), len(programs),
             top[["program", "signature", "overlap", "fold_enrichment", "fdr"]].to_string(index=False))

    # ---- 4. program activity vs signature score
    eig = pd.read_csv(os.path.join(mdir, sub, "eigengenes.csv"), index_col=0)
    eig.index = eig.index.astype(str)
    act = pd.DataFrame({k: eig.loc[[str(r) for r in v if str(r) in eig.index]].mean() for k, v in programs.items()})
    act = act.reindex(z.columns)
    score = pd.DataFrame({s: z.loc[sorted(g & set(z.index))].mean() for s, g in sigs.items() if len(g) >= mg})

    def corr(a, b):
        return pd.DataFrame(np.corrcoef(a.T.values, b.T.values)[:a.shape[1], a.shape[1]:], index=a.columns,
                            columns=b.columns).rename_axis("program")

    corr(act, score).to_csv(os.path.join(outdir, "program_signature_correlation_raw.tsv"), sep="\t", float_format="%.3f")
    # Each sample's mean z over all genes drives most program activities (median r 0.62) and signature
    # scores (median r 0.71), and tracks genes detected per sample (r 0.68; PROJECT_LOG). Regress it out of
    # both before correlating, so r reflects shared biology rather than the per-sample level.
    gm = z.mean(axis=0).reindex(act.index)
    gc = (gm - gm.mean()).values[:, None]

    def resid(X):
        X = X.loc[gm.index]
        return X - gc @ ((gc * (X - X.mean()).values).sum(0, keepdims=True) / (gc ** 2).sum())

    cor = corr(resid(act), resid(score))
    cor.to_csv(os.path.join(outdir, "program_signature_correlation.tsv"), sep="\t", float_format="%.3f")
    log.info("Program x signature r: raw median %.2f; global-mean-adjusted median %.2f, median best |r| %.2f",
             np.median(corr(act, score).values), np.median(cor.values), cor.abs().max(1).median())

    import qc_plots
    written = qc_plots.subtype_report(outdir, enr, ari, pt, cor, calls, samples, label_cols, states)
    log.info("QC figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
