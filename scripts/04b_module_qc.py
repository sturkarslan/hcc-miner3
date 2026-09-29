#!/usr/bin/env python
"""Step 04b: QC of MINER3 coexpression modules (and regulons, once mechinf has run).

Answers "is a very large module a problem?" with data. For every coexpression module:
  - coherence: PC1 and PC2 variance fraction, median gene-eigengene r, fraction of genes
    loosely attached (r < 0.3). A merged module of distinct signals has a weak PC1 and a strong PC2.
  - what the eigengene tracks: cohort (eta^2; should be ~0 after ComBat), LICA-FR labels,
    per-sample technical stats within cohort (TPM share in shared genes, genes detected),
    and marker programs (T cell, B/plasma, myeloid, stromal, endothelial, interferon,
    proliferation, hepatocyte, CTNNB1), plus ribosomal / mitochondrial / histone gene content.
If mechinf output exists, each regulon is assigned to the module holding most of its genes and
redundancy is reported per module (regulon count, regulators, mean pairwise Jaccard, mean
pairwise regulon-eigengene correlation).

Inputs: results/04_miner/<matrix>/{coexpr,mechinf}/, results/02_batch_corrected/expression_<matrix>_z.csv,
        results/01_harmonized/{samples,genes,sample_stats}.tsv
Outputs (results/04_miner/<matrix>/module_qc/): modules.tsv, regulons_by_module.tsv,
        module_<id>_genes.tsv for the largest modules, qc/m1-m4 PNGs
"""

import argparse
import json
import os

import numpy as np
import pandas as pd

from hcc_common import load_params, miner_id_backmap, p, setup_logging

MARKERS = {
    "T_cell": ["CD2", "CD3D", "CD3E", "CD8A", "GZMK", "LCK", "CCL5", "CXCL9", "CXCL10"],
    "B_plasma": ["CD79A", "MS4A1", "IGHM", "JCHAIN", "MZB1"],
    "myeloid": ["CD68", "CD163", "CSF1R", "C1QA", "C1QB", "LYZ", "MS4A7"],
    "stromal": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "THY1"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "KDR", "ENG"],
    "interferon": ["ISG15", "IFI6", "IFI44L", "MX1", "OAS1", "RSAD2", "IFIT1"],
}
GENE_CLASSES = {"ribosomal": r"^RP[LS]\d", "mitochondrial": r"^MT-", "histone": r"^H[1-4]"}


def eigengenes(z, modules):
    """PC1 per module (samples), sign-aligned to the module mean; plus PC1/PC2 variance and
    gene-eigengene correlations."""
    eig, stats = {}, {}
    zc = z.sub(z.mean(axis=1), axis=0)
    for k, genes in modules.items():
        x = zc.loc[genes].values
        u, s, vt = np.linalg.svd(x, full_matrices=False)
        var = s ** 2 / (s ** 2).sum()
        e = vt[0] * s[0]
        if np.corrcoef(e, x.mean(axis=0))[0, 1] < 0:
            e = -e
        xs = (x - x.mean(1, keepdims=True)) / (x.std(1, keepdims=True) + 1e-12)
        es = (e - e.mean()) / e.std()
        r = xs @ es / x.shape[1]
        eig[k] = e
        stats[k] = {"n_genes": len(genes), "pc1_var": var[0], "pc2_var": var[1] if len(var) > 1 else np.nan,
                    "median_r_to_eigengene": float(np.median(r)), "frac_genes_r_lt_0.3": float((r < 0.3).mean())}
    return pd.DataFrame(eig, index=z.columns), pd.DataFrame(stats).T, zc


def eta2(v, g):
    g = pd.Series(np.asarray(g))
    ok = g.notna().values
    v, g = np.asarray(v)[ok], g[ok].values
    tot = ((v - v.mean()) ** 2).sum()
    b = sum((g == lv).sum() * (v[g == lv].mean() - v.mean()) ** 2 for lv in pd.unique(g))
    return b / tot if tot > 0 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    ap.add_argument("--top", type=int, default=5, help="largest modules to profile in detail")
    args = ap.parse_args()
    P = load_params(args.params)
    matrix = args.matrix or P["miner"]["matrix"]
    res = p(P["paths"]["results"])
    mdir = os.path.join(res, "04_miner", matrix)
    outdir = os.path.join(mdir, "module_qc")
    log = setup_logging(outdir, "04b_module_qc")

    cdict = os.path.join(mdir, "coexpr", "coexpressionDictionary.json")
    z = pd.read_csv(os.path.join(res, "02_batch_corrected", f"expression_{matrix}_z.csv"), index_col=0)
    back = miner_id_backmap(p(P["miner"]["idmap"]), z.index)
    log.info("MINER renames %d of our gene IDs; mapping them back", len(back))
    modules = {k: [back.get(g, g) for g in v] for k, v in json.load(open(cdict)).items()}
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample").loc[z.columns]
    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    sstats = pd.read_csv(os.path.join(res, "01_harmonized", "sample_stats.tsv"), sep="\t", index_col="sample").loc[z.columns]
    sym = genes["symbol"].fillna("")
    missing = {k: [g for g in v if g not in z.index] for k, v in modules.items()}
    n_missing = sum(len(v) for v in missing.values())
    if n_missing:
        log.warning("%d module genes not in the expression matrix (ID conversion?); ignored", n_missing)
    modules = {k: [g for g in v if g in z.index] for k, v in modules.items()}
    modules = {k: v for k, v in modules.items() if len(v) >= 3}
    sizes = pd.Series({k: len(v) for k, v in modules.items()}).sort_values(ascending=False)
    log.info("%d modules, %d genes clustered of %d; size median %d, 90th pct %d, max %d (module %s)",
             len(sizes), len(set().union(*modules.values())), z.shape[0], sizes.median(),
             sizes.quantile(0.9), sizes.max(), sizes.index[0])

    E, M, zc = eigengenes(z, modules)
    batch = samples["cohort"]
    label_cols = [c for c in samples.columns if c not in ("cohort", "patient")]

    # marker program scores (mean z of the markers present in the matrix)
    sym2ens = sym.reset_index().query("symbol != ''").groupby("symbol")["ensembl"].first()
    programs = dict(MARKERS, **(P.get("batch", {}).get("qc", {}).get("programs", {}) or {}))
    prog = {}
    for name, syms in programs.items():
        ens = [sym2ens[s] for s in syms if s in sym2ens.index and sym2ens[s] in z.index]
        if len(ens) >= 3:
            prog[name] = z.loc[ens].mean()
    prog = pd.DataFrame(prog)
    log.info("Marker programs scored: %s", list(prog.columns))

    rows = []
    for k in M.index:
        e = E[k]
        r = {"module": k}
        r["cohort_eta2"] = eta2(e, batch)
        for lab in label_cols:
            r[f"{lab}_eta2"] = eta2(e, samples[lab])
        # technical stats: within-cohort correlation, worst cohort
        for st in sstats.columns:
            r[f"max_within_cohort_abs_r_{st}"] = max(
                abs(np.corrcoef(e[batch == c], sstats.loc[batch == c, st])[0, 1]) for c in batch.unique())
        for pr in prog.columns:
            r[f"r_{pr}"] = np.corrcoef(e, prog[pr])[0, 1]
        gsym = sym.reindex(modules[k]).fillna("")
        for cl, pat in GENE_CLASSES.items():
            r[f"frac_{cl}"] = gsym.str.match(pat).mean()
        r["top_genes"] = ",".join(gsym[gsym != ""].iloc[
            np.argsort(-np.abs(zc.loc[gsym[gsym != ""].index].values @ ((e - e.mean()) / e.std())))[:12]])
        rows.append(r)
    tab = M.join(pd.DataFrame(rows).set_index("module")).loc[sizes.index]
    rcols = [c for c in tab.columns if c.startswith("r_")]
    tab["best_marker_program"] = tab[rcols].abs().idxmax(axis=1).str[2:]
    tab["best_marker_r"] = [tab.loc[i, f"r_{b}"] for i, b in zip(tab.index, tab["best_marker_program"])]
    tab.loc[tab["best_marker_r"].abs() < 0.3, "best_marker_program"] = "none"
    tab.index.name = "module"

    # ---- regulons from mechinf, if available
    rdf_path = os.path.join(mdir, "mechinf", "regulonDf.csv")
    mech_dict = os.path.join(mdir, "mechinf", "coexpressionDictionary.json")
    if os.path.exists(rdf_path):
        if os.path.exists(mech_dict):
            mm = json.load(open(mech_dict))
            if {k: sorted(v) for k, v in mm.items()} != {k: sorted(v) for k, v in json.load(open(cdict)).items()}:  # both MINER IDs
                log.warning("mechinf re-clustered to different modules than coexpr; mapping regulons to coexpr modules")
        rdf = pd.read_csv(rdf_path)
        rdf["Gene"] = rdf["Gene"].map(lambda g: back.get(g, g))
        rdf["Regulator"] = rdf["Regulator"].map(lambda g: back.get(g, g))
        regs = rdf.groupby("Regulon_ID")["Gene"].apply(list)
        regulator = rdf.groupby("Regulon_ID")["Regulator"].first()
        gene2mod = {g: k for k, v in modules.items() for g in v}
        reg_mod = regs.apply(lambda gs: pd.Series([gene2mod.get(g) for g in gs]).mode().iloc[0]
                             if any(g in gene2mod for g in gs) else None)
        rtab = pd.DataFrame({"regulator": regulator, "module": reg_mod, "n_genes": regs.apply(len)})
        rtab.to_csv(os.path.join(outdir, "regulons_by_module.tsv"), sep="\t")
        per = rtab.groupby("module").agg(n_regulons=("regulator", "size"), n_regulators=("regulator", "nunique"),
                                         median_regulon_genes=("n_genes", "median"))
        red = {}
        for k in sizes.index[:args.top]:
            ids = rtab.index[rtab["module"] == k]
            if len(ids) < 2:
                continue
            sets = [set(regs[i]) for i in ids]
            jac = [len(a & b) / len(a | b) for i, a in enumerate(sets) for b in sets[i + 1:]]
            re_ = np.array([zc.loc[[g for g in regs[i] if g in zc.index]].mean().values for i in ids])
            c = np.corrcoef(re_)
            red[k] = {"regulon_mean_jaccard": float(np.mean(jac)),
                      "regulon_mean_eigengene_r": float(c[np.triu_indices_from(c, 1)].mean())}
        tab = tab.join(per).join(pd.DataFrame(red).T)
        log.info("Regulons: %d total from %d regulators; %d come from the %d largest modules",
                 len(rtab), rtab["regulator"].nunique(), rtab["module"].isin(sizes.index[:args.top]).sum(), args.top)
    else:
        log.info("No mechinf output yet (%s); regulon redundancy not computed", rdf_path)

    tab.to_csv(os.path.join(outdir, "modules.tsv"), sep="\t", float_format="%.4f")
    for k in sizes.index[:args.top]:
        e = E[k]
        es = (e - e.mean()) / e.std()
        g = modules[k]
        r = (zc.loc[g].values - zc.loc[g].values.mean(1, keepdims=True)) / (zc.loc[g].values.std(1, keepdims=True) + 1e-12)
        pd.DataFrame({"symbol": sym.reindex(g).values, "r_to_eigengene": r @ es / len(es)}, index=pd.Index(g, name="ensembl")) \
            .sort_values("r_to_eigengene", ascending=False).to_csv(os.path.join(outdir, f"module_{k}_genes.tsv"), sep="\t",
                                                                    float_format="%.3f")
    show = ["n_genes", "pc1_var", "pc2_var", "median_r_to_eigengene", "frac_genes_r_lt_0.3", "cohort_eta2",
            "best_marker_program", "best_marker_r"] + [c for c in ("n_regulons", "n_regulators", "regulon_mean_jaccard",
                                                                   "regulon_mean_eigengene_r") if c in tab.columns]
    log.info("Largest modules:\n%s", tab[show].head(args.top).to_string(float_format=lambda v: f"{v:.2f}"))
    log.info("All modules, medians: pc1_var %.2f, median_r_to_eigengene %.2f",
             tab["pc1_var"].median(), tab["median_r_to_eigengene"].median())
    for k in sizes.index[:args.top]:
        log.info("module %s top genes: %s", k, tab.loc[k, "top_genes"])

    import qc_plots
    written = qc_plots.module_report(outdir, tab, E, zc, modules, samples, sizes.index[:args.top], prog)
    log.info("QC figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
