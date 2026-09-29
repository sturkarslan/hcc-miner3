#!/usr/bin/env python
"""Step 03 (genomics part): binary genomic feature matrix for MINER3 causal inference (step 05).

One matrix, features x expression samples (all 929 network samples):
  1 = altered, 0 = profiled and not altered, empty = not profiled for that feature.
MINER's causalNetworkAnalysis() treats every non-1 column as wild type, so step 05 must run
each feature on its profiled (non-empty) samples only. Cohorts differ in what was profiled:

  feature type      TCGA                      CLCA                      LICA-FR
  gene mutation     WES MAF (cBioPortal)      WGS MAF (cBioPortal)      WES/WGS mutation table
  pathway mutation  OR of member gene mutations
  TERT promoter     not profiled (WES)        MAF, hotspots C228T/C250T driver table ("Promoter")
  arm-level CNA     cBioPortal arm calls      not available             segments (hg38), >= arm_min_fraction
  focal amp / HD    GISTIC +2 / -2            not available             driver table FA / HD

Gene and pathway features are limited to a curated HCC driver list (config genomics.driver_genes);
frequent long-gene passengers (TTN, MUC16, ...) are not tested. A feature is kept when altered in
>= min_altered samples and >= min_freq of its profiled samples.

Outputs (results/03_genomics_clinical/):
  genomic_features.csv        features x samples (1/0/empty), MINER mutation-matrix orientation
  genomic_features_info.tsv   per feature: type, per-cohort profiled / altered / frequency, kept
  qc/g1_feature_frequency.png
"""

import argparse
import os
import re

import numpy as np
import pandas as pd

from hcc_common import load_params, p, setup_logging

OUT = "03_genomics_clinical"


def read_case_list(path):
    for line in open(path):
        if line.startswith("case_list_ids:"):
            return [x for x in line.split(":", 1)[1].strip().split("\t") if x]
    raise ValueError(f"no case_list_ids in {path}")


def read_maf(path):
    return pd.read_csv(path, sep="\t", comment="#", dtype=str, low_memory=False,
                       usecols=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode",
                                "Chromosome", "Start_Position"])


def gene_calls(maf, sample_col, classes, genes):
    m = maf[maf["Variant_Classification"].isin(classes) & maf["Hugo_Symbol"].isin(genes)]
    return m.groupby("Hugo_Symbol")[sample_col].apply(set).to_dict()


def arm_table(cytoband_path):
    """Arm start/end from a UCSC cytoBand file (p arm ends at the start of the first acen band)."""
    b = pd.read_csv(cytoband_path, sep="\t", header=None, names=["chrom", "start", "end", "band", "stain"])
    b = b[b["chrom"].str.fullmatch(r"chr(\d+|X)")]
    rows = []
    for chrom, g in b.groupby("chrom"):
        cen = g.loc[g["stain"] == "acen", "start"].min()
        c = chrom[3:]
        rows += [(c, f"{c}p", 0, cen), (c, f"{c}q", cen, g["end"].max())]
    return pd.DataFrame(rows, columns=["chrom", "arm", "start", "end"])


def arm_calls_from_segments(seg, arms, min_frac, log):
    """Gain/Loss per sample and arm: fraction of the arm's segment-covered length with GNL >= 1 / <= -1."""
    out = {}
    seg = seg.copy()
    seg["chromosome"] = seg["chromosome"].astype(str)
    for a in arms.itertuples(index=False):
        s = seg[seg["chromosome"] == a.chrom]
        ov = (np.minimum(s["end"], a.end) - np.maximum(s["start"], a.start)).clip(lower=0)
        s = s.assign(ov=ov)[ov > 0]
        if s.empty:
            continue
        cov = s.groupby("CHCID")["ov"].sum()
        gain = s[s["GNL"] >= 1].groupby("CHCID")["ov"].sum().reindex(cov.index, fill_value=0) / cov
        loss = s[s["GNL"] <= -1].groupby("CHCID")["ov"].sum().reindex(cov.index, fill_value=0) / cov
        # arms with little segment coverage (acrocentric p arms) are skipped
        if cov.median() < 0.2 * (a.end - a.start):
            continue
        out[a.arm] = pd.DataFrame({"gain": gain >= min_frac, "loss": loss >= min_frac})
    log.info("LICA-FR arm calls for %d arms (min fraction %.2f)", len(out), min_frac)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    G = P["genomics"]
    outdir = p(os.path.join(P["paths"]["results"], OUT))
    log = setup_logging(outdir, "03b_genomic_features")

    samples = pd.read_csv(p(os.path.join(P["paths"]["results"], "01_harmonized", "samples.tsv")), sep="\t")
    cohort = samples.set_index("sample")["cohort"]
    S = list(samples["sample"])
    drivers = sorted(set(G["driver_genes"]))
    classes = set(G["mutation_classes"])
    feats, ftype = {}, {}   # feature -> pd.Series over S (1/0/NaN)

    def add(name, typ, altered, profiled):
        """altered, profiled: sets of expression sample IDs."""
        v = pd.Series(np.nan, index=S)
        prof = [s for s in S if s in profiled]
        v[prof] = 0.0
        v[[s for s in prof if s in altered]] = 1.0
        if name in feats:  # combine cohorts: fill this cohort's entries
            feats[name] = feats[name].where(feats[name].notna(), v)
        else:
            feats[name], ftype[name] = v, typ

    # ---------- TCGA
    T = G["TCGA"]
    tsm = samples[samples.cohort == "TCGA"].set_index("patient")["sample"]
    to_expr = lambda ids: {tsm[i[:12]] for i in ids if i[:12] in tsm.index}
    tseq = to_expr(read_case_list(p(T["sequenced_case_list"])))
    maf = read_maf(p(T["maf"]))
    calls = gene_calls(maf, "Tumor_Sample_Barcode", classes, drivers)
    for g in drivers:
        add(f"MUT_{g}", "gene_mutation", to_expr(calls.get(g, set())), tseq)
    log.info("TCGA: %d sequenced network samples; %d driver genes with calls", len(tseq), len(calls))
    arm = pd.read_csv(p(T["arm_cna"]), sep="\t", index_col=0).drop(columns=["NAME", "DESCRIPTION"])
    for a_id, row in arm.iterrows():
        a = a_id.replace("_status", "")
        if a in [str(x) for x in G.get("exclude_arms", [])]:
            continue
        prof = to_expr(row.index[row.notna()])
        add(f"ARM_{a}_gain", "arm_cna", to_expr(row.index[row == "Gain"]), prof)
        add(f"ARM_{a}_loss", "arm_cna", to_expr(row.index[row == "Loss"]), prof)
    gistic = pd.read_csv(p(T["gistic"]), sep="\t").drop(columns=["Entrez_Gene_Id"]).groupby("Hugo_Symbol").first()
    tcna = to_expr(read_case_list(p(T["cna_case_list"])))
    for name, spec in G["focal"].items():
        genes = [g for g in spec["genes"] if g in gistic.index]
        val = 2 if spec["type"] == "amp" else -2
        hit = (gistic.loc[genes] == val).any(axis=0) if genes else pd.Series(False, index=gistic.columns)
        add(f"{spec['type'].upper()}_{name}", "focal_cna", to_expr(hit.index[hit]), tcna)
    log.info("TCGA: %d CNA-profiled network samples", len(tcna))

    # ---------- CLCA
    C = G["CLCA"]
    cexpr = set(samples.loc[samples.cohort == "CLCA", "sample"])
    cseq = set(read_case_list(p(C["sequenced_case_list"]))) & cexpr
    maf = read_maf(p(C["maf"]))
    calls = gene_calls(maf, "Tumor_Sample_Barcode", classes, drivers)
    for g in drivers:
        add(f"MUT_{g}", "gene_mutation", calls.get(g, set()) & cseq, cseq)
    hs = {int(x) for x in C["tert_promoter_hotspots_hg19"]}
    tert = maf[(maf["Hugo_Symbol"] == "TERT") & (maf["Variant_Classification"] == "5'Flank")
               & maf["Start_Position"].astype(int).isin(hs)]
    add("TERT_promoter", "tert_promoter", set(tert["Tumor_Sample_Barcode"]) & cseq, cseq)
    log.info("CLCA: %d sequenced network samples; TERT promoter hotspot in %d",
             len(cseq), len(set(tert["Tumor_Sample_Barcode"]) & cseq))

    # ---------- LICA-FR
    L = G["LICA_FR"]
    import pyreadr
    lexpr = set(samples.loc[samples.cohort == "LICA_FR", "sample"])
    mt = pyreadr.read_r(p(L["mutations"]))[L["mutations_object"]]
    lseq = set(mt["Sample"]) & lexpr
    # The table keeps low-VAF calls (median coding VAF 0.04); VAF >= min_vaf reproduces the authors'
    # curated driver table (see PROJECT_LOG). TCGA/CLCA cBioPortal MAFs are already filtered calls.
    mt = mt[mt["Mutation_Class"].isin(L["mutation_class_keep"]) & mt["Hugo_Symbol"].isin(drivers)
            & (mt["TUMOR_AF.brc"] >= L["min_vaf"])]
    calls = mt.groupby("Hugo_Symbol")["Sample"].apply(set).to_dict()
    for g in drivers:
        add(f"MUT_{g}", "gene_mutation", calls.get(g, set()) & lseq, lseq)
    drv = pd.read_excel(p(L["driver_table"])).set_index("Sample")
    ldrv = set(drv.index) & lexpr
    add("TERT_promoter", "tert_promoter", set(drv.index[drv["TERT_Mut"] == "Promoter"]) & ldrv, ldrv)
    for name, spec in G["focal"].items():
        col, vals = spec.get("licafr_column"), spec.get("licafr_values")
        if col:
            add(f"{spec['type'].upper()}_{name}", "focal_cna", set(drv.index[drv[col].isin(vals)]) & ldrv, ldrv)
    seg = pyreadr.read_r(p(L["segments"]))[L["segments_object"]]
    arms = arm_table(p(L["cytoband"]))
    arms = arms[~arms["arm"].isin([str(a) for a in G.get("exclude_arms", [])])]
    ac = arm_calls_from_segments(seg, arms, G["arm_min_fraction"], log)
    lseg = set(seg["CHCID"]) & lexpr
    for a, df in ac.items():
        df = df.loc[df.index.isin(lexpr)]
        add(f"ARM_{a}_gain", "arm_cna", set(df.index[df["gain"]]), lseg)
        add(f"ARM_{a}_loss", "arm_cna", set(df.index[df["loss"]]), lseg)
    log.info("LICA-FR: %d mutation-profiled, %d driver-table, %d segment network samples",
             len(lseq), len(ldrv), len(lseg))

    # ---------- pathways (OR of member gene mutations; profiled where all members are)
    for pw, genes in G["pathways"].items():
        cols = [f"MUT_{g}" for g in genes if f"MUT_{g}" in feats]
        m = pd.concat([feats[c] for c in cols], axis=1)
        v = (m == 1).any(axis=1).astype(float).where(m.notna().all(axis=1))
        feats[f"PATH_{pw}"], ftype[f"PATH_{pw}"] = v, "pathway_mutation"

    X = pd.DataFrame(feats).T[S]
    info = pd.DataFrame({"type": pd.Series(ftype)})
    for c in ["TCGA", "CLCA", "LICA_FR"]:
        sub = X.loc[:, cohort[X.columns] == c]
        info[f"profiled_{c}"] = sub.notna().sum(axis=1)
        info[f"altered_{c}"] = (sub == 1).sum(axis=1)
        info[f"freq_{c}"] = (info[f"altered_{c}"] / info[f"profiled_{c}"]).round(3)
    info["profiled"] = X.notna().sum(axis=1)
    info["altered"] = (X == 1).sum(axis=1)
    info["freq"] = (info["altered"] / info["profiled"]).round(3)
    info["kept"] = (info["altered"] >= G["min_altered"]) & (info["freq"] >= G["min_freq"])
    info.index.name = "feature"
    info.sort_values(["type", "freq"], ascending=[True, False]).to_csv(
        os.path.join(outdir, "genomic_features_info.tsv"), sep="\t")
    Xk = X.loc[info.index[info["kept"]]]
    Xk.index.name = "feature"
    Xk.to_csv(os.path.join(outdir, "genomic_features.csv"), float_format="%.0f")
    log.info("Features: %d built, %d kept (>= %d altered and >= %.0f%% of profiled). By type: %s",
             len(info), int(info["kept"].sum()), G["min_altered"], 100 * G["min_freq"],
             info.loc[info["kept"], "type"].value_counts().to_dict())
    show = info[info["kept"]].sort_values("freq", ascending=False)
    log.info("Kept features:\n%s", show[["type", "freq_TCGA", "freq_CLCA", "freq_LICA_FR", "profiled", "altered"]]
             .to_string())

    import qc_plots
    written = qc_plots.genomic_report(outdir, info[info["kept"]], ["TCGA", "CLCA", "LICA_FR"])
    log.info("QC figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
