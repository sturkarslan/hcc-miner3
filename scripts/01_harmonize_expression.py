#!/usr/bin/env python
"""Step 01: harmonize expression across TCGA-LIHC, CLCA 2024 and LICA-FR.

TPM for every cohort -> Ensembl gene IDs -> primary HCC tumors only ->
genes shared by all cohorts -> optional TPM renormalization -> gene filtering ->
log2(TPM + 1).

Outputs (results/01_harmonized/):
  expression_log2tpm1.csv      genes x samples, all cohorts (input to step 02)
  samples.tsv                  sample, cohort, patient, LICA-FR labels
  genes.tsv                    ensembl, symbol, per-cohort detection, kept flag
  id_mapping_<cohort>.tsv      input identifier -> Ensembl, with mapping tier or failure reason
  summary.tsv                  counts at each step
"""

import argparse
import os
import re

import numpy as np
import pandas as pd

from hcc_common import load_params, p, require_columns, setup_logging

OUT = "01_harmonized"


# ------------------------------------------------------------------ loading

def load_xena_log2tpm(cfg, log):
    df = pd.read_csv(p(cfg["expression"]), sep="\t", index_col=0)
    df = df[~df.index.str.endswith("_PAR_Y")]
    df.index = df.index.str.replace(r"\.\d+$", "", regex=True)
    tpm = np.power(2.0, df) - cfg.get("pseudocount", 0.001)
    tpm = tpm.clip(lower=0)
    log.info("TCGA: %d genes x %d samples; back-transformed log2(TPM+%g); column-sum median %.3g",
             *tpm.shape, cfg.get("pseudocount", 0.001), tpm.sum().median())
    return tpm


def load_cbioportal_tpm(cfg, log):
    df = pd.read_csv(p(cfg["expression"]), sep="\t", comment="#")
    drop = [c for c in ("Entrez_Gene_Id",) if c in df.columns]
    df = df.drop(columns=drop).dropna(subset=["Hugo_Symbol"]).set_index("Hugo_Symbol")
    df = df.apply(pd.to_numeric, errors="coerce")
    log.info("%s: %d rows x %d samples; column-sum median %.3g",
             cfg["_name"], *df.shape, df.sum().median())
    return df


def load_rdata_tpm(cfg, log):
    import pyreadr
    res = pyreadr.read_r(p(cfg["expression"]))
    obj = cfg.get("rdata_object", "tpm_")
    if obj not in res:
        raise KeyError(f"RData object {obj!r} not found; objects: {list(res.keys())}")
    df = res[obj]
    if isinstance(df.index, pd.RangeIndex) or df.index.astype(str).str.fullmatch(r"\d+").all():
        raise ValueError(
            f"{cfg['expression']}: row names were not preserved by pyreadr. Export from R and use "
            "format: tsv_tpm, e.g.\n  Rscript -e 'load(\"RNAseq_tpm_457s.RData\"); "
            "write.table(tpm_, \"RNAseq_tpm_457s.tsv\", sep=\"\\t\", quote=FALSE, col.names=NA)'")
    log.info("%s: %d rows x %d samples; column-sum median %.3g",
             cfg["_name"], *df.shape, df.sum().median())
    return df


def load_tsv_tpm(cfg, log):
    df = pd.read_csv(p(cfg["expression"]), sep="\t", index_col=0)
    log.info("%s: %d rows x %d samples", cfg["_name"], *df.shape)
    return df


LOADERS = {"xena_log2tpm": load_xena_log2tpm, "cbioportal_tpm": load_cbioportal_tpm,
           "rdata_tpm": load_rdata_tpm, "tsv_tpm": load_tsv_tpm}


# ------------------------------------------------------------------ ID mapping

def _split_multi(s):
    if pd.isna(s) or s == "":
        return []
    return [x.strip() for x in str(s).strip('"').split("|") if x.strip()]


def build_symbol_maps(idcfg, log):
    """Ordered list of (tier, {symbol: set(ensembl)}) plus an Ensembl -> symbol table."""
    tiers = {t: {} for t in ("hgnc_symbol", "miner_gene_name", "hgnc_prev", "alias")}
    ens2sym = {}
    protein_coding = None

    def add(tier, sym, ens):
        if sym and isinstance(ens, str) and ens.startswith("ENSG"):
            tiers[tier].setdefault(sym, set()).add(ens)

    hgnc_path = p(idcfg.get("hgnc"))
    if hgnc_path and os.path.exists(hgnc_path):
        h = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
        require_columns(h, ["symbol", "ensembl_gene_id", "alias_symbol", "prev_symbol"], "HGNC table")
        h = h.dropna(subset=["ensembl_gene_id"])
        for r in h.itertuples(index=False):
            add("hgnc_symbol", r.symbol, r.ensembl_gene_id)
            ens2sym.setdefault(r.ensembl_gene_id, r.symbol)
            for s in _split_multi(r.prev_symbol):
                add("hgnc_prev", s, r.ensembl_gene_id)
            for s in _split_multi(r.alias_symbol):
                add("alias", s, r.ensembl_gene_id)
        if "locus_group" in h.columns:
            protein_coding = set(h.loc[h["locus_group"] == "protein-coding gene", "ensembl_gene_id"])
        log.info("HGNC: %d approved symbols with Ensembl IDs", len(tiers["hgnc_symbol"]))
    else:
        log.warning("HGNC table not found (%s); using MINER identifier_mappings only", hgnc_path)

    m = pd.read_csv(p(idcfg["miner_idmap"]), sep="\t", dtype=str)
    require_columns(m, ["Preferred_Name", "Name", "Source"], "MINER identifier_mappings")
    for r in m[m.Source == "Gene Name"].itertuples(index=False):
        add("miner_gene_name", r.Name, r.Preferred_Name)
        ens2sym.setdefault(r.Preferred_Name, r.Name)
    for r in m[m.Source == "Synonym"].itertuples(index=False):
        add("alias", r.Name, r.Preferred_Name)
    log.info("MINER idmap: %d gene names", len(tiers["miner_gene_name"]))
    return list(tiers.items()), ens2sym, protein_coding


def map_symbols(ids, tiers):
    """Return DataFrame(input_id, ensembl, status). First tier with any hit wins;
    a symbol with >1 Ensembl ID in that tier is ambiguous and dropped."""
    rows = []
    for s in ids:
        s_str = str(s)
        if re.match(r"^ENSG\d+", s_str):
            rows.append((s, s_str.split(".")[0], "ensembl"))
            continue
        for tier, mp in tiers:
            hit = mp.get(s_str)
            if hit:
                rows.append((s, next(iter(hit)), tier) if len(hit) == 1
                            else (s, None, f"ambiguous:{tier}:{'|'.join(sorted(hit))}"))
                break
        else:
            rows.append((s, None, "unmapped"))
    return pd.DataFrame(rows, columns=["input_id", "ensembl", "status"])


def to_ensembl(df, cfg, tiers, outdir, log):
    name = cfg["_name"]
    df = df.groupby(level=0).sum()  # duplicated input identifiers: TPM is additive
    if cfg["id_type"] == "ensembl":
        mp = pd.DataFrame({"input_id": df.index, "ensembl": df.index, "status": "ensembl"})
    else:
        mp = map_symbols(df.index, tiers)
    mp.to_csv(os.path.join(outdir, f"id_mapping_{name}.tsv"), sep="\t", index=False)
    counts = mp.status.str.split(":").str[0].value_counts()
    log.info("%s ID mapping: %s", name, counts.to_dict())
    if mp.ensembl.notna().mean() < 0.8:
        log.warning("%s: only %.0f%% of identifiers mapped to Ensembl; check id_type and the ID column",
                    name, 100 * mp.ensembl.notna().mean())
    lost_tpm = df.loc[mp.loc[mp.ensembl.isna(), "input_id"]].sum().median() if mp.ensembl.isna().any() else 0.0
    log.info("%s: median TPM per sample in unmapped/ambiguous rows = %.1f", name, lost_tpm)
    ok = mp.dropna(subset=["ensembl"])
    out = df.loc[ok.input_id].set_axis(ok.ensembl.values).groupby(level=0).sum()
    return out, counts


# ------------------------------------------------------------------ sample selection

def select_tcga(df, cfg, log):
    keep = [c for c in df.columns if re.match(cfg["sample_regex"], c)]
    df = df[keep]
    excl_path = p(cfg.get("exclude_list"))
    if excl_path and os.path.exists(excl_path):
        excl = {l.split("\t")[0].strip() for l in open(excl_path) if l.strip() and not l.startswith("#")}
        drop = [c for c in df.columns if c in excl or c[:12] in excl]
        (log.info if drop else log.warning)("TCGA: excluding %d samples from %s", len(drop), excl_path)
        df = df.drop(columns=drop)
    else:
        log.warning("TCGA: exclude list %s not found; fibrolamellar cases are NOT removed", excl_path)
    patient = pd.Series([c[:12] for c in df.columns], index=df.columns)
    dup = patient[patient.duplicated(keep="first")]
    if len(dup):
        log.info("TCGA: dropping %d duplicate primary aliquots (kept first per patient)", len(dup))
        df = df.drop(columns=dup.index)
    samples = pd.DataFrame({"sample": df.columns, "patient": [c[:12] for c in df.columns]})
    return df, samples


def select_licafr(df, cfg, log):
    ann = pd.read_excel(p(cfg["annotation"]))
    idc, typc = cfg["annotation_id_column"], cfg["sample_type_column"]
    labels = cfg.get("label_columns", {}) or {}
    require_columns(ann, [idc, typc] + list(labels.values()), "LICA-FR annotation")
    ann = ann.copy()
    ann["sample"] = ann[idc].astype(str)
    for a, b in (cfg.get("annotation_id_prefix_map") or {}).items():
        ann["sample"] = ann["sample"].str.replace(f"^{re.escape(a)}", b, regex=True)
    log.info("LICA-FR annotation %s counts: %s", typc, ann[typc].value_counts(dropna=False).to_dict())
    keep = ann[ann[typc].isin(cfg["keep_values"])]
    in_expr = keep[keep["sample"].isin(df.columns)]
    log.info("LICA-FR: %d annotated primary HCC, %d with expression; %d expression columns not annotated",
             len(keep), len(in_expr), len(set(df.columns) - set(ann["sample"])))
    df = df[in_expr["sample"].tolist()]
    samples = in_expr[["sample"]].copy()
    samples["patient"] = samples["sample"]
    for new, old in labels.items():
        samples[new] = in_expr[old].values
    return df, samples.reset_index(drop=True)


def select_default(df, cfg, log):
    return df, pd.DataFrame({"sample": df.columns, "patient": df.columns})


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    args = ap.parse_args()

    P = load_params(args.params)
    outdir = p(os.path.join(P["paths"]["results"], OUT))
    log = setup_logging(outdir, "01_harmonize_expression")
    tiers, ens2sym, protein_coding = build_symbol_maps(P["id_mapping"], log)
    H = P["harmonize"]

    mats, sample_tables, summary = {}, [], []
    for name, cfg in P["cohorts"].items():
        cfg = dict(cfg, _name=name)
        raw = LOADERS[cfg["format"]](cfg, log)
        if raw.isna().any().any():
            log.warning("%s: %d NaN values set to 0", name, int(raw.isna().sum().sum()))
            raw = raw.fillna(0)
        ens, _ = to_ensembl(raw, cfg, tiers, outdir, log)
        selector = {"xena_log2tpm": select_tcga}.get(cfg["format"],
                                                    select_licafr if "annotation" in cfg else select_default)
        ens, st = selector(ens, cfg, log)
        if st["sample"].duplicated().any():
            raise ValueError(f"{name}: duplicated sample IDs after selection")
        st.insert(1, "cohort", name)
        mats[name] = ens
        sample_tables.append(st)
        summary.append({"cohort": name, "input_rows": raw.shape[0], "input_samples": raw.shape[1],
                        "ensembl_genes": ens.shape[0], "samples_kept": ens.shape[1]})

    samples = pd.concat(sample_tables, ignore_index=True)
    if samples["sample"].duplicated().any():
        raise ValueError(f"sample IDs collide across cohorts: {samples.loc[samples['sample'].duplicated(), 'sample'].tolist()[:10]}")

    shared = sorted(set.intersection(*(set(m.index) for m in mats.values())))
    log.info("Genes shared by all cohorts: %d", len(shared))
    if len(shared) < 1000:
        raise ValueError(f"only {len(shared)} genes shared by all cohorts; per-cohort Ensembl gene counts: "
                         f"{ {n: m.shape[0] for n, m in mats.items()} }")
    if H.get("protein_coding_only"):
        if protein_coding is None:
            raise ValueError("protein_coding_only needs the HGNC table with locus_group")
        shared = [g for g in shared if g in protein_coding]
        log.info("Protein-coding shared genes: %d", len(shared))

    det, sstats = {}, []
    for name in mats:
        m = mats[name].loc[shared]
        # share of each sample's mapped TPM that falls in the shared genes (1 - renormalization loss)
        frac_shared = m.sum() / mats[name].sum()
        if H.get("renormalize_tpm", True):
            m = m / m.sum() * 1e6
        mats[name] = m
        det[name] = (m >= H["min_tpm"]).mean(axis=1)
        sstats.append(pd.DataFrame({"tpm_frac_shared_genes": frac_shared,
                                    "n_genes_detected": (m >= H["min_tpm"]).sum()}))
    det = pd.DataFrame(det)
    kept = (det >= H["min_frac"]).all(axis=1)
    log.info("Genes passing TPM >= %g in >= %.0f%% of samples in every cohort: %d / %d",
             H["min_tpm"], 100 * H["min_frac"], kept.sum(), len(kept))

    genes = det.add_prefix("frac_detected_")
    genes.insert(0, "symbol", [ens2sym.get(g, "") for g in genes.index])
    genes["kept"] = kept
    genes.index.name = "ensembl"
    genes.to_csv(os.path.join(outdir, "genes.tsv"), sep="\t")

    expr = pd.concat([np.log2(mats[n].loc[kept[kept].index] + 1) for n in mats], axis=1)
    expr = expr[samples["sample"]]
    expr.index.name = "ensembl"
    expr.to_csv(os.path.join(outdir, "expression_log2tpm1.csv"))
    samples.to_csv(os.path.join(outdir, "samples.tsv"), sep="\t", index=False)
    sstats = pd.concat(sstats).loc[samples["sample"]]
    sstats.index.name = "sample"
    sstats.to_csv(os.path.join(outdir, "sample_stats.tsv"), sep="\t", float_format="%.4f")

    summary = pd.DataFrame(summary)
    summary["shared_genes"] = len(shared)
    summary["genes_after_filter"] = int(kept.sum())
    summary.to_csv(os.path.join(outdir, "summary.tsv"), sep="\t", index=False)
    log.info("Wrote %d genes x %d samples to %s\n%s", *expr.shape, outdir, summary.to_string(index=False))

    import qc_plots
    written = qc_plots.harmonization_report(outdir, P["harmonize"], list(P["cohorts"]))
    log.info("QC figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
