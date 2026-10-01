#!/usr/bin/env python
"""Extract gene sets from Gao et al., Cell 2019;179:561 (HBV-HCC proteogenomics) supplementary tables
into config/subtype_signatures_custom.tsv (rows for sets named GAO2019_* are replaced, others kept).

Usage: python scripts/tools/extract_gao2019_signatures.py --mmc3 ... --mmc5 ... --mmc7 ...

Sets (all thresholds are the authors'; tables list only significant entries):
  GAO2019_CTNNB1MUT_PROTEIN (up / dn)   Table S7 sheet 5: proteins different in CTNNB1-mutant vs WT tumours
  GAO2019_TP53MUT_PROTEIN (up / dn)     Table S7 sheet 2: proteins different in TP53-mutant vs WT tumours
                                        (fold change mut/wt > 1 = up, < 1 = dn; BH adj. P < 0.05)
  GAO2019_ADH1A_ASSOC_MRNA (up / dn)    Table S5 sheet 4: mRNA panels associated with ADH1A protein
  GAO2019_PYCR2_ASSOC_MRNA (up / dn)    Table S5 sheet 2: mRNA panels associated with PYCR2 protein
                                        ("bottom panel" = up, "up panel" = dn: panel = heatmap position); ADH1A high = good, PYCR2 high = poor prognosis
  GAO2019_PROGNOSTIC_PROTEIN_PROTECTIVE Table S5 sheet 1: candidate prognostic proteins with HR < 1
                                        (only 3 have HR > 1: too few for a set)
  GAO2019_TUMOR_VS_NORMAL_PROTEIN (up / dn)  Table S3 sheet 1: tumour vs adjacent non-tumour proteins
Genes listed in both directions of one set are dropped. Protein-level sets transfer to mRNA only partly.
The proteomic subgroups 1/2/3 are not in Tables S3-S7 (expected in Table S1/S2).
"""

import argparse
import os

import pandas as pd

REF = "Gao 2019 Cell 179:561"


def rows(name, ref, up, dn):
    up, dn = set(map(str.strip, map(str, up))), set(map(str.strip, map(str, dn)))
    both = up & dn
    out = [(name, ref, "up", g) for g in sorted(up - both)] + [(name, ref, "dn", g) for g in sorted(dn - both)]
    return out, len(both)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mmc3", required=True)
    ap.add_argument("--mmc5", required=True)
    ap.add_argument("--mmc7", required=True)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "..", "config",
                                                  "subtype_signatures_custom.tsv"))
    a = ap.parse_args()
    new = []

    for sheet, name in (("5. Fig 6E-CTNNB1 protein list", "GAO2019_CTNNB1MUT_PROTEIN"),
                        ("2. Fig 6C-TP53 protein list", "GAO2019_TP53MUT_PROTEIN")):
        d = pd.read_excel(a.mmc7, sheet_name=sheet).dropna(subset=["Protein"])
        fc = d["Fold change (mut vs. wt)"]
        r, nb = rows(name, f"{REF} Table S7 '{sheet}'", d.loc[fc > 1, "Protein"], d.loc[fc < 1, "Protein"])
        new += r
        print(f"{name}: {sum(x[2] == 'up' for x in r)} up, {sum(x[2] == 'dn' for x in r)} dn, {nb} conflicting dropped")

    for sheet, name in (("4. ADH1A List", "GAO2019_ADH1A_ASSOC_MRNA"), ("2. PYCR2 List", "GAO2019_PYCR2_ASSOC_MRNA")):
        d = pd.read_excel(a.mmc5, sheet_name=sheet).dropna(subset=["mRNA"])
        pan = d["Panel (mRNA)"].astype(str)
        # "up panel" / "bottom panel" are positions in the paper's heatmap, not directions: the bottom panel
        # holds the mRNAs positively associated with the protein (it contains PYCR2 / ADH1A themselves and
        # its score correlates r +0.49 with the gene in our data; the up panel r -0.35 to -0.49). So
        # bottom = up, up = dn. Checked on the server 2026-10-01.
        r, nb = rows(name, f"{REF} Table S5 '{sheet}' mRNA panels", d.loc[pan.str.contains("bottom"), "mRNA"],
                     d.loc[pan.str.contains("up"), "mRNA"])
        new += r
        print(f"{name}: {sum(x[2] == 'up' for x in r)} up, {sum(x[2] == 'dn' for x in r)} dn, {nb} conflicting dropped")

    d = pd.read_excel(a.mmc5, sheet_name="1.Biomarker analysis").dropna(subset=["gene"])
    r, _ = rows("GAO2019_PROGNOSTIC_PROTEIN_PROTECTIVE", f"{REF} Table S5 '1.Biomarker analysis' HR < 1",
                d.loc[d["HR"] < 1, "gene"], [])
    new += r
    print(f"GAO2019_PROGNOSTIC_PROTEIN_PROTECTIVE: {len(r)} genes (HR > 1: {(d['HR'] > 1).sum()}, not used)")

    d = pd.read_excel(a.mmc3, sheet_name="1. 1274 DF proteins").dropna(subset=["Gene symbol"])
    r, nb = rows("GAO2019_TUMOR_VS_NORMAL_PROTEIN", f"{REF} Table S3 '1. 1274 DF proteins'",
                 d.loc[d["Log2 FC"] > 0, "Gene symbol"], d.loc[d["Log2 FC"] < 0, "Gene symbol"])
    new += r
    print(f"GAO2019_TUMOR_VS_NORMAL_PROTEIN: {sum(x[2] == 'up' for x in r)} up, {sum(x[2] == 'dn' for x in r)} dn")

    lines = open(a.out).read().splitlines()
    comments = [l for l in lines if l.startswith("#")]
    body = pd.read_csv(a.out, sep="\t", comment="#", dtype=str)
    body = body[~body["set"].fillna("").str.startswith("GAO2019_")]
    body = pd.concat([body, pd.DataFrame(new, columns=["set", "reference", "direction", "gene_symbol"])])
    with open(a.out, "w") as fh:
        fh.write("\n".join(comments) + "\n")
        body.to_csv(fh, sep="\t", index=False)
    print(f"wrote {len(new)} GAO2019 rows to {a.out} ({len(body)} rows total)")


if __name__ == "__main__":
    main()
