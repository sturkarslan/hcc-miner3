#!/usr/bin/env python
"""Extract gene sets from Zhu et al., Nat Med 2022;28:1599 (atezolizumab + bevacizumab in HCC; GO30140 /
IMbrave150) Supplementary Table 2 (41591_2022_1868_MOESM3_ESM.xlsx) into config/subtype_signatures_custom.tsv
(rows for sets named ZHU2022_* are replaced, others kept).

Usage: python scripts/tools/extract_zhu2022_signatures.py --moesm3 41591_2022_1868_MOESM3_ESM.xlsx

Sets:
  ZHU2022_ATEZOBEV_RESPONDER_UP  genes higher in responders (CR/PR) than non-responders (SD/PD), overlapping
                                 between GO30140 group A and the IMbrave150 atezo-bev arm (logFC > 0; the table
                                 has only 2 genes with logFC < 0, too few for a set)
  ZHU2022_ABRS                   the 10-gene atezolizumab-bevacizumab response signature (ABRS): the genes the
                                 authors highlight in red font in the same table (read from the cell font colour)
The Teff, Treg, myeloid-inflammation and angiogenesis signatures of the paper are defined in its Methods, not in
the supplementary files; add them from the main text.
"""

import argparse
import os

import openpyxl
import pandas as pd

REF = "Zhu 2022 Nat Med 28:1599 Supplementary Table 2"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--moesm3", required=True)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "..", "config",
                                                  "subtype_signatures_custom.tsv"))
    a = ap.parse_args()
    sheet = "SupplementaryTable 2"
    d = pd.read_excel(a.moesm3, sheet_name=sheet, header=1)
    d = d[d["GeneID"].notna()]                       # drop footnote rows
    up = d.loc[d["logFC"] > 0, "Gene Symbol"].astype(str).str.strip()
    print(f"responder genes: {len(up)} up, {(d['logFC'] < 0).sum()} down (down not used)")

    ws = openpyxl.load_workbook(a.moesm3)[sheet]
    genes = set(d["Gene Symbol"].astype(str).str.strip())
    red = []
    for (c,) in ws.iter_rows(min_row=3, max_col=1):
        col = c.font.color.rgb if c.font and c.font.color is not None and c.font.color.type == "rgb" else ""
        if isinstance(col, str) and col.upper().endswith("FF0000") and str(c.value).strip() in genes:
            red.append(str(c.value).strip())
    print(f"ABRS (red font): {len(red)} genes: {red}")
    if len(red) != 10:
        raise SystemExit("expected 10 red-highlighted ABRS genes; check the table formatting")

    new = [("ZHU2022_ATEZOBEV_RESPONDER_UP", f"{REF} (logFC > 0)", "up", g) for g in sorted(set(up))]
    new += [("ZHU2022_ABRS", f"{REF} (red font: ABRS 10-gene signature)", "up", g) for g in red]
    comments = [l for l in open(a.out).read().splitlines() if l.startswith("#")]
    body = pd.read_csv(a.out, sep="\t", comment="#", dtype=str)
    body = body[~body["set"].fillna("").str.startswith("ZHU2022_")]
    body = pd.concat([body, pd.DataFrame(new, columns=["set", "reference", "direction", "gene_symbol"])])
    with open(a.out, "w") as fh:
        fh.write("\n".join(comments) + "\n")
        body.to_csv(fh, sep="\t", index=False)
    print(f"wrote {len(new)} ZHU2022 rows to {a.out} ({len(body)} rows total)")


if __name__ == "__main__":
    main()
