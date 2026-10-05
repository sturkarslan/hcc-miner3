#!/usr/bin/env python
"""Add PROGENy pathway-response footprints (Schubert et al., Nat Commun 2018) to config/subtype_signatures_custom.tsv
as signed sets custom:PROGENY_<PATHWAY> (up = positive weight, _DN = negative weight), top 100 genes per pathway by
p value (the standard PROGENy model). Source: OmniPath annotations (resources=PROGENy), saved to
data/reference/progeny_omnipath.tsv. Used by the step-10 drug -> network map (MAPK, VEGF, hypoxia)."""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(ROOT, "data", "reference", "progeny_omnipath.tsv")
OUT = os.path.join(ROOT, "config", "subtype_signatures_custom.tsv")
PATHWAYS = ["MAPK", "VEGF", "Hypoxia", "EGFR", "TNFa", "NFkB", "JAK-STAT", "TGFb", "PI3K", "p53", "WNT", "Trail"]


def main():
    d = pd.read_csv(SRC, sep="\t")
    w = d.pivot_table(index=["record_id", "genesymbol"], columns="label", values="value", aggfunc="first").reset_index()
    w["weight"] = pd.to_numeric(w["weight"])
    w["p_value"] = pd.to_numeric(w["p_value"])
    lines = [l for l in open(OUT) if not l.startswith("PROGENY_")]
    new = []
    for pw in PATHWAYS:
        x = w[w["pathway"] == pw].sort_values("p_value").drop_duplicates("genesymbol").head(100)
        name = "PROGENY_" + pw.upper().replace("-", "_")
        for _, r in x.iterrows():
            new.append(f"{name}\tPROGENy (Schubert 2018) top 100 by p, via OmniPath\t{'up' if r['weight'] > 0 else 'dn'}\t{r['genesymbol']}\n")
        print(name, len(x), "genes:", int((x["weight"] > 0).sum()), "up /", int((x["weight"] < 0).sum()), "dn")
    open(OUT, "w").writelines(lines + new)


if __name__ == "__main__":
    main()
