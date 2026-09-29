#!/usr/bin/env python
"""Step 04c: drop regulons built on technical coexpression modules before miner3-subtypes.

Step 04b found coexpression modules that track RNA quality rather than biology (config
miner.exclude_modules; see PROJECT_LOG). A regulon is dropped when at least
miner.exclude_min_fraction of its genes lie in those modules. The full regulon set is kept
(mechinf/regulons.json); the filtered set is written next to it with MINER's own IDs, so it
can go straight to miner3-subtypes (04_run_miner.py --steps subtypes_filtered).

Outputs (results/04_miner/<matrix>/mechinf/):
  regulons_filtered.json        input to miner3-subtypes -> subtypes_filtered/
  regulon_filter.tsv            per regulon: regulator (symbol), genes, fraction in excluded modules, kept
"""

import argparse
import json
import os

import pandas as pd

from hcc_common import load_params, miner_id_backmap, p, setup_logging


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    M = P["miner"]
    matrix = args.matrix or M["matrix"]
    res = p(P["paths"]["results"])
    mdir = os.path.join(res, "04_miner", matrix)
    log = setup_logging(os.path.join(mdir, "mechinf"), "04c_filter_regulons")

    exclude = [str(m) for m in M["exclude_modules"]]
    min_frac = float(M.get("exclude_min_fraction", 0.5))
    # Module membership in MINER IDs, from the modules mechinf itself used
    mods = json.load(open(os.path.join(mdir, "mechinf", "coexpressionDictionary.json")))
    missing = [m for m in exclude if m not in mods]
    if missing:
        raise KeyError(f"exclude_modules {missing} not in mechinf coexpressionDictionary.json")
    # The exclusion list was chosen on the coexpr modules (step 04b); mechinf re-derives them.
    # Refuse to filter if the excluded modules are not the same gene sets in both.
    cmods = json.load(open(os.path.join(mdir, "coexpr", "coexpressionDictionary.json")))
    for m in exclude:
        a, b = set(mods[m]), set(cmods.get(m, []))
        jac = len(a & b) / len(a | b)
        log.info("module %s: %d genes in mechinf, %d in coexpr, Jaccard %.3f", m, len(a), len(b), jac)
        if jac < 0.9:
            raise ValueError(f"module {m} differs between coexpr and mechinf (Jaccard {jac:.2f}); "
                             "re-run step 04b on the mechinf modules and update miner.exclude_modules")
    excl_genes = set().union(*(mods[m] for m in exclude))
    log.info("Excluding regulons with >= %.0f%% of genes in modules %s (%d genes)",
             100 * min_frac, exclude, len(excl_genes))

    regulons = json.load(open(os.path.join(mdir, "mechinf", "regulons.json")))
    rdf = pd.read_csv(os.path.join(mdir, "mechinf", "regulonDf.csv"))
    regulator = rdf.groupby(rdf["Regulon_ID"].astype(str))["Regulator"].first()

    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    back = miner_id_backmap(p(M["idmap"]), genes.index)
    sym = genes["symbol"]

    rows = []
    for rid, gs in regulons.items():
        frac = sum(g in excl_genes for g in gs) / len(gs)
        reg = regulator.get(rid, "")
        rows.append({"regulon": rid, "regulator": reg, "regulator_symbol": sym.get(back.get(reg, reg), ""),
                     "n_genes": len(gs), "frac_in_excluded_modules": round(frac, 3), "kept": frac < min_frac})
    tab = pd.DataFrame(rows).set_index("regulon")
    tab.to_csv(os.path.join(mdir, "mechinf", "regulon_filter.tsv"), sep="\t")
    kept = {rid: regulons[rid] for rid in tab.index[tab["kept"]]}
    json.dump(kept, open(os.path.join(mdir, "mechinf", "regulons_filtered.json"), "w"))

    dropped = tab[~tab["kept"]]
    log.info("Regulons: %d total, %d dropped, %d kept; regulators %d -> %d",
             len(tab), len(dropped), len(kept), tab["regulator"].nunique(), tab.loc[tab["kept"], "regulator"].nunique())
    log.info("Dropped regulators (symbol, n regulons): %s",
             dropped["regulator_symbol"].replace("", "?").value_counts().head(40).to_dict())
    log.info("Fraction-in-excluded-modules distribution: %s",
             tab["frac_in_excluded_modules"].describe().round(3).to_dict())


if __name__ == "__main__":
    main()
