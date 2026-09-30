#!/usr/bin/env python
"""Step 07a: build the liver-cancer subtype signature library (config/subtype_signatures.yaml).

Downloads the MSigDB GMT files if missing, selects the library sets (name patterns + hallmark),
adds custom signatures (config/subtype_signatures_custom.tsv), resolves the NTP classifier
templates, and maps gene symbols to the project's Ensembl IDs:
  results/01_harmonized/genes.tsv symbol -> HGNC approved symbol -> HGNC previous/alias symbol
(old array signatures often use retired symbols). Ambiguous previous/alias symbols are dropped.

Outputs (results/07_post/signatures/):
  signatures.tsv      set, source, direction, symbol, ensembl, mapped_via  (one row per gene)
  sets_summary.tsv    per set: genes, mapped to our IDs, in the expression matrix
  classifiers.json    {classifier: {class: {up: [ensembl], dn: [ensembl], sets: [...]}}}
  07a_subtype_signatures.log  (lists patterns that matched nothing)
"""

import argparse
import json
import os
import re
import urllib.request

import pandas as pd
import yaml

from hcc_common import load_params, p, setup_logging

OUT = os.path.join("07_post", "signatures")


def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                sets[f[0]] = [g for g in f[2:] if g]
    return sets


def fetch(base_url, dest, log):
    if os.path.exists(dest):
        return
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    url = f"{base_url.rstrip('/')}/{os.path.basename(dest)}"
    log.info("Downloading %s", url)
    urllib.request.urlretrieve(url, dest)


def symbol_mapper(genes, hgnc_path, log):
    """symbol -> (ensembl, how). genes: our genes.tsv (index ensembl, column symbol)."""
    ours = set(genes.index)
    direct = genes["symbol"].dropna()
    direct = direct[direct != ""]
    dup = direct[direct.duplicated(keep=False)].unique()
    sym2ens = {s: e for e, s in direct.items() if s not in set(dup)}
    alias2approved = {}
    if hgnc_path and os.path.exists(hgnc_path):
        h = pd.read_csv(hgnc_path, sep="\t", dtype=str, low_memory=False)
        for r in h.itertuples(index=False):
            if isinstance(r.ensembl_gene_id, str) and r.ensembl_gene_id in ours and r.symbol not in sym2ens:
                sym2ens[r.symbol] = r.ensembl_gene_id
        pairs = []
        for col in ("prev_symbol", "alias_symbol"):
            for r in h[["symbol", col]].dropna().itertuples(index=False):
                pairs += [(a.strip(), r.symbol) for a in str(r[1]).strip('"').split("|") if a.strip()]
        a = pd.DataFrame(pairs, columns=["alias", "approved"]).drop_duplicates()
        uniq = a.groupby("alias")["approved"].nunique()
        a = a[a["alias"].isin(uniq[uniq == 1].index)]
        alias2approved = dict(zip(a["alias"], a["approved"]))
    else:
        log.warning("HGNC table not found (%s): only direct symbol matches", hgnc_path)

    def map_one(s):
        if s in sym2ens:
            return sym2ens[s], "symbol"
        ap = alias2approved.get(s)
        if ap and ap in sym2ens:
            return sym2ens[ap], "prev_or_alias"
        return None, "unmapped"
    return map_one


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    cfg = yaml.safe_load(open(p(P["post"]["signatures_config"])))
    res = p(P["paths"]["results"])
    outdir = os.path.join(res, OUT)
    log = setup_logging(outdir, "07a_subtype_signatures")

    # ---- collect sets
    M = cfg["msigdb"]
    gmt = {}
    for key, rel in M["files"].items():
        path = p(rel)
        try:
            fetch(M["base_url"], path, log)
        except Exception as e:  # noqa: BLE001 - report and continue with what is available
            log.warning("Could not download %s (%s). Download it manually to %s", key, e, path)
        if os.path.exists(path):
            s = read_gmt(path)
            log.info("%s: %d sets from %s", key, len(s), path)
            gmt[key] = s
    allsets = {}
    for key, s in gmt.items():
        for name, g in s.items():
            allsets[name] = (g, f"MSigDB {M['release']} {key}")
    cpath = p(cfg.get("custom_file"))
    if cpath and os.path.exists(cpath):
        c = pd.read_csv(cpath, sep="\t", comment="#", dtype=str).dropna(subset=["set", "gene_symbol"])
        for (name, direction), d in c.groupby(["set", "direction"]):
            key = f"custom:{name}" + ("_DN" if str(direction).lower() == "dn" else "")
            allsets[key] = (d["gene_symbol"].str.strip().tolist(), d["reference"].iloc[0])
        log.info("custom: %d sets from %s", c["set"].nunique(), cpath)

    library = set()
    for pat in cfg.get("library_patterns", []):
        hits = [n for n in allsets if re.search(pat, n)]
        if not hits:
            log.warning("library pattern %r matched no set", pat)
        library |= set(hits)
    if cfg.get("include_hallmark") and "hallmark" in gmt:
        library |= set(gmt["hallmark"])
    library |= {n for n in allsets if n.startswith("custom:")}

    # ---- classifiers
    classifiers, used = {}, set()
    for cl, classes in (cfg.get("classifiers") or {}).items():
        resolved = {}
        for cls, spec in classes.items():
            entry = {}
            for d in ("up", "dn"):
                names = sorted({n for pat in spec.get(d, []) for n in allsets if re.search(pat, n)})
                for pat in spec.get(d, []):
                    if not any(re.search(pat, n) for n in allsets):
                        log.warning("classifier %s class %s %s: pattern %r matched no set", cl, cls, d, pat)
                entry[d] = names
            if entry.get("up") or entry.get("dn"):   # a class may be defined by DN genes only
                resolved[cls] = entry
                used |= set(entry.get("up", [])) | set(entry.get("dn", []))
        if len(resolved) >= 2:
            classifiers[cl] = resolved
            log.info("classifier %s: classes %s", cl, list(resolved))
        else:
            log.warning("classifier %s skipped: only %d class(es) resolved", cl, len(resolved))
    library |= used
    if not library:
        raise RuntimeError("no signatures found: check the MSigDB download and patterns")

    # ---- map genes
    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    expressed = set(genes.index[genes["kept"].astype(bool)]) if "kept" in genes else set(genes.index)
    mapper = symbol_mapper(genes, p(P["id_mapping"].get("hgnc")), log)
    rows = []
    for name in sorted(library):
        syms, src = allsets[name]
        direction = "dn" if re.search(r"(_DN|_DOWN)$", name) else "up"
        for s in dict.fromkeys(syms):
            e, how = mapper(s)
            rows.append((name, src, direction, s, e, how))
    sig = pd.DataFrame(rows, columns=["set", "source", "direction", "symbol", "ensembl", "mapped_via"])
    sig["in_expression"] = sig["ensembl"].isin(expressed)
    sig.to_csv(os.path.join(outdir, "signatures.tsv"), sep="\t", index=False)
    summ = sig.groupby("set").agg(source=("source", "first"), direction=("direction", "first"),
                                  genes=("symbol", "size"), mapped=("ensembl", lambda x: x.notna().sum()),
                                  in_expression=("in_expression", "sum"))
    summ["frac_in_expression"] = summ["in_expression"] / summ["genes"]
    summ.to_csv(os.path.join(outdir, "sets_summary.tsv"), sep="\t", float_format="%.3f")
    log.info("Library: %d sets; median %d genes, %.0f%% of genes in the expression matrix (median over sets)",
             len(summ), summ["genes"].median(), 100 * summ["frac_in_expression"].median())
    low = summ[summ["frac_in_expression"] < 0.5]
    if len(low):
        log.warning("%d sets have < 50%% of genes in the expression matrix, e.g. %s", len(low), list(low.index[:5]))
    log.info("Mapping: %s", sig["mapped_via"].value_counts().to_dict())

    by_set = sig[sig["in_expression"]].groupby("set")["ensembl"].apply(lambda x: sorted(set(x)))
    out = {cl: {cls: {"up": sorted({g for n in e.get("up", []) for g in by_set.get(n, [])}),
                      "dn": sorted({g for n in e.get("dn", []) for g in by_set.get(n, [])}),
                      "sets": e.get("up", []) + e.get("dn", [])}
                for cls, e in classes.items()} for cl, classes in classifiers.items()}
    for cl, classes in out.items():
        log.info("classifier %s template genes: %s", cl,
                 {c: (len(v["up"]), len(v["dn"])) for c, v in classes.items()})
    json.dump(out, open(os.path.join(outdir, "classifiers.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
