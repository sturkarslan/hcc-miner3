#!/usr/bin/env python
"""Step 04: MINER3 network inference on the batch-corrected matrix.

Runs, as in the working GBM pipeline (GBM-15370004/analysis_stringent_sct):
  miner3-coexpr   -> coexpr/     coexpression modules (diagnostic plots, coexpressionDictionary.json)
  miner3-mechinf  -> mechinf/    regulons.json, regulonDf.csv, eigengenes.csv, coregulationModules.json,
                                 mechanisticOutput.json (re-derives the coexpression modules itself)
  miner3-subtypes -> subtypes/   coherent/over/under/dysregulatedMembers.csv, transcriptional_programs.json,
                                 transcriptional_states.json, programs_vs_states.csv
  --steps subtypes_filtered -> subtypes_filtered/  same, on mechinf/regulons_filtered.json (step 04c)
All with --skip_tpm: the input is already z-scored (step 02). MINER keeps only genes in its
identifier_mappings.txt.

Run in the miner3 conda env. Parameters: config/params.yaml -> miner.
Outputs: results/04_miner/<matrix>/{coexpr,mechinf,subtypes}/, 04_run_miner.log, summary.tsv
"""

import argparse
import json
import os
import subprocess
import time

import pandas as pd

from hcc_common import load_params, p, setup_logging

OUT = "04_miner"


def run(cmd, log):
    log.info("RUN %s", " ".join(cmd))
    t0 = time.time()
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in res.stdout.splitlines():
        log.info("  | %s", line)
    log.info("exit %d after %.1f min", res.returncode, (time.time() - t0) / 60)
    if res.returncode != 0:
        raise RuntimeError(f"{cmd[0]} failed with exit code {res.returncode}")


def summarize(outdir, log):
    s = {}
    reg = os.path.join(outdir, "mechinf", "regulons.json")
    if os.path.exists(reg):
        regulons = json.load(open(reg))
        s["regulons"] = len(regulons)
        s["regulon_genes_median"] = int(pd.Series([len(v) for v in regulons.values()]).median())
    rdf = os.path.join(outdir, "mechinf", "regulonDf.csv")
    if os.path.exists(rdf):
        r = pd.read_csv(rdf)
        s["regulators"] = r["Regulator"].nunique() if "Regulator" in r else None
        s["regulon_genes_unique"] = r["Gene"].nunique() if "Gene" in r else None
    for sub in ("subtypes", "subtypes_filtered"):
        for key in ("transcriptional_programs", "transcriptional_states"):
            f = os.path.join(outdir, sub, f"{key}.json")
            if os.path.exists(f):
                s[f"{key.split('_')[1]}{sub[8:]}"] = len(json.load(open(f)))
    f = os.path.join(outdir, "mechinf", "regulons_filtered.json")
    if os.path.exists(f):
        s["regulons_filtered"] = len(json.load(open(f)))
    log.info("Summary: %s", s)
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None, help="override miner.matrix (e.g. cohort_z)")
    ap.add_argument("--steps", default="coexpr,mechinf,subtypes")
    args = ap.parse_args()

    P = load_params(args.params)
    M = P["miner"]
    matrix = args.matrix or M["matrix"]
    expfile = p(os.path.join(P["paths"]["results"], "02_batch_corrected", f"expression_{matrix}_z.csv"))
    outdir = p(os.path.join(P["paths"]["results"], OUT, matrix))
    default_steps = "coexpr,mechinf,subtypes"
    # separate log per partial run, so e.g. --steps subtypes_filtered does not overwrite the main log
    log = setup_logging(outdir, "04_run_miner" if args.steps == default_steps
                        else "04_run_miner_" + args.steps.replace(",", "_"))
    idmap = p(M["idmap"])
    log.info("Input %s; idmap %s; params %s", expfile, idmap, M)
    for f in (expfile, idmap):
        if not os.path.exists(f):
            raise FileNotFoundError(f)

    common = ["-mg", str(M["min_genes"]), "-moxs", str(M["min_overexp_samples"]),
              "-mx", str(M["max_exclusion"]), "-rs", str(M["random_state"]),
              "-oxt", str(M["overexp_threshold"])]
    steps = args.steps.split(",")
    if "coexpr" in steps:
        run(["miner3-coexpr", *common, expfile, idmap, os.path.join(outdir, "coexpr"), "--skip_tpm"], log)
    if "mechinf" in steps:
        cmd = ["miner3-mechinf", *common, "-mc", str(M["min_correlation"]),
               expfile, idmap, os.path.join(outdir, "mechinf"), "--skip_tpm"]
        if M.get("tfs2genes"):
            cmd[1:1] = ["--tfs2genes", p(M["tfs2genes"])]
        run(cmd, log)
    if "subtypes" in steps:
        run(["miner3-subtypes", expfile, idmap, os.path.join(outdir, "mechinf", "regulons.json"),
             os.path.join(outdir, "subtypes"), "--skip_tpm"], log)
    if "subtypes_filtered" in steps:
        # regulons without those built on technical modules (step 04c)
        run(["miner3-subtypes", expfile, idmap, os.path.join(outdir, "mechinf", "regulons_filtered.json"),
             os.path.join(outdir, "subtypes_filtered"), "--skip_tpm"], log)

    s = summarize(outdir, log)
    pd.DataFrame([dict(matrix=matrix, **s)]).to_csv(os.path.join(outdir, "summary.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
