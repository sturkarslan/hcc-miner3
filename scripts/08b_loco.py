#!/usr/bin/env python
"""Step 08b: leave-one-cohort-out (LOCO) stability of the network, programs, causal flows and risk.

For each held-out discovery cohort H (TCGA, CLCA, LICA_FR), matrix name loco_no<H>:
  prepare   (hcc-prep env) ComBat + gene z-score on the two remaining cohorts only (same genes as step 01,
            same settings as step 02) -> results/02_batch_corrected/expression_loco_no<H>_z.csv.
            H contributes nothing to batch correction, network or causal inference.
  [04]      scripts/04_run_miner.py --matrix loco_no<H> --steps mechinf     (miner3 env)
  filter    (hcc-prep env) drop regulons with >= miner.exclude_min_fraction of genes in the discovery
            technical modules (miner.exclude_modules of the full run). Gene-based, because module IDs differ
            between runs -> mechinf/regulons_filtered.json
  [04]      --steps subtypes_filtered;  [05] scripts/05_causal_inference.py --matrix loco_no<H>
  compare   (miner3 env: needs miner + lifelines) vs the full network (miner.matrix):
            regulons: best Jaccard of each full regulon to any LOCO regulon, and with the same regulator;
            regulators recovered; programs: best gene-set Jaccard, and activity correlation of each full
            program with its best LOCO match in the held-out cohort (its own within-cohort z);
            causal: full high-confidence driver -> regulator edges (with direction) found among LOCO flows;
            risk (H with survival): MINER ridge on LOCO program activity trained in the other survival
            cohort, scored in H (never seen by network or model): C-index, HR per s.d.
Outputs: results/08_validation/loco/{loco_summary.tsv, per-cohort tables, figures/l1_loco.png}
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

from hcc_common import load_params, miner_id_backmap, p, setup_logging

COHORTS = ["TCGA", "CLCA", "LICA_FR"]


def zrows(x):
    return x.sub(x.mean(1), axis=0).div(x.std(1).replace(0, np.nan), axis=0)


# ---------------------------------------------------------------- prepare
def prepare(P, held, log):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import importlib.util
    spec = importlib.util.spec_from_file_location("s02", os.path.join(os.path.dirname(__file__), "02_batch_correct.py"))
    s02 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(s02)
    res = p(P["paths"]["results"])
    expr = pd.read_csv(os.path.join(res, "01_harmonized", "expression_log2tpm1.csv"), index_col=0)
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    keep = samples.index[samples["cohort"] != held]
    expr = expr[keep]
    batch = samples.loc[keep, "cohort"]
    zero = pd.concat([expr.loc[:, batch == b].var(axis=1) == 0 for b in batch.unique()], axis=1).any(axis=1)
    expr = expr.loc[~zero]
    B = P["batch"]
    cor = s02.combat(expr, batch, B.get("combat_parametric", True), log)
    z = s02.zscore_rows(cor)
    if B.get("z_clip_low") is not None:
        z = z.clip(lower=B["z_clip_low"])
    out = os.path.join(res, "02_batch_corrected", f"expression_loco_no{held}_z.csv")
    z.to_csv(out)
    log.info("loco_no%s: %d genes x %d samples (%s) -> %s", held, *z.shape, batch.value_counts().to_dict(), out)


# ---------------------------------------------------------------- filter
def filt(P, held, log):
    res = p(P["paths"]["results"])
    M = P["miner"]
    full = json.load(open(os.path.join(res, "04_miner", M["matrix"], "mechinf", "coexpressionDictionary.json")))
    tech = set().union(*(full[str(m)] for m in M["exclude_modules"]))
    mdir = os.path.join(res, "04_miner", f"loco_no{held}", "mechinf")
    reg = json.load(open(os.path.join(mdir, "regulons.json")))
    frac = {k: sum(g in tech for g in v) / len(v) for k, v in reg.items()}
    kept = {k: v for k, v in reg.items() if frac[k] < M["exclude_min_fraction"]}
    json.dump(kept, open(os.path.join(mdir, "regulons_filtered.json"), "w"))
    log.info("loco_no%s: %d regulons, %d dropped as technical (>= %.0f%% genes in full-run modules %s), %d kept",
             held, len(reg), len(reg) - len(kept), 100 * M["exclude_min_fraction"], M["exclude_modules"], len(kept))


# ---------------------------------------------------------------- compare
def load_net(res, matrix, back):
    mdir = os.path.join(res, "04_miner", matrix)
    reg = {k: {back.get(g, g) for g in v} for k, v in json.load(open(os.path.join(mdir, "mechinf", "regulons_filtered.json"))).items()}
    rdf = pd.read_csv(os.path.join(mdir, "mechinf", "regulonDf.csv"), index_col=0)
    rdf["Regulon_ID"] = rdf["Regulon_ID"].astype(str)
    regulator = rdf.groupby("Regulon_ID")["Regulator"].first().map(lambda g: back.get(g, g))
    progs = json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_programs.json")))
    pgenes = {k: set().union(*(reg[str(r)] for r in v if str(r) in reg)) for k, v in progs.items()}
    return reg, regulator, progs, pgenes


def best_jaccard(A, B):
    """For each set in A, the best Jaccard with any set in B (inverted index)."""
    idx = {}
    for k, s in B.items():
        for g in s:
            idx.setdefault(g, []).append(k)
    out = {}
    for k, s in A.items():
        cand = {}
        for g in s:
            for kb in idx.get(g, ()):
                cand[kb] = cand.get(kb, 0) + 1
        best, arg = 0.0, None
        for kb, inter in cand.items():
            j = inter / (len(s) + len(B[kb]) - inter)
            if j > best:
                best, arg = j, kb
        out[k] = (best, arg)
    return out


def compare(P, log):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    from miner import miner
    from sklearn.linear_model import Ridge
    res = p(P["paths"]["results"])
    M, R = P["miner"], P["risk"]
    outdir = os.path.join(res, "08_validation", "loco")
    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    back = miner_id_backmap(p(M["idmap"]), genes.index)
    sym = genes["symbol"]
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    expr = pd.read_csv(os.path.join(res, "01_harmonized", "expression_log2tpm1.csv"), index_col=0)
    reg0, regr0, progs0, pg0 = load_net(res, M["matrix"], back)
    hc0 = pd.read_csv(os.path.join(res, "05_causal", M["matrix"], "highConfidenceCausalResults.csv"), index_col=0)
    e0 = set(zip(hc0["Mutation"], hc0["regulator_symbol"], np.sign(hc0["MutationRegulatorEdge"])))
    w0 = pd.read_csv(os.path.join(res, "06_risk", M["matrix"], "predictor_ridge_programs_TCGA_RFS_h36m", "weights.tsv"),
                     sep="\t", index_col=0)["weight"].rename(index=str)
    rows, prog_rows = [], []
    for held in COHORTS:
        mx = f"loco_no{held}"
        if not os.path.exists(os.path.join(res, "04_miner", mx, "subtypes_filtered", "transcriptional_programs.json")):
            log.warning("%s: network not finished; skipped", mx)
            continue
        reg, regr, progs, pg = load_net(res, mx, back)
        r = {"held_out": held, "regulons_full": len(reg0), "regulons_loco": len(reg), "programs_loco": len(progs)}
        bj = best_jaccard(reg0, reg)
        jv = np.array([v[0] for v in bj.values()])
        same = np.array([bj[k][1] is not None and regr.get(bj[k][1]) == regr0.get(k) for k in reg0])
        r.update(regulon_median_best_jaccard=float(np.median(jv)), regulon_frac_jaccard_ge_0_5=float((jv >= 0.5).mean()),
                 regulon_frac_jaccard_ge_0_5_same_regulator=float(((jv >= 0.5) & same).mean()),
                 regulators_recovered=float(len(set(regr0.values) & set(regr.values)) / len(set(regr0.values))))
        # regulons weighted by the risk model's programs: are adverse/protective programs' regulons recovered?
        pj = best_jaccard(pg0, pg)
        # program activity in the held-out cohort (its own within-cohort z), full vs best LOCO match
        zh = zrows(expr.loc[:, samples.index[samples["cohort"] == held]])
        act = lambda gs: zh.loc[[g for g in gs if g in zh.index]].mean()  # noqa: E731
        for k, (j, kb) in pj.items():
            rr = np.corrcoef(act(pg0[k]), act(pg[kb]))[0, 1] if kb is not None else np.nan
            prog_rows.append({"held_out": held, "program": k, "n_genes": len(pg0[k]), "best_loco_program": kb,
                              "jaccard": j, "activity_r_heldout": rr, "risk_weight": w0.get(k, np.nan)})
        pr = pd.DataFrame([x for x in prog_rows if x["held_out"] == held])
        r.update(program_median_jaccard=float(pr["jaccard"].median()),
                 program_median_activity_r=float(pr["activity_r_heldout"].median()),
                 program_frac_activity_r_ge_0_8=float((pr["activity_r_heldout"] >= 0.8).mean()),
                 risk_programs_median_activity_r=float(pr.loc[pr["risk_weight"].abs() >= w0.abs().quantile(0.75),
                                                                "activity_r_heldout"].median()))
        # causal: full high-confidence driver -> regulator edges recovered in the LOCO flows (same direction)
        cf = os.path.join(res, "05_causal", mx, "filteredCausalResults.csv")
        if os.path.exists(cf):
            fl = pd.read_csv(cf, index_col=0)
            hl = pd.read_csv(os.path.join(res, "05_causal", mx, "highConfidenceCausalResults.csv"), index_col=0)
            drivers = set(fl["Mutation"])
            e0d = {e for e in e0 if e[0] in drivers}
            ef = set(zip(fl["Mutation"], fl["regulator_symbol"], np.sign(fl["MutationRegulatorEdge"])))
            eh = set(zip(hl["Mutation"], hl["regulator_symbol"], np.sign(hl["MutationRegulatorEdge"])))
            r.update(causal_full_edges_testable=len(e0d), causal_edges_recovered_filtered=len(e0d & ef) / max(len(e0d), 1),
                     causal_edges_recovered_highconf=len(e0d & eh) / max(len(e0d), 1))
            for d in ("MUT_CTNNB1", "MUT_TP53", "MUT_AXIN1"):
                ed = {e for e in e0d if e[0] == d}
                if ed:
                    r[f"causal_{d}_recovered_filtered"] = len(ed & ef) / len(ed)
        # risk: LOCO programs -> ridge trained in the other survival cohort -> scored in the held-out cohort
        if held in ("TCGA", "CLCA"):
            train = "CLCA" if held == "TCGA" else "TCGA"
            zl = pd.read_csv(os.path.join(res, "02_batch_corrected", f"expression_{mx}_z.csv"), index_col=0)
            tr_ids = samples.index[samples["cohort"] == train]
            Xtr = pd.DataFrame({k: zl.loc[[g for g in v if g in zl.index], tr_ids].mean() for k, v in pg.items()}).T
            Xte = pd.DataFrame({k: zh.loc[[g for g in v if g in zh.index]].mean() for k, v in pg.items()}).T
            for ep in ("RFS", "OS") if train == "TCGA" else ("RFS",):
                s = pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_{train}_{ep}_h36m_miner.csv"), index_col=0)
                s.columns = ["duration", "observed"]
                g = miner.guanRank(miner.kmAnalysis(s.copy(), "duration", "observed"))
                Xs = Xtr.sub(Xtr.mean(1), axis=0).div(Xtr.std(1).replace(0, 1), axis=0)
                np.random.seed(R["ridge"]["seed"])
                alpha, _, _, _ = miner.optimize_ridge_model(Xs[g.index], g, n_iter=R["ridge"]["n_iter"],
                                                            train_cut_high=R["class1_proportion"],
                                                            train_cut_low=R["ridge"]["train_cut_low"],
                                                            max_range=R["ridge"]["max_alpha"], range_step=R["ridge"]["alpha_step"])
                n = len(g)
                ids = list(g.index[:round(n * R["class1_proportion"])]) + list(g.index[round(n * R["ridge"]["train_cut_low"]):])
                y = np.r_[np.ones(round(n * R["class1_proportion"])), np.zeros(len(ids) - round(n * R["class1_proportion"]))]
                mdl = Ridge(alpha=alpha, random_state=0).fit(Xs[ids].T.values, y)
                Xt = zrows(Xte)
                sc = pd.Series(mdl.predict(Xt.fillna(0).T.values), index=Xt.columns)
                st = pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_{held}_{ep}_h36m_miner.csv"), index_col=0)
                st.columns = ["duration", "observed"]
                st = st.loc[st.index.intersection(sc.index)]
                zz = (sc[st.index] - sc[st.index].mean()) / sc[st.index].std()
                c = CoxPHFitter().fit(st.assign(x=zz), "duration", "observed").summary.loc["x"]
                r[f"risk_{ep}_train_{train}_c_index"] = concordance_index(st["duration"], -sc[st.index], st["observed"])
                r[f"risk_{ep}_train_{train}_hr_per_sd"] = c["exp(coef)"]
                r[f"risk_{ep}_train_{train}_hr_p"] = c["p"]
                plt_close()
        rows.append(r)
        log.info("%s held out:\n%s", held, pd.Series(r).to_string())
    os.makedirs(outdir, exist_ok=True)
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(outdir, "loco_summary.tsv"), sep="\t", index=False, float_format="%.4g")
    PR = pd.DataFrame(prog_rows)
    PR["label"] = PR["program"].map(pd.read_csv(p("config/program_labels.tsv"), sep="\t", comment="#",
                                                dtype={"program": str}).set_index("program")["label"])
    PR.to_csv(os.path.join(outdir, "loco_programs.tsv"), sep="\t", index=False, float_format="%.4g")
    figure(S, PR, os.path.join(outdir, "figures"))


def plt_close():
    import matplotlib.pyplot as plt
    plt.close("all")


def figure(S, PR, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import qc_plots as Q
    os.makedirs(outdir, exist_ok=True)
    col = Q.cohort_colors(COHORTS)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    ax = axes[0]
    for i, h in enumerate(S["held_out"]):
        v = PR.loc[PR["held_out"] == h, "activity_r_heldout"].dropna()
        ax.hist(v, bins=np.linspace(0, 1, 26), histtype="step", lw=1.6, color=col[h], label=f"{h} held out (median {v.median():.2f})")
    ax.set_xlabel("program activity r in the held-out cohort\n(full-network program vs best LOCO match)")
    ax.set_ylabel("programs")
    ax.legend(loc="upper left", fontsize=7)
    ax = axes[1]
    mets = [("regulon_frac_jaccard_ge_0_5", "regulons recovered\n(Jaccard ≥ 0.5)"),
            ("regulators_recovered", "regulators\nrecovered"), ("program_frac_activity_r_ge_0_8", "programs with\nactivity r ≥ 0.8"),
            ("causal_edges_recovered_filtered", "causal edges\nrecovered")]
    x = np.arange(len(mets))
    for j, (_, r) in enumerate(S.iterrows()):
        ax.bar(x + (j - 1) * 0.26, [r.get(m, np.nan) for m, _ in mets], width=0.24, color=col[r["held_out"]],
               label=f"{r['held_out']} held out")
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in mets], fontsize=7)
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of full network")
    ax.legend(fontsize=7)
    ax = axes[2]
    rr = []
    for _, r in S.iterrows():
        for k in r.index:
            if k.startswith("risk_") and k.endswith("_c_index") and pd.notna(r[k]):
                ep, tr = k.split("_")[1], k.split("_")[3]
                rr.append((f"{r['held_out']} {ep}\n(network w/o {r['held_out']}; model {tr})", r[k], col[r["held_out"]]))
    ax.barh(range(len(rr)), [v for _, v, _ in rr], color=[c for _, _, c in rr])
    ax.axvline(0.5, color=Q.AXIS, lw=0.8)
    ax.set_yticks(range(len(rr)))
    ax.set_yticklabels([l for l, _, _ in rr], fontsize=7)
    ax.set_xlim(0.45, 0.75)
    ax.set_xlabel("C-index in the held-out cohort")
    for i, (_, v, _) in enumerate(rr):
        ax.text(v + 0.003, i, f"{v:.3f}", va="center", fontsize=7)
    fig.suptitle("Leave-one-cohort-out stability", x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "l1_loco.png"), dpi=150, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["prepare", "filter", "compare"])
    ap.add_argument("--held", choices=COHORTS)
    ap.add_argument("--params", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    log = setup_logging(p(os.path.join(P["paths"]["results"], "08_validation", "loco")),
                        f"08b_loco_{args.action}{'_' + args.held if args.held else ''}")
    {"prepare": lambda: prepare(P, args.held, log), "filter": lambda: filt(P, args.held, log),
     "compare": lambda: compare(P, log)}[args.action]()


if __name__ == "__main__":
    main()
