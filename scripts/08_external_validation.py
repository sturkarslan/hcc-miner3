#!/usr/bin/env python
"""Step 08: external validation in cohorts never used for the network or the risk models.

For each cohort in config validation.cohorts (GSE14520 / LCI; ICGC LIRI-JP via HCCDB18):
  1. Expression -> project Ensembl IDs (symbol, then HGNC previous/alias), tumours only, one sample per
     patient; genes z-scored within the cohort (the same per-cohort scale as the discovery matrix).
  2. Regulon score = mean z of its genes present (>= validation.min_regulon_coverage of genes);
     program activity = mean of its regulon scores. On the discovery data this reproduces the MINER
     eigengene-based program activities (median r 1.000) and the step-06 risk score (r >= 0.997).
  3. Fixed step-06 ridge models (weights.tsv; never refit) applied to program activities standardized within
     the cohort -> risk score. Evaluated per endpoint on the 36-month horizon: C-index, Cox HR per SD, HR of the
     within-cohort top 20% vs rest, stage-adjusted HR per SD, Kaplan-Meier.
  4. Head-to-head: Cox models on known-subtype scores trained in TCGA (as in 07c), applied here; LR test of
     adding the MINER score to all known scores.
  5. Program replication: Cox z per program here vs the discovery meta-z (step 06).
  6. States: each sample assigned to the discovery state with the most correlated program-activity
     centroid; state mean risk vs median GuanRank.
Outputs: results/08_validation/<cohort>/{scores,evaluation,program_cox,head_to_head}.tsv;
         results/08_validation/{summary.tsv, figures/v1..v4}
"""

import argparse
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from hcc_common import load_params, miner_id_backmap, p, setup_logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import qc_plots as Q  # noqa: E402

OUT = "08_validation"
HERE = os.path.dirname(os.path.abspath(__file__))


def _module(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ------------------------------------------------------------------ loaders (-> genes x samples, clinical)

def load_geo_gse14520(cfg, log):
    x = pd.read_csv(p(cfg["expression"]), sep="\t", comment="!", index_col=0)
    ann = pd.read_csv(p(cfg["platform_annot"]), sep="\t", skiprows=cfg["annot_skiprows"], low_memory=False,
                      usecols=["ID", "Gene symbol"]).dropna()
    ann = ann[~ann["Gene symbol"].str.contains("///")]
    x = x.loc[x.index.intersection(ann["ID"])]
    x["symbol"] = ann.set_index("ID").loc[x.index, "Gene symbol"].values
    # one probe per gene: highest mean expression
    x = x.assign(m=x.drop(columns="symbol").mean(1)).sort_values("m", ascending=False).drop_duplicates("symbol")
    x = x.drop(columns="m").set_index("symbol")
    cl = pd.read_csv(p(cfg["clinical"]), sep="\t")
    cl = cl[(cl["Tissue Type"].str.strip() == "Tumor") & cl["Affy_GSM"].isin(x.columns)]
    cl = cl.drop_duplicates("ID").set_index("Affy_GSM")
    num = lambda s: pd.to_numeric(s, errors="coerce")  # noqa: E731
    clin = pd.DataFrame({"OS_time": num(cl["Survival months"]) * 30.44, "OS_event": num(cl["Survival status"]),
                         "RFS_time": num(cl["Recurr months"]) * 30.44, "RFS_event": num(cl["Recurr status"]),
                         "stage": cl["BCLC staging"].map({"0": 0, "A": 1, "B": 2, "C": 3}).astype(float)},
                        index=cl.index)
    log.info("GSE14520: %d genes (symbols), %d tumours with expression", x.shape[0], len(clin))
    return x[clin.index], clin


def load_hccdb(cfg, log):
    """HCCDB level-3 matrix (Entrez_ID, Symbol, samples; already log2) + sample/patient tables (transposed).
    Survival columns, units, event values and the stage column are set per cohort in the config."""
    x = pd.read_csv(p(cfg["expression"]), sep="\t")
    x = x.drop(columns=[c for c in ("Entrez_ID",) if c in x.columns]).dropna(subset=["Symbol"]).set_index("Symbol")
    sa = pd.read_csv(p(cfg["sample"]), sep="\t", index_col=0, header=None).T
    pa = pd.read_csv(p(cfg["patient"]), sep="\t", index_col=0, header=None).T.set_index("PATIENT_ID")
    sa["TYPE"] = sa["TYPE"].str.strip()
    sa = sa[(sa["TYPE"] == "HCC") & sa["SAMPLE_ID"].isin(x.columns)]
    n0 = len(sa)
    sa = sa.drop_duplicates("PATIENT_ID").set_index("SAMPLE_ID")      # one tumour per patient (first listed)
    x = x[sa.index].apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    if cfg.get("log2"):
        x = np.log2(x + 1)
    pt = pa.reindex(sa["PATIENT_ID"])
    num = lambda s: pd.to_numeric(s, errors="coerce").values  # noqa: E731
    clin = pd.DataFrame(index=sa.index)
    for ep in ("OS", "RFS"):
        c = cfg.get(ep.lower())
        if not c:
            continue
        clin[f"{ep}_time"] = num(pt[c["time"]]) * float(c["to_days"])
        ev = pt[c["status"]].astype(str).str.strip()
        clin[f"{ep}_event"] = np.where(pt[c["status"]].isna(), np.nan, ev.isin(c["event_values"]).astype(float).values)
    clin["stage"] = pt[cfg["stage_column"]].astype(str).str.strip().map(cfg["stage_map"]).astype(float).values
    log.info("HCCDB %s: %d genes; %d tumour samples, %d patients kept (one tumour each)", cfg["expression"],
             x.shape[0], n0, len(clin))
    return x, clin


LOADERS = {"geo_gse14520": load_geo_gse14520, "hccdb": load_hccdb}


# ------------------------------------------------------------------ scoring

def to_ensembl(x, mapper):
    ids = [mapper(s)[0] for s in x.index]
    x = x.assign(ens=ids).dropna(subset=["ens"])
    return x.groupby("ens").mean()


def zrows(x):
    return x.sub(x.mean(1), axis=0).div(x.std(1).replace(0, np.nan), axis=0)


def score_programs(z, regulons, programs, min_cov):
    reg, cov = {}, {}
    for k, genes in regulons.items():
        g = [x for x in genes if x in z.index]
        cov[k] = len(g) / len(genes)
        if cov[k] >= min_cov and len(g) >= 3:
            reg[k] = z.loc[g].mean()
    reg = pd.DataFrame(reg).T
    prog = pd.DataFrame({k: reg.loc[[r for r in v if r in reg.index]].mean() for k, v in programs.items()
                         if any(r in reg.index for r in v)}).T
    return reg, prog, pd.Series(cov)


def evaluate(score, srv, stage, top):
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    s = score[srv.index]
    zz = (s - s.mean()) / s.std()
    r = {"n": len(srv), "events": int(srv["observed"].sum()),
         "c_index": concordance_index(srv["duration"], -s, srv["observed"])}
    c = CoxPHFitter().fit(srv.assign(x=zz), "duration", "observed").summary.loc["x"]
    r.update(hr_per_sd=c["exp(coef)"], hr_per_sd_lo=c["exp(coef) lower 95%"], hr_per_sd_hi=c["exp(coef) upper 95%"],
             hr_per_sd_p=c["p"])
    high = (s >= s.quantile(1 - top)).astype(int)
    c = CoxPHFitter().fit(srv.assign(high=high), "duration", "observed").summary.loc["high"]
    r.update(hr_top=c["exp(coef)"], hr_top_lo=c["exp(coef) lower 95%"], hr_top_hi=c["exp(coef) upper 95%"], hr_top_p=c["p"])
    d = srv.assign(x=zz, stage=stage.reindex(srv.index)).dropna()
    if d["stage"].nunique() > 1:
        c = CoxPHFitter().fit(d, "duration", "observed").summary.loc["x"]
        r.update(hr_per_sd_stage_adj=c["exp(coef)"], hr_per_sd_stage_adj_p=c["p"], n_stage_adj=len(d))
    return r


def horizon(clin, ep, days):
    s = clin[[f"{ep}_time", f"{ep}_event"]].dropna()
    s.columns = ["duration", "observed"]
    over = s["duration"] > days
    s.loc[over, "duration"] = days
    s.loc[over, "observed"] = 0
    return s[s["duration"] > 0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--cohorts", default=None, help="comma list (default: all configured with existing files)")
    args = ap.parse_args()
    P = load_params(args.params)
    V = P["validation"]
    matrix = P["miner"]["matrix"]
    res = p(P["paths"]["results"])
    outroot = os.path.join(res, OUT)
    log = setup_logging(outroot, "08_external_validation")
    m07a = _module("m07a", "07a_subtype_signatures.py")
    m07c = _module("m07c", "07c_figures.py")
    Rf = P["post"]["figures"]

    genes = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    mapper = m07a.symbol_mapper(genes, p(P["id_mapping"].get("hgnc")), log)
    back = miner_id_backmap(p(P["miner"]["idmap"]), genes.index)
    mdir = os.path.join(res, "04_miner", matrix)
    regulons = {k: [back.get(g, g) for g in v]
                for k, v in json.load(open(os.path.join(mdir, "mechinf", "regulons_filtered.json"))).items()}
    programs = {k: [str(r) for r in v]
                for k, v in json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_programs.json"))).items()}
    states = json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_states.json")))
    rdir = os.path.join(res, "06_risk", matrix)
    models = {m: pd.read_csv(os.path.join(rdir, f"predictor_ridge_programs_{m}", "weights.tsv"), sep="\t",
                             index_col=0)["weight"].rename(lambda i: str(i)) for m in V["models"]}

    # discovery references: program activities (for state centroids), known-subtype models trained in TCGA
    zdisc = pd.read_csv(os.path.join(res, "02_batch_corrected", f"expression_{matrix}_z.csv"), index_col=0)
    samples = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    _, pdisc, _ = score_programs(zdisc, regulons, programs, 0.0)
    pdisc_std = pd.concat([zrows(pdisc[samples.index[samples["cohort"] == c]]) for c in samples["cohort"].unique()], axis=1)
    cent = pd.DataFrame({k: pdisc_std[[s for s in v if s in pdisc_std.columns]].mean(1) for k, v in states.items()})
    KSd = m07c.known_scores(P, None, Rf)
    sig = pd.read_csv(os.path.join(res, "07_post", "signatures", "signatures.tsv"), sep="\t").dropna(subset=["ensembl"])
    disc_meta = {ep: pd.read_csv(os.path.join(rdir, f"prognostic_programs_h36m_{ep}.tsv"), sep="\t", index_col=0)
                 for ep in ("RFS", "OS")}

    from lifelines import CoxPHFitter
    # known-signature models from the reference panel (config/reference_panel.yaml), as in 07c
    known_sets = {name: cols for name, cols, *_ in m07c.PN.known_models(KSd.columns)}
    tcga_srv = {ep: pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_TCGA_{ep}_h36m_miner.csv"),
                                index_col=0).set_axis(["duration", "observed"], axis=1) for ep in ("RFS", "OS")}

    summary, cohort_data = [], {}
    names = args.cohorts.split(",") if args.cohorts else list(V["cohorts"])
    for name in names:
        cfg = V["cohorts"][name]
        if not os.path.exists(p(cfg["expression"])):
            log.warning("%s: expression file %s not found; skipped", name, cfg["expression"])
            continue
        outdir = os.path.join(outroot, name)
        os.makedirs(outdir, exist_ok=True)
        x, clin = LOADERS[cfg["loader"]](cfg, log)
        xe = to_ensembl(x, mapper)
        z = zrows(xe).dropna(how="all")
        reg, prog, cov = score_programs(z, regulons, programs, V["min_regulon_coverage"])
        cov_net = len(set(z.index) & set().union(*map(set, regulons.values()))) / len(set().union(*map(set, regulons.values())))
        log.info("%s: %d genes mapped; %.0f%% of network genes present; regulons scored %d / %d (median coverage "
                 "%.0f%%); programs %d / %d", name, z.shape[0], 100 * cov_net, reg.shape[0], len(regulons),
                 100 * cov.median(), prog.shape[0], len(programs))
        ps = zrows(prog)
        scores = pd.DataFrame(index=z.columns)
        for m, w in models.items():
            ww = w.reindex(ps.index).fillna(0)
            scores[f"risk_{m}"] = ps.fillna(0).T @ ww
        # known-subtype scores (global mean regressed out), same construction as 07c
        KS = m07c.PN.known_scores(z, sig)
        # states: nearest discovery centroid by correlation
        common = ps.index.intersection(cent.index)
        cc = np.corrcoef(np.c_[ps.loc[common].fillna(0).values, cent.loc[common].values].T)[:ps.shape[1], ps.shape[1]:]
        scores["state"] = [cent.columns[i] for i in cc.argmax(1)]
        scores["state_r"] = cc.max(1)
        scores = scores.join(clin)
        scores.to_csv(os.path.join(outdir, "scores.tsv"), sep="\t", float_format="%.5g")

        ev, h2h, pc = [], [], []
        for ep in cfg["endpoints"]:
            srv = horizon(clin, ep, V["horizon_days"])
            for m in models:
                if not m.endswith(f"_{ep}_h36m") and not V.get("cross_endpoint", False):
                    continue
                r = evaluate(scores[f"risk_{m}"], srv, clin["stage"], V["top_fraction"])
                r.update(cohort=name, endpoint=ep, model=f"MINER programs ({m.split('_')[0]}-trained)")
                ev.append(r)
                log.info("%s %s | %s: n=%d ev=%d C=%.3f HR/SD %.2f (%.2f-%.2f) p=%.1e | top %.0f%% HR %.2f p=%.1e | "
                         "stage-adj HR/SD %.2f p=%.1e", name, ep, m, r["n"], r["events"], r["c_index"], r["hr_per_sd"],
                         r["hr_per_sd_lo"], r["hr_per_sd_hi"], r["hr_per_sd_p"], 100 * V["top_fraction"], r["hr_top"],
                         r["hr_top_p"], r.get("hr_per_sd_stage_adj", np.nan), r.get("hr_per_sd_stage_adj_p", np.nan))
            # head-to-head: known-score Cox models trained in TCGA on the same endpoint
            a = tcga_srv[ep]
            Xa = KSd.loc[a.index]
            b = srv
            Xb = KS.loc[b.index]
            for kname, cols in known_sets.items():
                cols = [c for c in cols if c in Xb.columns]
                mu, sd = Xa[cols].mean(), Xa[cols].std()
                mdl = CoxPHFitter(penalizer=0.1).fit(a.join((Xa[cols] - mu) / sd), "duration", "observed")
                sc = mdl.predict_partial_hazard((Xb[cols] - Xb[cols].mean()) / Xb[cols].std())
                from lifelines.utils import concordance_index
                h2h.append({"cohort": name, "endpoint": ep, "model": kname,
                            "c_index": concordance_index(b["duration"], -sc, b["observed"])})
            mkey = f"risk_TCGA_{ep}_h36m" if f"risk_TCGA_{ep}_h36m" in scores else scores.columns[0]
            d = b.join((Xb - Xb.mean()) / Xb.std())
            d["miner"] = (scores.loc[b.index, mkey] - scores.loc[b.index, mkey].mean()) / scores.loc[b.index, mkey].std()
            full = CoxPHFitter().fit(d, "duration", "observed")
            red = CoxPHFitter().fit(d.drop(columns="miner"), "duration", "observed")
            lr = 2 * (full.log_likelihood_ - red.log_likelihood_)
            h2h.append({"cohort": name, "endpoint": ep, "model": "MINER programs (TCGA-trained)",
                        "c_index": concordance_index(b["duration"], -scores.loc[b.index, mkey], b["observed"]),
                        "lr_added_to_known": lr, "lr_p": stats.chi2.sf(lr, 1),
                        "hr_per_sd_given_known": full.summary.loc["miner", "exp(coef)"]})
            # program replication
            for k in ps.index:
                dd = b.assign(x=ps.loc[k, b.index].values).dropna()
                try:
                    zc = CoxPHFitter().fit(dd, "duration", "observed").summary.loc["x", "z"]
                except Exception:  # noqa: BLE001
                    zc = np.nan
                pc.append({"cohort": name, "endpoint": ep, "program": k, "z_external": zc,
                           "meta_z_discovery": disc_meta[ep]["meta_z"].get(int(k), np.nan),
                           "discovery_prognostic": bool((disc_meta[ep]["meta_q"].get(int(k), 1) <= Rf["q_max"])
                                                        and disc_meta[ep]["consistent"].get(int(k), False))})
            # states: mean risk vs median GuanRank
            gr = m07c.guan_rank(srv)
            st = scores.loc[srv.index].assign(g=gr).groupby("state").agg(n=("g", "size"), g=("g", "median"),
                                                                         risk=(mkey, "mean"))
            st = st[st["n"] >= 5]
            rho, pr_ = stats.spearmanr(st["risk"], st["g"]) if len(st) >= 4 else (np.nan, np.nan)
            log.info("%s %s: %d states with >= 5 samples; state mean risk vs median GuanRank rho %.2f p %.2g",
                     name, ep, len(st), rho, pr_)
            summary.append({"cohort": name, "endpoint": ep, "states_ge5": len(st), "state_rho": rho, "state_p": pr_})
        ev, h2h, pc = pd.DataFrame(ev), pd.DataFrame(h2h), pd.DataFrame(pc)
        ev.to_csv(os.path.join(outdir, "evaluation.tsv"), sep="\t", index=False, float_format="%.4g")
        h2h.to_csv(os.path.join(outdir, "head_to_head.tsv"), sep="\t", index=False, float_format="%.4g")
        pc.to_csv(os.path.join(outdir, "program_cox.tsv"), sep="\t", index=False, float_format="%.4g")
        for ep in pc["endpoint"].unique():
            q = pc[(pc["endpoint"] == ep) & pc["discovery_prognostic"]].dropna()
            agree = (np.sign(q["z_external"]) == np.sign(q["meta_z_discovery"])).mean()
            rho, pr_ = stats.spearmanr(pc.loc[pc["endpoint"] == ep, "z_external"], pc.loc[pc["endpoint"] == ep, "meta_z_discovery"],
                                       nan_policy="omit")
            log.info("%s %s: discovery-prognostic programs with the same sign here: %d / %d (%.0f%%); "
                     "all programs z vs discovery meta-z rho %.2f p %.1e", name, ep, int((np.sign(q['z_external']) ==
                     np.sign(q['meta_z_discovery'])).sum()), len(q), 100 * agree, rho, pr_)
        log.info("%s head-to-head:\n%s", name, h2h.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        cohort_data[name] = (scores, clin, ev, h2h, pc)

    if summary:
        pd.DataFrame(summary).to_csv(os.path.join(outroot, "summary.tsv"), sep="\t", index=False, float_format="%.4g")
    if cohort_data:
        written = figures(cohort_data, V, os.path.join(outroot, "figures"))
        log.info("Figures: %s", ", ".join(written))


def figures(cohort_data, V, outdir):
    from lifelines import KaplanMeierFitter  # noqa: F401  (qc_plots.kaplan_meier is used)
    written = []
    panels = [(c, ep) for c, (s, cl, ev, h, pc) in cohort_data.items() for ep in ev["endpoint"].unique()]
    fig, axes = plt.subplots(1, len(panels), figsize=(4.9 * len(panels), 4.3), squeeze=False)
    for ax, (c, ep) in zip(axes[0], panels):
        s, cl, ev, h, pc = cohort_data[c]
        srv = horizon(cl, ep, V["horizon_days"])
        key = f"risk_TCGA_{ep}_h36m" if f"risk_TCGA_{ep}_h36m" in s else [k for k in s if k.startswith("risk_")][0]
        sc = s.loc[srv.index, key]
        high = sc >= sc.quantile(1 - V["top_fraction"])
        for grp, col, lab in ((True, Q.SLOTS[7], "predicted high risk (top 20%)"), (False, Q.SLOTS[0], "rest")):
            ids = srv.index[high == grp]
            ts, ss, t, e = Q.kaplan_meier(srv.loc[ids, "duration"] / 30.44, srv.loc[ids, "observed"])
            ax.step(np.append(ts, t.max()), np.append(ss, ss[-1]), where="post", color=col,
                    label=f"{lab} (n={len(t)}, events={int(e.sum())})")
        r = ev[(ev["endpoint"] == ep) & ev["model"].str.contains("TCGA")].iloc[0]
        ax.set_title(f"{c} — {ep} (external)\nC = {r['c_index']:.2f}; HR top vs rest {r['hr_top']:.2f} "
                     f"(p = {r['hr_top_p']:.1e}); stage-adj HR/SD {r.get('hr_per_sd_stage_adj', np.nan):.2f}",
                     fontsize=8.5, loc="left")
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("months")
        ax.set_ylabel("survival probability")
        ax.legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    Q._save(fig, outdir, "v1_km_external.png", written)
    # v2: head-to-head C-index
    H = pd.concat([v[3] for v in cohort_data.values()])
    order = [m for m in dict.fromkeys(H["model"]) if not m.startswith("MINER")] + ["MINER programs (TCGA-trained)"]
    groups = H[["cohort", "endpoint"]].drop_duplicates().values.tolist()
    fig, ax = plt.subplots(figsize=(8.5, 1.2 + 0.9 * len(order)))
    y = np.arange(len(order))
    for j, (c, ep) in enumerate(groups):
        h = H[(H["cohort"] == c) & (H["endpoint"] == ep)].set_index("model").reindex(order)
        off = (j - (len(groups) - 1) / 2) * (0.8 / len(groups))
        ax.barh(y + off, h["c_index"], height=0.75 / len(groups), color=Q.SLOTS[j], label=f"{c} {ep}")
        for yi, v in zip(y + off, h["c_index"]):
            ax.text(v + 0.003, yi, f"{v:.3f}", va="center", fontsize=6.5, color=Q.INK2)
    ax.axvline(0.5, color=Q.AXIS, lw=0.8)
    ax.set_xlim(0.45, 0.72)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("C-index (36 months; all models trained in TCGA, never refit)")
    m = H[H["model"].str.startswith("MINER")]
    ax.set_title("External cohorts: known subtype models vs MINER\nadding MINER to all known scores (LR): " +
                 "; ".join(f"{a} {b} p={c_:.2g}" for a, b, c_ in m[["cohort", "endpoint", "lr_p"]].values),
                 loc="left", fontsize=8.5)
    ax.legend(loc="lower right", fontsize=7)
    fig.tight_layout()
    Q._save(fig, outdir, "v2_head_to_head_external.png", written)
    # v3: program replication
    PC = pd.concat([v[4] for v in cohort_data.values()])
    combos = PC[["cohort", "endpoint"]].drop_duplicates().values.tolist()
    fig, axes = plt.subplots(1, len(combos), figsize=(4.6 * len(combos), 4.4), squeeze=False)
    for ax, (c, ep) in zip(axes[0], combos):
        q = PC[(PC["cohort"] == c) & (PC["endpoint"] == ep)].dropna(subset=["z_external", "meta_z_discovery"])
        ax.axhline(0, color=Q.AXIS, lw=0.8)
        ax.axvline(0, color=Q.AXIS, lw=0.8)
        ax.scatter(q.loc[~q["discovery_prognostic"], "meta_z_discovery"], q.loc[~q["discovery_prognostic"], "z_external"],
                   s=14, color=Q.BACKGROUND, edgecolor=Q.MUTED, linewidth=0.5, label="not prognostic in discovery")
        qq = q[q["discovery_prognostic"]]
        ax.scatter(qq["meta_z_discovery"], qq["z_external"], s=20,
                   color=[Q.SLOTS[7] if v > 0 else Q.SLOTS[0] for v in qq["meta_z_discovery"]],
                   label="prognostic in discovery")
        rho, pv = stats.spearmanr(q["meta_z_discovery"], q["z_external"])
        agree = (np.sign(qq["z_external"]) == np.sign(qq["meta_z_discovery"])).mean()
        ax.set_title(f"{c} {ep}: ρ = {rho:.2f} (p = {pv:.1e})\nsame sign: {agree:.0%} of {len(qq)} discovery-prognostic",
                     fontsize=8.5, loc="left")
        ax.set_xlabel("discovery meta-z (TCGA + CLCA)")
        ax.set_ylabel(f"Cox z in {c}")
        ax.legend(loc="upper left", fontsize=6.5)
    fig.tight_layout()
    Q._save(fig, outdir, "v3_program_replication.png", written)
    # v4: forest of HR per SD (TCGA-trained models) with a fixed-effect pooled estimate per endpoint
    E = pd.concat([v[2] for v in cohort_data.values()])
    E = E[E["model"].str.contains("TCGA")].copy()
    E["b"] = np.log(E["hr_per_sd"])
    E["se"] = (np.log(E["hr_per_sd_hi"]) - np.log(E["hr_per_sd_lo"])) / (2 * 1.96)
    rows = []
    for ep, g in E.groupby("endpoint"):
        for _, r in g.iterrows():
            rows.append((f"{r['cohort']} {ep} (n={int(r['n'])}, events={int(r['events'])})", r["b"], r["se"], False, r["c_index"]))
        w = 1 / g["se"] ** 2
        b = (w * g["b"]).sum() / w.sum()
        se = np.sqrt(1 / w.sum())
        qh = (w * (g["b"] - b) ** 2).sum()
        k = len(g)
        rows.append((f"pooled {ep}, fixed effect (heterogeneity p = {stats.chi2.sf(qh, k - 1):.3f})", b, se, True, np.nan))
        # DerSimonian-Laird random effects
        tau2 = max(0.0, (qh - (k - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if k > 1 else 0.0
        wr = 1 / (g["se"] ** 2 + tau2)
        br = (wr * g["b"]).sum() / wr.sum()
        rows.append((f"pooled {ep}, random effects (τ² = {tau2:.3f})", br, np.sqrt(1 / wr.sum()), True, np.nan))
    fig, ax = plt.subplots(figsize=(8.8, 0.9 + 0.42 * len(rows)))
    for i, (lab, b, se, pooled, c) in enumerate(rows):
        lo, hi = np.exp(b - 1.96 * se), np.exp(b + 1.96 * se)
        pz = 2 * stats.norm.sf(abs(b / se))
        ax.plot([lo, hi], [i, i], color=Q.INK if pooled else Q.MUTED, lw=2 if pooled else 1.4)
        ax.scatter(np.exp(b), i, marker="D" if pooled else "s", s=60 if pooled else 36,
                   color=Q.SLOTS[7] if pooled else Q.INK2, zorder=3)
        txt = f"HR/SD {np.exp(b):.2f} ({lo:.2f}-{hi:.2f}), p = {pz:.1e}" + ("" if pooled else f", C = {c:.2f}")
        ax.text(3.9, i, txt, va="center", fontsize=7, color=Q.INK2)
    ax.axvline(1, color=Q.AXIS, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(0.4, 3.8)
    ax.set_xticks([0.5, 1, 2, 3])
    ax.set_xticklabels(["0.5", "1", "2", "3"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows], fontsize=7.5)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Cox HR per SD of the program risk score (36 months; TCGA-trained, never refit)")
    ax.set_title("External validation of the MINER program risk score", loc="left", fontsize=9.5)
    fig.subplots_adjust(right=0.62)
    Q._save(fig, outdir, "v4_forest_external.png", written)
    return written


if __name__ == "__main__":
    main()
