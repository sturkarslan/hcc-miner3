#!/usr/bin/env python
"""Step 07c: figures relating MINER programs and states to published HCC subtypes and to risk.

f1_integrated_map      programs x states (mean regulon dysregulation, over - under), both clustered
                       (average linkage, correlation distance) with dendrograms. Tracks above: state size, mean
                       risk, cohort mix, published-class composition (reference panel f07c_tracks: Montironi 2023,
                       Sia 2017, WNT activation, Hoshida/Boyault/Chiang as references, LICA-FR author labels),
                       driver alteration frequencies. Left: risk weight and
                       prognostic meta-z per program. Right: top regulators and best known signature.
f2_known_vs_new        per program: best gene-set overlap (-log10 FDR) vs best activity correlation with
                       any known signature; colour = prognostic meta-z. "New" = prognostic, no overlap at
                       FDR < 0.05 with >= min_fold, and |r| < known_r.
f3_risk_by_subtype     risk score (within-cohort z) by published class, small multiples, per sample.
f4_program_signatures  programs x reference-panel signatures (config/reference_panel.yaml), both clustered: activity r;
                       dot = direct gene overlap (FDR < 0.05, >= min_fold).
f5_state_risk          state-level Cox z per cohort + meta-z, with each state's dominant subtype calls.
f7_state_risk_survival states ordered by mean risk score: subtype/mutation enrichment per state (fraction,
                       dot = one-sided Fisher FDR < 0.05), risk score and within-cohort GuanRank (RFS) per state.
f6_beyond_known_subtypes  per program: prognostic meta-z before vs after adjusting (within cohort) for known
                       subtype scores (reference panel adjust_signatures) and the per-sample global mean z;
                       the risk score is tested the same way (risk_score_beyond_known.tsv). Panel b:
                       cross-cohort C-index of Cox models on known-subtype scores vs the MINER program
                       model, and the LR test of adding the MINER score (risk_head_to_head.tsv).

Risk score = MINER ridge on program activity trained in TCGA (RFS, 36-month horizon; step 06). It is
in-sample for TCGA and out-of-sample for CLCA and LICA-FR; shown as a within-cohort z.
Tables: program_annotation.tsv, state_annotation_risk.tsv, program_labels_proposed.tsv (labels from the
reference panel, same format as config/program_labels.tsv, for curation) (results/07_post/figures/).
"""

import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from hcc_common import load_params, p, setup_logging
from hcc_panel import Panel

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402

import qc_plots as Q  # noqa: E402  (palette, rcParams, _save, cohort_colors)

OUT = os.path.join("07_post", "figures")
MUTS = ["MUT_CTNNB1", "MUT_AXIN1", "MUT_TP53", "TERT_promoter", "MUT_ARID1A", "MUT_ALB", "AMP_CCND1_FGF19",
        "ARM_8p_loss", "ARM_1q_gain"]
PN = Panel()   # reference signature panel (config/reference_panel.yaml): signatures, classes, known models


def short(sig):
    return PN.label(sig)


def load(P, matrix, log):
    res = p(P["paths"]["results"])
    sub = P["post"]["subtypes_dir"]
    mdir = os.path.join(res, "04_miner", matrix, sub)
    D = {}
    D["samples"] = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    E = pd.read_csv(os.path.join(mdir, "eigengenes.csv"), index_col=0)
    E.index = E.index.astype(str)
    over = pd.read_csv(os.path.join(mdir, "overExpressedMembers.csv"), index_col=0)
    under = pd.read_csv(os.path.join(mdir, "underExpressedMembers.csv"), index_col=0)
    over.index, under.index = over.index.astype(str), under.index.astype(str)
    D["diff"] = over - under
    D["programs"] = {k: [str(r) for r in v] for k, v in json.load(open(os.path.join(mdir, "transcriptional_programs.json"))).items()}
    D["states"] = json.load(open(os.path.join(mdir, "transcriptional_states.json")))
    Ez = E.sub(E.mean(1), axis=0)
    D["act"] = pd.DataFrame({k: Ez.loc[[r for r in v if r in Ez.index]].mean() for k, v in D["programs"].items()}).T
    rdir = os.path.join(res, "06_risk", matrix)
    R = P["post"]["figures"]
    pred = pd.read_csv(os.path.join(rdir, f"predictor_{R['risk_model']}", "predictions.tsv"), sep="\t", index_col=0)
    D["risk"] = pred.groupby("cohort")["risk_score"].transform(lambda v: (v - v.mean()) / v.std())
    D["weights"] = pd.read_csv(os.path.join(rdir, f"predictor_{R['risk_model']}", "weights.tsv"), sep="\t",
                               index_col=0)["weight"]
    D["weights"].index = D["weights"].index.astype(str)
    ep = R["endpoint"]
    for unit in ("programs", "states"):
        t = pd.read_csv(os.path.join(rdir, f"prognostic_{unit}_h36m_{ep}.tsv"), sep="\t", index_col=0)
        t.index = t.index.astype(str)
        D[f"prog_{unit}"] = t
    smap = os.path.join(res, "07_post", "subtype_mapping", sub)
    D["cor"] = pd.read_csv(os.path.join(smap, "program_signature_correlation.tsv"), sep="\t", index_col=0)
    D["cor"].index = D["cor"].index.astype(str)
    ov = pd.read_csv(os.path.join(smap, "program_signature_overlap.tsv"), sep="\t")
    ov["program"] = ov["program"].astype(str)
    D["overlap"] = ov
    D["ntp"] = PN.class_calls(smap, D["samples"])
    fam = pd.read_csv(os.path.join(res, "05_causal", matrix, "regulon_families.tsv"), sep="\t", index_col=0)
    fam["program"] = fam["program"].astype(str)
    D["fam"] = fam
    D["genomic"] = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
    log.info("Loaded %d programs, %d states, %d samples", len(D["programs"]), len(D["states"]), len(D["samples"]))
    return D


def known_scores(P, D, R):
    """Per-sample scores of the reference-panel adjust signatures (mean z of their genes; signed entries
    "A - B" = up minus down), each with the sample's global mean z regressed out; plus the global mean itself.
    Columns are the panel entries (display labels: PN.label)."""
    res = p(P["paths"]["results"])
    z = pd.read_csv(os.path.join(res, "02_batch_corrected", f"expression_{P['miner']['matrix']}_z.csv"), index_col=0)
    sig = pd.read_csv(os.path.join(res, "07_post", "signatures", "signatures.tsv"), sep="\t").dropna(subset=["ensembl"])
    return PN.known_scores(z, sig)


def cox_adjusted(values, surv, covars, R):
    """Cox z of each row of `values` (units x samples) within each cohort, alone and adjusted for covars;
    Stouffer meta-z over cohorts (weights sqrt(events)); BH; consistency."""
    from lifelines import CoxPHFitter
    zs = {"unadjusted": {}, "adjusted": {}}
    w = {}
    for c, srv in surv.items():
        ids = srv.index.intersection(values.columns).intersection(covars.index)
        srv = srv.loc[ids]
        w[c] = np.sqrt(srv["observed"].sum())
        cv = covars.loc[ids]
        cv = (cv - cv.mean()) / cv.std()
        for mode in zs:
            col = {}
            for k in values.index:
                x = values.loc[k, ids].astype(float)
                if x.std() == 0:
                    col[k] = np.nan
                    continue
                d = srv.assign(x=(x - x.mean()) / x.std())
                if mode == "adjusted":
                    d = d.join(cv)
                try:
                    col[k] = CoxPHFitter(penalizer=0.0).fit(d, "duration", "observed").summary.loc["x", "z"]
                except Exception:  # noqa: BLE001
                    col[k] = np.nan
            zs[mode][c] = pd.Series(col)
    out = pd.DataFrame(index=values.index)
    for mode, per in zs.items():
        Z = pd.DataFrame(per)
        ww = np.array([w[c] for c in Z.columns])
        for c in Z.columns:
            out[f"z_{mode}_{c}"] = Z[c]
        meta = (Z.values * ww).sum(1) / np.sqrt((ww ** 2).sum())
        out[f"meta_z_{mode}"] = meta
        pv = 2 * stats.norm.sf(np.abs(meta))
        out[f"meta_q_{mode}"] = _bh(pv)
        out[f"consistent_{mode}"] = np.sign(Z.values).min(1) == np.sign(Z.values).max(1)
    return out


def _bh(pv):
    pv = np.asarray(pv, float)
    q = np.full_like(pv, np.nan)
    ok = ~np.isnan(pv)
    n = ok.sum()
    o = np.argsort(pv[ok])
    r = pv[ok][o] * n / (np.arange(n) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(n)
    qq[o] = np.minimum(r, 1)
    q[ok] = qq
    return q


def load_surv(P, R):
    res = p(P["paths"]["results"])
    out = {}
    for c in ("TCGA", "CLCA"):
        s_ = pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_{c}_{R['endpoint']}_h36m_miner.csv"), index_col=0)
        s_.columns = ["duration", "observed"]
        out[c] = s_
    return out


def head_to_head(P, KS, R):
    """Cross-cohort C-index of Cox models on known-subtype scores (trained in one cohort, applied to the
    other) vs the step-06 MINER ridge program models; and the likelihood-ratio test of adding the MINER
    score to the known scores within the test cohort."""
    from lifelines import CoxPHFitter
    from lifelines.utils import concordance_index
    res = p(P["paths"]["results"])
    surv = load_surv(P, R)
    sets = {name: cols for name, cols, *_ in PN.known_models(KS.columns)}
    rows = []
    for tr, te in (("TCGA", "CLCA"), ("CLCA", "TCGA")):
        a, b = surv[tr], surv[te]
        Xa, Xb = KS.loc[a.index], KS.loc[b.index]
        for name, cols in sets.items():
            mu, sd = Xa[cols].mean(), Xa[cols].std()
            m = CoxPHFitter(penalizer=0.1).fit(a.join((Xa[cols] - mu) / sd), "duration", "observed")
            sc = m.predict_partial_hazard((Xb[cols] - mu) / sd)
            rows.append({"train": tr, "test": te, "model": name, "c_index": concordance_index(b["duration"], -sc, b["observed"])})
        mp = pd.read_csv(os.path.join(res, "06_risk", P["miner"]["matrix"], f"predictor_ridge_programs_{tr}_{R['endpoint']}_h36m",
                                      "predictions.tsv"), sep="\t", index_col=0)["risk_score"]
        rows.append({"train": tr, "test": te, "model": "MINER programs (ridge)",
                     "c_index": concordance_index(b["duration"], -mp[b.index], b["observed"])})
        d = b.join((Xb - Xb.mean()) / Xb.std())
        d["miner"] = (mp[b.index] - mp[b.index].mean()) / mp[b.index].std()
        full = CoxPHFitter().fit(d, "duration", "observed")
        red = CoxPHFitter().fit(d.drop(columns="miner"), "duration", "observed")
        lr = 2 * (full.log_likelihood_ - red.log_likelihood_)
        rows[-1].update(lr_added_to_known=lr, lr_p=stats.chi2.sf(lr, 1),
                        hr_per_sd_given_known=full.summary.loc["miner", "exp(coef)"])
    return pd.DataFrame(rows)


def fig_beyond_known(prog, adj, risk_adj, h2h, outdir, R, written):
    """Prognostic meta-z of each program before vs after adjusting for known-subtype scores."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13, 6.2), gridspec_kw={"width_ratios": [1.15, 1]})
    lim = np.nanmax(np.abs(adj[["meta_z_unadjusted", "meta_z_adjusted"]].values)) * 1.08
    new = (adj["meta_q_adjusted"] <= R["q_max"]) & adj["consistent_adjusted"]
    ax.axhspan(-3, 3, color=Q.GRID, alpha=0.5, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color=Q.MUTED, ls="--", lw=1)
    ax.axhline(0, color=Q.AXIS, lw=0.8)
    ax.axvline(0, color=Q.AXIS, lw=0.8)
    ax.scatter(adj.loc[~new, "meta_z_unadjusted"], adj.loc[~new, "meta_z_adjusted"], s=18, color=Q.BACKGROUND,
               edgecolor=Q.MUTED, linewidth=0.5, label="not prognostic after adjustment", zorder=2)
    if new.any():
        ax.scatter(adj.loc[new, "meta_z_unadjusted"], adj.loc[new, "meta_z_adjusted"], s=34,
                   color=[Q.SLOTS[7] if v > 0 else Q.SLOTS[0] for v in adj.loc[new, "meta_z_adjusted"]],
                   edgecolor=Q.SURFACE, linewidth=1,
                   label="prognostic beyond known subtypes (meta q ≤ %.2f, same sign)" % R["q_max"], zorder=3)
    for k in adj.index[new]:
        ax.annotate(f"P{k}", (adj.loc[k, "meta_z_unadjusted"], adj.loc[k, "meta_z_adjusted"]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points", color=Q.INK)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("prognostic meta-z, unadjusted (> 0 adverse)")
    ax.set_ylabel("meta-z adjusted for known subtype scores + global mean")
    ax.set_title(f"a. Programs: prognosis before vs after adjustment ({int(new.sum())} of {len(adj)} remain)",
                 loc="left", fontsize=9.5)
    import textwrap
    txt = "; ".join(f"{c} HR/SD {r['hr']:.2f} (p={r['p']:.1e})" for c, r in risk_adj.iterrows())
    ax.text(0.01, 0.01, "\n".join(textwrap.wrap("adjusted within cohort for: " + ", ".join(R["adjust_labels"]), 80))
            + "\nMINER risk score, adjusted: " + txt, transform=ax.transAxes, fontsize=6.5, color=Q.INK2, va="bottom")
    ax.legend(loc="upper left", fontsize=7)
    # b. head-to-head C-index in the external cohort
    order = [m for m in dict.fromkeys(h2h["model"]) if m != "MINER programs (ridge)"] + ["MINER programs (ridge)"]
    y = np.arange(len(order))
    for j, (tr, te) in enumerate((("TCGA", "CLCA"), ("CLCA", "TCGA"))):
        h = h2h[(h2h["train"] == tr)].set_index("model").reindex(order)
        ax2.barh(y + (j - 0.5) * 0.36, h["c_index"], height=0.34, color=Q.cohort_colors(["TCGA", "CLCA"])[te],
                 label=f"trained in {tr}, tested in {te}")
        for yi, v in zip(y + (j - 0.5) * 0.36, h["c_index"]):
            ax2.text(v + 0.003, yi, f"{v:.3f}", va="center", fontsize=7, color=Q.INK2)
    ax2.axvline(0.5, color=Q.AXIS, lw=0.8)
    ax2.set_xlim(0.45, 0.70)
    ax2.set_yticks(y)
    ax2.set_yticklabels(order, fontsize=8)
    ax2.invert_yaxis()
    ax2.grid(axis="y", visible=False)
    ax2.set_xlabel("C-index in the external cohort (RFS, 36 months)")
    m = h2h[h2h["model"] == "MINER programs (ridge)"].set_index("test")
    ax2.set_title("b. Known subtype models vs MINER, cross-cohort", loc="left", fontsize=9.5)
    ax2.text(0.01, 0.01, "adding the MINER score to all known scores (LR test): " +
             "; ".join(f"in {t} p={r['lr_p']:.2g}" for t, r in m.iterrows()), transform=ax2.transAxes,
             fontsize=6.8, color=Q.INK2, va="bottom")
    ax2.legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    Q._save(fig, outdir, "f6_beyond_known_subtypes.png", written)


def annotate_programs(D, R):
    """Per program: weight, prognostic meta-z, top regulators, best signature by activity and by overlap,
    known/new call."""
    rows = []
    ov = D["overlap"]
    for k in D["programs"]:
        o = ov[(ov["program"] == k) & (ov["fdr"] < 0.05) & (ov["fold_enrichment"] >= R["min_fold"])].sort_values("fdr")
        c = D["cor"].loc[k]
        pr = D["prog_programs"].loc[k]
        regs = D["fam"].loc[D["fam"]["program"] == k, "regulator_symbol"].value_counts().index[:5]
        cp = c[[s_["set"] for s_ in PN.signatures() if s_["set"] in c.index]]
        bp = cp.abs().idxmax() if len(cp) else ""
        rows.append({"program": k, "n_regulons": len(D["programs"][k]), "risk_weight": D["weights"].get(k, np.nan),
                     "meta_z": pr["meta_z"], "meta_q": pr["meta_q"], "consistent": bool(pr["consistent"]),
                     "top_regulators": ",".join(regs),
                     "best_activity_signature": c.abs().idxmax(), "best_activity_r": c[c.abs().idxmax()],
                     "best_panel_signature": bp, "best_panel_r": cp[bp] if bp else np.nan,
                     "best_overlap_signature": o["signature"].iloc[0] if len(o) else "",
                     "best_overlap_fold": o["fold_enrichment"].iloc[0] if len(o) else np.nan,
                     "best_overlap_neglog10_fdr": -np.log10(o["fdr"].iloc[0]) if len(o) else 0.0})
    t = pd.DataFrame(rows).set_index("program")
    t["prognostic"] = (t["meta_q"] <= R["q_max"]) & t["consistent"]
    t["known"] = (t["best_overlap_signature"] != "") | (t["best_activity_r"].abs() >= R["known_r"])
    t["class"] = np.select([t["prognostic"] & ~t["known"], t["prognostic"], t["known"]],
                           ["prognostic, new", "prognostic, known", "known, not prognostic"], "other")
    return t


def proposed_labels(prog, outdir, R, log):
    """Program labels from the reference panel, in config/program_labels.tsv format, for curation:
    best panel signature when |r| >= R['label_r'] (default 0.6; "low ..." when r < 0), else named by
    the top regulators. Evidence lists the panel hit, the best signature overall and the regulators."""
    try:
        cur = pd.read_csv(p("config/program_labels.tsv"), sep="\t", comment="#", dtype={"program": str}) \
            .set_index("program")["label"]
    except FileNotFoundError:
        cur = pd.Series(dtype=str)
    thr = R.get("label_r", 0.6)
    rows = []
    for k, r in prog.iterrows():
        bp, rp = r["best_panel_signature"], r["best_panel_r"]
        regs = [x for x in str(r["top_regulators"]).split(",") if x][:3]
        if bp and abs(rp) >= thr:
            lab = PN.label(bp) if rp > 0 else "low " + PN.label(bp)
        else:
            lab = "/".join(regs) + " program" if regs else f"program {k}"
        ev = (f"{PN.label(bp)} ({PN.tag(bp)}) r {rp:.2f}; " if bp else "") + \
             f"best overall {r['best_activity_signature']} r {r['best_activity_r']:.2f}; " + "/".join(regs)
        rows.append({"program": k, "label": lab, "evidence": ev, "current_label": cur.get(str(k), "")})
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(outdir, "program_labels_proposed.tsv"), sep="\t", index=False)
    changed = (t["current_label"] != "") & (t["label"] != t["current_label"])
    log.info("Proposed program labels (reference panel): %d programs; %d curated labels differ -> review "
             "program_labels_proposed.tsv and update config/program_labels.tsv", len(t), int(changed.sum()))


def state_table(D):
    S = D["samples"]
    rows = []
    for k, v in D["states"].items():
        v = [s for s in v if s in S.index]
        r = {"state": k, "n": len(v), "risk_mean_z": D["risk"].reindex(v).mean()}
        for c in ("TCGA", "CLCA", "LICA_FR"):
            r[f"frac_{c}"] = (S.loc[v, "cohort"] == c).mean()
        rows.append(r)
    t = pd.DataFrame(rows).set_index("state")
    ps = D["prog_states"]
    t["meta_z"] = ps["meta_z"].reindex(t.index)
    t["meta_q"] = ps["meta_q"].reindex(t.index)
    return t


def frac_table(states, labels, levels):
    out = pd.DataFrame(index=levels, columns=list(states), dtype=float)
    for k, v in states.items():
        x = labels.reindex(v).dropna()
        x = x[x != "unassigned"]
        for lv in levels:
            out.loc[lv, k] = (x == lv).mean() if len(x) else np.nan
    return out


def cluster(X):
    """Average-linkage clustering of the rows of X on correlation distance; returns (linkage, order)."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist
    A = np.nan_to_num(np.asarray(X, float))
    d = np.nan_to_num(pdist(A, "correlation"), nan=1.0)
    Z = linkage(d, "average")
    return Z, leaves_list(Z)


def draw_dendro(ax, Z, orientation):
    from scipy.cluster.hierarchy import dendrogram
    dendrogram(Z, ax=ax, orientation=orientation, no_labels=True, color_threshold=0,
               above_threshold_color=Q.MUTED, link_color_func=lambda k: Q.MUTED)
    ax.axis("off")


def fig_integrated(D, prog, st, outdir, R, written):
    """Programs x states, both clustered (average linkage, correlation distance), with dendrograms."""
    progs0, states0 = list(D["programs"]), list(D["states"])
    diff = D["diff"]
    M0 = pd.DataFrame({s: [diff.loc[D["programs"][k], [x for x in D["states"][s] if x in diff.columns]].values.mean()
                           for k in progs0] for s in states0}, index=progs0)
    Zr, ro = cluster(M0.values)
    Zc, co = cluster(M0.values.T)
    p_order = [progs0[i] for i in ro]
    s_order = [states0[i] for i in co]
    M = M0.loc[p_order, s_order]
    tracks = [("Cohort", pd.DataFrame({s: [st.loc[s, f"frac_{c}"] for c in ("TCGA", "CLCA", "LICA_FR")]
                                       for s in s_order}, index=["TCGA", "CLCA", "LICA-FR"]))]
    for key in PN.figures.get("f07c_tracks", []):
        if key not in D["ntp"]:
            continue
        lv = list(PN.classifiers[key]["levels"])
        T = frac_table(D["states"], D["ntp"][key], lv)[s_order]
        T.index = [PN.class_label(key, x) for x in lv]
        tracks.append((PN.class_title(key), T))
    G = D["genomic"]
    tracks.append(("Alteration frequency", pd.DataFrame(
        {s: [G.loc[m, [x for x in D["states"][s] if x in G.columns]].mean() for m in MUTS] for s in s_order},
        index=[m.replace("MUT_", "").replace("_", " ") for m in MUTS])))

    h_rows = [2.4, 2.2, 2.2] + [len(t[1]) for t in tracks] + [len(M)]
    ncol = len(s_order)
    cell = 0.15
    fig = plt.figure(figsize=(2.6 + 0.7 + 0.30 * ncol + 5.6, 2.0 + cell * (sum(h_rows) + 1.2 * len(h_rows))))
    gs = gridspec.GridSpec(len(h_rows) + 1, 4, height_ratios=h_rows + [3.0],
                           width_ratios=[2.6, 0.7, 0.30 * ncol, 5.6], hspace=0.28, wspace=0.02)
    HM = 2
    x = np.arange(ncol)

    def track_title(row, label):
        axt = fig.add_subplot(gs[row, 0])
        axt.axis("off")
        axt.text(0.0, 0.5, label, ha="left", va="center", fontsize=7.5, color=Q.INK2, transform=axt.transAxes)

    axd = fig.add_subplot(gs[0, HM])
    draw_dendro(axd, Zc, "top")
    axd.set_title("MINER programs × states, both clustered (average linkage, correlation distance)",
                  loc="left", fontsize=9.5)
    for row, vals, label, colored in ((1, st.loc[s_order, "n"], "samples per state", False),
                                      (2, st.loc[s_order, "risk_mean_z"], "mean risk score (z)", True)):
        ax = fig.add_subplot(gs[row, HM])
        cols = [Q.SLOTS[7] if v > 0 else Q.SLOTS[0] for v in vals] if colored else Q.MUTED
        ax.bar(x, vals, color=cols, width=0.7)
        if colored:
            ax.axhline(0, color=Q.AXIS, lw=0.8)
        ax.set_xlim(-0.5, ncol - 0.5)
        ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=6)
        track_title(row, label)
    for i, (title, T) in enumerate(tracks):
        ax = fig.add_subplot(gs[3 + i, HM])
        ax.imshow(T.values.astype(float), aspect="auto", cmap=Q.SEQ, vmin=0, vmax=1, interpolation="none")
        ax.set_yticks(range(len(T)))
        ax.set_yticklabels(T.index, fontsize=6.5)
        ax.set_xticks([])
        ax.grid(False)
        track_title(3 + i, title)
    last = len(h_rows) - 1
    ax = fig.add_subplot(gs[last, HM])
    vmax = np.nanpercentile(np.abs(M.values), 98)
    im = ax.imshow(M.values, aspect="auto", cmap=Q.DIV, vmin=-vmax, vmax=vmax, interpolation="none")
    ax.set_xticks(x)
    ax.set_xticklabels([f"S{s}" for s in s_order], fontsize=6.5, rotation=90)
    ax.set_yticks([])
    ax.grid(False)
    axd2 = fig.add_subplot(gs[last, 1])
    draw_dendro(axd2, Zr, "left")
    axd2.set_ylim(len(p_order) * 10, 0)   # leaves at 5, 15, 25, ... top to bottom, matching the imshow rows
    axl = fig.add_subplot(gs[last, 0])
    y = np.arange(len(p_order))
    w = prog.loc[p_order, "risk_weight"] / prog["risk_weight"].abs().max()
    mz = prog.loc[p_order, "meta_z"] / prog["meta_z"].abs().max()
    axl.barh(y - 0.2, w, height=0.4, color=[Q.SLOTS[7] if v > 0 else Q.SLOTS[0] for v in w])
    axl.barh(y + 0.2, mz, height=0.4, color=Q.BACKGROUND)
    axl.axvline(0, color=Q.AXIS, lw=0.8)
    axl.set_xlim(-1.1, 1.1)
    axl.set_ylim(len(p_order) - 0.5, -0.5)
    axl.set_xticks([-1, 0, 1])
    axl.set_xticklabels(["protective", "0", "adverse"], fontsize=6.5)
    axl.set_yticks(y)
    axl.set_yticklabels([f"P{k}" for k in p_order], fontsize=5.8)
    axl.grid(False)
    axl.set_xlabel("colored: risk-model weight\ngray: prognostic meta-z (RFS)\n(each scaled to max |value|)",
                   fontsize=6.5)
    axr = fig.add_subplot(gs[last, 3])
    axr.set_ylim(len(p_order) - 0.5, -0.5)
    axr.axis("off")
    for i, k in enumerate(p_order):
        r = prog.loc[k]
        sig = short(r["best_overlap_signature"]) if r["best_overlap_signature"] else short(r["best_activity_signature"])
        mark = {"prognostic, new": "★", "prognostic, known": "•"}.get(r["class"], " ")
        axr.text(0.01, i, f"{mark} P{k:<3} {r['top_regulators'][:32]:<32} {sig[:46]}", fontsize=5.6, va="center",
                 family="monospace", color=Q.INK if r["prognostic"] else Q.MUTED)
    axr.text(0.01, -1.6, "  prog top regulators                    best known signature", fontsize=6.2,
             family="monospace", color=Q.INK2, va="bottom")
    axc = fig.add_subplot(gs[-1, HM])
    axc.axis("off")
    cax = axc.inset_axes([0.0, 0.25, 0.45, 0.22])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("mean regulon dysregulation (over − under)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    cax2 = axc.inset_axes([0.55, 0.25, 0.45, 0.22])
    cb2 = fig.colorbar(plt.cm.ScalarMappable(cmap=Q.SEQ, norm=plt.Normalize(0, 1)), cax=cax2, orientation="horizontal")
    cb2.set_label("fraction of state (tracks above)", fontsize=7)
    cb2.ax.tick_params(labelsize=6)
    axn = fig.add_subplot(gs[-1, 3])
    axn.axis("off")
    axn.text(0.01, 0.8, "★ prognostic, not matched by known signatures (f2)\n• prognostic (meta q ≤ %.2f, same sign "
             "in TCGA and CLCA), matched\nblack = prognostic; gray = not prognostic" % R["q_max"],
             fontsize=6.8, va="top", color=Q.INK2)
    Q._save(fig, outdir, "f1_integrated_map.png", written)
    return M


def fig_program_signatures(D, prog, outdir, R, written):
    """Programs x curated signatures (activity r, global mean regressed out), both clustered."""
    cols, groups = [], []
    for s_ in PN.signatures(available=set(D["cor"].columns)):
        cols.append(s_["set"])
        groups.append(s_["group"])
    C0 = D["cor"].loc[list(D["programs"]), cols]
    Zr, ro = cluster(C0.values)
    Zc, co = cluster(C0.values.T)
    C = C0.iloc[ro, co]
    grp = [groups[i] for i in co]
    p_order = list(C.index)
    ov = D["overlap"]
    hit = ov[(ov["fdr"] < 0.05) & (ov["fold_enrichment"] >= R["min_fold"])]
    hit = set(zip(hit["program"], hit["signature"]))
    nr, nc = C.shape
    fig = plt.figure(figsize=(4.6 + 0.2 * nc, 3.4 + 0.15 * nr))
    gs = gridspec.GridSpec(3, 4, height_ratios=[2.0, 0.35, 0.15 * nr], width_ratios=[1.4, 0.2 * nc, 0.75, 0.3],
                           hspace=0.02, wspace=0.02)
    axd = fig.add_subplot(gs[0, 1])
    draw_dendro(axd, Zc, "top")
    axd.set_title("Program activity vs known signatures, clustered (r, per-sample global mean regressed out);\n"
                  "dot = direct gene overlap FDR < 0.05, ≥ %d-fold" % R["min_fold"], loc="left", fontsize=9)
    axg = fig.add_subplot(gs[1, 1])
    gnames = list(dict.fromkeys(groups))
    gcol = {g: Q.SLOTS[i] for i, g in enumerate(gnames)}
    for j, g in enumerate(grp):
        axg.add_patch(plt.Rectangle((j - 0.5, 0), 1, 1, color=gcol[g], lw=0))
    axg.set_xlim(-0.5, nc - 0.5)
    axg.set_ylim(0, 1)
    axg.axis("off")
    ax = fig.add_subplot(gs[2, 1])
    im = ax.imshow(C.values, aspect="auto", cmap=Q.DIV, vmin=-1, vmax=1, interpolation="none")
    for i, k in enumerate(p_order):
        for j, s_ in enumerate(C.columns):
            if (k, s_) in hit:
                ax.scatter(j, i, s=6, color=Q.INK, zorder=3)
    ax.yaxis.tick_right()
    ax.set_yticks(range(nr))
    ax.set_yticklabels([f"P{k}" + (" ↑" if prog.loc[k, "risk_weight"] > 0 else " ↓") for k in p_order], fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks(range(nc))
    ax.set_xticklabels([short(s_) for s_ in C.columns], fontsize=6, rotation=90)
    ax.grid(False)
    axl = fig.add_subplot(gs[2, 0])
    draw_dendro(axl, Zr, "left")
    axl.set_ylim(nr * 10, 0)
    cax = fig.add_subplot(gs[2, 3])
    cax.axis("off")
    cb = fig.colorbar(im, ax=cax, fraction=0.9, aspect=30)
    cb.set_label("Pearson r across samples")
    handles = [plt.Rectangle((0, 0), 1, 1, color=gcol[g]) for g in gnames]
    fig.legend(handles, gnames, loc="upper right", fontsize=7, title="signature group", title_fontsize=7,
               bbox_to_anchor=(0.99, 0.99))
    fig.text(0.01, 0.005, "program labels: ↑ adverse / ↓ protective weight in the risk model", fontsize=7, color=Q.INK2)
    Q._save(fig, outdir, "f4_program_signatures.png", written)


def guan_rank(srv):
    """MINER's GuanRank (miner.guanRank on miner.kmAnalysis), reimplemented without the miner package:
    per-sample risk rank from censored survival, scaled to max 1 (higher = earlier event = higher risk)."""
    from lifelines import KaplanMeierFitter
    srv = srv.sort_values("duration")
    km = KaplanMeierFitter().fit(srv["duration"], srv["observed"]).survival_function_.iloc[:, 0]
    times = km.index.values
    def s_at(t):
        if t in km.index:
            return km.loc[t]
        i = np.where(times < t)[0][-1]
        return 0.5 * (km.iloc[i] + km.iloc[i + 1])
    T = srv["duration"].values
    E = srv["observed"].values
    Sp = np.array([s_at(t) for t in T])
    score = np.zeros(len(T))
    for a in range(len(T)):
        tot = 0.0
        for b in range(len(T)):
            if a == b:
                continue
            if E[a] == 1:
                if T[b] > T[a]:
                    tot += 1
                if T[b] <= T[a] and E[b] == 0:
                    tot += Sp[a] / Sp[b]
                if T[b] == T[a] and E[b] == 1:
                    tot += 0.5
            else:
                if T[b] >= T[a]:
                    tot += (1 - 0.5 * Sp[b] / Sp[a]) if E[b] == 0 else (1 - Sp[b] / Sp[a])
                if T[b] < T[a] and E[b] == 0:
                    tot += 0.5 * Sp[a] / Sp[b]
        score[a] = tot
    return pd.Series(score / score.max(), index=srv.index)


def fig_state_survival(P, D, st, outdir, R, log, written):
    """States ordered by mean risk score: subtype enrichment (top), risk score and GuanRank per state."""
    res = p(P["paths"]["results"])
    S = D["samples"]
    s_order = st.sort_values("risk_mean_z").index.tolist()
    # GuanRank within each cohort (RFS, 36-month horizon), comparable across cohorts
    guan = pd.concat([guan_rank(s_) for s_ in load_surv(P, R).values()])
    # enrichment rows: fraction of the state, dot = one-sided Fisher FDR < 0.05 (07b; mutations here)
    enr = pd.read_csv(os.path.join(res, "07_post", "subtype_mapping", P["post"]["subtypes_dir"], "state_enrichment.tsv"),
                      sep="\t")
    enr["state"] = enr["state"].astype(str)
    blocks = []
    for key in PN.figures.get("f07c_tracks", []):
        spec = PN.classifiers[key]
        ann = spec.get("column", key) if spec.get("source") == "label" else f"ntp_{key}"
        blocks.append((PN.class_title(key), ann, list(spec["levels"])))
    blocks += [("Molecular group (LICA-FR)", "molecular_group", None),
              ("Immune class (LICA-FR)", "immune_class", ["hot", "cold"]),
              ("Cohort", "cohort", ["TCGA", "CLCA", "LICA_FR"])]
    rows, sig, labels, starts = [], [], [], []
    for title, ann, levels in blocks:
        e = enr[enr["annotation"] == ann]
        if levels is None:
            levels = sorted(e["level"].unique())
        starts.append((len(rows), title, len(levels)))
        for lv in levels:
            x = e[e["level"] == lv].set_index("state")
            rows.append([x["frac_of_state"].get(s_, np.nan) for s_ in s_order])
            sig.append([(x["fdr"].get(s_, 1) < 0.05) and (x["odds_ratio"].get(s_, 0) > 1) for s_ in s_order])
            labels.append(str(lv))
    # driver alterations: fraction among profiled, one-sided Fisher per state vs rest, BH within feature
    G = D["genomic"]
    starts.append((len(rows), "Alterations", len(MUTS)))
    for m in MUTS:
        g = G.loc[m].dropna()
        fr, pv = [], []
        for s_ in s_order:
            ins = g.index.isin(D["states"][s_])
            a_, b_ = int((g[ins] == 1).sum()), int((g[ins] == 0).sum())
            c_, d_ = int((g[~ins] == 1).sum()), int((g[~ins] == 0).sum())
            fr.append(a_ / (a_ + b_) if a_ + b_ else np.nan)
            pv.append(stats.fisher_exact([[a_, b_], [c_, d_]], alternative="greater")[1] if a_ + b_ else 1.0)
        rows.append(fr)
        sig.append(list(_bh(pv) < 0.05))
        labels.append(m.replace("MUT_", "").replace("_", " "))
    F = np.array(rows, float)
    Sg = np.array(sig, bool)

    nrow = len(rows)
    ncol = len(s_order)
    fig = plt.figure(figsize=(3.4 + 0.34 * ncol, 3.0 + 0.13 * nrow + 5.2))
    gs = gridspec.GridSpec(3, 2, height_ratios=[0.13 * nrow, 2.4, 2.8], width_ratios=[2.9, 0.34 * ncol],
                           hspace=0.08, wspace=0.02)
    x = np.arange(ncol)
    ax = fig.add_subplot(gs[0, 1])
    im = ax.imshow(F, aspect="auto", cmap=Q.SEQ, vmin=0, vmax=1, interpolation="none")
    yy, xx = np.where(Sg)
    ax.scatter(xx, yy, s=7, color=Q.INK, zorder=3)
    for st0, _, n in starts[1:]:
        ax.axhline(st0 - 0.5, color=Q.SURFACE, lw=2)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(labels, fontsize=6.3)
    ax.set_xticks([])
    ax.grid(False)
    ax.set_title("States ordered by mean program risk score: subtype enrichment, predicted risk and observed "
                 "outcome (RFS, 36 months)", loc="left", fontsize=9.5)
    axt = fig.add_subplot(gs[0, 0], sharey=ax)
    axt.axis("off")
    for st0, title, n in starts:
        axt.text(0.0, st0 + (n - 1) / 2, title, fontsize=7.2, va="center", ha="left", color=Q.INK2)
    colors = Q.cohort_colors(["TCGA", "CLCA", "LICA_FR"])
    rng = np.random.default_rng(0)

    def box_panel(axp, values, ylabel, cohorts):
        data = [values.reindex(D["states"][s_]).dropna() for s_ in s_order]
        bp = axp.boxplot(data, positions=x, widths=0.6, showfliers=False, patch_artist=True,
                         medianprops=dict(color=Q.INK, lw=1.6), whiskerprops=dict(color=Q.MUTED),
                         capprops=dict(color=Q.MUTED), boxprops=dict(facecolor=Q.SURFACE, edgecolor=Q.MUTED))
        for i, dd in enumerate(data):
            xs = i + rng.uniform(-0.22, 0.22, len(dd))
            axp.scatter(xs, dd, s=5, alpha=0.6, linewidth=0,
                        color=[colors.get(S.loc[s_, "cohort"], Q.MUTED) for s_ in dd.index], zorder=3)
        axp.set_xlim(-0.6, ncol - 0.4)
        axp.set_ylabel(ylabel, fontsize=8)
        axp.grid(axis="x", visible=False)
        return [len(dd) for dd in data], [dd.median() if len(dd) else np.nan for dd in data]

    ax.tick_params(labelbottom=False)
    ax1 = fig.add_subplot(gs[1, 1], sharex=ax)
    box_panel(ax1, D["risk"], "risk score\n(within-cohort z)", None)
    ax1.axhline(0, color=Q.AXIS, lw=0.8)
    ax1.tick_params(labelbottom=False)
    ax2 = fig.add_subplot(gs[2, 1], sharex=ax)
    n2, med2 = box_panel(ax2, guan, "GuanRank (RFS)\nwithin cohort; 1 = earliest event", None)
    ax2.set_xticks(x)
    ax2.set_xticklabels([f"S{s_}\n({n})" for s_, n in zip(s_order, n2)], fontsize=6.3)
    ax2.set_xlabel("state (samples with survival: TCGA + CLCA)", fontsize=8)
    ok = [i for i, n in enumerate(n2) if n >= 5]
    rho, pr = stats.spearmanr(st.loc[[s_order[i] for i in ok], "risk_mean_z"], [med2[i] for i in ok])
    ax2.text(0.01, 0.97, f"state mean risk vs median GuanRank (states with ≥ 5 samples): Spearman ρ = {rho:.2f}, "
             f"p = {pr:.1e}", transform=ax2.transAxes, fontsize=7, va="top", color=Q.INK2)
    handles = [plt.Line2D([], [], marker="o", ls="", color=colors[c], label=c) for c in ("TCGA", "CLCA", "LICA_FR")]
    ax1.legend(handles=handles, loc="upper left", fontsize=7, ncol=3)
    cax = fig.add_subplot(gs[1, 0])
    cax.axis("off")
    cb = fig.colorbar(im, ax=cax, fraction=0.5, location="left", aspect=12)
    cb.set_label("fraction of state\n(dot: enriched, FDR < 0.05)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    Q._save(fig, outdir, "f7_state_risk_survival.png", written)
    log.info("f7: state mean risk vs median GuanRank Spearman rho %.2f (p %.1e, %d states)", rho, pr, len(ok))


def fig_known_new(prog, outdir, R, written):
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    lim = max(3, np.nanmax(prog["meta_z"].abs()))
    sc = ax.scatter(prog["best_activity_r"].abs(), prog["best_overlap_neglog10_fdr"].clip(upper=60),
                    c=prog["meta_z"], cmap=Q.DIV, vmin=-lim, vmax=lim, s=12 + 1.2 * prog["n_regulons"],
                    edgecolor=Q.SURFACE, linewidth=1, zorder=3)
    ax.axvline(R["known_r"], color=Q.MUTED, ls="--", lw=1)
    ax.axhline(-np.log10(0.05), color=Q.MUTED, ls="--", lw=1)
    ax.text(0.02, 0.97, "lower left: no gene overlap and |r| < %.1f (not matched by known signatures)" % R["known_r"],
            transform=ax.transAxes, fontsize=7.5,
            color=Q.INK2, va="top")
    lab = prog[prog["prognostic"] & ((prog["class"] == "prognostic, new") | (prog["meta_z"].abs() >= 4.5))]
    lab = lab.assign(x=lab["best_activity_r"].abs()).sort_values("x")
    for i, (k, r) in enumerate(lab.iterrows()):
        dy = [10, -12, 22, -24][i % 4] if r["best_overlap_neglog10_fdr"] == 0 else 4
        ax.annotate(f"P{k}", (abs(r["best_activity_r"]), min(r["best_overlap_neglog10_fdr"], 60)), fontsize=7,
                    xytext=(3, dy), textcoords="offset points", color=Q.INK,
                    arrowprops=dict(arrowstyle="-", color=Q.MUTED, lw=0.5) if dy != 4 else None)
    ax.set_xlabel("best activity correlation with a known signature (|r|, per-sample global mean regressed out)")
    ax.set_ylabel("best gene-set overlap, −log10 FDR (capped at 60)")
    ax.set_title("Programs: similarity to known signatures (size = regulons, colour = prognostic meta-z, RFS)\n"
                 "similarity only; for prognosis beyond known subtypes see f6", fontsize=9.5)
    cb = fig.colorbar(sc, ax=ax, fraction=0.04)
    cb.set_label("meta-z (> 0 adverse)")
    fig.tight_layout()
    Q._save(fig, outdir, "f2_known_vs_new.png", written)


def fig_risk_by_subtype(D, outdir, written):
    S = D["samples"]
    risk = D["risk"]
    panels = [(PN.class_title(k), D["ntp"][k], list(PN.classifiers[k]["levels"]))
              for k in PN.figures.get("f07c_tracks", []) if k in D["ntp"]]
    panels += [("Molecular group (LICA-FR authors)", S["molecular_group"], None),
               ("Immune class (LICA-FR authors)", S["immune_class"], ["hot", "cold"])]
    colors = Q.cohort_colors(["TCGA", "CLCA", "LICA_FR"])
    nr = int(np.ceil(len(panels) / 3))
    fig, axes = plt.subplots(nr, 3, figsize=(15, 4.1 * nr))
    for ax in axes.ravel()[len(panels):]:
        ax.set_visible(False)
    rng = np.random.default_rng(0)
    for ax, (title, lab, levels) in zip(axes.ravel(), panels):
        lab = lab.reindex(risk.index)
        if levels is None:
            med = risk.groupby(lab).median().sort_values()
            levels = [lv for lv in med.index if (lab == lv).sum() >= 5]
        groups = [risk[lab == lv].dropna() for lv in levels]
        for i, g in enumerate(groups):
            xs = i + rng.uniform(-0.28, 0.28, len(g))
            ax.scatter(xs, g, s=8, color=[colors.get(S.loc[s, "cohort"], Q.MUTED) for s in g.index], alpha=0.7,
                       linewidth=0, zorder=2)
            if len(g):
                ax.plot([i - 0.34, i + 0.34], [g.median()] * 2, color=Q.INK, lw=2, zorder=3)
        ok = [g for g in groups if len(g) >= 3]
        pk = stats.kruskal(*ok).pvalue if len(ok) >= 2 else np.nan
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels([f"{lv}\n(n={len(g)})" for lv, g in zip(levels, groups)], fontsize=6.5,
                           rotation=30 if len(levels) > 6 else 0)
        ax.axhline(0, color=Q.AXIS, lw=0.8)
        ax.set_title(f"{title}   Kruskal p = {pk:.1e}", fontsize=9)
        ax.set_ylabel("risk score (within-cohort z)")
        ax.grid(axis="x", visible=False)
    handles = [plt.Line2D([], [], marker="o", ls="", color=colors[c], label=c) for c in ("TCGA", "CLCA", "LICA_FR")]
    handles.append(plt.Line2D([], [], color=Q.INK, lw=2, label="median"))
    fig.legend(handles=handles, loc="upper right", ncol=4, fontsize=8)
    fig.suptitle("Program-based risk score across published HCC subtypes (model trained in TCGA, RFS)", x=0.01,
                 ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    Q._save(fig, outdir, "f3_risk_by_subtype.png", written)


def fig_state_risk(D, st, outdir, written):
    ps = D["prog_states"].copy()
    ps = ps.loc[st.index]
    order = ps.sort_values("meta_z").index
    zc = [c for c in ps.columns if c.startswith("z_")]
    colors = Q.cohort_colors([c[2:] for c in zc])
    lab = {}
    for k in order:
        v = D["states"][k]
        parts = []
        for cl in [k for k, _ in PN.figures.get("fig1c_tracks", []) if k in D["ntp"]]:
            c = D["ntp"][cl].reindex(v)
            c = c[c != "unassigned"].dropna()
            if len(c):
                top = c.value_counts()
                parts.append(f"{top.index[0]} {top.iloc[0] / len(c):.0%}")
        dom = max(("TCGA", "CLCA", "LICA_FR"), key=lambda c: st.loc[k, f"frac_{c}"])
        lab[k] = f"S{k} (n={st.loc[k, 'n']}; {dom} {st.loc[k, f'frac_{dom}']:.0%})  " + " · ".join(parts)
    fig, ax = plt.subplots(figsize=(8.8, 0.6 + 0.24 * len(order)))
    y = np.arange(len(order))
    for j, c in enumerate(zc):
        ax.scatter(ps.loc[order, c], y + (j - (len(zc) - 1) / 2) * 0.22, s=16, color=colors[c[2:]], label=c[2:],
                   zorder=3)
    sig = (ps.loc[order, "meta_q"] <= 0.1) & ps.loc[order, "consistent"]
    ax.scatter(ps.loc[order, "meta_z"], y, marker="|", s=140, color=Q.INK, label="meta-z", zorder=4)
    ax.axvline(0, color=Q.AXIS, lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([("★ " if s else "") + lab[k] for k, s in zip(order, sig)], fontsize=6.5)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Cox z for state membership (> 0: worse RFS); ★ meta q ≤ 0.1, same sign in both cohorts")
    ax.legend(loc="lower right", fontsize=7)
    ax.set_title("State-level risk (RFS, 36 months) with dominant NTP subtype calls", loc="left")
    fig.tight_layout()
    Q._save(fig, outdir, "f5_state_risk.png", written)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    matrix = args.matrix or P["miner"]["matrix"]
    R = P["post"]["figures"]
    outdir = p(os.path.join(P["paths"]["results"], OUT))
    log = setup_logging(outdir, "07c_figures")
    PN.log = log
    D = load(P, matrix, log)
    prog = annotate_programs(D, R)
    st = state_table(D)
    prog.to_csv(os.path.join(outdir, "program_annotation.tsv"), sep="\t", float_format="%.4g")
    st.to_csv(os.path.join(outdir, "state_annotation_risk.tsv"), sep="\t", float_format="%.4g")
    log.info("Program classes: %s", prog["class"].value_counts().to_dict())
    log.info("Prognostic programs not explained by known signatures:\n%s",
             prog[prog["class"] == "prognostic, new"][["n_regulons", "risk_weight", "meta_z", "top_regulators",
                                                        "best_activity_signature", "best_activity_r"]].to_string())
    # prognostic value beyond known subtypes
    KS = known_scores(P, D, R)
    R["adjust_labels"] = [PN.label(c) for c in KS.columns]
    surv = load_surv(P, R)
    adj = cox_adjusted(D["act"], surv, KS, R)
    adj.index.name = "program"
    adj.to_csv(os.path.join(outdir, "program_prognosis_beyond_known.tsv"), sep="\t", float_format="%.4g")
    from lifelines import CoxPHFitter
    rows = {}
    for c, srv in surv.items():
        ids = srv.index.intersection(KS.index)
        cv = KS.loc[ids]
        cv = (cv - cv.mean()) / cv.std()
        rk = D["risk"].reindex(ids)
        d = srv.loc[ids].assign(x=(rk - rk.mean()) / rk.std()).join(cv)
        r_ = CoxPHFitter().fit(d, "duration", "observed").summary.loc["x"]
        rows[c] = {"hr": r_["exp(coef)"], "lo": r_["exp(coef) lower 95%"], "hi": r_["exp(coef) upper 95%"], "p": r_["p"],
                   "in_sample": c == R["risk_model"].split("_")[2]}
    risk_adj = pd.DataFrame(rows).T
    risk_adj.to_csv(os.path.join(outdir, "risk_score_beyond_known.tsv"), sep="\t", float_format="%.4g")
    new = adj[(adj["meta_q_adjusted"] <= R["q_max"]) & adj["consistent_adjusted"]]
    prog = prog.join(adj[["meta_z_adjusted", "meta_q_adjusted", "consistent_adjusted"]])
    prog["beyond_known"] = prog.index.isin(new.index)
    prog.to_csv(os.path.join(outdir, "program_annotation.tsv"), sep="\t", float_format="%.4g")
    log.info("Adjusting for: %s", R["adjust_labels"])
    log.info("Programs prognostic unadjusted: %d; still prognostic after adjustment: %d -> %s",
             int(((adj["meta_q_unadjusted"] <= R["q_max"]) & adj["consistent_unadjusted"]).sum()), len(new),
             prog.loc[prog["beyond_known"], ["meta_z", "meta_z_adjusted", "top_regulators", "best_activity_signature",
                                             "best_activity_r"]].to_string())
    log.info("Risk score adjusted for known subtype scores:\n%s", risk_adj.to_string())
    h2h = head_to_head(P, KS, R)
    h2h.to_csv(os.path.join(outdir, "risk_head_to_head.tsv"), sep="\t", index=False, float_format="%.4g")
    log.info("Head-to-head (external C-index):\n%s", h2h.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    written = []
    fig_beyond_known(prog, adj, risk_adj, h2h, outdir, R, written)
    fig_integrated(D, prog, st, outdir, R, written)
    fig_known_new(prog, outdir, R, written)
    fig_risk_by_subtype(D, outdir, written)
    fig_program_signatures(D, prog, outdir, R, written)
    fig_state_risk(D, st, outdir, written)
    fig_state_survival(P, D, st, outdir, R, log, written)
    proposed_labels(prog, outdir, R, log)
    log.info("Figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
