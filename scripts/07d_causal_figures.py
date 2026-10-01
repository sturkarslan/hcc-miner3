#!/usr/bin/env python
"""Step 07d: causal-flow figures (genomic driver -> regulator -> regulon/program), from step 05's
high-confidence flows, linked to the step-06 risk programs and the step-07 program annotation.

f8_causal_driver_program   a. drivers x programs: dot size = regulon families with a high-confidence flow,
                              colour = mean Cohen's d (up/down in altered tumours); programs ordered by risk
                              weight (protective -> adverse), weight bar on top.
                           b. per driver, "net risk push" = sum over its causal families of sign(d) x risk
                              weight of the family's program.
                           c. net risk push vs the driver's observed association with recurrence (Cox z per
                              cohort on profiled samples, Stouffer meta-z over TCGA and CLCA).
f9_causal_flows            one panel per selected driver: driver -> top regulators (best flow per regulon
                           family, ranked by |d|) -> programs. Edge colour: red = up / activates, blue = down /
                           represses; width ~ |d|. Program boxes: border = risk-model direction; label = best
                           HCC / hallmark signature by activity r (hepatoblastoma sets excluded).

A family is a group of regulons with near-identical activity (step 05), so each counts once.
Inputs: results/05_causal/<matrix>/highConfidenceCausalResults.csv, results/06_risk/<matrix>/,
        results/07_post/figures/program_annotation.tsv, results/03_genomics_clinical/
Outputs: results/07_post/figures/f8_*.png, f9_*.png, causal_driver_summary.tsv, causal_flow_edges.tsv
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats

from hcc_common import load_params, p, setup_logging

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.patches import FancyBboxPatch, PathPatch  # noqa: E402

import qc_plots as Q  # noqa: E402

OUT = os.path.join("07_post", "figures")


def label(f):
    return (f.replace("MUT_", "").replace("PATH_", "pathway: ").replace("AMP_", "amp ").replace("HD_", "del ")
            .replace("ARM_", "").replace("_", " "))


def short_sig(s):
    if not isinstance(s, str) or not s:
        return ""
    s = s.replace("HALLMARK_", "H: ").replace("_LIVER_CANCER", "").replace("_HEPATOCELLULAR_CARCINOMA", "")
    s = s.replace("_SUBCLASS", "").replace("LIVER_CANCER_", "").replace("_", " ")
    return s.title().replace("Ctnnb1", "CTNNB1").replace("Epcam", "EPCAM").replace("Krt19", "KRT19") \
        .replace("Myc", "MYC").replace("Mtorc1", "mTORC1").replace(" Up", " up").replace(" Dn", " dn")


def driver_survival(F, surv, cohort):
    """Cox z of each binary feature within each cohort (profiled samples), Stouffer meta-z."""
    from lifelines import CoxPHFitter
    rows = {}
    for f in F.index:
        zs, ws = [], []
        for c, s in surv.items():
            x = F.loc[f].reindex(s.index).dropna()
            if x.sum() < 5 or (x == 0).sum() < 5:
                continue
            d = s.loc[x.index].assign(x=x.values)
            try:
                z = CoxPHFitter().fit(d, "duration", "observed").summary.loc["x", "z"]
            except Exception:  # noqa: BLE001
                continue
            zs.append(z)
            ws.append(np.sqrt(d["observed"].sum()))
        if zs:
            ws = np.array(ws)
            rows[f] = {"surv_meta_z": float((np.array(zs) * ws).sum() / np.sqrt((ws ** 2).sum())),
                       "surv_cohorts": len(zs)}
    return pd.DataFrame(rows).T


def fig_driver_program(hc, fams, prog, drv, outdir, R, written):
    feats = drv.sort_values("net_risk_push").index.tolist()
    progs = prog.sort_values("risk_weight").index.tolist()
    progs = [k for k in progs if k in set(fams["program"])]
    n = fams.groupby(["Mutation", "program"]).agg(n=("family", "size"), d=("cohen_d", "mean")).reset_index()
    fig = plt.figure(figsize=(3.2 + 0.19 * len(progs) + 6.0, 1.8 + 0.28 * len(feats)))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1.0, 0.28 * len(feats)], width_ratios=[0.19 * len(progs), 2.0, 3.6],
                           hspace=0.05, wspace=0.25)
    # a. weight bar
    axw = fig.add_subplot(gs[0, 0])
    w = prog.loc[progs, "risk_weight"]
    axw.bar(range(len(progs)), w, color=[Q.SLOTS[7] if v > 0 else Q.SLOTS[0] for v in w], width=0.75)
    axw.axhline(0, color=Q.AXIS, lw=0.8)
    axw.set_xlim(-0.6, len(progs) - 0.4)
    axw.tick_params(axis="x", labelbottom=False, length=0)
    axw.set_ylabel("risk\nweight", fontsize=7)
    axw.tick_params(axis="y", labelsize=6)
    axw.set_title("a. Causal effects of drivers on programs (high-confidence flows; one dot per driver × program)",
                  loc="left", fontsize=9)
    ax = fig.add_subplot(gs[1, 0], sharex=axw)
    xi = {k: i for i, k in enumerate(progs)}
    yi = {f: i for i, f in enumerate(feats)}
    m = n[n["Mutation"].isin(feats) & n["program"].isin(progs)]
    sc = ax.scatter([xi[k] for k in m["program"]], [yi[f] for f in m["Mutation"]], s=10 + 9 * m["n"], c=m["d"],
                    cmap=Q.DIV, vmin=-1.2, vmax=1.2, edgecolor=Q.SURFACE, linewidth=0.6, zorder=3)
    ax.set_yticks(range(len(feats)))
    ax.set_yticklabels([label(f) for f in feats], fontsize=7)
    ax.set_xticks(range(len(progs)))
    ax.set_xticklabels([f"P{k}" for k in progs], fontsize=5.8, rotation=90)
    ax.set_ylim(-0.7, len(feats) - 0.3)
    ax.set_xlabel("program (ordered by risk-model weight: protective → adverse)", fontsize=7.5)
    ax.grid(True, color=Q.GRID, lw=0.4)
    cax = ax.inset_axes([1.01, 0.0, 0.012, 0.35])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("mean Cohen's d\n(> 0 up in altered)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    for nn in (1, 5, 15):
        ax.scatter([], [], s=10 + 9 * nn, color=Q.MUTED, label=f"{nn} families")
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=6.5, frameon=False, labelspacing=1.2)
    # b. net push
    axb = fig.add_subplot(gs[1, 1], sharey=ax)
    v = drv.loc[feats, "net_risk_push"]
    axb.barh(range(len(feats)), v, color=[Q.SLOTS[7] if x > 0 else Q.SLOTS[0] for x in v], height=0.6)
    axb.axvline(0, color=Q.AXIS, lw=0.8)
    axb.tick_params(axis="y", labelleft=False)
    axb.set_xlabel("net risk push\nΣ sign(d) × program weight", fontsize=7.5)
    axb.grid(axis="y", visible=False)
    axb.set_title("b.", loc="left", fontsize=9)
    # c. push vs observed survival association
    axc = fig.add_subplot(gs[:, 2])
    d = drv.dropna(subset=["surv_meta_z"])
    axc.axhline(0, color=Q.AXIS, lw=0.8)
    axc.axvline(0, color=Q.AXIS, lw=0.8)
    axc.scatter(d["net_risk_push"], d["surv_meta_z"], s=16 + 0.15 * d["n_altered"],
                color=[Q.SLOTS[1] if t == "arm_cna" else Q.SLOTS[6] for t in d["type"]],
                edgecolor=Q.SURFACE, linewidth=0.8, zorder=3)
    for f, r in d.iterrows():
        if r["type"] != "arm_cna" or abs(r["surv_meta_z"]) > 2 or abs(r["net_risk_push"]) > 0.3:
            axc.annotate(label(f), (r["net_risk_push"], r["surv_meta_z"]), fontsize=6.5, xytext=(3, 2),
                         textcoords="offset points", color=Q.INK)
    rho, pv = stats.spearmanr(d["net_risk_push"], d["surv_meta_z"])
    axc.set_xlabel("net risk push from causal flows (panel b)")
    axc.set_ylabel("observed association with recurrence\n(Cox meta-z, TCGA + CLCA; > 0 worse)")
    axc.set_title(f"c. Does the causal layer explain driver prognosis?  Spearman ρ = {rho:.2f} (p = {pv:.1e}, "
                  f"{len(d)} drivers)", loc="left", fontsize=9)
    axc.scatter([], [], color=Q.SLOTS[6], label="mutation / focal CNA / pathway")
    axc.scatter([], [], color=Q.SLOTS[1], label="arm-level CNA")
    axc.legend(loc="upper left", fontsize=7)
    Q._save(fig, outdir, "f8_causal_driver_program.png", written)
    return rho, pv


def curve(ax, x0, y0, x1, y1, color, lw, ls="-", alpha=0.85):
    xm = (x0 + x1) / 2
    path = Path([(x0, y0), (xm, y0), (xm, y1), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, ls=ls, alpha=alpha, zorder=1))


def fig_flows(fams, prog, prog_label, drivers, outdir, R, written, log):
    n_top = R["n_regulators"]
    ncols = 2
    nrows = int(np.ceil(len(drivers) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 5.4 * nrows))
    axes = np.atleast_1d(axes).ravel()
    edges = []
    for ax, f in zip(axes, drivers):
        x = fams[fams["Mutation"] == f].copy()
        if x.empty:
            ax.axis("off")
            continue
        x = x.assign(ad=x["cohen_d"].abs()).sort_values("ad", ascending=False).head(n_top)
        # regulators (one row per family's best flow) and the programs they reach
        regs = x.drop_duplicates("regulator_symbol")
        pg = x.groupby("program").agg(n=("family", "size"), d=("cohen_d", "mean")).reset_index()
        pg = pg.assign(w=pg["program"].map(prog["risk_weight"])).sort_values("w")
        yr = np.linspace(0.95, 0.05, len(regs)) if len(regs) > 1 else np.array([0.5])
        yp = np.linspace(0.9, 0.1, len(pg)) if len(pg) > 1 else np.array([0.5])
        ry = dict(zip(regs["regulator_symbol"], yr))
        py = dict(zip(pg["program"], yp))
        X0, X1, X2 = 0.08, 0.42, 0.66
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        for _, r in x.iterrows():
            ed = r["MutationRegulatorEdge"]
            curve(ax, X0 + 0.07, 0.5, X1 - 0.055, ry[r["regulator_symbol"]], Q.SLOTS[7] if ed > 0 else Q.SLOTS[0],
                  0.6 + 1.6 * min(r["ad"], 2) / 2)
            act = r["RegulatorRegulon_Spearman_R"] > 0
            curve(ax, X1 + 0.055, ry[r["regulator_symbol"]], X2 - 0.005, py[r["program"]],
                  Q.SLOTS[7] if r["cohen_d"] > 0 else Q.SLOTS[0], 0.6 + 2.4 * min(r["ad"], 2) / 2,
                  ls="-" if act else (0, (3, 2)))
            edges.append({"driver": f, "regulator": r["regulator_symbol"], "regulator_up_in_altered": ed > 0,
                          "regulator_activates_regulon": act, "regulon": r["Regulon"], "family": r["family"],
                          "program": r["program"], "cohen_d": r["cohen_d"]})
        # nodes
        ax.add_patch(FancyBboxPatch((X0 - 0.07, 0.44), 0.14, 0.12, boxstyle="round,pad=0.01", fc=Q.INK, ec=Q.INK))
        ax.text(X0, 0.5, label(f).replace("amp ", "amp\n").replace("pathway: ", "pathway:\n"), color=Q.SURFACE,
                ha="center", va="center", fontsize=8.5, fontweight="bold")
        ax.text(X0, 0.40, f"{int(x['n_altered'].iloc[0])} altered", ha="center", va="top", fontsize=6.5, color=Q.INK2)
        for _, r in regs.iterrows():
            up = r["MutationRegulatorEdge"] > 0
            ax.text(X1, ry[r["regulator_symbol"]], f"{r['regulator_symbol']} {'↑' if up else '↓'}", ha="center",
                    va="center", fontsize=7.2, color=Q.INK,
                    bbox=dict(boxstyle="round,pad=0.25", fc=Q.SURFACE, ec=Q.MUTED, lw=0.6))
        for _, r in pg.iterrows():
            k = r["program"]
            ann = prog.loc[k]
            sig = prog_label.get(k, "")
            adverse = ann["risk_weight"] > 0
            ax.text(X2, py[k], f"P{k} {'↑' if r['d'] > 0 else '↓'}  ({int(r['n'])})", ha="left", va="center",
                    fontsize=7.5, color=Q.INK, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.3", fc=Q.SURFACE, ec=Q.SLOTS[7] if adverse else Q.SLOTS[0], lw=1.6))
            ax.text(X2 + 0.105, py[k], sig[:44], ha="left", va="center", fontsize=6.5, color=Q.INK2)
        net = np.sum(np.sign(x["cohen_d"]) * x["program"].map(prog["risk_weight"]))
        ax.set_title(f"{label(f)}: top {len(x)} regulon families by |d| → {len(pg)} programs", loc="left", fontsize=9.5)
    for ax in axes[len(drivers):]:
        ax.axis("off")
    fig.text(0.01, 0.005, "Regulator ↑/↓: regulator up/down in altered tumours (red/blue edge from driver). Edge to program: "
             "red/blue = regulon up/down in altered tumours, solid = regulator activates, dashed = represses; width ~ |d|. "
             "Program box border: red = adverse, blue = protective weight in the risk model; (n) = families.",
             fontsize=7, color=Q.INK2)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    Q._save(fig, outdir, "f9_causal_flows.png", written)
    return pd.DataFrame(edges)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--matrix", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    matrix = args.matrix or P["miner"]["matrix"]
    R = P["post"]["causal_figures"]
    res = p(P["paths"]["results"])
    outdir = os.path.join(res, OUT)
    log = setup_logging(outdir, "07d_causal_figures")

    hc = pd.read_csv(os.path.join(res, "05_causal", matrix, "highConfidenceCausalResults.csv"), index_col=0)
    hc["program"] = hc["program"].astype(int).astype(str)
    prog = pd.read_csv(os.path.join(outdir, "program_annotation.tsv"), sep="\t", index_col=0)
    prog.index = prog.index.astype(str)
    # one row per driver x family: the flow with the largest |d| (families count once)
    fams = (hc.assign(ad=hc["cohen_d"].abs()).sort_values("ad", ascending=False)
            .drop_duplicates(["Mutation", "family"]).drop(columns="ad"))
    info = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features_info.tsv"), sep="\t", index_col=0)
    # drivers shown: all non-arm features with high-confidence flows + the most frequent arm events
    nonarm = [f for f in fams["Mutation"].unique() if info.loc[f, "type"] != "arm_cna"]
    arms = [a for a in R["arm_events"] if a in set(fams["Mutation"])]
    feats = nonarm + arms
    drv = pd.DataFrame(index=feats)
    drv["type"] = info.loc[feats, "type"]
    drv["n_altered"] = info.loc[feats, "altered"]
    x = fams[fams["Mutation"].isin(feats)]
    drv["families"] = x.groupby("Mutation")["family"].nunique()
    drv["net_risk_push"] = x.assign(v=np.sign(x["cohen_d"]) * x["program"].map(prog["risk_weight"])).groupby("Mutation")["v"].sum()
    F = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
    surv = {}
    for c in ("TCGA", "CLCA"):
        s = pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_{c}_{R['endpoint']}_h36m_miner.csv"),
                        index_col=0)
        s.columns = ["duration", "observed"]
        surv[c] = s
    drv = drv.join(driver_survival(F.loc[feats], surv, None))
    drv.index.name = "feature"
    drv.to_csv(os.path.join(outdir, "causal_driver_summary.tsv"), sep="\t", float_format="%.4g")
    log.info("Drivers:\n%s", drv.sort_values("net_risk_push").to_string(float_format=lambda v: f"{v:.3f}"))

    written = []
    rho, pv = fig_driver_program(hc, fams[fams["Mutation"].isin(feats)], prog, drv, outdir, R, written)
    log.info("Net risk push vs observed recurrence association: Spearman rho %.2f, p %.1e", rho, pv)
    # program labels: signature by activity r (07b, global mean regressed
    # out); hepatoblastoma sets are excluded because they mostly restate the proliferation axis
    cor = pd.read_csv(os.path.join(res, "07_post", "subtype_mapping", P["post"]["subtypes_dir"],
                                   "program_signature_correlation.tsv"), sep="\t", index_col=0)
    cor.index = cor.index.astype(str)
    cor = cor[[c for c in cor.columns if "HEPATOBLAST" not in c]]
    # prefer the reference panel (config/reference_panel.yaml; recent classifications first) when it matches
    # with |r| >= 0.5, otherwise the best signature overall
    from hcc_panel import Panel
    PN = Panel(log=log)
    pset = [x["set"] for x in PN.signatures(available=set(cor.columns))]
    prog_label = {}
    for k in cor.index:
        bp = cor.loc[k, pset].abs().idxmax() if pset else None
        if bp is not None and abs(cor.loc[k, bp]) >= 0.5:
            prog_label[k] = f"{PN.label(bp)} (r {cor.loc[k, bp]:.2f})"
        else:
            b_ = cor.loc[k].abs().idxmax()
            prog_label[k] = f"{short_sig(b_)} (r {cor.loc[k, b_]:.2f})"
    edges = fig_flows(fams, prog, prog_label, [d for d in R["flow_drivers"] if d in set(fams["Mutation"])], outdir,
                      R, written, log)
    edges.to_csv(os.path.join(outdir, "causal_flow_edges.tsv"), sep="\t", index=False, float_format="%.4g")
    log.info("Figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
