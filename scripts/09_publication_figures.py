#!/usr/bin/env python
"""Step 09: two publication figures (Nature double column, 183 mm), assembled from saved results only.

Figure 1  Data, network and biology
  a study design (cohorts, data types, workflow, external cohorts)   b network by the numbers
  c programs x states (clustered) with published-class biology, drivers and risk
  d what the programs are: program activity vs curated signatures (biology labels)
  e CTNNB1 causal flow (driver -> regulators -> programs)             f causal net risk push vs observed prognosis
Figure 2  Risk and validation
  a risk-model program weights         b risk score across published classes (biology labels)
  c states ordered by risk: class composition, risk score, observed GuanRank
  d cross-cohort KM (discovery)         e external KM
  f HR per SD, all cohorts + pooled    g C-index vs known signatures     h leave-one-cohort-out (placeholder)

Labels: programs from config/program_labels.tsv (curated, with evidence); published classes, signatures, biology
blocks and the known-signature benchmark from the reference panel (config/reference_panel.yaml: recent
classifications first - Montironi 2023, Haber 2023, Zhu 2022, Gao 2019, Sia 2017, Désert 2017 - with the
2007-2009 classes as references), shown by their biology (class code / source as a small grey tag). Outputs: results/09_figures/figure{1,2}.{pdf,png}, source tables.
"""

import argparse
import importlib.util
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from hcc_common import load_params, miner_id_backmap, p, setup_logging
from hcc_panel import Panel

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import gridspec  # noqa: E402
from matplotlib.patches import FancyBboxPatch, PathPatch  # noqa: E402
from matplotlib.path import Path  # noqa: E402

import qc_plots as Q  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
MM = 1 / 25.4
W = 183 * MM

# ---------------------------------------------------------------- style (Nature)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6, "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
    "legend.fontsize": 5.5, "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 0.9, "axes.titleweight": "normal",
    "axes.titlelocation": "left", "axes.titlepad": 3, "axes.labelpad": 2, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 300, "axes.grid": False, "legend.frameon": False, "legend.handlelength": 1.2,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
})
ADV, PROT = Q.SLOTS[7], Q.SLOTS[0]          # adverse red / protective blue, everywhere
COH = Q.cohort_colors(["TCGA", "CLCA", "LICA_FR"])
EXT_COL = {"GSE14520": Q.SLOTS[3], "LIRI_JP": Q.SLOTS[4], "GSE76427": Q.SLOTS[6]}
PN = Panel()   # reference signature panel (config/reference_panel.yaml)
_BCOL = {"adverse": ADV, "protective": PROT, "immune": Q.SLOTS[2]}
BLOCK_COL = {b: _BCOL.get(spec.get("colour"), Q.MUTED) for b, spec in PN.blocks.items()}
AXIS_COL = {"proliferation": ADV, "wnt": PROT, "immune": Q.SLOTS[2], "other": Q.MUTED}
MODEL_COL = {"background": Q.BACKGROUND, "muted": Q.MUTED, "slot2": Q.SLOTS[2], "slot3": Q.SLOTS[3],
             "slot4": Q.SLOTS[4], "slot5": Q.SLOTS[5], "slot6": Q.SLOTS[6]}


# Red-blue is reserved for risk and regulon dysregulation; fractions use viridis and signed
# correlations a purple-green diverging scale.
FRAC_CMAP = "viridis"
COR_CMAP = "PRGn"  # signed correlations: purple negative, white zero, green positive


def class_tag(key, level):
    """Small grey tag next to a class label: the class code, or the source when the code is the label."""
    lab = PN.class_label(key, level)
    return str(level) if str(level).lower() not in lab.lower() else PN.class_title(key).split(" (")[-1].rstrip(")")


def _module(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def letter(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="right")


def despine(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def km(ax, dur_days, ev, color, label, lw=0.9):
    ts, ss, t, e = Q.kaplan_meier(dur_days / 30.44, ev)
    ax.step(np.append(ts, t.max()), np.append(ss, ss[-1]), where="post", color=color, lw=lw, label=label)


# ================================================================ data
def load_all(P, log):
    res = p(P["paths"]["results"])
    mx = P["miner"]["matrix"]
    D = {"res": res}
    mdir = os.path.join(res, "04_miner", mx)
    D["samples"] = pd.read_csv(os.path.join(res, "01_harmonized", "samples.tsv"), sep="\t", index_col="sample")
    D["genes"] = pd.read_csv(os.path.join(res, "01_harmonized", "genes.tsv"), sep="\t", index_col="ensembl")
    D["programs"] = {k: [str(r) for r in v] for k, v in json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_programs.json"))).items()}
    D["states"] = json.load(open(os.path.join(mdir, "subtypes_filtered", "transcriptional_states.json")))
    over = pd.read_csv(os.path.join(mdir, "subtypes_filtered", "overExpressedMembers.csv"), index_col=0)
    under = pd.read_csv(os.path.join(mdir, "subtypes_filtered", "underExpressedMembers.csv"), index_col=0)
    over.index, under.index = over.index.astype(str), under.index.astype(str)
    D["diff"] = over - under
    D["cor"] = pd.read_csv(os.path.join(res, "07_post", "subtype_mapping", "subtypes_filtered",
                                        "program_signature_correlation.tsv"), sep="\t", index_col=0).rename(index=str)
    D["ntp"] = PN.class_calls(os.path.join(res, "07_post", "subtype_mapping", "subtypes_filtered"), D["samples"])
    D["genomic"] = pd.read_csv(os.path.join(res, "03_genomics_clinical", "genomic_features.csv"), index_col=0)
    rdir = os.path.join(res, "06_risk", mx)
    D["weights"] = pd.read_csv(os.path.join(rdir, "predictor_ridge_programs_TCGA_RFS_h36m", "weights.tsv"), sep="\t",
                               index_col=0)["weight"].rename(index=str)
    pred = pd.read_csv(os.path.join(rdir, "predictor_ridge_programs_TCGA_RFS_h36m", "predictions.tsv"), sep="\t", index_col=0)
    D["risk"] = pred.groupby("cohort")["risk_score"].transform(lambda v: (v - v.mean()) / v.std())
    D["prog_meta"] = pd.read_csv(os.path.join(rdir, "prognostic_programs_h36m_RFS.tsv"), sep="\t", index_col=0).rename(index=str)
    D["surv"] = {(c, ep): pd.read_csv(os.path.join(res, "03_genomics_clinical", f"survival_{c}_{ep}_h36m_miner.csv"),
                                      index_col=0).set_axis(["duration", "observed"], axis=1)
                 for c in ("TCGA", "CLCA") for ep in ("RFS", "OS")}
    D["rdir"] = rdir
    lab = pd.read_csv(p("config/program_labels.tsv"), sep="\t", comment="#", dtype={"program": str}).set_index("program")
    D["plabel"] = lab["label"].to_dict()
    # biology block per program: anchor signature group with the highest mean activity r (>= 0.3), else other
    D["block"] = pd.Series({k: PN.block_of_program(D["cor"].loc[k]) for k in D["cor"].index})
    log.info("Program blocks: %s", D["block"].value_counts().to_dict())
    return D


def plab(D, k, short=False):
    return f"P{k} {D['plabel'][k]}" if k in D["plabel"] else f"P{k}"


# ================================================================ Figure 1
def f1a_design(ax, D, P):
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 40)
    ax.axis("off")

    def box(x, y, w, h, title, lines, ec, fc="white", tc=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.8", fc=fc, ec=ec, lw=0.7))
        ax.text(x + 0.8, y + h - 1.0, title, fontsize=5.8, fontweight="bold", va="top", color=tc or Q.INK)
        ax.text(x + 0.8, y + h - 3.6, "\n".join(lines), fontsize=4.9, va="top", color=Q.INK2, linespacing=1.25)

    def arrow(x0, y0, x1, y1):
        ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=Q.MUTED, shrinkA=0, shrinkB=0,
                                                              mutation_scale=6))
    ax.text(0, 39.5, "Discovery (929 primary HCC)", fontsize=5.8, fontweight="bold", va="top")
    coh = [("TCGA-LIHC", "TCGA", ["n = 366 · RNA-seq (TPM)", "WES · CNA · OS / PFS", "mixed etiology, US"]),
           ("CLCA", "CLCA", ["n = 239 · RNA-seq (TPM)", "WGS · OS / RFS", "HBV-dominant, China"]),
           ("LICA-FR", "LICA_FR", ["n = 324 · RNA-seq (TPM)", "WES/WGS · CNA · no survival", "alcohol/NASH/viral, France"])]
    for i, (t, c, lines) in enumerate(coh):
        box(0, 26 - i * 12.5, 20.3, 11, t, lines, COH[c])
    steps = [("Harmonize", ["13,866 genes", "log2 TPM", "ComBat by cohort"]),
             ("MINER network", ["regulons (186", "technical excluded)", "→ programs", "→ states"]),
             ("Causal inference", ["113 genomic", "features: driver", "→ regulator", "→ regulon"]),
             ("Risk model", ["ridge on", "76 programs;", "TCGA ↔ CLCA", "36 months"])]
    for i, (t, lines) in enumerate(steps):
        x = 22.6 + i * 15.1
        box(x, 13.5, 14.1, 16, t, lines, Q.INK)
        arrow(x - 1.0 if i else 20.7, 21.5, x - 0.15, 21.5)
    ax.text(22.6, 11.2, "Annotation: immune classes (Montironi 2023), reference classes (Hoshida, Boyault), hallmarks",
            fontsize=4.9, color=Q.INK2, va="top")
    ax.text(84.5, 39.5, "External validation", fontsize=5.8, fontweight="bold", va="top")
    ext = [("GSE14520", ["n = 221 · Affymetrix", "OS / RFS · HBV, China"]),
           ("LIRI-JP", ["n = 203 · RNA-seq", "OS · HCV-dominant, Japan"]),
           ("GSE76427", ["n = 115 · Illumina", "OS / RFS · Singapore"])]
    for i, (t, lines) in enumerate(ext):
        box(84.5, 27.5 - i * 11.5, 15.5, 9.8, t, lines, EXT_COL[t.replace("LIRI-JP", "LIRI_JP")])
    arrow(83.2, 21.5, 84.3, 21.5)
    ax.set_ylim(4, 40)
    ax.text(84.5, 2.2, "fixed weights, never refit", fontsize=4.9, color=Q.INK2)


def f1b_numbers(ax, D, P):
    res = D["res"]
    mdir = os.path.join(res, "04_miner", P["miner"]["matrix"])
    rdf = pd.read_csv(os.path.join(mdir, "mechinf", "regulonDf.csv"), index_col=0)
    keep = set(json.load(open(os.path.join(mdir, "mechinf", "regulons_filtered.json"))))
    rk = rdf[rdf["Regulon_ID"].astype(str).isin(keep)]
    hc = pd.read_csv(os.path.join(res, "05_causal", P["miner"]["matrix"], "highConfidenceCausalResults.csv"), index_col=0)
    fam = pd.read_csv(os.path.join(res, "05_causal", P["miner"]["matrix"], "regulon_families.tsv"), sep="\t")
    rows = [("genes", int(D["genes"]["kept"].sum())),
            ("co-expression modules", len(json.load(open(os.path.join(mdir, "coexpr", "coexpressionDictionary.json"))))),
            ("regulons", len(keep)), ("regulators", rk["Regulator"].nunique()), ("regulon genes", rk["Gene"].nunique()),
            ("regulon families", fam["family"].nunique()), ("programs", len(D["programs"])), ("states", len(D["states"])),
            ("genomic features", D["genomic"].shape[0]), ("causal drivers", hc["Mutation"].nunique()),
            ("causal flows", len(hc))]
    y = np.arange(len(rows))
    v = np.array([r[1] for r in rows], float)
    ax.barh(y, np.log10(v), color=[Q.MUTED] * 6 + [Q.INK] * 2 + [Q.MUTED] * 3, height=0.62)
    for yi, (lab, n) in zip(y, rows):
        ax.text(np.log10(n) + 0.06, yi, f"{n:,}", va="center", fontsize=5.2, color=Q.INK)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=5.3)
    ax.invert_yaxis()
    ax.set_xlim(0, 5.9)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xticklabels(["1", "10", "10²", "10³", "10⁴"])
    ax.set_xlabel("count (log scale)")
    despine(ax)
    ax.set_title("Network by the numbers")
    return pd.DataFrame(rows, columns=["item", "n"])


def cluster(X):
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist
    Z = linkage(np.nan_to_num(pdist(np.nan_to_num(np.asarray(X, float)), "correlation"), nan=1.0), "average")
    return Z, leaves_list(Z)


def frac_table(states, labels, levels):
    out = pd.DataFrame(index=levels, columns=list(states), dtype=float)
    for k, v in states.items():
        x = labels.reindex(v).dropna()
        x = x[x != "unassigned"]
        for lv in levels:
            out.loc[lv, k] = (x == lv).mean() if len(x) else np.nan
    return out


def f1c_map(fig, sub, D):
    progs0, states0 = list(D["programs"]), list(D["states"])
    diff = D["diff"]
    M0 = pd.DataFrame({s: [diff.loc[D["programs"][k], [x for x in D["states"][s] if x in diff.columns]].values.mean()
                           for k in progs0] for s in states0}, index=progs0)
    from scipy.cluster.hierarchy import dendrogram
    Zc, co = cluster(M0.values.T)
    s_order = [states0[i] for i in co]
    # programs grouped by biology block, clustered within block
    p_order = []
    for b in list(BLOCK_COL) + ["Other"]:
        ks = [k for k in progs0 if D["block"].get(k) == b]
        if len(ks) > 2:
            _, o = cluster(M0.loc[ks].values)
            ks = [ks[i] for i in o]
        p_order += ks
    M = M0.loc[p_order, s_order]
    rz = pd.Series({s: D["risk"].reindex(D["states"][s]).mean() for s in s_order})
    tracks = []
    for cl, lv in PN.figures["fig1c_tracks"]:
        if cl not in D["ntp"]:
            continue
        T = frac_table(D["states"], D["ntp"][cl], lv)[s_order]
        T.index = [f"{PN.class_label(cl, l)}|{class_tag(cl, l)}" for l in lv]
        tracks.append(T)
    G = D["genomic"]
    mut = pd.DataFrame({s: [G.loc[m, [x for x in D["states"][s] if x in G.columns]].mean()
                            for m in ("MUT_CTNNB1", "MUT_TP53", "MUT_AXIN1", "TERT_promoter")] for s in s_order},
                       index=["CTNNB1 mut|", "TP53 mut|", "AXIN1 mut|", "TERT promoter|"])
    tracks.append(mut)
    nrow = [1.6, 2.2] + [len(t) for t in tracks] + [len(M) * 0.42]
    gs = gridspec.GridSpecFromSubplotSpec(len(nrow), 4, subplot_spec=sub, height_ratios=nrow,
                                          width_ratios=[0.36, 0.05, 1, 0.34], hspace=0.12, wspace=0.02)
    x = np.arange(len(s_order))
    axd = fig.add_subplot(gs[0, 2])
    dendrogram(Zc, ax=axd, no_labels=True, color_threshold=0, link_color_func=lambda k: Q.MUTED)
    for c in axd.collections:
        c.set_linewidth(0.4)
    axd.axis("off")
    letter(axd, "c", x=-0.42, y=0.5)
    axd.set_title("Programs × states (clustered): published-class biology, drivers and risk", fontsize=6.5, pad=1)
    axr = fig.add_subplot(gs[1, 2])
    axr.bar(x, rz.values, color=[ADV if v > 0 else PROT for v in rz.values], width=0.75)
    axr.axhline(0, color=Q.AXIS, lw=0.4)
    axr.set_xlim(-0.5, len(x) - 0.5)
    axr.set_xticks([])
    axr.set_yticks([])
    for sp in axr.spines.values():
        sp.set_visible(False)
    axr.text(-0.6, 0, "risk", ha="right", va="center", fontsize=5.2)
    for i, T in enumerate(tracks):
        ax = fig.add_subplot(gs[2 + i, 2])
        imt = ax.imshow(T.values.astype(float), aspect="auto", cmap=FRAC_CMAP, vmin=0, vmax=1, interpolation="none")
        ax.set_yticks(range(len(T)))
        ax.set_yticklabels([t.split("|")[0] for t in T.index], fontsize=4.8)
        ax.set_xticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.3)
        axt = fig.add_subplot(gs[2 + i, 3])
        axt.set_ylim(len(T) - 0.5, -0.5)
        axt.axis("off")
        for j, t in enumerate(T.index):
            axt.text(0.03, j, t.split("|")[1], fontsize=4.4, color=Q.MUTED, va="center")
    ax = fig.add_subplot(gs[-1, 2])
    v = np.nanpercentile(np.abs(M.values), 98)
    im = ax.imshow(M.values, aspect="auto", cmap=Q.DIV, vmin=-v, vmax=v, interpolation="none")
    ax.set_yticks([])
    ax.set_xticks(x)
    ax.set_xticklabels([f"S{s}" for s in s_order], fontsize=4.2, rotation=90)
    ax.set_xlabel("MINER transcriptional state", fontsize=5.5, labelpad=1)
    for sp in ax.spines.values():
        sp.set_linewidth(0.3)
    # block strip + labels
    axb = fig.add_subplot(gs[-1, 1])
    blocks = [D["block"][k] for k in p_order]
    for i, b in enumerate(blocks):
        axb.add_patch(plt.Rectangle((0, i - 0.5), 1, 1, color=BLOCK_COL.get(b, Q.BACKGROUND), lw=0))
    axb.set_xlim(0, 1)
    axb.set_ylim(len(blocks) - 0.5, -0.5)
    axb.axis("off")
    axl = fig.add_subplot(gs[-1, 3])
    axl.set_ylim(len(blocks) - 0.5, -0.5)
    axl.axis("off")
    start = 0
    for b in list(BLOCK_COL) + ["Other"]:
        n = blocks.count(b)
        if n:
            axl.text(0.03, start + n / 2 - 0.5, f"{b}\n({n} programs)", fontsize=4.8, va="center",
                     color=BLOCK_COL.get(b, Q.MUTED), fontweight="bold" if b != "Other" else "normal")
            start += n
    cax = axl.inset_axes([0.05, -0.075, 0.85, 0.018])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("regulon dysregulation", fontsize=4.6, labelpad=1)
    cb.ax.tick_params(labelsize=4.2, length=1.5, width=0.3)
    cb.outline.set_linewidth(0.3)
    cax = axl.inset_axes([0.05, -0.2, 0.85, 0.018])
    cb = fig.colorbar(imt, cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
    cb.set_label("fraction of tumours (tracks)", fontsize=4.6, labelpad=1)
    cb.ax.tick_params(labelsize=4.2, length=1.5, width=0.3)
    cb.outline.set_linewidth(0.3)
    return M, pd.DataFrame(tracks[0])


def f1d_signatures(ax, D):
    ks = [k for k in D["plabel"] if k in D["cor"].index]
    ks = sorted(ks, key=lambda k: (list(BLOCK_COL).index(D["block"][k]) if D["block"][k] in BLOCK_COL else 9,
                                   -D["weights"].get(k, 0)))
    panel = PN.signatures(available=set(D["cor"].columns))
    sigs = [x["set"] for x in panel]
    C = D["cor"].loc[ks, sigs]
    im = ax.imshow(C.values, aspect="auto", cmap=COR_CMAP, vmin=-1, vmax=1, interpolation="none")
    ax.set_yticks(range(len(ks)))
    ax.set_yticklabels([plab(D, k) for k in ks], fontsize=4.6)
    for t, k in zip(ax.get_yticklabels(), ks):
        t.set_color(ADV if D["weights"].get(k, 0) > 0 else PROT)
    info = {x["set"]: (x["label"], x["tag"], x["block"]) for x in panel}
    ax.set_xticks(range(len(sigs)))
    ax.set_xticklabels([f"{info[s][0]} ({info[s][1]})" for s in sigs], fontsize=4.1, rotation=60, ha="right",
                       rotation_mode="anchor")
    for t, s in zip(ax.get_xticklabels(), sigs):
        t.set_color(BLOCK_COL[info[s][2]])
    for sp in ax.spines.values():
        sp.set_linewidth(0.3)
    # group separators
    g = [info[s][2] for s in sigs]
    for j in range(1, len(g)):
        if g[j] != g[j - 1]:
            ax.axvline(j - 0.5, color="white", lw=1.2)
    bk = [D["block"][k] for k in ks]
    for i in range(1, len(bk)):
        if bk[i] != bk[i - 1]:
            ax.axhline(i - 0.5, color="white", lw=1.2)
    ax.set_title("What the programs are (activity r with reference signatures)", pad=13)
    cax = ax.inset_axes([0.78, 1.02, 0.22, 0.018])
    cb = plt.colorbar(im, cax=cax, orientation="horizontal", ticks=[-1, 0, 1])
    cb.ax.xaxis.set_ticks_position("top")
    cb.ax.tick_params(labelsize=4.2, length=1.5, width=0.3, pad=1)
    cb.outline.set_linewidth(0.3)
    ax.text(0.77, 1.029, "r", transform=ax.transAxes, fontsize=5, ha="right", va="center")
    ax.text(0.0, 1.035, "row colour: risk weight (red adverse, blue protective)", transform=ax.transAxes, fontsize=4.4,
            color=Q.INK2, va="bottom", ha="left")
    return C


def curve(ax, x0, y0, x1, y1, color, lw, ls="-"):
    xm = (x0 + x1) / 2
    path = Path([(x0, y0), (xm, y0), (xm, y1), (x1, y1)], [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
    ax.add_patch(PathPatch(path, facecolor="none", edgecolor=color, lw=lw, ls=ls, alpha=0.85, zorder=1))


def f1e_ctnnb1(ax, D, P, n_top=12, force=("LEF1", "TCF7", "TCF7L1")):
    res = D["res"]
    hc = pd.read_csv(os.path.join(res, "05_causal", P["miner"]["matrix"], "highConfidenceCausalResults.csv"), index_col=0)
    hc["program"] = hc["program"].astype(int).astype(str)
    x = hc[hc["Mutation"] == "MUT_CTNNB1"].assign(ad=lambda d: d["cohen_d"].abs())
    x = x.sort_values("ad", ascending=False).drop_duplicates("family")
    sel = pd.concat([x.head(n_top), x[x["regulator_symbol"].isin(force)]]).drop_duplicates("family")
    sel = sel.drop_duplicates("regulator_symbol")
    sel = sel.sort_values(["MutationRegulatorEdge", "ad"], ascending=[False, False])
    pg = sel.groupby("program").agg(n=("family", "size"), d=("cohen_d", "mean")).reset_index()
    pg = pg.assign(w=pg["program"].map(D["weights"])).sort_values("w")
    ry = dict(zip(sel["regulator_symbol"], np.linspace(0.95, 0.05, len(sel))))
    py = dict(zip(pg["program"], np.linspace(0.92, 0.08, len(pg))))
    X0, X1, X2 = 0.07, 0.36, 0.55
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    for _, r in sel.iterrows():
        curve(ax, X0 + 0.055, 0.5, X1 - 0.05, ry[r["regulator_symbol"]], ADV if r["MutationRegulatorEdge"] > 0 else PROT,
              0.35 + 0.9 * min(r["ad"], 2) / 2)
        curve(ax, X1 + 0.05, ry[r["regulator_symbol"]], X2 - 0.005, py[r["program"]], ADV if r["cohen_d"] > 0 else PROT,
              0.35 + 1.2 * min(r["ad"], 2) / 2, ls="-" if r["RegulatorRegulon_Spearman_R"] > 0 else (0, (2.5, 1.5)))
    ax.add_patch(FancyBboxPatch((X0 - 0.055, 0.45), 0.11, 0.1, boxstyle="round,pad=0.01", fc=Q.INK, ec=Q.INK))
    ax.text(X0, 0.5, "CTNNB1\nmutation", color="white", ha="center", va="center", fontsize=5.5, fontweight="bold")
    ax.text(X0, 0.42, f"n = {int(x['n_altered'].iloc[0])}", ha="center", va="top", fontsize=4.6, color=Q.INK2)
    for _, r in sel.iterrows():
        up = r["MutationRegulatorEdge"] > 0
        ax.text(X1, ry[r["regulator_symbol"]], f"{r['regulator_symbol']} {'↑' if up else '↓'}", ha="center", va="center",
                fontsize=5, fontweight="bold" if r["regulator_symbol"] in force else "normal",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=Q.MUTED, lw=0.4))
    for _, r in pg.iterrows():
        k = r["program"]
        adverse = D["weights"].get(k, 0) > 0
        ax.text(X2, py[k], f"{plab(D, k)} {'↑' if r['d'] > 0 else '↓'}", ha="left", va="center", fontsize=4.9,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=ADV if adverse else PROT, lw=0.8))
    ax.set_title("CTNNB1 causal flow: driver → regulators → programs", fontsize=6.5)
    ax.text(0.0, -0.04, "edges: red up / blue down in mutant tumours; dashed = repressor; width ~ |d|; "
            "box border: red adverse / blue protective program", transform=ax.transAxes, fontsize=4.4, color=Q.INK2,
            va="top")
    return sel


def f1f_axin1(fig, sub, D, P):
    """AXIN1-only vs CTNNB1-only mutants against double wild type (step 07f): every program, then marker genes."""
    adir = os.path.join(D["res"], "07_post", "axin1_ctnnb1")
    T = pd.read_csv(os.path.join(adir, "program_effects.tsv"), sep="\t", dtype={"program": str}).set_index("program")
    G = pd.read_csv(os.path.join(adir, "gene_effects.tsv"), sep="\t", index_col=0)
    grp = pd.read_csv(os.path.join(adir, "groups.tsv"), sep="\t", index_col=0)["group"].value_counts()
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=sub, width_ratios=[1.5, 1], wspace=0.5)
    ax = fig.add_subplot(gs[0])
    col = {"CTNNB1 only": Q.SLOTS[0], "shared": Q.INK, "opposite": Q.SLOTS[1], "AXIN1 only": Q.SLOTS[2], "neither": "#c4c4c0"}
    ax.axhline(0, color=Q.AXIS, lw=0.4)
    ax.axvline(0, color=Q.AXIS, lw=0.4)
    for k in ("neither", "CTNNB1 only", "shared", "AXIN1 only", "opposite"):
        d = T[T["pattern"] == k]
        ax.scatter(d["beta_CTNNB1"], d["beta_AXIN1"], s=7, color=col[k], edgecolor="white", linewidth=0.25, zorder=3,
                   label=f"{k} ({len(d)})")
    names = {"2": ("WNT/β-catenin", -2, -7, "right"), "7": ("CTNNB1 proteome", -2, 3, "right"),
             "47": ("NFE2L2/JUN", 3, -5, "left"), "64": ("KLF16/NFAT5", 3, -2, "left"), "21": ("Progenitor", 3, 2, "left"),
             "60": ("Interferon-γ", -4, -4, "right"), "11": ("Immune class", -4, 0, "right"),
             "1": ("FOXJ3/NR3C1", 3, -5, "left")}
    for k, (nm, dx, dy, ha) in names.items():
        ax.annotate(nm, (T.loc[k, "beta_CTNNB1"], T.loc[k, "beta_AXIN1"]), fontsize=4.1, xytext=(dx, dy),
                    textcoords="offset points", ha=ha)
    ax.set_xlabel("CTNNB1-mutant effect (s.d.)")
    ax.set_ylabel("AXIN1-mutant effect (s.d.)")
    ax.margins(x=0.1, y=0.12)
    ax.set_xlim(left=-1.9)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.14), ncol=3, fontsize=4.1, handletextpad=0.0, columnspacing=0.4,
              borderaxespad=0)
    despine(ax)
    fig.text(ax.get_position().x0, ax.get_position().y1 + 0.036,
             f"AXIN1 (n = {grp['AXIN1']}) and CTNNB1 (n = {grp['CTNNB1']}) mutations: program effects", fontsize=6.5, va="bottom")

    ax2 = fig.add_subplot(gs[1])
    names = {"wnt_targets": "WNT targets", "antigen_presentation": "MHC class I", "t_cell": "T cell"}
    y, yt, yl = 0, [], []
    for st, nm in names.items():
        g = G[G["set"] == st]
        ax2.text(-0.02, y - 0.75, nm, fontsize=4.4, color=Q.MUTED, transform=ax2.get_yaxis_transform(), ha="right", va="center")
        for gene, r in g.iterrows():
            ax2.plot([r["beta_AXIN1"], r["beta_CTNNB1"]], [y, y], color=Q.AXIS, lw=0.5, zorder=1)
            ax2.scatter(r["beta_AXIN1"], y, s=9, color=Q.SLOTS[2], zorder=3, label="AXIN1" if y == 0 else None)
            ax2.scatter(r["beta_CTNNB1"], y, s=9, color=Q.SLOTS[0], zorder=3, label="CTNNB1" if y == 0 else None)
            yt.append(y)
            yl.append(gene)
            y += 1
        y += 1.3
    ax2.axvline(0, color=Q.AXIS, lw=0.4)
    ax2.set_yticks(yt)
    ax2.set_yticklabels(yl, fontsize=4.4, style="italic")
    ax2.set_ylim(y - 1.6, -1.6)
    ax2.set_xlabel("effect on expression (s.d.)")
    ax2.legend(loc="upper left", bbox_to_anchor=(0.0, 1.14), ncol=2, fontsize=4.1, handletextpad=0.0, columnspacing=0.4,
               borderaxespad=0)
    despine(ax2)
    return ax, T, G


def f1f_push(fig, sub, D, P):
    """Driver risk push vs driver prognosis, out of sample only (step 07e): program weights trained in one
    cohort, driver Cox z measured in the other. Pre-specified predictor and endpoint (post.driver_push.primary)."""
    ddir = os.path.join(D["res"], "07_post", "driver_push")
    prim = P["post"]["driver_push"]["primary"]
    pred, ep = prim["predictor"], prim["endpoint"]
    T = pd.read_csv(os.path.join(ddir, "driver_table.tsv"), sep="\t")
    S = pd.read_csv(os.path.join(ddir, "summary.tsv"), sep="\t")
    T = T[(T["kind"] == "out_of_sample") & (T["endpoint"] == ep)]
    S = S[(S["kind"] == "out_of_sample") & (S["endpoint"] == ep) & (S["adjustment"] == prim["adjustment"])
          & (S["driver_set"] == prim["driver_set"]) & (S["predictor"] == pred)].set_index("pairing")
    order = ["CLCA->TCGA", "TCGA->CLCA"]
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=sub, wspace=0.34)
    lab = lambda f: (f.replace("MUT_", "").replace("PATH_", "").replace("AMP_", "amp ").replace("ARM_", "")  # noqa: E731
                     .replace("_", " ").replace("p53 cell cycle", "p53 pathway").replace("NRF2 oxidative stress", "NRF2")
                     .replace("CCND1 FGF19", "CCND1"))
    show = {"MUT_TP53", "MUT_CTNNB1", "MUT_AXIN1", "MUT_NFE2L2", "MUT_RB1", "MUT_TSC2", "MUT_BAP1", "MUT_APC",
            "PATH_WNT", "PATH_p53_cell_cycle", "ARM_17p_loss", "ARM_13q_loss", "MUT_ARID1A"}
    axes = []
    for j, pr in enumerate(order):
        ax = fig.add_subplot(gs[j])
        axes.append(ax)
        d = T[T["pairing"] == pr].dropna(subset=[pred, "z"])
        ax.axhline(0, color=Q.AXIS, lw=0.4)
        ax.axvline(0, color=Q.AXIS, lw=0.4)
        x = 1000 * d[pred]
        ax.scatter(x, d["z"], s=5 + 0.04 * d["n_altered_test"], color=[Q.SLOTS[1] if t == "arm_cna" else Q.INK for t in d["type"]],
                   edgecolor="white", linewidth=0.3, zorder=3)
        for _, r in d.iterrows():
            if r["driver"] in show:
                right = 1000 * r[pred] > x.min() + 0.7 * (x.max() - x.min())
                ax.annotate(lab(r["driver"]), (1000 * r[pred], r["z"]), fontsize=4.1, textcoords="offset points",
                            xytext=(-2, 1.5) if right else (2, 1.5), ha="right" if right else "left")
        r = S.loc[pr]
        tr, te = pr.split("->")
        ax.text(0.0, 1.015, f"ρ = {r['rho']:.2f} ({r['boot_lo']:.2f} to {r['boot_hi']:.2f}); "
                f"P = {r['p_perm_weights']:.2f}\n{int(r['n'])} drivers, {int(r['n_clusters'])} clusters",
                transform=ax.transAxes, fontsize=4.6, va="bottom")
        ax.set_title(f"weights {tr} → outcome {te}", fontsize=5.6, pad=14, loc="left")
        ax.set_xlabel("risk push ×10³ (mean d × weight)", fontsize=5.2)
        ax.margins(x=0.12, y=0.12)
        despine(ax)
    axes[0].set_ylabel("driver recurrence association (Cox z)")
    axes[0].scatter([], [], color=Q.INK, s=7, label="mutation / focal / pathway")
    axes[0].scatter([], [], color=Q.SLOTS[1], s=7, label="arm-level CNA")
    axes[0].legend(loc="lower right", fontsize=4.2, handletextpad=0.1, borderaxespad=0.1)
    fig.text(axes[0].get_position().x0, axes[0].get_position().y1 + 0.036,
             "Driver risk push vs driver prognosis, out of sample", fontsize=6.5, va="bottom")
    return axes[0], T


def figure1(D, P, outdir, log):
    fig = plt.figure(figsize=(W, 238 * MM))
    gs = gridspec.GridSpec(3, 1, height_ratios=[40, 120, 62], hspace=0.36, left=0.03, right=0.985, top=0.985, bottom=0.03)
    top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[3.25, 1], wspace=0.28)
    ax = fig.add_subplot(top[0])
    f1a_design(ax, D, P)
    letter(ax, "a", x=0.0, y=0.99)
    ax = fig.add_subplot(top[1])
    tb = f1b_numbers(ax, D, P)
    letter(ax, "b", x=-0.33)
    mid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1], width_ratios=[1.25, 1], wspace=0.55)
    M, _ = f1c_map(fig, mid[0], D)
    dsub = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=mid[1], height_ratios=[1, 0.4], hspace=0.0)
    axd = fig.add_subplot(dsub[0])
    C = f1d_signatures(axd, D)
    letter(axd, "d", x=-0.42, y=1.07)
    bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], width_ratios=[1.38, 1], wspace=0.16)
    ax = fig.add_subplot(bot[0])
    e = f1e_ctnnb1(ax, D, P)
    letter(ax, "e", x=0.0, y=1.0)
    ax, f, f2 = f1f_axin1(fig, bot[1], D, P)
    letter(ax, "f", x=-0.22, y=1.3)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"figure1.{ext}"))
    plt.close(fig)
    with pd.ExcelWriter(os.path.join(outdir, "figure1_source_data.xlsx")) as xw:
        tb.to_excel(xw, "1b_numbers", index=False)
        M.to_excel(xw, "1c_programs_x_states")
        C.to_excel(xw, "1d_program_signature_r")
        e.to_excel(xw, "1e_ctnnb1_flows")
        f.to_excel(xw, sheet_name="1f_programs")
        f2.to_excel(xw, sheet_name="1f_genes")
    log.info("Figure 1 written")


# ================================================================ Figure 2
def f2a_weights(ax, D, n=8):
    w = D["weights"].sort_values()
    sel = pd.concat([w.head(n), w.tail(n)])
    y = np.r_[np.arange(n), np.arange(n) + n + 1]
    ax.barh(y, 1000 * sel.values, color=[ADV if v > 0 else PROT for v in sel], height=0.7)
    ax.axvline(0, color=Q.AXIS, lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels([plab(D, k) for k in sel.index], fontsize=4.8)
    for t, v in zip(ax.get_yticklabels(), sel.values):
        t.set_color(ADV if v > 0 else PROT)
    ax.text(0, n, f"… {len(w) - 2 * n} programs with smaller weights …", ha="center", va="center", fontsize=4.5,
            color=Q.MUTED)
    ax.set_ylim(-0.7, 2 * n + 0.7)
    ax.set_xlabel("ridge weight ×10³ (TCGA, RFS)\n← protective   adverse →")
    despine(ax)
    ax.set_title(f"Risk model: largest of {len(w)} program weights")
    return w


def f2b_classes(fig, sub, D):
    risk = D["risk"]
    panels = [(cl, D["ntp"][cl], lv, PN.class_title(cl)) for cl, lv in PN.figures["fig2b_panels"] if cl in D["ntp"]]
    gs = gridspec.GridSpecFromSubplotSpec(1, len(panels), subplot_spec=sub, wspace=0.12,
                                          width_ratios=[max(3, len(lv)) for _, _, lv, _ in panels])
    rng = np.random.default_rng(0)
    axes = []
    for j, (cl, lab, levels, name) in enumerate(panels):
        ax = fig.add_subplot(gs[j], sharey=axes[0] if axes else None)
        axes.append(ax)
        lab = lab.reindex(risk.index)
        groups = [risk[lab == lv].dropna() for lv in levels]
        for i, g in enumerate(groups):
            col = AXIS_COL.get(PN.class_axis(cl, levels[i]), Q.MUTED)
            ax.scatter(i + rng.uniform(-0.25, 0.25, len(g)), g, s=1.2, color=col, alpha=0.45, lw=0)
            if len(g):
                q1, md, q3 = np.percentile(g, [25, 50, 75])
                ax.plot([i, i], [q1, q3], color=Q.INK, lw=1.2)
                ax.plot([i - 0.28, i + 0.28], [md, md], color=Q.INK, lw=1.0)
        ok = [g for g in groups if len(g) >= 3]
        pk = stats.kruskal(*ok).pvalue if len(ok) >= 2 else np.nan
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels([PN.class_label(cl, lv) for lv in levels], fontsize=4.4, rotation=45, ha="right",
                           rotation_mode="anchor")
        ax.set_title(f"{name}\nP = {pk:.0e}".replace("e-0", "e-"), fontsize=4.8)
        ax.axhline(0, color=Q.AXIS, lw=0.4)
        despine(ax)
        if j:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel("risk score (z)")
    return axes[0]


def f2c_states(fig, sub, D, m07c, P):
    S = D["samples"]
    st = pd.Series({k: D["risk"].reindex(v).mean() for k, v in D["states"].items()})
    s_order = st.sort_values().index.tolist()
    guan = pd.concat([m07c.guan_rank(D["surv"][(c, "RFS")]) for c in ("TCGA", "CLCA")])
    rows, labs = [], []
    for cl, lv in PN.figures["fig2c_rows"]:
        if cl not in D["ntp"]:
            continue
        rows.append(frac_table(D["states"], D["ntp"][cl], [lv]).loc[lv, s_order].values.astype(float))
        labs.append(f"{PN.class_label(cl, lv)} ({class_tag(cl, lv)})")
    G = D["genomic"]
    for m, lab in (("MUT_TP53", "TP53 mut"), ("MUT_CTNNB1", "CTNNB1 mut")):
        rows.append([G.loc[m, [x for x in D["states"][s] if x in G.columns]].mean() for s in s_order])
        labs.append(lab)
    gs = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=sub, height_ratios=[len(rows), 5, 5], hspace=0.1,
                                          width_ratios=[0.2, 1], wspace=0.0)
    ax0 = fig.add_subplot(gs[0, 1])
    im0 = ax0.imshow(np.array(rows), aspect="auto", cmap=FRAC_CMAP, vmin=0, vmax=1, interpolation="none")
    cax = ax0.inset_axes([0.86, 1.04, 0.14, 0.06])
    cb = fig.colorbar(im0, cax=cax, orientation="horizontal", ticks=[0, 0.5, 1])
    cb.ax.xaxis.set_ticks_position("top")
    cb.ax.tick_params(labelsize=4.2, length=1.5, width=0.3, pad=1)
    cb.outline.set_linewidth(0.3)
    ax0.text(0.855, 1.07, "fraction of tumours", transform=ax0.transAxes, fontsize=4.4, ha="right", va="center", color=Q.INK2)
    ax0.set_yticks(range(len(labs)))
    ax0.set_yticklabels(labs, fontsize=4.5)
    ax0.set_xticks([])
    for sp in ax0.spines.values():
        sp.set_linewidth(0.3)
    ax0.set_title("States ordered by mean risk: class composition, predicted risk, observed outcome", loc="left", y=1.02)
    x = np.arange(len(s_order))
    for axi, vals, lab in ((1, D["risk"], "risk score (z)"), (2, guan, "GuanRank (RFS)")):
        ax = fig.add_subplot(gs[axi, 1], sharex=ax0)
        data = [vals.reindex(D["states"][s]).dropna() for s in s_order]
        bp = ax.boxplot(data, positions=x, widths=0.55, showfliers=False, patch_artist=True,
                        medianprops=dict(color=Q.INK, lw=0.8), whiskerprops=dict(color=Q.MUTED, lw=0.4),
                        capprops=dict(color=Q.MUTED, lw=0.4), boxprops=dict(lw=0.4))
        med = np.array([np.median(d) if len(d) else np.nan for d in data])
        cmap = plt.get_cmap(Q.DIV)
        for patch, s in zip(bp["boxes"], s_order):
            patch.set_facecolor(cmap(0.5 + np.clip(st[s], -1.2, 1.2) / 2.4))
            patch.set_edgecolor(Q.MUTED)
        ax.set_ylabel(lab, fontsize=5)
        ax.set_xlim(-0.6, len(x) - 0.4)
        despine(ax)
        if axi == 1:
            ax.axhline(0, color=Q.AXIS, lw=0.4)
            plt.setp(ax.get_xticklabels(), visible=False)
        else:
            n = [len(d) for d in data]
            ok = [i for i, k in enumerate(n) if k >= 5]
            rho, pv = stats.spearmanr(st[[s_order[i] for i in ok]], med[ok])
            ax.text(0.01, 0.97, f"state risk vs median GuanRank: ρ = {rho:.2f}, P = {pv:.1e} ({len(ok)} states)",
                    transform=ax.transAxes, fontsize=4.8, va="top")
            ax.set_xticks(x)
            ax.set_xticklabels([f"S{s}" for s in s_order], fontsize=4.2, rotation=90)
    return ax0


def f2d_e_km(fig, sub, D, P):
    res = D["res"]
    gs = gridspec.GridSpecFromSubplotSpec(1, 5, subplot_spec=sub, wspace=0.18)
    panels = []
    for tr, te in (("TCGA", "CLCA"), ("CLCA", "TCGA")):
        pdir = os.path.join(D["rdir"], f"predictor_ridge_programs_{tr}_RFS_h36m")
        pr = pd.read_csv(os.path.join(pdir, "predictions.tsv"), sep="\t", index_col=0)
        ev = pd.read_csv(os.path.join(pdir, "evaluation.tsv"), sep="\t").set_index("cohort").loc[te]
        s = D["surv"][(te, "RFS")]
        panels.append((f"{te} RFS", f"trained in {tr}", s, pr.loc[s.index, "predicted_high_risk"].astype(bool), ev, COH[te]))
    for c, ep in (("GSE14520", "OS"), ("LIRI_JP", "OS"), ("GSE76427", "RFS")):
        sc = pd.read_csv(os.path.join(res, "08_validation", c, "scores.tsv"), sep="\t", index_col=0)
        ev = pd.read_csv(os.path.join(res, "08_validation", c, "evaluation.tsv"), sep="\t")
        m = "TCGA"   # pre-specified primary model for every external cohort
        e = ev[(ev["endpoint"] == ep) & ev["model"].str.contains(m)].iloc[0]
        s = sc[[f"{ep}_time", f"{ep}_event"]].dropna().set_axis(["duration", "observed"], axis=1)
        over = s["duration"] > 1096
        s.loc[over, "duration"], s.loc[over, "observed"] = 1096, 0
        r = sc.loc[s.index, f"risk_{m}_{ep}_h36m"]
        panels.append((f"{c.replace('_', '-')} {ep}", f"external, {m}-trained", s, r >= r.quantile(0.8),
                       e.rename({"hr_top": "hr", "hr_top_p": "hr_p"}), EXT_COL[c]))
    axes = []
    for j, (title, sub_, s, high, ev, col) in enumerate(panels):
        ax = fig.add_subplot(gs[j], sharey=axes[0] if axes else None)
        axes.append(ax)
        km(ax, s.loc[high[high].index, "duration"], s.loc[high[high].index, "observed"], ADV, "top 20%")
        km(ax, s.loc[high[~high].index, "duration"], s.loc[high[~high].index, "observed"], PROT, "rest")
        ax.set_ylim(0, 1.03)
        ax.set_xlim(0, 36.5)
        ax.set_xticks([0, 12, 24, 36])
        ax.set_title(f"{title}\n{sub_}", fontsize=5.2, color=col if j < 2 else Q.INK)
        ax.text(0.04, 0.05, f"C = {ev['c_index']:.2f}\nHR = {ev['hr']:.2f}, P = {ev['hr_p']:.0e}".replace("e-0", "e-"),
                transform=ax.transAxes, fontsize=4.6)
        ax.set_xlabel("months")
        despine(ax)
        if j:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel("event-free probability")
            ax.legend(loc="upper right", fontsize=4.5)
    return axes[0], axes[2]


def f2f_forest(fig, sub, D):
    res = D["res"]
    rows = []
    for tr, ep, te in (("TCGA", "RFS", "CLCA"), ("TCGA", "OS", "CLCA"), ("CLCA", "RFS", "TCGA")):
        e = pd.read_csv(os.path.join(D["rdir"], f"predictor_ridge_programs_{tr}_{ep}_h36m", "evaluation.tsv"), sep="\t")
        e = e.set_index("cohort").loc[te]
        rows.append((f"{te} {ep} (from {tr})", "discovery", e["hr_per_sd"], e["hr_per_sd_lo"], e["hr_per_sd_hi"], COH[te]))
    ext = []
    for c in ("GSE14520", "LIRI_JP", "GSE76427"):
        e = pd.read_csv(os.path.join(res, "08_validation", c, "evaluation.tsv"), sep="\t")
        for _, r in e[e["model"].str.contains("TCGA")].iterrows():
            rows.append((f"{c.replace('_', '-')} {r['endpoint']}", "external", r["hr_per_sd"], r["hr_per_sd_lo"],
                         r["hr_per_sd_hi"], EXT_COL[c]))
            ext.append((r["endpoint"], np.log(r["hr_per_sd"]), (np.log(r["hr_per_sd_hi"]) - np.log(r["hr_per_sd_lo"])) / 3.92))
    E = pd.DataFrame(ext, columns=["ep", "b", "se"])
    for ep, g in E.groupby("ep"):
        w = 1 / g["se"] ** 2
        q = (w * (g["b"] - (w * g["b"]).sum() / w.sum()) ** 2).sum()
        tau2 = max(0, (q - (len(g) - 1)) / (w.sum() - (w ** 2).sum() / w.sum())) if len(g) > 1 else 0
        wr = 1 / (g["se"] ** 2 + tau2)
        b, se = (wr * g["b"]).sum() / wr.sum(), np.sqrt(1 / wr.sum())
        rows.append((f"pooled external {ep} (RE)", "pooled", np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), Q.INK))
    g3 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=sub, width_ratios=[0.62, 1, 0.5], wspace=0.02)
    ax = fig.add_subplot(g3[1])
    axl = fig.add_subplot(g3[0], sharey=ax)
    axn = fig.add_subplot(g3[2], sharey=ax)
    y = np.arange(len(rows))
    for i, (lab, kind, hr, lo, hi, col) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=col, lw=0.8)
        ax.scatter(hr, i, marker="D" if kind == "pooled" else "s", s=10 if kind == "pooled" else 6, color=col, zorder=3)
        axn.text(0.05, i, f"{hr:.2f} ({lo:.2f}–{hi:.2f})", fontsize=4.4, va="center")
        axl.text(0.98, i, lab, fontsize=4.7, va="center", ha="right", color=Q.INK)
    ax.axvline(1, color=Q.AXIS, lw=0.4)
    ax.set_xscale("log")
    ax.set_xlim(0.45, 3.4)
    ax.set_xticks([0.5, 1, 2, 3])
    ax.set_xticklabels(["0.5", "1", "2", "3"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_yticks([])
    ax.set_ylim(len(rows) - 0.5, -0.5)
    for a_ in (axl, axn):
        a_.axis("off")
    for i in range(1, len(rows)):
        if rows[i][1] != rows[i - 1][1]:
            ax.axhline(i - 0.5, color=Q.GRID, lw=0.5)
    ax.set_xlabel("HR per s.d. of risk score")
    despine(ax)
    ax.set_title("HR per s.d., all test cohorts")
    ax.spines["left"].set_visible(False)
    return ax, pd.DataFrame(rows, columns=["cohort", "kind", "hr", "lo", "hi", "color"]).drop(columns="color")


def f2g_benchmark(ax, D):
    res = D["res"]
    h = pd.read_csv(os.path.join(res, "07_post", "figures", "risk_head_to_head.tsv"), sep="\t")
    h["set"] = h["test"] + " RFS (from " + h["train"] + ")"
    hs = [h]
    for c in ("GSE14520", "LIRI_JP", "GSE76427"):
        x = pd.read_csv(os.path.join(res, "08_validation", c, "head_to_head.tsv"), sep="\t")
        x["set"] = x["cohort"].str.replace("_", "-") + " " + x["endpoint"]
        hs.append(x)
    H = pd.concat(hs)
    H["model"] = H["model"].replace({"MINER programs (ridge)": "MINER programs", "MINER programs (TCGA-trained)": "MINER programs"})
    sets = list(dict.fromkeys(H["set"]))
    models = [(m["name"], m["name"], MODEL_COL.get(m.get("colour"), Q.MUTED), m.get("marker", "o"))
              for m in PN.cfg.get("known_models", [])] + [("MINER programs", "MINER programs", ADV, "D")]
    y = np.arange(len(sets))
    for m, lab, col, mk in models:
        v = H[H["model"] == m].set_index("set").reindex(sets)["c_index"]
        ax.scatter(v, y, s=10 if m != "MINER programs" else 12, color=col, marker=mk, label=lab, zorder=3,
                   edgecolor=Q.INK if m == "MINER programs" else "none", linewidth=0.3)
    for i in y:
        ax.axhline(i, color=Q.GRID, lw=0.3, zorder=0)
    ax.axvline(0.5, color=Q.AXIS, lw=0.4)
    ax.set_yticks(y)
    ax.set_yticklabels(sets, fontsize=4.7)
    ax.invert_yaxis()
    ax.set_xlim(0.42, 0.78)
    ax.set_xlabel("C-index")
    ax.legend(loc="lower left", bbox_to_anchor=(-0.5, 1.0), ncol=2, fontsize=4.2, handletextpad=0.2, columnspacing=0.6,
              borderaxespad=0.2)
    despine(ax)
    lr = H[(H["model"] == "MINER programs") & H["lr_p"].notna()]
    ax.set_title("vs published signatures", pad=24, loc="left", x=-0.5)
    ax.text(0.99, 0.875, f"MINER added to all published: P ≥ {lr['lr_p'].min():.2f}", transform=ax.transAxes, fontsize=4.3,
            ha="right", va="center", color=Q.INK2)
    return H


def f2h_loco(ax, D):
    """Leave-one-cohort-out stability (step 08b); placeholder if not run yet."""
    f = os.path.join(D["res"], "08_validation", "loco", "loco_summary.tsv")
    if not os.path.exists(f):
        ax.axis("off")
        ax.text(0.5, 0.5, "Leave-one-cohort-out\n(pending)", ha="center", va="center", fontsize=5.5, color=Q.MUTED)
        return None
    S = pd.read_csv(f, sep="\t").set_index("held_out")
    mets = [("regulators_recovered", "regulators\nrecovered"),
            ("program_median_activity_r", "program\nactivity r"),
            ("risk_programs_median_activity_r", "risk-program\nactivity r"),
            ("causal_MUT_CTNNB1_recovered_filtered", "CTNNB1\nedges"),
            ("causal_edges_recovered_filtered", "all causal\nedges")]
    x = np.arange(len(mets))
    hs = [h for h in ("TCGA", "CLCA", "LICA_FR") if h in S.index]
    for j, h in enumerate(hs):
        ax.bar(x + (j - (len(hs) - 1) / 2) * 0.26, [S.loc[h, m] for m, _ in mets], width=0.24, color=COH[h],
               label=f"{h.replace('_', '-')} held out")
    ax.set_xticks(x)
    ax.set_xticklabels([l.replace("\n", " ") for _, l in mets], fontsize=4.4, rotation=40, ha="right",
                       rotation_mode="anchor")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("recovered / median r", fontsize=5)
    despine(ax)
    ax.legend(loc="upper left", fontsize=4.3, handlelength=0.8, bbox_to_anchor=(0.0, -0.42), ncol=1)
    txt = []
    for h, keys in (("TCGA", ["risk_RFS_train_CLCA_c_index"]), ("CLCA", ["risk_RFS_train_TCGA_c_index", "risk_OS_train_TCGA_c_index"])):
        if h in S.index:
            vals = [f"{k.split('_')[1]} {S.loc[h, k]:.2f}" for k in keys if k in S.columns and pd.notna(S.loc[h, k])]
            txt.append(f"{h}: {', '.join(vals)}")
    ax.text(0.0, -0.78, "C-index, network + model\nwithout the test cohort:\n" + "\n".join(txt), transform=ax.transAxes,
            fontsize=4.4, va="top", color=Q.INK)
    ax.set_title("Leave-one-out", fontsize=6.2)
    return S


def figure2(D, P, outdir, log, m07c):
    fig = plt.figure(figsize=(W, 240 * MM))
    gs = gridspec.GridSpec(4, 1, height_ratios=[50, 62, 42, 62], hspace=0.5, left=0.035, right=0.985, top=0.975,
                           bottom=0.035)
    r1 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0], width_ratios=[1, 1.45], wspace=0.12)
    r1a = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=r1[0], width_ratios=[0.62, 1], wspace=0.0)
    ax = fig.add_subplot(r1a[1])
    w = f2a_weights(ax, D)
    fig.text(0.005, 0.978, "a", fontsize=8, fontweight="bold", va="bottom")
    r1b = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=r1[1], width_ratios=[0.08, 1], wspace=0.0)
    axb = f2b_classes(fig, r1b[1], D)
    letter(axb, "b", x=-0.35)
    axc = f2c_states(fig, gs[1], D, m07c, P)
    letter(axc, "c", x=-0.1)
    r3 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], width_ratios=[0.045, 1], wspace=0.0)
    axd, axe = f2d_e_km(fig, r3[1], D, P)
    letter(axd, "d", x=-0.3)
    letter(axe, "e", x=-0.1)
    r4 = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[3], width_ratios=[1.45, 1.12, 0.78], wspace=0.3)
    axf, fo = f2f_forest(fig, r4[0], D)
    letter(axf, "f", x=-0.6)
    r4g = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=r4[1], width_ratios=[0.5, 1], wspace=0.0)
    ax = fig.add_subplot(r4g[1])
    H = f2g_benchmark(ax, D)
    letter(ax, "g", x=-0.55)
    r4h = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=r4[2], height_ratios=[1, 0.9], hspace=0.0)
    ax = fig.add_subplot(r4h[0])
    loco = f2h_loco(ax, D)
    letter(ax, "h", x=-0.08, y=1.02)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"figure2.{ext}"))
    plt.close(fig)
    with pd.ExcelWriter(os.path.join(outdir, "figure2_source_data.xlsx")) as xw:
        w.to_frame("weight").to_excel(xw, "2a_weights")
        fo.to_excel(xw, "2f_forest", index=False)
        H.to_excel(xw, "2g_benchmark", index=False)
        if loco is not None:
            loco.to_excel(xw, "2h_loco")
    log.info("Figure 2 written")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    ap.add_argument("--only", default="1,2")
    args = ap.parse_args()
    P = load_params(args.params)
    outdir = p(os.path.join(P["paths"]["results"], "09_figures"))
    log = setup_logging(outdir, "09_publication_figures")
    D = load_all(P, log)
    m07c = _module("m07c", "07c_figures.py")
    if "1" in args.only:
        figure1(D, P, outdir, log)
    if "2" in args.only:
        figure2(D, P, outdir, log, m07c)


if __name__ == "__main__":
    main()
