"""QC figures for steps 01 (harmonization), 02 (batch correction) and 03 (survival).

Static PNGs written to <step results>/qc/. Colour follows the entity, in fixed order:
cohorts take categorical slots 1-3 (the only slots validated for all-pairs use in
scatter plots); correction methods take gray (uncorrected) + violet + green so they
never borrow a cohort colour; magnitudes use one blue ramp; correlations use
blue <-> gray <-> red. Labels with more than three levels are shown as small
multiples (one level highlighted per panel), never as more colours.
"""

import os

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# ------------------------------------------------------------------ palette (reference instance, light mode)
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
METHOD_COLORS = {"uncorrected": MUTED, "combat": "#4a3aa7", "cohort_z": "#008300"}
BACKGROUND = "#d9d8d2"  # de-emphasized points in highlight panels
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#104281"])
DIV = LinearSegmentedColormap.from_list("div_blue_red", ["#104281", "#3987e5", "#f0efec", "#e34948", "#9e2a2a"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.family": "sans-serif", "font.size": 9,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 2,
})


def cohort_colors(cohorts):
    """Fixed colour per cohort across all figures: slot order follows `cohorts:` in params.yaml,
    whatever subset or order a figure passes in."""
    try:
        from hcc_common import load_params
        order = list(load_params()["cohorts"])
    except Exception:
        order = []
    order += [c for c in cohorts if c not in order]
    if len(order) > 3:
        raise ValueError("more than 3 cohorts: facet instead of adding scatter colours")
    return {c: SLOTS[order.index(c)] for c in cohorts}


def _save(fig, outdir, name, written):
    os.makedirs(outdir, exist_ok=True)
    fig.savefig(os.path.join(outdir, name), dpi=150, bbox_inches="tight")
    plt.close(fig)
    written.append(name)


def _cohort_legend(ax, colors, counts=None, **kw):
    h = [Line2D([], [], marker="o", ls="", ms=6, mfc=c, mec=SURFACE,
                label=f"{k} (n={counts[k]})" if counts is not None else k) for k, c in colors.items()]
    ax.legend(handles=h, **kw)


def _ordered_samples(samples, cohorts, key=None):
    """Sample order: by cohort, then by `key` (a Series) within cohort."""
    s = samples.copy()
    s["_c"] = s["cohort"].map({c: i for i, c in enumerate(cohorts)})
    if key is not None:
        s["_k"] = key.reindex(s.index)
        s = s.sort_values(["_c", "_k"])
    else:
        s = s.sort_values("_c", kind="stable")
    return s.index


def _cohort_blocks(ax, samples_sorted_cohort, colors, axis="x"):
    """Thin cohort colour band + boundaries along one axis of a sample-ordered plot."""
    edges = np.flatnonzero(samples_sorted_cohort.values[1:] != samples_sorted_cohort.values[:-1]) + 1
    for e in edges:
        (ax.axvline if axis == "x" else ax.axhline)(e - 0.5, color=AXIS, lw=0.8)


def eta_squared(values, groups):
    """Fraction of variance explained by a categorical grouping (one-way ANOVA R^2).
    values: samples x k array; groups: length-n labels. Returns length-k array."""
    v = np.asarray(values, float)
    g = pd.Series(np.asarray(groups))
    ok = g.notna().values
    v, g = v[ok], g[ok].values
    tot = ((v - v.mean(0)) ** 2).sum(0)
    between = np.zeros(v.shape[1])
    for lev in pd.unique(g):
        sel = g == lev
        between += sel.sum() * (v[sel].mean(0) - v.mean(0)) ** 2
    return np.where(tot > 0, between / tot, np.nan)


def pca_scores(z, n=10):
    from sklearn.decomposition import PCA
    x = z.sub(z.mean(axis=1), axis=0).T.values
    n = min(n, min(x.shape) - 1)
    pca = PCA(n_components=n, random_state=0).fit(x)
    return pd.DataFrame(pca.transform(x), index=z.columns), pca.explained_variance_ratio_


def correlation_heatmap(ax, expr, samples, cohorts, colors, title, n_genes=2000):
    """Sample-sample Pearson correlation on the most variable genes (genes centered)."""
    top = expr.var(axis=1).sort_values(ascending=False).index[:n_genes]
    x = expr.loc[top]
    x = x.sub(x.mean(axis=1), axis=0)
    pc1 = pd.Series(np.linalg.svd(x.values, full_matrices=False)[2][0], index=x.columns)
    order = _ordered_samples(samples, cohorts, key=pc1)
    c = np.corrcoef(x[order].T.values)
    im = ax.imshow(c, cmap=DIV, vmin=-1, vmax=1, interpolation="nearest", aspect="auto")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    coh = samples.loc[order, "cohort"]
    # cohort band along the top and left edges
    band = np.array([[matplotlib.colors.to_rgb(colors[k]) for k in coh]])
    ax.imshow(band, extent=(-0.5, len(order) - 0.5, -0.02 * len(order) - 0.5, -0.5), aspect="auto", clip_on=False)
    ax.imshow(band.transpose(1, 0, 2), extent=(-0.02 * len(order) - 0.5, -0.5, len(order) - 0.5, -0.5),
              aspect="auto", clip_on=False)
    ax.set_xlim(-0.02 * len(order) - 0.5, len(order) - 0.5)
    ax.set_ylim(len(order) - 0.5, -0.02 * len(order) - 0.5)
    _cohort_blocks(ax, coh, colors, "x")
    _cohort_blocks(ax, coh, colors, "y")
    ax.set_title(title)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


# ================================================================== step 01

def harmonization_report(outdir, H, cohorts):
    qc = os.path.join(outdir, "qc")
    written = []
    colors = cohort_colors(cohorts)
    summary = pd.read_csv(os.path.join(outdir, "summary.tsv"), sep="\t").set_index("cohort").loc[cohorts]
    samples = pd.read_csv(os.path.join(outdir, "samples.tsv"), sep="\t", index_col="sample")
    sstats = pd.read_csv(os.path.join(outdir, "sample_stats.tsv"), sep="\t", index_col="sample")
    genes = pd.read_csv(os.path.join(outdir, "genes.tsv"), sep="\t", index_col="ensembl")
    expr = pd.read_csv(os.path.join(outdir, "expression_log2tpm1.csv"), index_col=0)
    counts = samples["cohort"].value_counts()

    # ---- h1: sample and gene flow
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 3.2), gridspec_kw={"width_ratios": [1, 1.3]}, layout="constrained")
    y = np.arange(len(cohorts))[::-1]
    for yi, c in zip(y, cohorts):
        n0, n1 = summary.loc[c, "input_samples"], summary.loc[c, "samples_kept"]
        a.plot([n1, n0], [yi, yi], color=AXIS, lw=2, zorder=1)
        a.scatter([n0], [yi], s=60, color=SURFACE, edgecolor=colors[c], lw=2, zorder=2)
        a.scatter([n1], [yi], s=60, color=colors[c], edgecolor=SURFACE, lw=1.5, zorder=3)
        a.annotate(f"{n1} kept of {n0}", (max(n0, n1), yi), xytext=(8, 0), textcoords="offset points",
                   va="center", color=INK2, fontsize=8)
    a.set_yticks(y, cohorts)
    a.set_xlim(0, summary["input_samples"].max() * 1.45)
    a.set_xlabel("samples")
    a.grid(axis="y", visible=False)
    a.set_title("Samples: input (open) → primary HCC kept (filled)")
    stages = ["input_rows", "ensembl_genes", "shared_genes", "genes_after_filter"]
    labels = ["input rows", "mapped to\nEnsembl", "shared by\nall cohorts", "pass detection\nfilter"]
    for c in cohorts:
        v = summary.loc[c, stages].values.astype(float)
        b.plot(range(4), v, color=colors[c], marker="o", ms=6, mec=SURFACE,
               label=f"{c}: {int(v[0]):,} input rows, {int(v[1]):,} Ensembl genes")
    for i, s in enumerate(stages[2:], start=2):
        b.annotate(f"{int(summary[s].iloc[0]):,}", (i, summary[s].iloc[0]), xytext=(0, -14),
                   textcoords="offset points", ha="center", color=INK2, fontsize=8)
    b.set_xticks(range(4), labels)
    b.set_xlim(-0.3, 3.3)
    b.set_ylabel("genes")
    b.yaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter("{x:,.0f}"))
    b.grid(axis="x", visible=False)
    b.legend(loc="lower left")
    b.set_title("Genes at each harmonization step")
    _save(fig, qc, "h1_sample_gene_flow.png", written)

    # ---- h2: identifier mapping outcome
    tiers = ["ensembl", "hgnc_symbol", "miner_gene_name", "hgnc_prev", "alias", "ambiguous", "unmapped"]
    tier_lab = ["Ensembl ID", "HGNC symbol", "MINER gene name", "HGNC previous", "alias",
                "ambiguous (dropped)", "unmapped (dropped)"]
    tier_col = SLOTS[:5] + ["#b5b3ab", "#6f6d68"]
    frac = {}
    for c in cohorts:
        mp = pd.read_csv(os.path.join(outdir, f"id_mapping_{c}.tsv"), sep="\t")
        frac[c] = mp["status"].str.split(":").str[0].value_counts().reindex(tiers, fill_value=0)
    frac = pd.DataFrame(frac).T
    fig, ax = plt.subplots(figsize=(9, 0.6 * len(cohorts) + 1.6))
    left = np.zeros(len(cohorts))
    pct = frac.div(frac.sum(axis=1), axis=0) * 100
    for t, lab, col in zip(tiers, tier_lab, tier_col):
        ax.barh(y, pct[t].values, left=left, color=col, edgecolor=SURFACE, lw=2, height=0.6, label=lab)
        for yi, v, l0, n in zip(y, pct[t].values, left, frac[t].values):
            if v >= 6:
                ax.text(l0 + v / 2, yi, f"{n:,}", ha="center", va="center", fontsize=7.5, color=SURFACE)
        left += pct[t].values
    ax.set_yticks(y, cohorts)
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of input identifiers (numbers are identifier counts)")
    ax.grid(axis="y", visible=False)
    ax.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.02), fontsize=7.5)
    ax.set_title("Identifier → Ensembl mapping, by first matching tier", pad=38)
    _save(fig, qc, "h2_id_mapping.png", written)

    # ---- h3: per-sample distributions
    med, q1, q3 = expr.median(), expr.quantile(0.25), expr.quantile(0.75)
    order = _ordered_samples(samples, cohorts, key=med)
    xs = np.arange(len(order))
    coh = samples.loc[order, "cohort"]
    fig, axes = plt.subplots(3, 1, figsize=(11, 7), sharex=True)
    for c in cohorts:
        sel = (coh == c).values
        axes[0].vlines(xs[sel], q1[order][sel], q3[order][sel], color=colors[c], lw=1, alpha=0.6)
        axes[0].scatter(xs[sel], med[order][sel], s=6, color=colors[c], lw=0, label=c)
        axes[1].scatter(xs[sel], sstats.loc[order, "n_genes_detected"][sel], s=6, color=colors[c], lw=0)
        axes[2].scatter(xs[sel], 100 * sstats.loc[order, "tpm_frac_shared_genes"][sel], s=6, color=colors[c], lw=0)
    axes[0].set_ylabel("log2(TPM+1)\nmedian and IQR")
    axes[1].set_ylabel(f"genes with\nTPM ≥ {H['min_tpm']:g}")
    axes[2].set_ylabel("% of mapped TPM\nin shared genes")
    axes[2].set_xlabel("samples (ordered by cohort, then median)")
    for ax in axes:
        _cohort_blocks(ax, coh, colors)
        ax.grid(axis="x", visible=False)
    _cohort_legend(axes[0], colors, counts, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.0))
    axes[0].set_title("Per-sample expression, detection and renormalization loss", pad=24)
    axes[2].set_xticks([])
    fig.tight_layout()
    _save(fig, qc, "h3_per_sample.png", written)

    # ---- h4: detection fraction per gene, with the filter threshold
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    bins = np.linspace(0, 1, 41)
    for c in cohorts:
        h, _ = np.histogram(genes[f"frac_detected_{c}"], bins=bins)
        ax.stairs(h, bins, color=colors[c], lw=2, label=c)
    ax.axvline(H["min_frac"], color=INK2, lw=1, ls="--")
    ax.annotate(f"keep if ≥ {H['min_frac']:.0%} in every cohort\n{int(genes['kept'].sum()):,} of {len(genes):,} genes",
                (H["min_frac"], 0.55), xycoords=("data", "axes fraction"), xytext=(6, 0),
                textcoords="offset points", va="center", color=INK2, fontsize=8)
    ax.set_xlabel(f"fraction of samples with TPM ≥ {H['min_tpm']:g}")
    ax.set_ylabel("genes (shared by all cohorts)")
    ax.set_yscale("log")
    ax.legend(loc="upper center", ncol=3)
    ax.set_title("Gene detection by cohort (log scale)")
    _save(fig, qc, "h4_gene_detection.png", written)

    # ---- h5: per-gene mean expression, cohort vs cohort (gene-specific cohort effects)
    means = pd.DataFrame({c: expr.loc[:, samples.loc[expr.columns, "cohort"] == c].mean(axis=1) for c in cohorts})
    pairs = [(cohorts[i], cohorts[j]) for i in range(len(cohorts)) for j in range(i + 1, len(cohorts))]
    fig, axes = plt.subplots(1, len(pairs), figsize=(3.6 * len(pairs) + 0.8, 3.9), squeeze=False,
                             layout="constrained")
    lim = (0, float(np.ceil(means.values.max())))
    for ax, (u, v) in zip(axes[0], pairs):
        hb = ax.hexbin(means[u], means[v], gridsize=45, cmap=SEQ, mincnt=1, bins="log", extent=lim + lim, lw=0)
        ax.plot(lim, lim, color=INK2, lw=1, ls="--")
        d = means[v] - means[u]
        ax.text(0.03, 0.97, f"r = {np.corrcoef(means[u], means[v])[0, 1]:.3f}\n"
                            f"|Δ| > 1: {(d.abs() > 1).mean():.1%} of genes",
                transform=ax.transAxes, va="top", fontsize=8, color=INK2)
        ax.set_xlabel(f"{u} mean log2(TPM+1)")
        ax.set_ylabel(f"{v} mean log2(TPM+1)")
        ax.set_aspect("equal")
        ax.grid(False)
    fig.colorbar(hb, ax=axes[0].tolist(), shrink=0.8, label="genes per bin (log)")
    fig.suptitle("Per-gene mean expression between cohorts (before batch correction)", x=0.01, ha="left",
                 fontweight="bold", color=INK, fontsize=10)
    _save(fig, qc, "h5_gene_means_between_cohorts.png", written)

    # ---- h6: sample-sample correlation, uncorrected
    fig, ax = plt.subplots(figsize=(6.6, 6))
    im = correlation_heatmap(ax, expr, samples, cohorts, colors,
                             "Sample–sample correlation, 2,000 most variable genes (uncorrected)")
    fig.colorbar(im, ax=ax, shrink=0.7, label="Pearson r")
    _cohort_legend(ax, colors, counts, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.01))
    _save(fig, qc, "h6_sample_correlation.png", written)
    return written


# ================================================================== step 02

def batch_report(outdir, finals, corrected, samples, genes, qc_metrics, prog_scores, label_cols):
    """finals: {method: z matrix}; corrected: {method: matrix before final z};
    prog_scores: {method: samples x programs DataFrame}."""
    qc = os.path.join(outdir, "qc")
    written = []
    cohorts = list(pd.unique(samples["cohort"]))
    colors = cohort_colors(cohorts)
    counts = samples["cohort"].value_counts()
    methods = list(finals)
    batch = samples.loc[finals[methods[0]].columns, "cohort"]
    pcs = {m: pca_scores(z, 10) for m, z in finals.items()}

    # ---- b1: PCA coloured by cohort
    fig, axes = plt.subplots(1, len(methods), figsize=(4.2 * len(methods), 4.2), squeeze=False)
    for ax, m in zip(axes[0], methods):
        X, ev = pcs[m]
        for c in cohorts:
            sel = (batch == c).values
            ax.scatter(X.values[sel, 0], X.values[sel, 1], s=12, color=colors[c], edgecolor=SURFACE, lw=0.4,
                       alpha=0.85)
        r = qc_metrics.set_index("matrix").loc[m]
        ax.set_title(m)
        ax.text(0.02, 0.02, f"cohort silhouette {r['silhouette_cohort']:.2f}\nkNN mixing {r['knn_mixing_cohort']:.2f}",
                transform=ax.transAxes, fontsize=8, color=INK2,
                bbox=dict(boxstyle="round,pad=0.3", fc=SURFACE, ec=GRID, alpha=0.9))
        ax.set_xlabel(f"PC1 ({ev[0]:.1%})")
        ax.set_ylabel(f"PC2 ({ev[1]:.1%})")
    _cohort_legend(axes[0][0], colors, counts, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.08))
    fig.tight_layout()
    _save(fig, qc, "b1_pca_by_cohort.png", written)

    # ---- b2: PCA small multiples per label level (labelled samples only)
    for lab in label_cols:
        has = samples[lab].notna()
        levels = samples.loc[has, lab].astype(str).value_counts()
        if len(levels) < 2:
            continue
        if len(levels) > 8:
            keep = levels.index[:7]
            lab_vals = samples.loc[has, lab].astype(str).where(lambda s: s.isin(keep), "Other")
        else:
            lab_vals = samples.loc[has, lab].astype(str)
        levs = sorted(lab_vals.unique(), key=lambda v: (v == "Other", v))
        fig, axes = plt.subplots(len(methods), len(levs), figsize=(2.3 * len(levs), 2.3 * len(methods)),
                                 squeeze=False)
        for i, m in enumerate(methods):
            X = pcs[m][0].loc[lab_vals.index]
            for j, lv in enumerate(levs):
                ax = axes[i][j]
                sel = (lab_vals == lv).values
                ax.scatter(X.values[~sel, 0], X.values[~sel, 1], s=6, color=BACKGROUND, lw=0)
                ax.scatter(X.values[sel, 0], X.values[sel, 1], s=12, color=SLOTS[0], edgecolor=SURFACE, lw=0.4)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.grid(False)
                if i == 0:
                    ax.set_title(f"{lv} (n={sel.sum()})", fontsize=9)
                if j == 0:
                    ax.set_ylabel(m, color=INK, fontweight="bold")
        fig.suptitle(f"{lab}: PC1 vs PC2, labelled samples only; one level highlighted per panel",
                     x=0.01, ha="left", fontweight="bold", color=INK, fontsize=10)
        fig.tight_layout()
        _save(fig, qc, f"b2_pca_by_{lab}.png", written)

    # ---- b3: variance of each PC explained by cohort and by labels
    vars_ = ["cohort"] + [lab for lab in label_cols if samples[lab].notna().sum() > 0]
    fig, axes = plt.subplots(1, len(methods), figsize=(1.0 + 0.75 * len(vars_) * len(methods), 4.2),
                             squeeze=False, sharey=True)
    for ax, m in zip(axes[0], methods):
        X, ev = pcs[m]
        R = np.array([eta_squared(X.values, samples.loc[X.index, v].values) for v in vars_]).T
        im = ax.imshow(R, cmap=SEQ, vmin=0, vmax=1, aspect="auto")
        for (r, c), val in np.ndenumerate(R):
            if not np.isnan(val):
                ax.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color=SURFACE if val > 0.5 else INK2)
        ax.set_xticks(range(len(vars_)), vars_, rotation=40, ha="right")
        ax.set_yticks(range(R.shape[0]), [f"PC{k + 1} ({e:.0%})" for k, e in enumerate(ev)])
        ax.grid(False)
        ax.set_title(m)
    fig.colorbar(im, ax=axes[0].tolist(), shrink=0.8, label="variance explained (η²)")
    fig.suptitle("Principal components explained by cohort and by LICA-FR labels (labels: labelled samples only)",
                 x=0.01, ha="left", fontweight="bold", color=INK, fontsize=10)
    _save(fig, qc, "b3_pc_association.png", written)

    # ---- b4: per-gene variance explained by cohort (ECDF)
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    for m in methods:
        r2 = np.sort(eta_squared(finals[m].T.values, batch.values))
        r2 = r2[~np.isnan(r2)]
        yv = np.arange(1, len(r2) + 1) / len(r2)
        ax.plot(r2, yv, color=METHOD_COLORS.get(m, INK2),
                label=f"{m}: median {np.median(r2):.3f}, 90th pct {np.quantile(r2, 0.9):.3f}")
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_xlim(0, 1)
    ax.set_xlabel("per-gene variance explained by cohort (η²)")
    ax.set_ylabel("cumulative fraction of genes")
    ax.legend(loc="lower right")
    ax.set_title("How much of each gene's variance is cohort")
    _save(fig, qc, "b4_gene_cohort_variance.png", written)

    # ---- b5: summary metrics, one panel per metric
    qm = qc_metrics.set_index("matrix")
    metric_cols = [c for c in qm.columns if qm[c].notna().any()]
    ncol = 3
    nrow = int(np.ceil(len(metric_cols) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, (0.35 * len(methods) + 1.0) * nrow),
                             squeeze=False)
    yy = np.arange(len(methods))[::-1]
    for ax, col in zip(axes.flat, metric_cols):
        vals = qm.loc[methods, col].values
        for yi, m, v in zip(yy, methods, vals):
            if np.isnan(v):
                continue
            ax.scatter([v], [yi], s=50, color=METHOD_COLORS.get(m, INK2), edgecolor=SURFACE, zorder=3)
            ax.annotate(f"{v:.2f}", (v, yi), xytext=(6, 0), textcoords="offset points", va="center",
                        fontsize=7.5, color=INK2)
        ax.set_yticks(yy, methods)
        ax.set_title(col.replace("_", " "), fontsize=8.5)
        ax.axvline(0, color=AXIS, lw=0.8)
        ax.grid(axis="y", visible=False)
        # fixed scales so panels are comparable and near-zero values are not blown up
        if col.startswith("silhouette"):
            ax.set_xlim(min(-0.25, np.nanmin(vals) - 0.05), max(0.75, np.nanmax(vals) + 0.25))
        else:
            ax.set_xlim(0, max(1.25, np.nanmax(vals) + 0.25))
    for ax in list(axes.flat)[len(metric_cols):]:
        ax.set_visible(False)
    fig.suptitle("Batch-correction QC metrics. Cohort silhouette: lower is better; kNN mixing: 1 = well mixed; "
                 "label silhouettes and within-cohort ρ: should stay near uncorrected",
                 x=0.01, ha="left", fontweight="bold", color=INK, fontsize=9)
    fig.tight_layout()
    _save(fig, qc, "b5_qc_metrics.png", written)

    # ---- b6: marker-program scores per cohort
    progs = list(prog_scores[methods[0]].columns)
    if progs:
        fig, axes = plt.subplots(len(progs), len(methods), figsize=(3.4 * len(methods), 2.3 * len(progs)),
                                 squeeze=False, sharey="row")
        for i, pr in enumerate(progs):
            for j, m in enumerate(methods):
                ax = axes[i][j]
                data = [prog_scores[m].loc[batch.index[batch == c], pr].values for c in cohorts]
                bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                                medianprops=dict(color=SURFACE, lw=1.5), whiskerprops=dict(color=AXIS),
                                capprops=dict(color=AXIS))
                for patch, c in zip(bp["boxes"], cohorts):
                    patch.set(facecolor=colors[c], edgecolor=colors[c])
                ax.set_xticks(range(1, len(cohorts) + 1), cohorts)
                ax.grid(axis="x", visible=False)
                if i == 0:
                    ax.set_title(m)
                if j == 0:
                    ax.set_ylabel(f"{pr}\nmean z", color=INK, fontweight="bold")
        fig.suptitle("Marker-program scores by cohort", x=0.01, ha="left", fontweight="bold", color=INK, fontsize=10)
        fig.tight_layout()
        _save(fig, qc, "b6_programs_by_cohort.png", written)

        # ---- b7: programs by LICA-FR label (should be unchanged by correction)
        for lab in label_cols:
            lv = samples[lab].dropna().astype(str)
            groups = [g for g, n in lv.value_counts().sort_index().items() if n >= 3]
            if len(groups) < 2:
                continue
            fig, axes = plt.subplots(len(progs), len(methods),
                                     figsize=(max(3.4, 0.45 * len(groups) + 1.2) * len(methods), 2.3 * len(progs)),
                                     squeeze=False, sharey="row")
            for i, pr in enumerate(progs):
                for j, m in enumerate(methods):
                    ax = axes[i][j]
                    data = [prog_scores[m].loc[lv.index[lv == g], pr].values for g in groups]
                    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                                    medianprops=dict(color=SURFACE, lw=1.5), whiskerprops=dict(color=AXIS),
                                    capprops=dict(color=AXIS))
                    for patch in bp["boxes"]:
                        patch.set(facecolor=SLOTS[0], edgecolor=SLOTS[0])
                    ax.set_xticks(range(1, len(groups) + 1), groups, rotation=30 if len(groups) > 4 else 0)
                    ax.grid(axis="x", visible=False)
                    if i == 0:
                        ax.set_title(m)
                    if j == 0:
                        ax.set_ylabel(f"{pr}\nmean z", color=INK, fontweight="bold")
            fig.suptitle(f"Marker-program scores by {lab} (labelled samples only)", x=0.01, ha="left",
                         fontweight="bold", color=INK, fontsize=10)
            fig.tight_layout()
            _save(fig, qc, f"b7_programs_by_{lab}.png", written)

    # ---- b8: relative log expression per sample, before vs after ComBat
    rle_methods = [m for m in ("uncorrected", "combat") if m in corrected]
    fig, axes = plt.subplots(len(rle_methods), 1, figsize=(11, 2.4 * len(rle_methods)), sharex=True, sharey=True,
                             squeeze=False)
    order = _ordered_samples(samples.loc[batch.index], cohorts)
    coh = samples.loc[order, "cohort"]
    xs = np.arange(len(order))
    for ax, m in zip(axes[:, 0], rle_methods):
        rle = corrected[m].sub(corrected[m].median(axis=1), axis=0)[order]
        q1, md, q3 = rle.quantile(0.25), rle.median(), rle.quantile(0.75)
        for c in cohorts:
            sel = (coh == c).values
            ax.vlines(xs[sel], q1.values[sel], q3.values[sel], color=colors[c], lw=1, alpha=0.6)
            ax.scatter(xs[sel], md.values[sel], s=5, color=colors[c], lw=0)
        ax.axhline(0, color=INK2, lw=0.8)
        _cohort_blocks(ax, coh, colors)
        ax.set_ylabel(f"{m}\nRLE", color=INK, fontweight="bold")
        ax.grid(axis="x", visible=False)
    axes[-1, 0].set_xticks([])
    axes[-1, 0].set_xlabel("samples (ordered by cohort)")
    _cohort_legend(axes[0, 0], colors, counts, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.0))
    axes[0, 0].set_title("Relative log expression (gene-median-centered log2): median and IQR per sample", pad=24)
    fig.tight_layout()
    _save(fig, qc, "b8_rle.png", written)

    # ---- b9: sample correlation after each correction
    show = [m for m in methods if m != "uncorrected"]
    if show:
        fig, axes = plt.subplots(1, len(show), figsize=(6.2 * len(show), 6), squeeze=False)
        for ax, m in zip(axes[0], show):
            im = correlation_heatmap(ax, corrected[m], samples.loc[batch.index], cohorts, colors, m)
        fig.colorbar(im, ax=axes[0].tolist(), shrink=0.7, label="Pearson r")
        _cohort_legend(axes[0][0], colors, counts, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.01))
        fig.suptitle("Sample–sample correlation after correction (2,000 most variable genes)", x=0.01, ha="left",
                     fontweight="bold", color=INK, fontsize=10)
        _save(fig, qc, "b9_sample_correlation_corrected.png", written)
    return written


# ================================================================== step 03

def kaplan_meier(time, event):
    """Return step-function arrays (t, S) and the at-risk function."""
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    ok = ~np.isnan(t)
    t, e = t[ok], e[ok]
    times = np.unique(t[e == 1])
    s, ts, ss = 1.0, [0.0], [1.0]
    for u in times:
        n = (t >= u).sum()
        d = ((t == u) & (e == 1)).sum()
        s *= 1 - d / n
        ts.append(u)
        ss.append(s)
    return np.array(ts), np.array(ss), t, e


def survival_report(outdir, surv, cohorts, endpoints=(("OS", "Overall survival"), ("RFS", "Relapse-free survival")),
                    horizon_months=None):
    """surv: DataFrame with cohort, <EP>_time (days), <EP>_event columns.
    horizon_months: draw the common analysis horizon as a dashed line."""
    qc = os.path.join(outdir, "qc")
    written = []
    colors = cohort_colors(cohorts)
    fig, axes = plt.subplots(1, len(endpoints), figsize=(5.2 * len(endpoints), 4.4), squeeze=False)
    for ax, (ep, title) in zip(axes[0], endpoints):
        tcol, ecol = f"{ep}_time", f"{ep}_event"
        if tcol not in surv.columns:
            ax.set_visible(False)
            continue
        tmax = 0
        rows = []
        for c in cohorts:
            d = surv[(surv["cohort"] == c)].dropna(subset=[tcol, ecol])
            if d.empty:
                continue
            ts, ss, t, e = kaplan_meier(d[tcol] / 30.44, d[ecol])
            tmax = max(tmax, t.max())
            ax.step(np.append(ts, t.max()), np.append(ss, ss[-1]), where="post", color=colors[c],
                    label=f"{c} (n={len(t)}, events={int(e.sum())})")
            cens = t[e == 0]
            cs = [ss[np.searchsorted(ts, x, side="right") - 1] for x in cens]
            ax.scatter(cens, cs, marker="|", s=30, color=colors[c], lw=1)
            rows.append((c, t))
        ax.set_ylim(0, 1.02)
        if horizon_months:
            ax.axvline(horizon_months, color=MUTED, ls="--", lw=1)
            ax.text(horizon_months, 1.0, f" {horizon_months}-month horizon", color=MUTED, fontsize=7.5, va="top")
        ax.set_xlabel("months from surgery / diagnosis")
        ax.set_ylabel("survival probability")
        ax.set_title(title)
        ax.legend(loc="lower left")
        # number at risk
        ticks = np.arange(0, tmax + 1, 12 if tmax <= 72 else 24)
        ax.set_xticks(ticks)
        ax.text(0, -0.2, "number at risk", transform=ax.transAxes, fontsize=7.5, color=MUTED)
        for k, (c, t) in enumerate(rows):
            yk = -0.27 - 0.07 * k
            ax.text(-0.02, yk, c, transform=ax.transAxes, ha="right", va="top", fontsize=7.5, color=INK2)
            for x in ticks:
                ax.text(x, yk, str(int((t >= x).sum())), transform=ax.get_xaxis_transform(),
                        ha="center", va="top", fontsize=7.5, color=INK2)
    fig.tight_layout()
    _save(fig, qc, "s1_kaplan_meier.png", written)
    return written
