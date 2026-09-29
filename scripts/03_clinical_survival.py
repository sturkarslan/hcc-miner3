#!/usr/bin/env python
"""Step 03 (clinical part): survival tables.

CLCA 2024: built from Nature 2024 Supplementary Table 1a (41586_2024_7054_MOESM3_ESM.csv).
  The table gives operation date (day), recurrence and death status, recurrence and
  death dates (month only) and RFS days for recurrences. It has NO last-follow-up date,
  so patients without an event are censored at an assumed administrative cutoff
  (config clinical.CLCA.cutoff_date; last recorded events are in 2021-09).
  If the cBioPortal patient file is present, the cutoff implied by its censored RFS
  times is reported, and with use_cbioportal_censoring: true those censoring times
  replace the cutoff-based ones.

Outputs (results/03_genomics_clinical/):
  survival_<cohort>.tsv            patient, OS/RFS time (days) and event, time source, clinical covariates
  survival_<cohort>_<EP>_miner.csv        duration, observed (MINER3 survival format), indexed by
                                          expression sample ID, network samples with data only
  survival_<cohort>_<EP>_h<N>m_miner.csv  same, administratively censored at survival.horizon_months
  qc/s1_kaplan_meier.png
"""

import argparse
import os

import numpy as np
import pandas as pd

from hcc_common import load_params, p, setup_logging

OUT = "03_genomics_clinical"


def _col(df, prefix):
    hits = [c for c in df.columns if str(c).strip().startswith(prefix)]
    if len(hits) != 1:
        raise KeyError(f"expected one column starting with {prefix!r}, found {hits}. Columns: {list(df.columns)}")
    return hits[0]


def _num(s):
    return pd.to_numeric(s.astype(str).str.strip().replace({"/": np.nan, "NA": np.nan, "": np.nan}),
                         errors="coerce")


def _month_mid(s):
    d = pd.to_datetime(s.astype(str).str.strip(), format="%Y.%m", errors="coerce")
    return d + pd.Timedelta(days=14)


def read_clca_supp(path):
    d = pd.read_csv(path, skiprows=2, encoding="utf-8-sig", dtype=str)
    d = d.loc[:, ~d.columns.str.startswith("Unnamed")]
    d = d[d["CLCA_ID"].str.match(r"CLCA_\d+", na=False)].set_index("CLCA_ID")
    return d


def build_clca(cfg, log):
    d = read_clca_supp(p(cfg["supplementary"]))
    log.info("CLCA supplementary table: %d patients", len(d))
    op = pd.to_datetime(d[_col(d, "Date of Operation")].str.strip(), format="%Y.%m.%d", errors="coerce")
    rec = _num(d[_col(d, "Recurrence (")])
    rec_date = _month_mid(d[_col(d, "Date of Recurrence")])
    dead = _num(d[_col(d, "Status (")])
    death_date = _month_mid(d[_col(d, "Date of Death")])
    rfs_days = _num(d[_col(d, "Relapse-Free Survival")])
    cutoff = pd.Timestamp(cfg["cutoff_date"])
    last_event = pd.concat([rec_date, death_date]).max()
    log.info("Latest recorded recurrence/death month: %s; censoring cutoff: %s", last_event.date(), cutoff.date())
    if last_event > cutoff:
        raise ValueError("cutoff_date is before the latest recorded event")

    out = pd.DataFrame(index=d.index)
    out.index.name = "patient"
    out["cohort"] = "CLCA"
    # ---- OS
    out["OS_event"] = dead
    os_t = pd.Series(np.nan, index=d.index)
    ev = dead == 1
    os_t[ev] = (death_date[ev] - op[ev]).dt.days
    os_t[dead == 0] = (cutoff - op[dead == 0]).dt.days
    out["OS_time"] = os_t.clip(lower=1)
    out["OS_time_source"] = np.select([ev & os_t.notna(), dead == 0], ["death_month_mid", "censored_at_cutoff"], "missing")
    # ---- RFS
    out["RFS_event"] = rec
    rfs_t = pd.Series(np.nan, index=d.index)
    er = rec == 1
    rfs_t[er] = rfs_days[er].fillna((rec_date[er] - op[er]).dt.days)
    rfs_t[rec == 0] = (cutoff - op[rec == 0]).dt.days
    out["RFS_time"] = rfs_t.clip(lower=1)
    out["RFS_time_source"] = np.select([er & rfs_days.notna(), er & rfs_t.notna(), rec == 0],
                                       ["table_days", "recurrence_month_mid", "censored_at_cutoff"], "missing")

    for ep in ("OS", "RFS"):
        log.info("CLCA %s: %s", ep, out[f"{ep}_time_source"].value_counts().to_dict())
    bad = (out["OS_event"] == 1) & (out["RFS_event"] == 1) & (out["OS_time"] < out["RFS_time"])
    if bad.any():
        log.warning("CLCA: %d patients with death before recurrence (month rounding); OS set to RFS time", bad.sum())
        out.loc[bad, "OS_time"] = out.loc[bad, "RFS_time"]

    # ---- optional check against cBioPortal censoring times
    cb = p(cfg.get("cbioportal_patient"))
    if cb and os.path.exists(cb):
        c = pd.read_csv(cb, sep="\t", comment="#", dtype=str)
        idc, tc, sc = cfg["cbioportal_id_column"], cfg["cbioportal_rfs_time_column"], cfg["cbioportal_rfs_status_column"]
        missing = [x for x in (idc, tc, sc) if x not in c.columns]
        if missing:
            raise KeyError(f"cBioPortal patient file: missing {missing}. Columns: {list(c.columns)}")
        c = c.set_index(idc)
        common = out.index.intersection(c.index)
        log.info("cBioPortal patient file: %d patients, %d match the supplementary IDs", len(c), len(common))
        t = pd.to_numeric(c.loc[common, tc], errors="coerce") * float(cfg.get("cbioportal_rfs_time_to_days", 1))
        censored = c.loc[common, sc].astype(str).str.startswith("0")
        implied = (op[common] + pd.to_timedelta(t, unit="D"))[censored]
        if len(implied.dropna()):
            log.info("Cutoff implied by cBioPortal censored RFS times: median %s, range %s to %s (n=%d)",
                     implied.median().date(), implied.min().date(), implied.max().date(), implied.notna().sum())
        if cfg.get("use_cbioportal_censoring"):
            idx = censored[censored].index.intersection(out.index[out["RFS_event"] == 0])
            out.loc[idx, "RFS_time"] = t[idx]
            out.loc[idx, "RFS_time_source"] = "cbioportal_censored"
            log.info("Replaced %d RFS censoring times with cBioPortal values", len(idx))
    else:
        log.warning("cBioPortal patient file not found (%s); censoring uses the assumed cutoff only", cb)

    # clinical covariates for later stratification
    cov = {"age": "Age", "sex_male": "Gender", "HBV": "HBV (", "HCV": "HCV (", "Edmondson": "Edmondson",
           "MVI": "MVI", "BCLC": "BCLC", "AFP_ge20": "AFP (>=", "cirrhosis_fibrosis": "Cirrhosis/Fibrosis"}
    for new, pre in cov.items():
        try:
            out[new] = d[_col(d, pre)].str.strip()
        except KeyError as e:
            log.warning("CLCA covariate %s not found (%s)", new, str(e).split(".")[0])
    return out


def build_tcga(cfg, log):
    """TCGA-LIHC from the cBioPortal PanCan patient file (times from TCGA-CDR, in months).
    RFS endpoint = PFS (CDR progression-free interval, recommended for LIHC)."""
    c = pd.read_csv(p(cfg["cbioportal_patient"]), sep="\t", comment="#", dtype=str).set_index("PATIENT_ID")
    out = pd.DataFrame(index=c.index)
    out.index.name = "patient"
    out["cohort"] = "TCGA"
    m2d = float(cfg.get("months_to_days", 30.44))
    for ep, (tcol, scol) in {"OS": cfg["os_columns"], "RFS": cfg["rfs_columns"]}.items():
        missing = [x for x in (tcol, scol) if x not in c.columns]
        if missing:
            raise KeyError(f"TCGA patient file: missing {missing}. Columns: {list(c.columns)}")
        out[f"{ep}_event"] = pd.to_numeric(c[scol].str.split(":").str[0], errors="coerce")
        out[f"{ep}_time"] = (pd.to_numeric(c[tcol], errors="coerce") * m2d).clip(lower=1)
        out[f"{ep}_time_source"] = np.where(out[f"{ep}_time"].notna(), f"cbioportal_{tcol}", "missing")
        log.info("TCGA %s (%s/%s): %d with time, %d events", ep, tcol, scol,
                 out[f"{ep}_time"].notna().sum(), int(out.loc[out[f"{ep}_time"].notna(), f"{ep}_event"].sum()))
    for new, col in (cfg.get("covariates") or {}).items():
        if col in c.columns:
            out[new] = c[col]
    return out


BUILDERS = {"CLCA": build_clca, "TCGA": build_tcga}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=None)
    args = ap.parse_args()
    P = load_params(args.params)
    outdir = p(os.path.join(P["paths"]["results"], OUT))
    log = setup_logging(outdir, "03_clinical_survival")

    expr_samples = None
    sp = p(os.path.join(P["paths"]["results"], "01_harmonized", "samples.tsv"))
    if os.path.exists(sp):
        expr_samples = pd.read_csv(sp, sep="\t")

    horizon_months = (P.get("survival") or {}).get("horizon_months")
    horizon = round(horizon_months * 30.44) if horizon_months else None
    if horizon:
        log.info("Common horizon: %s months (%d days); events after it are censored in the _h files",
                 horizon_months, horizon)

    tables = []
    for name, cfg in (P.get("clinical") or {}).items():
        if name not in BUILDERS:
            log.warning("no survival builder for %s yet; skipped", name)
            continue
        try:
            s = BUILDERS[name](cfg, log)
        except FileNotFoundError as e:
            log.warning("%s: input file missing (%s); skipped", name, e)
            continue
        s.to_csv(os.path.join(outdir, f"survival_{name}.tsv"), sep="\t")
        # MINER3 survival files: one per cohort and endpoint (GuanRank must be computed within a
        # cohort), indexed by expression sample ID, network samples only. Full length and truncated
        # at the common horizon (administrative censoring), which is the version for cross-cohort work.
        es = expr_samples[expr_samples["cohort"] == name] if expr_samples is not None else None
        for ep in ("OS", "RFS"):
            m = s[[f"{ep}_time", f"{ep}_event"]].dropna().astype({f"{ep}_event": int})
            m.columns = ["duration", "observed"]
            if es is not None:
                m = es.set_index("sample")[["patient"]].join(m, on="patient", how="inner").drop(columns="patient")
                m.index.name = "sample"
            versions = {"": m}
            if horizon:
                h = m.copy()
                over = h["duration"] > horizon
                h.loc[over, "duration"] = horizon
                h.loc[over, "observed"] = 0
                versions[f"_h{horizon_months}m"] = h
            for suffix, mm in versions.items():
                mm.to_csv(os.path.join(outdir, f"survival_{name}_{ep}{suffix}_miner.csv"))
                log.info("%s %s%s: %d %s, %d events, median follow-up of censored %.0f days",
                         name, ep, suffix, len(mm), "network samples" if es is not None else "patients",
                         mm["observed"].sum(), mm.loc[mm.observed == 0, "duration"].median())
        if expr_samples is not None:
            es = expr_samples[expr_samples["cohort"] == name]
            n = es["patient"].isin(s.index).sum()
            log.info("%s: %d of %d expression samples have a survival record", name, n, len(es))
            if n < len(es):
                log.warning("%s: unmatched expression IDs, e.g. %s", name,
                            es.loc[~es["patient"].isin(s.index), "patient"].head(5).tolist())
        tables.append(s)

    if tables:
        import qc_plots
        allsurv = pd.concat(tables)
        written = qc_plots.survival_report(outdir, allsurv, list(pd.unique(allsurv["cohort"])),
                                           endpoints=(("OS", "Overall survival"),
                                                      ("RFS", "Relapse-free (CLCA) / progression-free (TCGA) survival")),
                                           horizon_months=horizon_months)
        log.info("QC figures: %s", ", ".join(written))


if __name__ == "__main__":
    main()
