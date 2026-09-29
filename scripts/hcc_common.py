"""Shared helpers for the HCC MINER3 pipeline scripts."""

import logging
import os
import sys

import yaml


def project_root():
    """Project root = parent of scripts/, unless HCC_ROOT is set."""
    return os.environ.get("HCC_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load_params(path=None):
    path = path or os.path.join(project_root(), "config", "params.yaml")
    with open(path) as fh:
        return yaml.safe_load(fh)


def p(rel):
    """Resolve a project-relative path; absolute paths pass through."""
    if rel is None:
        return None
    return rel if os.path.isabs(rel) else os.path.join(project_root(), rel)


def setup_logging(outdir, name):
    os.makedirs(outdir, exist_ok=True)
    log = logging.getLogger(name)
    log.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S")
    for h in (logging.StreamHandler(sys.stdout), logging.FileHandler(os.path.join(outdir, f"{name}.log"), mode="w")):
        h.setFormatter(fmt)
        log.addHandler(h)
    return log


def require_columns(df, cols, what):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{what}: missing column(s) {missing}. Available: {list(df.columns)}. "
                       "Fix the column names in config/params.yaml.")


def miner_id_backmap(idmap_path, our_ids):
    """MINER's identifier conversion renames some Ensembl IDs to a different Preferred_Name
    (often an alternate-locus ID; ~500 of our GENCODE v36 genes). Return {miner_id: our_id}
    for those, so MINER outputs (regulons, modules, regulators) can be mapped back to the
    IDs used in results/01_harmonized/genes.tsv. IDs MINER leaves unchanged are not listed."""
    import pandas as pd
    m = pd.read_csv(idmap_path, sep="\t", dtype=str)
    e = m[(m["Source"] == "Ensembl Gene ID") & m["Name"].isin(set(our_ids)) & (m["Name"] != m["Preferred_Name"])]
    return dict(zip(e["Preferred_Name"], e["Name"]))
