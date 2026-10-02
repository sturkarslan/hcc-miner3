"""Download and inventory the iCCA cohorts listed in config/icca_cohorts.yaml.

    python scripts/tools/fetch_icca_cohorts.py                    # everything that can be fetched automatically
    python scripts/tools/fetch_icca_cohorts.py --only GSE244807 TCGA_CHOL
    python scripts/tools/fetch_icca_cohorts.py --dry-run          # list what would be downloaded
    python scripts/tools/fetch_icca_cohorts.py --inventory-only   # re-inventory files already on disk
    python scripts/tools/fetch_icca_cohorts.py --fastq            # also download FASTQ from ENA (large)

Per source:
  geo   series matrix + every supplementary file (files above --max-gb are skipped with a warning),
        a parsed sample table (samples_geo.tsv: title, source, every "key: value" characteristic)
        and a FASTQ manifest from ENA (fastq_manifest.tsv). FASTQ is downloaded only with --fastq.
  tcga  UCSC Xena GDC hub files + a GDC API case table (gdc_cases.tsv) with diagnosis site, stage
        and OS fields, and an is_icca flag (tissue_or_organ_of_origin contains 'intrahepatic').
  node  cannot be fetched without a NODE login. Writes manual/README_DOWNLOAD.txt with what to get;
        files placed in manual/ are inventoried.

Outputs: <data>/iCCA/<key>/... and <results>/00_icca/{inventory.tsv, files.tsv, fetch_icca.log}.
Needs internet access: run on a login/transfer node, not a compute node. Re-running skips files
already downloaded with the expected size.
"""

import argparse
import gzip
import html.parser
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request

import pandas as pd
import yaml

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from hcc_common import p, setup_logging  # noqa: E402

GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
GDC_API = "https://api.gdc.cancer.gov/cases"
UA = {"User-Agent": "hcc-miner3-fetch/1.0"}
TABLE_EXT = (".tsv", ".txt", ".csv", ".tsv.gz", ".txt.gz", ".csv.gz")

log = None


# ---------------------------------------------------------------- http helpers
def _open(url, data=None, headers=None, timeout=120):
    req = urllib.request.Request(url, data=data, headers={**UA, **(headers or {})})
    return urllib.request.urlopen(req, timeout=timeout)


def get_text(url, data=None, headers=None, retries=4):
    for i in range(retries):
        try:
            with _open(url, data, headers) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise
            err = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            err = e
        time.sleep(2 ** (i + 1))
    raise err


def remote_size(url):
    try:
        req = urllib.request.Request(url, method="HEAD", headers=UA)
        with urllib.request.urlopen(req, timeout=60) as r:
            n = r.headers.get("Content-Length")
            return int(n) if n else None
    except Exception:
        return None


def download(url, dest, max_bytes=None, dry=False, retries=4):
    """Stream url to dest (via dest.part). Returns a status string."""
    size = remote_size(url)
    if os.path.exists(dest) and (size is None or os.path.getsize(dest) == size):
        return "present"
    if max_bytes and size and size > max_bytes:
        log.warning("skip %s: %.1f GB > --max-gb", os.path.basename(dest), size / 1e9)
        return f"skipped_size_{size / 1e9:.1f}GB"
    if dry:
        log.info("[dry-run] %s -> %s (%s)", url, dest, f"{size / 1e6:.1f} MB" if size else "size unknown")
        return "dry_run"
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + ".part"
    for i in range(retries):
        try:
            with _open(url, timeout=300) as r, open(part, "wb") as fh:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    fh.write(chunk)
            if size is not None and os.path.getsize(part) != size:
                raise IOError(f"size mismatch {os.path.getsize(part)} != {size}")
            os.replace(part, dest)
            log.info("downloaded %s (%.1f MB)", os.path.basename(dest), os.path.getsize(dest) / 1e6)
            return "downloaded"
        except urllib.error.HTTPError as e:
            if e.code == 404:
                log.warning("not found: %s", url)
                return "not_found"
            log.warning("attempt %d failed for %s: %s", i + 1, url, e)
        except Exception as e:
            log.warning("attempt %d failed for %s: %s", i + 1, url, e)
        time.sleep(2 ** (i + 1))
    return "failed"


class _Links(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            h = dict(attrs).get("href")
            if h and not h.startswith(("?", "/", "http")) and not h.endswith("/"):
                self.links.append(h)


def list_dir(url):
    try:
        lp = _Links()
        lp.feed(get_text(url))
        return sorted(set(lp.links))
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise


# ---------------------------------------------------------------- GEO
def geo_base(acc):
    return f"{GEO_FTP}/{acc[:-3]}nnn/{acc}"


def parse_series_matrix(path):
    """(series fields dict, samples DataFrame) from a GEO series_matrix.txt.gz."""
    series, rows = {}, {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!"):
                continue
            key, *vals = line.rstrip("\n").split("\t")
            vals = [v.strip('"') for v in vals]
            if key.startswith("!Series_"):
                series.setdefault(key[8:], []).extend(vals)
            elif key.startswith("!Sample_"):
                rows.setdefault(key[8:], []).append(vals)
    acc = rows.get("geo_accession", [[]])[0]
    df = pd.DataFrame(index=acc)
    for k in ("title", "source_name_ch1", "description", "molecule_ch1", "platform_id"):
        if k in rows:
            df[k] = rows[k][0]
    for vals in rows.get("characteristics_ch1", []):
        for s, v in zip(acc, vals):
            if ":" in v:
                ck, cv = v.split(":", 1)
                df.loc[s, "ch:" + ck.strip()] = cv.strip()
            elif v:
                df.loc[s, "ch:unlabelled"] = v
    for vals in rows.get("relation", []):
        for s, v in zip(acc, vals):
            m = re.search(r"(SRX\d+)", v)
            if m:
                df.loc[s, "srx"] = m.group(1)
    df.index.name = "gsm"
    return series, df


def ena_manifest(project):
    q = (f"{ENA_API}?accession={project}&result=read_run&format=tsv&fields="
         "run_accession,experiment_accession,sample_alias,sample_title,library_layout,"
         "read_count,fastq_ftp,fastq_bytes,fastq_md5")
    txt = get_text(q)
    from io import StringIO
    return pd.read_csv(StringIO(txt), sep="\t") if txt.strip() else pd.DataFrame()


def fetch_geo(key, spec, out, args, files):
    acc = spec["accession"]
    base = geo_base(acc)
    for sub in ("matrix", "suppl"):
        names = list_dir(f"{base}/{sub}/")
        if not names:
            log.warning("%s: no files listed under %s/", acc, sub)
        for n in names:
            st = download(f"{base}/{sub}/{n}", os.path.join(out, sub, n), args.max_bytes, args.dry_run)
            files.append({"cohort": key, "file": f"{sub}/{n}", "url": f"{base}/{sub}/{n}", "status": st})
    if args.dry_run:
        return
    mats = sorted(f for f in os.listdir(os.path.join(out, "matrix")) if f.endswith("series_matrix.txt.gz")) \
        if os.path.isdir(os.path.join(out, "matrix")) else []
    tabs, series = [], {}
    for f in mats:
        series, df = parse_series_matrix(os.path.join(out, "matrix", f))
        tabs.append(df)
    if tabs:
        samples = pd.concat(tabs)
        samples.to_csv(os.path.join(out, "samples_geo.tsv"), sep="\t")
        log.info("%s: %d samples; characteristics: %s", acc, len(samples),
                 ", ".join(c[3:] for c in samples.columns if c.startswith("ch:")))
    # FASTQ manifest via the BioProject / SRA study named in the series relations
    rel = " ".join(series.get("relation", []))
    proj = re.search(r"(PRJNA\d+)", rel) or re.search(r"(SRP\d+)", rel)
    if not proj:
        log.warning("%s: no BioProject/SRA study in series relations; no FASTQ manifest", acc)
        return
    try:
        man = ena_manifest(proj.group(1))
    except Exception as e:
        log.warning("%s: ENA manifest for %s failed: %s", acc, proj.group(1), e)
        return
    man.to_csv(os.path.join(out, "fastq_manifest.tsv"), sep="\t", index=False)
    gb = pd.to_numeric(man.get("fastq_bytes", pd.Series(dtype=str)).astype(str).str.split(";").explode(),
                       errors="coerce").sum() / 1e9
    log.info("%s: FASTQ manifest %s, %d runs, %.1f GB", acc, proj.group(1), len(man), gb)
    if args.fastq:
        for _, r in man.iterrows():
            for u in str(r.get("fastq_ftp", "")).split(";"):
                if u and u != "nan":
                    st = download("https://" + u, os.path.join(out, "fastq", os.path.basename(u)), None)
                    files.append({"cohort": key, "file": "fastq/" + os.path.basename(u), "url": u, "status": st})


# ---------------------------------------------------------------- TCGA
GDC_FIELDS = [
    "submitter_id", "project.project_id",
    "demographic.gender", "demographic.vital_status", "demographic.days_to_death", "demographic.age_at_index",
    "diagnoses.primary_diagnosis", "diagnoses.tissue_or_organ_of_origin", "diagnoses.site_of_resection_or_biopsy",
    "diagnoses.morphology", "diagnoses.ajcc_pathologic_stage", "diagnoses.days_to_last_follow_up",
    "diagnoses.classification_of_tumor",
]


def gdc_cases(project):
    body = json.dumps({
        "filters": {"op": "in", "content": {"field": "project.project_id", "value": [project]}},
        "fields": ",".join(GDC_FIELDS), "format": "JSON", "size": 2000,
    }).encode()
    hits = json.loads(get_text(GDC_API, data=body, headers={"Content-Type": "application/json"}))["data"]["hits"]
    rows = []
    for h in hits:
        dem = h.get("demographic") or {}
        diags = h.get("diagnoses") or [{}]
        # primary diagnosis record first when there are several
        diags = sorted(diags, key=lambda d: d.get("classification_of_tumor") != "primary")
        d = diags[0]
        r = {"patient": h.get("submitter_id"), "n_diagnoses": len(diags)}
        r.update({f"demographic.{k}": dem.get(k) for k in
                  ("gender", "vital_status", "days_to_death", "age_at_index")})
        r.update({f"diagnoses.{k}": d.get(k) for k in
                  ("primary_diagnosis", "tissue_or_organ_of_origin", "site_of_resection_or_biopsy",
                   "morphology", "ajcc_pathologic_stage", "days_to_last_follow_up")})
        rows.append(r)
    df = pd.DataFrame(rows).set_index("patient").sort_index()
    site = (df["diagnoses.tissue_or_organ_of_origin"].fillna("") + " " +
            df["diagnoses.site_of_resection_or_biopsy"].fillna("")).str.lower()
    df["is_icca"] = site.str.contains("intrahepatic")
    dead = df["demographic.vital_status"].eq("Dead")
    df["os_event"] = dead.astype(int)
    df["os_days"] = pd.to_numeric(df["demographic.days_to_death"].where(dead,
                                                                       df["diagnoses.days_to_last_follow_up"]),
                                  errors="coerce")
    return df


def fetch_tcga(key, spec, out, args, files):
    for n in spec["xena_files"]:
        url = f"{spec['xena_hub']}/{n}"
        files.append({"cohort": key, "file": n, "url": url,
                      "status": download(url, os.path.join(out, n), args.max_bytes, args.dry_run)})
    if args.dry_run:
        return
    try:
        cases = gdc_cases(spec["project"])
    except Exception as e:
        log.warning("%s: GDC API query failed: %s", key, e)
        return
    cases.to_csv(os.path.join(out, "gdc_cases.tsv"), sep="\t")
    log.info("%s: %d cases, %d intrahepatic; sites: %s", key, len(cases), int(cases["is_icca"].sum()),
             cases["diagnoses.tissue_or_organ_of_origin"].value_counts(dropna=False).to_dict())


# ---------------------------------------------------------------- NODE
def fetch_node(key, spec, out, args, files):
    man = os.path.join(out, "manual")
    os.makedirs(man, exist_ok=True)
    readme = os.path.join(man, "README_DOWNLOAD.txt")
    with open(readme, "w") as fh:
        fh.write(f"{key} ({spec['accession']}): download by hand from {spec['url']}\n")
        fh.write("Log in to NODE (free account), open the project, and download the processed files below\n")
        fh.write("into this folder. Keep the original file names. Raw FASTQ needs a data-use application.\n\n")
        for m in spec.get("manual_files", []):
            fh.write(f"  - {m}\n")
        if spec.get("paper"):
            fh.write(f"\nPaper: {spec['paper']} (supplementary tables from the journal site)\n")
        if spec.get("notes"):
            fh.write(f"Notes: {spec['notes']}\n")
    got = [f for f in os.listdir(man) if f != "README_DOWNLOAD.txt"]
    if got:
        log.info("%s: %d manual file(s) present: %s", key, len(got), ", ".join(sorted(got)))
    else:
        log.warning("%s: nothing downloaded yet - see %s", key, readme)


# ---------------------------------------------------------------- inventory
def sniff(path):
    """Short description of a tabular file: rows x columns and the first column names."""
    low = path.lower()
    try:
        if low.endswith((".xlsx", ".xls")):
            x = pd.ExcelFile(path)
            return "; ".join(f"{s}: {x.parse(s, header=None).shape}" for s in x.sheet_names[:6])
        if low.endswith(TABLE_EXT) and "series_matrix" not in low:
            sep = "," if ".csv" in low else "\t"
            head = pd.read_csv(path, sep=sep, nrows=5, comment="!" if "matrix" in low else None)
            opener = gzip.open if low.endswith(".gz") else open
            with opener(path, "rt", errors="replace") as fh:
                nrow = sum(1 for _ in fh) - 1
            return f"{nrow} rows x {head.shape[1]} cols; cols: {', '.join(map(str, head.columns[:6]))}"
        if low.endswith(".tar"):
            import tarfile
            with tarfile.open(path) as t:
                names = t.getnames()
            return f"tar with {len(names)} members, e.g. {', '.join(names[:3])}"
    except Exception as e:
        return f"unreadable: {e}"
    return ""


def inventory(cfg, root, keys):
    rows, frows = [], []
    for key in keys:
        spec = cfg["cohorts"][key]
        out = os.path.join(root, key)
        n_samples, n_icca = "", ""
        if os.path.exists(os.path.join(out, "samples_geo.tsv")):
            n_samples = len(pd.read_csv(os.path.join(out, "samples_geo.tsv"), sep="\t"))
        if os.path.exists(os.path.join(out, "gdc_cases.tsv")):
            c = pd.read_csv(os.path.join(out, "gdc_cases.tsv"), sep="\t")
            n_samples, n_icca = len(c), int(c["is_icca"].sum())
        nfiles = 0
        for dp, _, fs in os.walk(out):
            for f in fs:
                if f.endswith(".part") or f == "README_DOWNLOAD.txt":
                    continue
                fp = os.path.join(dp, f)
                nfiles += 1
                frows.append({"cohort": key, "file": os.path.relpath(fp, out),
                              "mb": round(os.path.getsize(fp) / 1e6, 2),
                              "summary": "" if "/fastq/" in fp else sniff(fp)})
        rows.append({"cohort": key, "source": spec["source"],
                     "accession": spec.get("accession", spec.get("project")), "role": spec.get("role"),
                     "expected_n": spec.get("expected_n"), "samples_found": n_samples,
                     "icca_flagged": n_icca, "files": nfiles,
                     "status": "ok" if nfiles else ("manual download needed" if spec["source"] == "node"
                                                    else "missing")})
    return pd.DataFrame(rows), pd.DataFrame(frows)


# ---------------------------------------------------------------- main
FETCH = {"geo": fetch_geo, "tcga": fetch_tcga, "node": fetch_node}


def main():
    global log
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/icca_cohorts.yaml")
    ap.add_argument("--only", nargs="+", help="cohort keys to process (default: all)")
    ap.add_argument("--max-gb", type=float, default=20.0, help="skip single files larger than this (GEO/Xena)")
    ap.add_argument("--fastq", action="store_true", help="also download FASTQ listed in the ENA manifest")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--inventory-only", action="store_true")
    args = ap.parse_args()
    args.max_bytes = int(args.max_gb * 1e9) if args.max_gb else None

    cfg = yaml.safe_load(open(p(args.config)))
    root, inv_dir = p(cfg["data_dir"]), p(cfg["inventory_dir"])
    log = setup_logging(inv_dir, "fetch_icca")
    keys = args.only or list(cfg["cohorts"])
    bad = [k for k in keys if k not in cfg["cohorts"]]
    if bad:
        sys.exit(f"unknown cohort(s) {bad}; choose from {list(cfg['cohorts'])}")

    files = []
    if not args.inventory_only:
        for key in keys:
            spec = cfg["cohorts"][key]
            out = os.path.join(root, key)
            os.makedirs(out, exist_ok=True)
            log.info("==== %s (%s %s)", key, spec["source"], spec.get("accession", spec.get("project")))
            try:
                FETCH[spec["source"]](key, spec, out, args, files)
            except Exception as e:
                log.error("%s failed: %s", key, e)
                files.append({"cohort": key, "file": "", "url": "", "status": f"error: {e}"})
        if files:
            pd.DataFrame(files).to_csv(os.path.join(inv_dir, "fetch_status.tsv"), sep="\t", index=False)
    if args.dry_run:
        return
    inv, ftab = inventory(cfg, root, keys)
    inv.to_csv(os.path.join(inv_dir, "inventory.tsv"), sep="\t", index=False)
    ftab.to_csv(os.path.join(inv_dir, "files.tsv"), sep="\t", index=False)
    log.info("inventory:\n%s", inv.drop(columns=["expected_n"]).to_string(index=False))
    failed = [f for f in files if f["status"] in ("failed", "not_found") or f["status"].startswith("error")]
    if failed:
        log.warning("%d download(s) failed or not found - see fetch_status.tsv", len(failed))


if __name__ == "__main__":
    main()
