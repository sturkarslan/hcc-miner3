"""Reference signature panel (config/reference_panel.yaml): which published signatures and classes the
downstream analyses (07c, 08) and publication figures (09) use. Keeps those choices in one config file
instead of hard-coded lists.

    from hcc_panel import Panel
    PN = Panel()                       # or Panel(path, log)
    PN.signatures(available)           # [{set, label, tag, block, group}], only sets in `available`
    PN.blocks                          # {block: {"colour": ..., "anchors": [...]}}
    PN.class_calls(smap_dir, samples)  # {classifier key: Series of calls}
    PN.known_scores(z, sig)            # samples x (adjust entries + "global mean z")
    PN.known_models(columns)           # [(name, [columns], marker, colour)]
    PN.label(entry)                    # display label for a set or signed entry
"""

import os

import numpy as np
import pandas as pd
import yaml

from hcc_common import project_root


def short(sig):
    """Fallback display label for a gene-set name."""
    s = sig.replace("custom:", "").replace("HALLMARK_", "H: ").replace("_LIVER_CANCER", "")
    s = s.replace("_HEPATOCELLULAR_CARCINOMA", "").replace("_SUBCLASS", "").replace("LIVER_CANCER_", "")
    s = s.replace("_", " ")
    return s.title().replace("Ctnnb1", "CTNNB1").replace("Epcam", "EPCAM").replace("Krt19", "KRT19") \
        .replace("Myc", "MYC").replace("Mtorc1", "mTORC1").replace("Tgf", "TGF").replace(" Up", " up") \
        .replace(" Dn", " dn").replace("Wnt Beta Catenin", "WNT/β-catenin").replace("Ifnap", "IFNAP") \
        .replace("Abrs", "ABRS").replace("Tp53", "TP53").replace("Adh1A", "ADH1A").replace("Pycr2", "PYCR2")


class Panel:
    def __init__(self, path=None, log=None):
        path = path or os.path.join(project_root(), "config", "reference_panel.yaml")
        self.cfg = yaml.safe_load(open(path))
        self.log = log
        self.blocks = self.cfg["blocks"]
        self.block_min_r = float(self.cfg.get("block_min_r", 0.3))
        self.classifiers = self.cfg["classifiers"]
        self.figures = self.cfg.get("figures", {})
        self._info = {s["set"]: s for s in self.cfg.get("signatures", [])}

    def _warn(self, msg, *a):
        if self.log:
            self.log.warning(msg, *a)

    # ---------------------------------------------------------------- signatures
    def signatures(self, available=None):
        out = []
        for s in self.cfg.get("signatures", []):
            if available is not None and s["set"] not in available:
                self._warn("panel signature %s not available; skipped", s["set"])
                continue
            out.append(s)
        return out

    def label(self, entry):
        if entry in self._info:
            return self._info[entry]["label"]
        if " - " in entry:
            a = entry.split(" - ")[0].strip()
            return (self._info[a]["label"] if a in self._info else short(a)) + " (signed)"
        return short(entry)

    def tag(self, entry):
        return self._info.get(entry, {}).get("tag", "")

    def block_of_program(self, cor_row):
        """Biology block for one program from its activity correlations (Series indexed by set)."""
        best, val = "Other", -np.inf
        for b, spec in self.blocks.items():
            a = [s for s in spec["anchors"] if s in cor_row.index]
            if not a:
                continue
            v = cor_row[a].mean()
            if v > val:
                best, val = b, v
        return best if val >= self.block_min_r else "Other"

    # ---------------------------------------------------------------- classes
    def class_label(self, key, level):
        lv = self.classifiers.get(key, {}).get("levels", {}).get(level)
        return lv[0] if lv else str(level)

    def class_axis(self, key, level):
        lv = self.classifiers.get(key, {}).get("levels", {}).get(level)
        return lv[1] if lv else "other"

    def class_title(self, key):
        return self.classifiers.get(key, {}).get("title", key)

    def class_calls(self, smap_dir, samples):
        """{key: Series of calls indexed by sample} for every classifier with available calls."""
        out = {}
        for key, spec in self.classifiers.items():
            if spec.get("source") == "label":
                col = spec.get("column", key)
                if col in samples.columns:
                    out[key] = samples[col]
                else:
                    self._warn("class labels %s: column %s not in samples.tsv", key, col)
                continue
            f = os.path.join(smap_dir, f"ntp_calls_{key}.tsv")
            if os.path.exists(f):
                out[key] = pd.read_csv(f, sep="\t", index_col=0)["call"]
            else:
                self._warn("class calls %s not found (%s); rerun 07b", key, f)
        return out

    # ---------------------------------------------------------------- scores
    @staticmethod
    def _set_mean(z, sig, name, min_genes=5):
        g = sorted(set(sig.loc[sig["set"] == name, "ensembl"]) & set(z.index))
        return z.loc[g].mean() if len(g) >= min_genes else None

    def score(self, z, sig, entry, min_genes=5):
        """Per-sample score of a set ("A") or signed pair ("A - B"); None if too few genes."""
        if " - " in entry:
            a, b = [x.strip() for x in entry.split(" - ")]
            va, vb = self._set_mean(z, sig, a, min_genes), self._set_mean(z, sig, b, min_genes)
            if va is None:
                return None
            return va - vb if vb is not None else va
        return self._set_mean(z, sig, entry, min_genes)

    def known_scores(self, z, sig, entries=None):
        """samples x (one column per adjust entry, per-sample global mean regressed out) + 'global mean z'.
        Columns are the entry strings; use label() for display."""
        entries = entries or self.cfg["adjust_signatures"]
        gm = z.mean(axis=0)
        gc = gm - gm.mean()
        out = {}
        for e in entries:
            v = self.score(z, sig, e)
            if v is None:
                self._warn("adjust signature %s: fewer than 5 genes available; skipped", e)
                continue
            out[e] = v - gc * ((gc * (v - v.mean())).sum() / (gc ** 2).sum())
        X = pd.DataFrame(out)
        X["global mean z"] = gm
        return X

    def known_models(self, columns):
        """[(name, columns present, marker, colour key)] for the head-to-head prognostic comparison."""
        out = []
        for m in self.cfg.get("known_models", []):
            cols = list(columns) if m["sets"] == "all" else [c for c in m["sets"] if c in columns]
            if not cols:
                self._warn("known model %s: none of its signatures available; skipped", m["name"])
                continue
            out.append((m["name"], cols, m.get("marker", "o"), m.get("colour", "muted")))
        return out
