#!/usr/bin/env python
"""Extract immune-class and immunotherapy-response gene sets from the supplements of
  Montironi et al., Gut 2023;72:129   (Supp1 PDF: Suppl Table 9; Supp2 xlsx: Suppl Table 3)
  Sia et al., Gastroenterology 2017;153:812  (Supplementary Materials PDF: Suppl Tables 3 and 8)
  Haber et al., Gastroenterology 2023   (Supplementary Tables PDF: Tables S6 and S8)
into config/subtype_signatures_custom.tsv (rows for sets with these prefixes are replaced, others kept).

Requires pdfplumber (pip install pdfplumber). Every list is parsed from the table itself and checked
against the size the paper states, so a layout change stops the script instead of giving a wrong set.

Sets:
  MONTIRONI2023_INFLAMED            Inflamed 20-gene signature (Suppl Table 9)
  MONTIRONI2023_WNT_ACTIVATION (up/dn)  Wnt-beta-catenin activation signature (Suppl Table 3; Lachenmayer 2012)
  SIA2017_IMMUNE_CLASS_UP           genes over-expressed in the Immune class (Suppl Table 3; multi-gene
                                    probes "A /// B" contribute every symbol)
  SIA2017_CLASSIFIER_IMMUNE / SIA2017_CLASSIFIER_REST   Immune class gene classifier (Suppl Table 8)
  HABER2023_IFNAP                   IFNAP 11-gene signature (Table S8)
  HABER2023_S6_<NAME> (up/dn)       Table S6: previously reported response signatures and immune-cell
                                    subsets (Bindea 2013 cell types, POPLAR Teff, Sangro inflammatory,
                                    cytolytic activity, Ayers IFN, TLS, IMPRES, MDSC, WNT/TGFb, Jerby-Arnon,
                                    melanoma ICI/MAPKi resistance...). Parsed by word position, because the
                                    rows wrap; labels ending in up/down become the two directions of one set.
                                    Numeric tokens (Excel-corrupted symbols, e.g. 44261) are dropped.
"""

import argparse
import os
import re

import pandas as pd
import pdfplumber

PREFIXES = ("MONTIRONI2023_", "SIA2017_", "HABER2023_")
GENE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9.\-@_]*$")


def page_texts(path):
    with pdfplumber.open(path) as pdf:
        return [(p.extract_text() or "") for p in pdf.pages]


def lines_after(texts, start_pat, stop_pat, occurrence=-1):
    """Lines between the chosen occurrence of start_pat and the next stop_pat (across pages)."""
    lines = [l for t in texts for l in t.splitlines()]
    starts = [i for i, l in enumerate(lines) if re.search(start_pat, l)]
    if not starts:
        raise SystemExit(f"table header not found: {start_pat}")
    s = starts[occurrence]
    e = next((i for i in range(s + 1, len(lines)) if re.search(stop_pat, lines[i])), len(lines))
    return lines[s + 1:e]


def montironi(supp1, supp2):
    out = []
    ls = lines_after(page_texts(supp1), r"Supplementary Table 9: Genes comprising the Inflamed signature",
                     r"Supplementary Table 10")
    genes = [l.split()[0] for l in ls if l.strip() and l.split()[0] not in ("Gene",) and GENE.match(l.split()[0])
             and not l.strip()[0].isdigit()]
    if len(genes) != 20:
        raise SystemExit(f"Montironi Inflamed signature: expected 20 genes, got {len(genes)}: {genes}")
    out += [("MONTIRONI2023_INFLAMED", "Montironi 2023 Gut 72:129 Suppl Table 9", "up", g) for g in genes]
    d = pd.read_excel(supp2, sheet_name="Suppl Table 3")
    ref = "Montironi 2023 Gut 72:129 Suppl Table 3 (Wnt-bcatenin activation; Lachenmayer 2012)"
    up, dn = d["Upregulated Gene"].dropna().astype(str).str.strip(), d["Downregulated Gene"].dropna().astype(str).str.strip()
    out += [("MONTIRONI2023_WNT_ACTIVATION", ref, "up", g) for g in sorted(set(up))]
    out += [("MONTIRONI2023_WNT_ACTIVATION", ref, "dn", g) for g in sorted(set(dn))]
    print(f"Montironi: Inflamed {len(genes)} genes; Wnt activation {up.nunique()} up / {dn.nunique()} dn")
    return out


def sia(pdf):
    texts = page_texts(pdf)
    out = []
    ls = lines_after(texts, r"Supplementary Table 3\. List of genes", r"Supplementary Table 4\.")
    up = []
    for l in ls:
        m = re.match(r"^\s*((?:[A-Za-z0-9.\-]+)(?:\s+///\s+[A-Za-z0-9.\-]+)*)\s", l)
        if m and m.group(1) != "Feature" and re.search(r"-?\d+\.\d\d\s+\d\.\d\d\s+\d\.\d\d\s+\d", l):
            up += [g.strip() for g in m.group(1).split("///")]
    up = list(dict.fromkeys(up))
    if len(up) < 50:
        raise SystemExit(f"Sia Table 3: only {len(up)} genes parsed")
    out += [("SIA2017_IMMUNE_CLASS_UP", "Sia 2017 Gastroenterology 153:812 Suppl Table 3", "up", g) for g in up]
    ls = lines_after(texts, r"Supplementary Table 8\. Immune subclass gene classifier", r"Supplementary Table 9\.")
    cls = {"Immune class": [], "Rest": []}
    for l in ls:
        m = re.match(r"^\s*(\S+)\s+(\S+)\s+(Immune class|Rest)\s*$", l)
        if m:
            cls[m.group(3)].append(m.group(2))
    ref = "Sia 2017 Gastroenterology 153:812 Suppl Table 8 (Immune class classifier)"
    for c, name in (("Immune class", "SIA2017_CLASSIFIER_IMMUNE"), ("Rest", "SIA2017_CLASSIFIER_REST")):
        g = list(dict.fromkeys(cls[c]))
        out += [(name, ref, "up", x) for x in g]
    print(f"Sia: Immune class up {len(up)} genes; classifier {len(set(cls['Immune class']))} immune / "
          f"{len(set(cls['Rest']))} rest")
    if not cls["Immune class"] or not cls["Rest"]:
        raise SystemExit("Sia Table 8 classifier not parsed")
    return out


def haber(pdf_path):
    out = []
    texts = page_texts(pdf_path)
    ls = lines_after(texts, r"Table S8: Genes incorporated in the IFNAP signature", r"Table S9")
    genes = [l.strip() for l in ls if l.strip() and l.strip() != "Genename" and GENE.match(l.strip())]
    if len(genes) != 11:
        raise SystemExit(f"Haber IFNAP: expected 11 genes, got {len(genes)}: {genes}")
    out += [("HABER2023_IFNAP", "Haber 2023 Gastroenterology Table S8 (IFNAP)", "up", g) for g in genes]
    n_ifnap = len(genes)

    sets, pending = [], []
    with pdfplumber.open(pdf_path) as pdf:
        s = next(i for i, t in enumerate(texts) if "Table S6" in t)
        e = next(i for i, t in enumerate(texts) if "Table S7" in t)
        for pi in range(s, e + 1):
            ws = pdf.pages[pi].extract_words()
            if pi == s:
                first = next(w for w in ws if w["text"] == "aDC")
                G = next(w for w in ws if w["text"] == "CCL1")["x0"] - 5
                ws = [w for w in ws if w["top"] >= first["top"] - 2]
            if pi == e:
                stop = min(w["top"] for w in ws if w["text"] == "Table")
                ws = [w for w in ws if w["top"] < stop - 2]
            rows = {}
            for w in ws:
                k = next((t for t in rows if abs(t - w["top"]) < 3), w["top"])
                rows.setdefault(k, []).append(w)
            for t in sorted(rows):
                r = sorted(rows[t], key=lambda w: w["x0"])
                lab = [w["text"] for w in r if w["x0"] < G]
                gen = [w["text"] for w in r if w["x0"] >= G]
                if lab and not gen:
                    pending += lab
                elif lab:
                    sets.append([" ".join(pending + lab), gen])
                    pending = []
                elif gen and sets:
                    sets[-1][1] += gen
    if len(sets) != 46:
        raise SystemExit(f"Haber Table S6: expected 46 sets, parsed {len(sets)}")
    ref = "Haber 2023 Gastroenterology Table S6"
    for label, gs in sets:
        m = re.search(r"[\s_](up|down)$", label, flags=re.I)
        direction = "dn" if m and m.group(1).lower() == "down" else "up"
        base = label[:m.start()] if m else label
        name = "HABER2023_S6_" + re.sub(r"[^A-Za-z0-9]+", "_", base).strip("_").upper()
        g = [x for x in dict.fromkeys(gs) if GENE.match(x) and not x.isdigit()]
        out += [(name, f"{ref} '{label}'", direction, x) for x in g]
    print(f"Haber: IFNAP {n_ifnap} genes; Table S6 {len(sets)} sets")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--montironi-supp1", required=True)
    ap.add_argument("--montironi-supp2", required=True)
    ap.add_argument("--sia", required=True)
    ap.add_argument("--haber", required=True)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "..", "config",
                                                  "subtype_signatures_custom.tsv"))
    a = ap.parse_args()
    new = montironi(a.montironi_supp1, a.montironi_supp2) + sia(a.sia) + haber(a.haber)
    comments = [l for l in open(a.out).read().splitlines() if l.startswith("#")]
    body = pd.read_csv(a.out, sep="\t", comment="#", dtype=str)
    body = body[~body["set"].fillna("").str.startswith(PREFIXES)]
    body = pd.concat([body, pd.DataFrame(new, columns=["set", "reference", "direction", "gene_symbol"])])
    with open(a.out, "w") as fh:
        fh.write("\n".join(comments) + "\n")
        body.to_csv(fh, sep="\t", index=False)
    print(f"wrote {len(new)} rows ({len({(r[0], r[2]) for r in new})} set-directions) to {a.out}; {len(body)} rows total")


if __name__ == "__main__":
    main()
