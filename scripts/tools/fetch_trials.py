#!/usr/bin/env python
"""Fetch HCC clinical trials from ClinicalTrials.gov (API v2) and extract, per arm, what the trial-emulation step needs:
objective response rate (ORR) and its denominator, assessment criteria and reviewer, and baseline characteristics
(mean / median age, % female, % Asian); plus the eligibility criteria text. Values come from the posted results only.

    python scripts/tools/fetch_trials.py            # trials in config/clinical_trials.yaml
Outputs: data/trials/<NCT>.json (raw), results/10_response/trials/{trial_arms.tsv, orr_measures.tsv, eligibility/<NCT>.txt}
"""

import json
import os
import re
import time
import urllib.request

import pandas as pd
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
API = "https://clinicaltrials.gov/api/v2/studies/{}"


def get(nct):
    f = os.path.join(ROOT, "data", "trials", f"{nct}.json")
    if not os.path.exists(f):
        for i in range(4):
            try:
                with urllib.request.urlopen(API.format(nct), timeout=60) as r:
                    open(f, "wb").write(r.read())
                break
            except Exception:  # noqa: BLE001
                time.sleep(2 ** i)
    return json.load(open(f))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config", "clinical_trials.yaml")))
    out = os.path.join(ROOT, "results", "10_response", "trials")
    os.makedirs(os.path.join(out, "eligibility"), exist_ok=True)
    arms, orrs = [], []
    for t in cfg["trials"]:
        d = get(t["nct"])
        ps, rs = d["protocolSection"], d.get("resultsSection")
        open(os.path.join(out, "eligibility", f"{t['nct']}.txt"), "w").write(ps["eligibilityModule"].get("eligibilityCriteria", ""))
        if not rs:
            print(f"{t['nct']} {t['name']}: no posted results")
            continue
        bl = rs["baselineCharacteristicsModule"]
        gid = {g["id"]: g["title"] for g in bl["groups"]}
        denom = {c["groupId"]: num(c["value"]) for c in bl["denoms"][0]["counts"]}
        base = {g: {"baseline_n": denom.get(g)} for g in gid}
        for m in bl["measures"]:
            title = m["title"].lower()
            for c in m.get("classes", []):
                if c.get("title") and c["title"].lower() not in ("global", "overall", ""):
                    continue                                   # first (global) class only
                for cc in c.get("categories", []):
                    for x in cc["measurements"]:
                        g, v = x["groupId"], num(x["value"])
                        if g not in base or v is None:
                            continue
                        if title.startswith("age") and ("continuous" in title or m.get("paramType", "") in ("MEAN", "MEDIAN")):
                            base[g].setdefault("age", v)
                            base[g].setdefault("age_stat", m.get("paramType"))
                        elif "sex" in title and cc.get("title", "").lower() == "female":
                            base[g].setdefault("female_n", v)
                        elif "race" in title and cc.get("title", "").lower() == "asian":
                            base[g].setdefault("asian_n", v)
        for g, b in base.items():
            if gid[g].lower() == "total":
                continue
            n = b.get("baseline_n")
            arms.append({"nct": t["nct"], "trial": t["name"], "group": gid[g], "baseline_n": n, "age": b.get("age"), "age_stat": b.get("age_stat"),
                         "female_frac": b["female_n"] / n if n and b.get("female_n") is not None else None,
                         "asian_frac": b["asian_n"] / n if n and b.get("asian_n") is not None else None})
        for o in rs["outcomeMeasuresModule"]["outcomeMeasures"]:
            ti = o["title"]
            if not re.search(r"objective response|overall response|\bORR\b", ti, re.I) or "duration" in ti.lower():
                continue
            og = {g["id"]: g["title"] for g in o["groups"]}
            den = {c["groupId"]: num(c["value"]) for c in (o.get("denoms") or [{}])[0].get("counts", [])}
            for c in o.get("classes", []):
                for cc in c.get("categories", []):
                    for x in cc["measurements"]:
                        orrs.append({"nct": t["nct"], "trial": t["name"], "measure": ti, "type": o.get("type"), "unit": o.get("unitOfMeasure"),
                                     "group": og.get(x["groupId"]), "class": c.get("title", ""), "category": cc.get("title", ""),
                                     "value": num(x["value"]), "n": den.get(x["groupId"])})
        print(f"{t['nct']} {t['name']}: {sum(1 for a in arms if a['nct'] == t['nct'])} baseline groups, "
              f"{sum(1 for r in orrs if r['nct'] == t['nct'])} ORR values")
    pd.DataFrame(arms).to_csv(os.path.join(out, "trial_arms.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(orrs).to_csv(os.path.join(out, "orr_measures.tsv"), sep="\t", index=False, float_format="%.4g")


if __name__ == "__main__":
    main()
