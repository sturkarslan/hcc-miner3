#!/usr/bin/env python
"""Drug -> targets and action type from the Open Targets Platform (GraphQL), for drug-constrained network activity
(step 10b). Same source as the GBM DCNA work (drug_repurposing/scripts/03_opentargets_query.py), queried by drug name.

    python scripts/tools/fetch_drug_targets.py sorafenib lenvatinib ...   (default: HCC standard-of-care + TACE agents)
Output: data/reference/ot_drug_targets.tsv  (drug, chembl_id, mechanism, action_type, target_ensembl, target_symbol)
"""

import json
import os
import sys
import time
import urllib.request

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
URL = "https://api.platform.opentargets.org/api/v4/graphql"
DEFAULT = ["sorafenib", "lenvatinib", "regorafenib", "cabozantinib", "ramucirumab", "bevacizumab", "atezolizumab",
           "durvalumab", "tremelimumab", "nivolumab", "ipilimumab", "pembrolizumab", "tislelizumab", "camrelizumab",
           "apatinib", "doxorubicin", "epirubicin", "cisplatin", "mitomycin", "fluorouracil", "oxaliplatin"]
Q_SEARCH = """query s($q:String!){search(queryString:$q, entityNames:["drug"], page:{index:0,size:3}){hits{id name entity}}}"""
Q_DRUG = """query d($id:String!){drug(chemblId:$id){id name mechanismsOfAction{rows{mechanismOfAction actionType
           targets{id approvedSymbol}}}}}"""


def gql(q, v):
    for i in range(4):
        try:
            req = urllib.request.Request(URL, data=json.dumps({"query": q, "variables": v}).encode(),
                                         headers={"Content-Type": "application/json"})
            return json.load(urllib.request.urlopen(req, timeout=60))["data"]
        except Exception as e:  # noqa: BLE001
            err = e
            time.sleep(2 ** i)
    raise err


def main():
    if sys.argv[1:2] == ["--file"]:
        names = [l.strip() for l in open(sys.argv[2]) if l.strip()]
    else:
        names = sys.argv[1:] or DEFAULT
    rows = []
    for n in names:
        try:
            hits = [h for h in gql(Q_SEARCH, {"q": n})["search"]["hits"] if h["entity"] == "drug"]
        except Exception as e:  # noqa: BLE001
            print(f"{n}: query failed ({e})", flush=True)
            continue
        if not hits:
            print(f"{n}: not found")
            continue
        h = next((h for h in hits if h["name"].lower() == n.lower()), hits[0])
        try:
            d = gql(Q_DRUG, {"id": h["id"]})["drug"]
        except Exception as e:  # noqa: BLE001
            print(f"{n}: drug query failed ({e})", flush=True)
            continue
        if not d or not d.get("mechanismsOfAction"):
            print(f"{n}: no mechanism of action")
            continue
        for m in d["mechanismsOfAction"]["rows"]:
            for t in m["targets"] or []:
                rows.append({"drug": n.lower(), "chembl_id": d["id"], "ot_name": d["name"], "mechanism": m["mechanismOfAction"],
                             "action_type": m["actionType"], "target_ensembl": t["id"], "target_symbol": t["approvedSymbol"]})
        print(f"{n}: {d['name']} ({d['id']}), {sum(1 for r in rows if r['drug'] == n.lower())} targets", flush=True)
    out = os.path.join(ROOT, "data", "reference", "ot_drug_targets.tsv")
    old = pd.read_csv(out, sep="\t") if os.path.exists(out) else pd.DataFrame()
    new = pd.DataFrame(rows)
    if len(old):
        new = pd.concat([old[~old["drug"].isin(new["drug"])], new])
    new.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main()
