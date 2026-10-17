"""Grok S63: collect dhritz.py (FORM=DH) runs into dh_scan_summary.json (sorted by L, K)."""
import json, glob
from fractions import Fraction
rows=[]
for f in glob.glob('dh_L*.jsonl'):
    for l in open(f):
        d=json.loads(l)
        for s in ('even','odd'):
            if s in d: rows.append({"L":d["L"],"Lf":float(Fraction(d["L"])),"K":d["K"],"sector":s,"ritz1":d[s]["ritz1"],"ritz2":d[s]["ritz2"],
                                    "RQ_arb":d[s]["RQ_arb"],"RQ_upper":d[s]["RQ_upper"],"negative_certified":d[s]["RQ_upper"]<0})
rows.sort(key=lambda r:(r["sector"],r["Lf"],r["K"]))
json.dump(rows,open('dh_scan_summary.json','w'),indent=1)
for r in rows: print(f'{r["sector"]:4s} L={r["L"]:>9s} ({r["Lf"]:.4f}) K={r["K"]:3d}  RQ={r["RQ_arb"]:34s} {"NEG-CERT" if r["negative_certified"] else ""}')
