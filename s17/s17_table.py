"""Collect ritz_M*.json into a markdown table (Ritz min per M; certified upper (largest M); certified lower)."""
import json, glob, re, os
d = os.path.dirname(os.path.abspath(__file__))
Ms = sorted(int(re.findall(r"M(\d+)", f)[0]) for f in glob.glob(os.path.join(d, "ritz_M*.json")))
R = {M: {(r["h"], r["eta"], r["sector"]): r for r in json.load(open(os.path.join(d, f"ritz_M{M}.json")))} for M in Ms}
keys = list(R[Ms[-1]].keys())
short = lambda e: e.split(" (")[0]
lines = ["| h | eta | g-sector | " + " | ".join(f"Ritz M={M}" for M in Ms) + f" | CERT upper (M={Ms[-1]}) | CERT lower |", "|"+"---|"*(5+len(Ms))]
for k in keys:
    h, e, s = k
    lines.append(f"| {h} | {short(e)} | {s} | " + " | ".join(f'{float(R[M][k]["ritz_min"]):.6g}' if k in R[M] else "" for M in Ms)
                 + f" | {R[Ms[-1]][k]['cert_upper']} | {R[Ms[-1]][k]['cert_lower']} |")
open(os.path.join(d, "TABLE.md"), "w").write("\n".join(lines)+"\n"); print("\n".join(lines))
