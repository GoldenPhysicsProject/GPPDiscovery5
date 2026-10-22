"""Grok S83 H80c step 3: compare a rebuilt T dump with the S71 reference dump entrywise (exact decimal via flint arb at 2048 bits).
Pass iff max|mid_new - mid_ref| <= radT_ref + radT_new (row-sum radius bounds from the cert jsonl, each >= every entry radius)
+ 2 print units. CHECKED-level comparison only (the PSD certificates are rerun separately). usage: cmp_dumps.py new ref radnew radref"""
import sys, json
from flint import arb, ctx
ctx.prec = 2048
new, ref = json.load(open(sys.argv[1])), json.load(open(sys.argv[2])); rn, rr = float(sys.argv[3]), float(sys.argv[4])
n = len(ref); assert len(new) == n
mx = arb(0); arg = None; same = 0
for i in range(n):
    for j in range(n):
        if new[i][j] == ref[i][j]: same += 1; continue
        d = abs(arb(new[i][j]) - arb(ref[i][j]))
        if d > mx: mx = d; arg = (i, j)
tol = rn + rr
print(json.dumps(dict(K=n, identical_strings=same, total=n * n, max_abs_diff=float(mx.upper()), argmax=arg, tol=tol, ratio=float(mx.upper()) / tol, PASS=bool(float(mx.upper()) <= tol))))
