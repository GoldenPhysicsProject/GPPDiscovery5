import json, sys, mpmath as mp
mp.mp.dps=50
rows=[]
for f in sys.argv[1:]:
    rows+=json.load(open(f))
by={}
for r in rows:
    if r["dps"]!=120: continue
    by.setdefault(r["L"],{})[r["N"]]=r
TOL=mp.mpf("1e-3")
def fmt(x): return mp.nstr(mp.mpf(x),6)
out=[]
for q in ["lam1e","lam2e","lam1o","lam2o"]:
    out.append(f"\n### {q}\n\n| L | "+" | ".join(f"N={n}" for n in sorted({n for d in by.values() for n in d}))+" | rel. change (last two N) | converged (tol 1e-3)? |")
    Ns=sorted({n for d in by.values() for n in d})
    out.append("|"+"---|"*(len(Ns)+3))
    for L in sorted(by):
        d=by[L]; vals=[d[n][q] if n in d else None for n in Ns]
        have=[n for n in Ns if n in d]
        if len(have)>=2:
            a=mp.mpf(d[have[-1]][q]); b=mp.mpf(d[have[-2]][q]); rc=abs(a-b)/abs(a)
            conv="yes" if rc<TOL else "NO"
            rcs=mp.nstr(rc,2)
        else: rcs,conv="–","–"
        out.append(f"| {L} | "+" | ".join(fmt(v) if v else "–" for v in vals)+f" | {rcs} | {conv} |")
print("\n".join(out))
