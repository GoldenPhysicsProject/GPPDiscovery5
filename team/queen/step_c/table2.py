import json, mpmath as mp
mp.mp.dps=50
rows=json.load(open("grid.json"))+json.load(open("grid_hiN.json"))+json.load(open("grid_hiN2.json"))
by={}
for r in rows:
    if r["dps"]==120: by.setdefault(r["L"],{})[r["N"]]=r
Ns=[24,32,40,48,56,64]; TOL=mp.mpf("0.02")
lines=[]; summary={}
for q,name in [("lam1e","λ₁(R_e)"),("lam2e","λ₂(R_e)"),("lam1o","λ₁(R_o)"),("lam2o","λ₂(R_o)")]:
    lines.append(f"\n#### {name}\n\n| L | "+" | ".join(f"N={n}" for n in Ns)+" | rel. change, last two N | converged (2%)? |\n|"+"---|"*(len(Ns)+3))
    for L in sorted(by):
        d=by[L]; have=[n for n in Ns if n in d]
        a=mp.mpf(d[have[-1]][q]); b=mp.mpf(d[have[-2]][q]); rc=abs(a-b)/abs(a)
        ok = rc<TOL and have[-1]>=56
        conv = "**yes**" if ok else ("yes at N=48 only (no N≥56 run)" if rc<TOL else "NO")
        lines.append(f"| {L} | "+" | ".join(mp.nstr(mp.mpf(d[n][q]),5) if n in d else "–" for n in Ns)+f" | {mp.nstr(rc,2)} (N{have[-2]}→{have[-1]}) | {conv} |")
        summary.setdefault(q,{})[L]=dict(value=d[have[-1]][q],N=have[-1],relchange=mp.nstr(rc,3),converged_2pct=bool(rc<TOL))
open("tables_c.md","w").write("\n".join(lines)); json.dump(summary,open("summary_c.json","w"),indent=1)
print("\n".join(lines))
# eigenvector export: highest N per L, plus overlap with previous N (nested Legendre basis -> zero-pad)
ev=[]
for L in sorted(by):
    d=by[L]; have=sorted(d); n1,n0=have[-1],have[-2]
    v1=[mp.mpf(x) for x in d[n1]["v2e"]]; v0=[mp.mpf(x) for x in d[n0]["v2e"]]+[mp.mpf(0)]*(n1-n0)
    ov=abs(mp.fsum(a*b for a,b in zip(v1,v0)))/mp.sqrt(mp.fsum(a*a for a in v1)*mp.fsum(b*b for b in v0))
    ev.append(dict(L=L,N=n1,dps=120,eigenvalue_lam2_Re=d[n1]["lam2e"],eigenvalue_converged_2pct=summary["lam2e"][L]["converged_2pct"],
                   abs_cos_with_N_prev=mp.nstr(ov,12),N_prev=n0,coefficients=d[n1]["v2e"]))
    print(L,n1,mp.nstr(1-ov,3))
json.dump(dict(description="Queen S1c step (c): second even eigenvector of R_e = A_e - P_e (zeta window form minus pole; S87 hp_forms.py, ZETA), "
  "basis b_j(s) = sqrt((2j+1)/L) P_j(2s/L - 1), s = log x in [0,L], extended evenly to [-L,L] (even sector); coefficients are in this orthonormal basis "
  "with metric L^2(ds) on [0,L] (half the L^2[-L,L] norm: multiply by 1/sqrt(2) for unit norm on [-L,L]). Half-width L, primes n <= e^{2L}. "
  "Sign of eigenvector arbitrary. HIGH-PRECISION NUMERICAL (mpmath dps 120; dps 170 check identical to 29+ digits at (L,N)=(1.0,64),(1.2,48),(2.0,48)). "
  "Eigenvalues for L >= 1.1 are NOT converged in N (upper bounds only); eigenvectors there are Ritz vectors of unconverged eigenvalues. RH OPEN.",
  vectors=ev),open("lambda2_Re_eigenvectors.json","w"),indent=1)
