"""Grok S75: resource extrapolation (NUMERICAL ESTIMATE, not a certificate) for the Mori-route chain at the DH window.
Astra necessary bound: log Omega > log 2pi + s + m, s = Schur upper bound of C_L (s/b_L ~ 1.0505 at L=7/5 and 71/50).
Empirical cure factor F = Omega_cert/Omega_bound: 2.12 (7/5, both parities), 2.13 even / 2.45 odd (71/50).
N = 1.671*Omega*L; K = N/2 per parity. Scaling anchored on S73a odd run (N 2732, K 1366): assembly 17,868 worker-s,
3.5 GB/worker, LDLt 503 s (640 bits), dump size measured."""
import mpmath as mp, os, json
mp.mp.dps=30
def mang(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return mp.log(p) if m==1 else 0
def bL(L): return 2*sum(mang(n)/mp.sqrt(n)*(1-mp.log(n)/(2*L)) for n in range(2,int(mp.exp(2*L))+1))
dump=os.path.getsize('../s71/run/T_L7150_m12_O1150_p1.json')
N0,K0,Om0,asm0,mem0,ldl0=2732,1366,1150,17868.,3.5,503.
rows=[]
for Ls in ['7/5','71/50','3/2','8/5','17/10','137/80']:
    p,q=map(int,Ls.split('/')); L=mp.mpf(p)/q; b=bL(L); s=1.0505*b; Ob=2*mp.pi*mp.exp(s+0.5)
    for F in (2.12,2.45):
        Om=F*Ob; N=1.671*Om*L; K=N/2; r=K/K0
        rows.append(dict(L=Ls,bL=float(b),Om_bound=float(Ob),F=F,Omega=float(Om),N=float(N),
          asm_worker_h_N2=float(asm0*r**2/3600), asm_worker_h_N2Om=float(asm0*r**2*Om/Om0/3600),
          mem_GB_per_worker_linN=float(mem0*r), dump_GB_N2=float(dump*r**2/1e9), ldlt_h_K3=float(ldl0*r**3/3600),
          mlog2_dc=float(2*(Om*L+1)/mp.log(2))))
for r in rows: print(json.dumps({k:(round(v,3) if isinstance(v,float) else v) for k,v in r.items()}))
json.dump(rows,open('mori_cost_dh.json','w'),indent=1)
