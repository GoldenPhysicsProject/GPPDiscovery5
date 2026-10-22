"""Queen S1b: zero-side decomposition of Q_zeta(f)/||f||^2, f = e^{-i g0 x} sinh(b x) on [-L,L] (half-width L),
b = 0.308517, g0 = 85.699348 (DH off-line zero). Explicit formula: Q(f) = sum_{gamma} |F(gamma)|^2 (zeros on the line,
RH verified far beyond this height), F(t) = int f e^{itx} = sinh(a+L)/a+ - sinh(a-L)/a-, a± = ±b + i(t-g0).
Compares with the prime-side value from dh_sinh.json. NUMERICAL."""
import mpmath as mp, json, math
mp.mp.dps=20
b,g0=0.308517,85.699348
NZ=300
zs=[float(mp.zetazero(k).imag) for k in range(1,NZ+1)]
json.dump(zs,open("zeta_zeros_300.json","w"))
def F(t,L):
    out=0
    for s in (1,-1):
        a=complex(s*b,t-g0); out+= s*mp.sinh(a*L)/a
    return complex(out)
nrm=lambda L: math.sinh(2*b*L)/(2*b)-L
ref={r["L"]:r["Q_over_norm"] for r in json.load(open("dh_sinh.json"))["ZETA"]}
res=[]
for L in [1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0]:
    n=nrm(L); terms=[(abs(F(g,L))**2+abs(F(-g,L))**2)/n for g in zs]
    tot=sum(terms); near=sum(t for g,t in zip(zs,terms) if abs(g-g0)<5)
    # tail estimate beyond last zero: density log(t/2pi)/2pi times kernel bound ~ beta/(beta^2+D^2) averaged
    T=zs[-1]; tail=(math.log(T/(2*math.pi))/(2*math.pi))*b/(T-g0)  # mean kernel b/(b^2+D^2) integrated
    mean_k=sum(b/(b*b+(g-g0)**2) for g in zs)
    r=dict(L=L,zero_side=tot,tail_est=tail,near_5=near,prime_side=ref.get(L),mean_kernel_pred=mean_k,
           top=[(round(g,3),round(t,3)) for g,t in sorted(zip(zs,terms),key=lambda z:-z[1])[:4]])
    print(r,flush=True); res.append(r)
json.dump(res,open("zeta_zero_side.json","w"),indent=1)
