"""Queen S1b: exact quartet term on the matched vector f = e^{-i g0 x} sinh(b x), [-L,L] half-width.
For the zero quartet {1/2 +- b +- i g0}: ghat(z) = F(z) conj(F(conj z)), F(t) = int f e^{itx}.
At z = g0 -+ i b: ghat = -(S-L)^2 with S = sinh(2bL)/(2b), and ||f||^2 = S-L, so these two give exactly -2||f||^4.
Residual R = Q/||f||^2 - quartet/||f||^2 = contribution of all other zeros. NUMERICAL."""
import mpmath as mp, json, math
def F(t,L,b,g0):
    o=0
    for s in (1,-1):
        a=s*b+1j*(t-g0); o+= s*mp.sinh(a*L)/a if abs(a)>1e-30 else s*L
    return complex(o)
def quartet(L,b,g0):
    tot=0
    for z in (g0-1j*b, g0+1j*b, -g0-1j*b, -g0+1j*b):
        tot+=F(z,L,b,g0)*F(z.conjugate(),L,b,g0).conjugate()
    return tot.real
out={}
dh=json.load(open("../GPPDiscovery5/team/queen/dh_sinh.json"))["DH"] if False else json.load(open("/workspace/queen/GPPDiscovery5/team/queen/dh_sinh.json"))["DH"]
ep=json.load(open("epstein_growth.json"))["2x2+2xy+3y2"]
for name,rows,b,g0 in [("DH",dh,0.308517,85.699348),("Epstein 2x2+2xy+3y2",ep,1.57941410481,7.62886758778)]:
    out[name]=[]
    for r in rows:
        L=r["L"]; n=r["norm2"]; q=quartet(L,b,g0)/n
        d=dict(L=L,Q_over_norm=r["Q_over_norm"],minus2norm=-2*n,quartet=q,residual=r["Q_over_norm"]-q)
        out[name].append(d); print(name,d)
json.dump(out,open("quartet_identity.json","w"),indent=1)
