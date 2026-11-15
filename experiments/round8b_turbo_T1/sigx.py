"""Band-limited effective symbol sigma_x(t): Q(g) = (1/pi) int_0^inf sigma_x(t) |G(t)|^2 dt for g in L^2[0,1/4] (any g; exact identity
because |G|^2 has spectrum in |x| <= 2l = 1/2). Float NUMERICAL study of where sigma_x > 0 and the low-frequency mass of H_J."""
import numpy as np, json
from scipy import special
from scipy.integrate import quad
l = 0.25
NSn = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
h = lambda u: 2*np.exp(-u/2)/(1-np.exp(-2*u))
w0 = special.digamma(0.25) - np.log(np.pi)
H5 = quad(h, 0.5, np.inf, limit=200)[0]
C0 = 3*w0 + 3*H5
k = lambda s: h(1+s) - 0.25*h(2+s)
kap = lambda s: k(s) + k(-s)
mus = []
for n,p in NSn:
    lam = np.log(n); ln = np.log(p)/np.sqrt(n)
    if abs(lam-1) < 0.5: mus.append((abs(lam-1), 2*ln, n))
    if abs(lam-2) < 0.5: mus.append((abs(lam-2), -0.5*ln, n))
xg, wg = np.polynomial.legendre.leggauss(40)
def gl(a,b,n):
    e = np.linspace(a,b,n+1); return (np.concatenate([(p+q)/2+(q-p)/2*xg for p,q in zip(e[:-1],e[1:])]), np.concatenate([(q-p)/2*wg for p,q in zip(e[:-1],e[1:])]))
u, wu = gl(0, 0.5, 400)
def sigx(t):
    t = np.atleast_1d(t)
    a = 6*(np.cos(np.outer(t,u))*-1 + 1) @ (wu*h(u))
    b = 4*np.cos(np.outer(t,u)) @ (wu*kap(u))
    c = 4*sum(cn*np.cos(m*t) for m,cn,_ in mus)
    return 2*C0 + a + b + c
wsym = lambda t: special.digamma(0.25+0.5j*t).real - np.log(np.pi)
tt = np.linspace(0, 3000, 300001)
sx = np.concatenate([sigx(tt[i:i+20000]) for i in range(0, len(tt), 20000)])
neg = tt[sx < 0]
out = dict(C0=C0, H5=H5, w0=w0, mus=[(m,c,n) for m,c,n in mus], sum_abs_c4=4*sum(abs(c) for _,c,_ in mus),
           last_negative_t=float(neg.max()) if len(neg) else None, sigx_at=[[float(a), float(sigx(a)[0])] for a in [0,10,50,100,163.36,500,917.3,3000]])
for T in [163.36, 917.3]:
    m = tt >= T; out[f'min_sigx_over_6w_t>={T}'] = float(np.min(sx[m]/(6*wsym(tt[m])))); out[f'min_sigx_t>={T}'] = float(sx[m].min())
print(json.dumps(out, indent=1)); json.dump(out, open('sigx_results.json','w'), indent=1)
