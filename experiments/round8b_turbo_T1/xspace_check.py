"""Round 8 (Turbo 6890033593) (b): x-space form of B on Psi = span{psi_j}, psi_j = phi_j - w_j delta = -(phi_j/omega_j^2)''.
Check: B(psi_j,psi_k) (spectral, Round 7 code) == Q(g_j,g_k), g_j = phi_j/omega_j^2, with Q in autocorrelation form."""
import numpy as np, json
from scipy import special
exec(open('p2_boundary_float.py').read().split("out = dict")[0])
l = 0.25
NSn = [(2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),(11,11)]
h = lambda u: 2*np.exp(-u/2)/(1-np.exp(-2*u))
w0 = special.digamma(0.25) - np.log(np.pi)
# Gauss-Legendre on [0,0.5] in many panels
xg, wg = np.polynomial.legendre.leggauss(30)
def gl(a, b, n):
    e = np.linspace(a, b, n+1)
    return (np.concatenate([(p+q)/2+(q-p)/2*xg for p,q in zip(e[:-1],e[1:])]), np.concatenate([(q-p)/2*wg for p,q in zip(e[:-1],e[1:])]))
u, wu = gl(0, 0.5, 2000)
from scipy.integrate import quad
H5 = quad(h, 0.5, np.inf, limit=200)[0]
def Ajk(j, k, s):   # (2/l) int_{-l}^{l-s} sin(a x) sin(b(x+s)) dx
    a, b = j*np.pi/l, k*np.pi/l
    def I(c, x): return np.where(np.abs(c) < 1e-14, x, np.sin(c*x)/np.where(np.abs(c)<1e-14,1,c))
    # sin(ax)sin(bx+bs) = 1/2[cos((a-b)x - bs) - cos((a+b)x + bs)]
    def prim(x):
        t1 = np.where(abs(a-b)<1e-14, x*np.cos(b*s), (np.sin((a-b)*x - b*s))/np.where(abs(a-b)<1e-14,1,(a-b)))
        t2 = np.sin((a+b)*x + b*s)/(a+b)
        return 0.5*(t1 - t2)
    return (2/l)*(prim(l - s) - prim(-l))
def Asym(j,k,s): return 0.5*(Ajk(j,k,s)+Ajk(k,j,s))
mus = []
for n,p in NSn:
    lam = np.log(n); ln = np.log(p)/np.sqrt(n)
    if abs(lam-1) < 0.5: mus.append((abs(lam-1), 2*ln))
    if abs(lam-2) < 0.5: mus.append((abs(lam-2), -0.5*ln))
s_, ws_ = gl(0, 0.5, 2000)
ka = lambda s: h(1+s) - 0.25*h(2+s)
def Qx(j, k):
    om_j, om_k = j*np.pi/l, k*np.pi/l
    d = 1.0 if j == k else 0.0          # <phi_j, phi_k>; A(0) = 2 delta
    main = 3*w0*d + 3*H5*d + 1.5*np.sum(wu*h(u)*(2*d - Asym(j,k,u)))
    arch_shift = np.sum(ws_*(ka(s_)+ka(-s_))*Asym(j,k,s_))
    prime = sum(c*Asym(j,k,np.array([m]))[0] for m,c in mus)
    return (main + arch_shift + prime)/(om_j**2*om_k**2), dict(main=main, ashift=arch_shift, prime=prime)
def Bpsi(j, k):
    Fj = S(j) - w_(j)*Fd; sj = sig(j) - w_(j)*sd
    Fk = S(k) - w_(k)*Fd; sk = sig(k) - w_(k)*sd
    return B(Fj, sj, Fk, sk), -2*bI2**2*sj*sk
rows = []
for (j,k) in [(13,13),(13,14),(13,20),(30,30),(30,31),(72,73),(73,73),(100,100),(100,103),(200,200)]:
    q, parts = Qx(j,k); b, pole = Bpsi(j,k)
    rows.append(dict(j=j,k=k,B_spectral=b, Q_x_noPole=q, Q_x=q+pole, rel=(q+pole-b)/abs(b), parts=parts))
    print(rows[-1], flush=True)
json.dump(dict(mus=mus, w0=w0, H5=H5, rows=rows), open('xspace_check_results.json','w'), indent=1, default=float)
