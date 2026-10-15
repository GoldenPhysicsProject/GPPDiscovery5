"""S22 [design, UNTRUSTED float LP]: positive odd-admissible symbol for zeta's pole-free form B on (-h,h).
m(xi) = sigma_A(xi) + C*wDC(xi) - sum_k a_k cos(u_k xi) >= c, u_k = 2h + k D,
sigma_A = Re psi(1/4 + i xi/2) - log pi - 2 kappa cos(xi log2),  kappa = log2/sqrt2,
wDC = FT of the trapezoid W (W=1 on |t|<=2h, linear to 0 at 2h+Ddc): wDC = 2 sin(a xi)/xi * sinc(Ddc xi/2), a=2h+Ddc/2.
For ODD g supported in (-h,h): int |ghat|^2 cos(u xi) = 0 (u>=2h) and int |ghat|^2 wDC = 2pi (int g)^2 = 0.
So B(g) = (1/2pi) int |ghat|^2 m  on odd window functions.  Rigour later (Arb)."""
import numpy as np, scipy.special as sp, json, sys, time
from scipy.optimize import linprog
l2 = np.log(2); kap = l2/np.sqrt(2)
def sigA(xi): return sp.psi(0.25+0.5j*xi).real-np.log(np.pi)-2*kap*np.cos(xi*l2)
def wDC(xi, h, Dd):
    a = 2*h+Dd/2; xi = np.asarray(xi, float)
    return 2*a*np.sinc(a*xi/np.pi)*np.sinc(Dd*xi/2/np.pi)
def lp(h, U=12.0, D=0.05, Dd=0.5, B=3.0, useDC=True, Xi=None):
    K = int(round((U-2*h)/D))+1; u = 2*h+D*np.arange(K)
    nph = 8
    if Xi is None:
        xs = np.arange(0, 5000, 0.01); low = sp.psi(0.25+0.5j*xs).real-np.log(np.pi)-2*kap
        Xi = xs[np.argmax(low > B/np.cos(np.pi/nph)+0.05)]
    xi = np.arange(0, Xi, 0.005); s = sigA(xi)
    cols = [np.ones((xi.size, 1))] + ([-wDC(xi, h, Dd)[:, None]] if useDC else []) + [np.cos(np.outer(xi, u))]
    A1 = np.hstack(cols); nd = 1 if useDC else 0
    th = np.linspace(0, 2*np.pi, 1500, endpoint=False)
    A2 = np.hstack([np.zeros((1500*nph, 1+nd)), np.vstack([np.cos(np.outer(th, np.arange(K))+ph) for ph in 2*np.pi*np.arange(nph)/nph])])
    bnds = [(None, None)]+[(0, None)]*nd+[(None, None)]*K
    r = linprog(np.r_[-1, np.zeros(nd+K)], A_ub=np.vstack([A1, A2]), b_ub=np.r_[s, np.full(1500*nph, B)], bounds=bnds, method="highs-ds")
    return r, u, Xi
if __name__ == "__main__" and sys.argv[1] != "save":
    for hs in sys.argv[1:]:
        h = np.log(3)/2 if hs == "hstar" else float(hs)
        for useDC in (False, True):
            t0 = time.time(); r, u, Xi = lp(h, useDC=useDC)
            print(f"h={hs} DC={useDC}: status={r.status} c_LP={-r.fun if r.status==0 else None} C={r.x[1] if (useDC and r.status==0) else '-'} Xi={Xi:.1f} ({time.time()-t0:.0f}s)", flush=True)
def save(h, hs, D=0.05, Dd=0.5, U=12.0, B=3.0, c_frac=0.5):
    from fractions import Fraction
    r, u, Xi = lp(h, U=U, D=D, Dd=Dd, B=B, useDC=True)
    cLP, C, a = -r.fun, r.x[1], r.x[2:]
    out = dict(h=hs, D=D, Dd=Dd, U=U, K=len(a), c_LP=cLP, C_num=int(round(C*10**15)), a_num=[int(round(x*10**15)) for x in a], den=10**15,
               c_target=float(cLP*c_frac), Xi_LP=float(Xi))
    json.dump(out, open(f"s22_cert_h{hs}.json", "w")); print("saved", hs, cLP, C, flush=True)
if __name__ == "__main__" and sys.argv[1] == "save":
    for hs in sys.argv[2:]:
        save(np.log(3)/2 if hs == "hstar" else float(hs), hs)
