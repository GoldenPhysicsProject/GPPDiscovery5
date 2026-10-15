"""S20 [E]: (UNTRUSTED, floating LP) dual cosine+sine window certificate for a COMPLEX character, complex g.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH OPEN, GRH OPEN.
sigma_chi(xi) = Re psi((1+2eps)/4 + i xi/2) + log(f/pi) - 2 kappa Re(chi(2) e^{i xi log2})   (not even in xi).
Find c, a_k, b_k (u_k = log3 + kD) with sigma - c - sum a_k cos(u_k xi) + b_k sin(u_k xi) >= 0 on all real xi, and
Re(e^{i phi} P(e^{i th})) <= B (P(z) = sum (a_k - i b_k) z^k).  sin/cos at u >= 2h annihilate |ghat|^2 for supp g in (-h,h).
Covering all real xi also covers the conjugate character (xi -> -xi)."""
import numpy as np, scipy.special as sp, json, time, sys
from fractions import Fraction
from scipy.optimize import linprog
t = np.log(2); kap = t/np.sqrt(2); L = np.log(3)
CFG = {"chi7cubic": dict(f=7, eps=0, chi2=(-0.5, np.sqrt(3)/2), chi2_exact="exp(2 pi i/3)", U=12, D=Fraction(1, 20), B=3.0, c_frac=0.5, dxi=0.005),
       "chi7cubic_small": dict(f=7, eps=0, chi2=(-0.5, np.sqrt(3)/2), chi2_exact="exp(2 pi i/3)", U=7, D=Fraction(1, 10), B=2.0, c_frac=0.5, dxi=0.01),
       "chi7cubic_med": dict(f=7, eps=0, chi2=(-0.5, np.sqrt(3)/2), chi2_exact="exp(2 pi i/3)", U=10, D=Fraction(1, 10), B=3.6, c_frac=0.5, dxi=0.01)}

for _nm, _f, _eps, _p, _q in [("chi5quartic", 5, 1, 1, 4), ("chi9cubic", 9, 0, 1, 3), ("chi11quintic1", 11, 0, 1, 5), ("chi11quintic2", 11, 0, 2, 5), ("chi7cubic_med2", 7, 0, 1, 3)]:
    CFG[_nm] = dict(f=_f, eps=_eps, chi2=(np.cos(2*np.pi*_p/_q), np.sin(2*np.pi*_p/_q)), chi2_exact=f"exp(2 pi i {_p}/{_q})", chi2_angle=[_p, _q],
                    U=10, D=Fraction(1, 10), B=3.6 if _eps == 0 else 3.0, c_frac=0.5, dxi=0.01)
nm = sys.argv[1] if len(sys.argv) > 1 else "chi7cubic"; cf = CFG[nm]
f, eps, D = cf["f"], cf["eps"], float(cf["D"]); cr, ci = cf["chi2"]
def sigma(xi): return sp.psi((1+2*eps)/4+0.5j*xi).real+np.log(f/np.pi)-2*kap*(cr*np.cos(xi*t)-ci*np.sin(xi*t))
def siglow(xi): return sp.psi((1+2*eps)/4+0.5j*xi).real+np.log(f/np.pi)-2*kap
K = int(round((cf["U"]-L)/D))+1; u = L+D*np.arange(K)
nph = 8; Bp = cf["B"]/np.cos(np.pi/nph)
xs = np.arange(0, 3000, 0.01); Xi = xs[np.argmax(siglow(xs) > Bp+0.05)]
xi = np.arange(-Xi, Xi+cf["dxi"], cf["dxi"]); s = sigma(xi)
A1 = np.hstack([np.ones((xi.size, 1)), np.cos(np.outer(xi, u)), np.sin(np.outer(xi, u))])
th = np.linspace(0, 2*np.pi, 800, endpoint=False); kk = np.arange(K)
# Re(e^{i ph} sum (a_k - i b_k) e^{i k th}) = sum a_k cos(k th+ph) + b_k sin(k th+ph)
A2 = np.vstack([np.hstack([np.zeros((th.size, 1)), np.cos(np.outer(th, kk)+ph), np.sin(np.outer(th, kk)+ph)]) for ph in 2*np.pi*np.arange(nph)/nph])
t0 = time.time()
r = linprog(np.r_[-1, np.zeros(2*K)], A_ub=np.vstack([A1, A2]), b_ub=np.r_[s, np.full(A2.shape[0], cf["B"])],
            bounds=[(None, None)]*(2*K+1), method="highs-ds")
cLP = -r.fun; a, b = r.x[1:K+1], r.x[K+1:]
c = Fraction(cLP*cf["c_frac"]).limit_denominator(10**4)
aq = [int(round(x*10**15)) for x in a]; bq = [int(round(x*10**15)) for x in b]
xx = np.arange(-Xi, Xi, 0.0005); ff = sigma(xx)-float(c)-np.cos(np.outer(xx, u))@(np.array(aq)/1e15)-np.sin(np.outer(xx, u))@(np.array(bq)/1e15)
out = dict(name=nm, f=f, eps=eps, chi2_exact=cf["chi2_exact"], chi2_angle=cf.get("chi2_angle", [1, 3]), L="log(3)", D=[cf["D"].numerator, cf["D"].denominator], K=K,
           a_num=aq, b_num=bq, a_den=10**15, c=[c.numerator, c.denominator], c_LP=cLP, B_LP=cf["B"], Xi_LP=float(Xi),
           float_min_f=float(ff.min()), float_argmin=float(xx[np.argmin(ff)]), lp_time=time.time()-t0)
json.dump(out, open(f"s20_cert_{nm}.json", "w"))
print(f"{nm}: K={K} c_LP={cLP:.6f} c_cert={c} ({float(c):.6f}) Xi={Xi:.1f} float min f={ff.min():.6f} @ {xx[np.argmin(ff)]:.3f} ({time.time()-t0:.0f}s)", flush=True)
