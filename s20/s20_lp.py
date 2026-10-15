"""S20 [B]: generate (UNTRUSTED, floating LP) dual 'Krein/Beurling-Selberg-type' window certificates for chi3, chi5q.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH OPEN, GRH OPEN.
Find c and a_k (u_k = log3 + k D, k=0..K-1) with
   f(xi) = sigma_chi(xi) - c - sum_k a_k cos(u_k xi) >= 0   for all real xi,
plus |P(e^{i th})| <= B for P(z) = sum a_k z^k (periodic sup control of the tail).
Since every u_k >= log3 = 2h*, the cosine terms annihilate |ghat|^2 for supp g in (-h*,h*) => Q_chi >= c ||g||^2.
Output s20_cert_<name>.json (rational a_k); ALL rigour is in s20_verify_arb.py."""
import numpy as np, scipy.special as sp, json, time, sys
from fractions import Fraction
from scipy.optimize import linprog
t = np.log(2); kap = t/np.sqrt(2); L = np.log(3)
def sigma(xi, f, eps, chi2): return sp.psi((1+2*eps)/4+0.5j*xi).real+np.log(f/np.pi)-2*kap*chi2*np.cos(xi*t)
def siglow(xi, f, eps): return sp.psi((1+2*eps)/4+0.5j*xi).real+np.log(f/np.pi)-2*kap
CFG = {"chi3": dict(f=3, eps=1, chi2=-1, U=12, D=Fraction(1, 20), B=3.0, c_frac=0.5),
       "chi5q": dict(f=5, eps=0, chi2=-1, U=12, D=Fraction(1, 20), B=3.0, c_frac=0.5)}
for nm in (sys.argv[1:] or CFG):
    cf = CFG[nm]; f, eps, chi2, D = cf["f"], cf["eps"], cf["chi2"], float(cf["D"])
    K = int(round((cf["U"]-L)/D))+1; u = L+D*np.arange(K)
    nph = 8; Bp = cf["B"]/np.cos(np.pi/nph)
    xs = np.arange(0, 3000, 0.01); Xi = xs[np.argmax(siglow(xs, f, eps) > Bp+0.05)]
    xi = np.arange(0, Xi+0.005, 0.005); s = sigma(xi, f, eps, chi2)
    A1 = np.hstack([np.ones((xi.size, 1)), np.cos(np.outer(xi, u))])
    th = np.linspace(0, 2*np.pi, 1500, endpoint=False)
    A2 = np.hstack([np.zeros((1500*nph, 1)), np.vstack([np.cos(np.outer(th, np.arange(K))+ph) for ph in 2*np.pi*np.arange(nph)/nph])])
    t0 = time.time()
    r = linprog(np.r_[-1, np.zeros(K)], A_ub=np.vstack([A1, A2]), b_ub=np.r_[s, np.full(1500*nph, cf["B"])],
                bounds=[(None, None)]*(K+1), method="highs-ds")
    cLP, a = -r.fun, r.x[1:]
    c = Fraction(cLP*cf["c_frac"]).limit_denominator(10**4)
    aq = [int(round(x*10**15)) for x in a]
    # float diagnostics
    xx = np.arange(0, Xi, 0.0005); ff = sigma(xx, f, eps, chi2)-float(c)-np.cos(np.outer(xx, u))@(np.array(aq)/1e15)
    tt = np.linspace(0, np.pi, 200001); Pm = np.abs(np.polyval((np.array(aq)/1e15)[::-1], np.exp(1j*tt))).max()
    out = dict(name=nm, f=f, eps=eps, chi2=chi2, L="log(3)", D=[cf["D"].numerator, cf["D"].denominator], K=K,
               a_num=aq, a_den=10**15, c=[c.numerator, c.denominator], c_LP=cLP, B_LP=cf["B"], Xi_LP=float(Xi),
               float_min_f=float(ff.min()), float_argmin=float(xx[np.argmin(ff)]), float_supP=float(Pm), lp_time=time.time()-t0)
    json.dump(out, open(f"s20_cert_{nm}.json", "w"))
    print(f"{nm}: K={K} c_LP={cLP:.6f} c_cert={c} ({float(c):.6f}) Xi={Xi:.1f} sum|a|={np.abs(a).sum():.2f} "
          f"float min f={ff.min():.6f} @ {xx[np.argmin(ff)]:.3f}  sup|P|~{Pm:.4f}  ({time.time()-t0:.0f}s)", flush=True)
