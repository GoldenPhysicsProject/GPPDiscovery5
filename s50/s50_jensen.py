# S50 (Grok) H50b: GORZ shifted Jensen polynomials J^{d,n}(X)=sum_j C(d,j) gam(n+j) X^j, gam(n)=n! [u^n] f(u),
# f(u)=xi(1/2+sqrt u)/xi(1/2)=prod(1+u/g_j^2) (genus 0), with a planted mid-gap off-line quartet (u-pair a=(d+ig)^2, abar).
# Log-domain elementary symmetric functions in tau=u/T, T=gp^2; tail beyond the zero list as exp(c_tail u).
import json, math, numpy as np, time
from scipy.special import gammaln, logsumexp
G = np.loadtxt('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt', usecols=1); G2 = G**2
Gmax = G[-1]; c_tail = (math.log(Gmax/(2*math.pi)) + 1)/(2*math.pi*Gmax)
def logE(T, K):
    lE = np.full(K + 1, -np.inf); lE[0] = 0.0
    for lw in np.log(T/G2):
        lE[1:] = np.logaddexp(lE[1:], lw + lE[:-1])
    lc = math.log(c_tail*T); j = np.arange(K + 1); lp = j*lc - gammaln(j + 1)
    out = np.array([logsumexp(lE[:n + 1][::-1] + lp[:n + 1]) for n in range(K + 1)])
    return out
def plant(lE, gp, d, T):
    A2 = (gp*gp + d*d)**2
    s1 = 2*(gp*gp - d*d)*T/A2; s2 = T*T/A2
    K = len(lE) - 1; out = np.empty(K + 1)
    for n in range(K + 1):
        terms = [lE[n]] + ([math.log(s1) + lE[n - 1]] if n >= 1 else []) + ([math.log(s2) + lE[n - 2]] if n >= 2 else [])
        out[n] = logsumexp(terms)
    return out
def plant_pair(lE, a, T):  # general complex u-pair a, abar: (1-u/a)(1-u/abar) = 1 + s1 tau + s2 tau^2 (s1 may be negative)
    s1 = -2*(1/a).real*T; s2 = T*T/abs(a)**2
    E = np.exp(lE - lE.max()); K = len(E) - 1
    out = E.copy(); out[1:] += s1*E[:-1]; out[2:] += s2*E[:-2]
    return out, lE.max()
def jensen_tests(lC):
    K = len(lC) - 1; lg = lC + gammaln(np.arange(K + 1) + 1)
    d2 = 1 - np.exp(lg[:-2] + lg[2:] - 2*lg[1:-1])  # >0 <=> J^{2,n} hyperbolic
    d3 = []
    for n in range(K - 2):
        A, B, C, D = np.exp(lg[n:n + 4] - lg[n + 1])
        a, b, c, d_ = D, 3*C, 3*B, A
        disc = 18*a*b*c*d_ - 4*b**3*d_ + b*b*c*c - 4*a*c**3 - 27*a*a*d_*d_
        d3.append(disc/(b*b*c*c))  # scale-free relative discriminant
    return d2, np.array(d3)
res = {'c_tail': c_tail}
for target in [85.7, 1000]:
    j = np.searchsorted(G, target); gp = 0.5*(G[j - 1] + G[j]); T = gp*gp; N = j
    K = int(5*N) + 60
    t0 = time.time(); lE = logE(T, K)
    # sanity: p2 = e1^2 - 2 e2 = sum g^-4 (Arb value 3.71725992852696861648662624717e-5) -- unscaled
    e1 = math.exp(lE[1])/T; e2 = math.exp(lE[2])/T**2
    d2, d3 = jensen_tests(lE)
    rec = dict(gamma=float(gp), N_below=int(N), K=K, sum_g4_check=e1*e1 - 2*e2, sum_g2=e1,
               xi_min_d2=float(d2.min()), xi_argmin_d2=int(d2.argmin()), xi_min_d3=float(d3.min()), xi_argmin_d3=int(d3.argmin()),
               secs=time.time() - t0, planted={})
    q2, q3 = jensen_tests(plant(lE, gp, 0.0, T))  # on-line double zero at mid-gap (delta=0) reference
    rec['online_double_min_d2'] = float(q2.min())
    for dl in [0.01, 0.1, 0.3, 0.5, 2.0, 10.0]:
        lP = plant(lE, gp, dl, T); p2, p3 = jensen_tests(lP)
        bad2 = np.where(p2 <= 0)[0]; bad3 = np.where(p3 <= 0)[0]
        rec['planted'][str(dl)] = dict(min_d2=float(p2.min()), n_fail_d2=bad2.tolist()[:20], min_d3=float(p3.min()), n_fail_d3=bad3.tolist()[:20],
                                       max_rel_change_d2=float(np.max(np.abs(p2 - d2)/np.abs(d2))),
                                       argmax_change_d2=int(np.argmax(np.abs(p2 - d2)/np.abs(d2))),
                                       max_rel_offline_vs_online_d2=float(np.max(np.abs(p2 - q2)/np.abs(q2))),
                                       max_rel_offline_vs_online_d3=float(np.max(np.abs(p3 - q3)/np.abs(q3))))
    res[str(target)] = rec
    print(json.dumps(rec), flush=True)
# control: a u-pair OUTSIDE the cone |Im a|<|Re a| (not a zeta-admissible location) must be detected at some shift
j = np.searchsorted(G, 85.7); gp = 0.5*(G[j - 1] + G[j]); T = gp*gp; lE = logE(T, 300)
ctrl = {}
for ang in [30, 44, 46, 60, 80]:
    a = 100*np.exp(1j*math.radians(180 - ang))  # arg(-a) = ang degrees, |a|=100 < g1^2
    E, sh = plant_pair(lE, a, T)
    if (E <= 0).any():
        ctrl[str(ang)] = 'negative Taylor coefficient'; continue
    p2, p3 = jensen_tests(np.log(E))
    ctrl[str(ang)] = dict(n_fail_d2=np.where(p2 <= 0)[0].tolist()[:10], min_d2=float(p2.min()), n_fail_d3=np.where(p3 <= 0)[0].tolist()[:10])
print('control', ctrl)
res['control_cone'] = ctrl
json.dump(res, open('s50_jensen.json', 'w'), indent=1, default=str)
