"""Round 31, task 1: the prime torus T = prod_p S^1, its winding t -> (p^{-it})_p, and the Bohr/Bagchi picture.
No zeta zeros are inputs. Values of zeta / Hurwitz zeta are computed with python-flint (arb), 53-bit.
 T1  winding: continued fractions of log q/log p (how 'golden' the prime winding is), equidistribution discrepancy.
 T2  Bagchi/Bohr-Jessen limit law: time samples of log|Z(sigma+it)| vs the random Euler product on T.
 T3  'turns vs stretches': mean of log|Z(sigma+it)| at increasing heights (stationary = rotation; drift = stretch).
 T4  variance of the random model on T: sum_p |a_p|^2 p^{-2 sigma} -> finite iff sigma > 1/2 (+theta for F_theta).
 T5  strong recurrence (Bagchi 1981/1987) illustration: density of tau with max_K |Z(s+i tau)-Z(s)| < eps.
     K = circle of radius 0.05 (max principle) around s0; eps = 0.25 (0.5 at the base point s0 = 0.75).
     The DH disk around 0.808+85.70i is centred on a known off-line DH zero (Spira; Round 17): this is a control
     location only, no zero enters any argument.
Models on T: ZETA a_p = 1 (one circle, one species); FTH = zeta(s+th)zeta(s-th) (a_p = p^th + p^-th: two radii);
DH = Davenport-Heilbronn (no Euler product: a function on T but not a product over circles)."""
import json, math, sys, time, numpy as np
from flint import acb, arb, ctx
ctx.prec = 53
TH = 0.25
KAP = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
rng = np.random.default_rng(31)
out = {}

def Z(name, s):
    s = acb(s.real, s.imag)
    if name == 'ZETA': v = s.zeta()
    elif name == 'FTH': v = (s + TH).zeta()*(s - TH).zeta()
    elif name == 'DH':
        b = [1, KAP, -KAP, -1]
        v = sum((acb(5)**(-s))*bb*s.zeta(acb(a/5)) for a, bb in zip(range(1, 5), b))
    return complex(float(v.real.mid()), float(v.imag.mid()))

def primes(N):
    s = np.ones(N + 1, bool); s[:2] = False
    for p in range(2, int(N**0.5) + 1):
        if s[p]: s[p*p::p] = False
    return np.nonzero(s)[0]

# ---------- T1: the winding ----------
def cf(x, n=16):
    a = []
    for _ in range(n):
        q = math.floor(x); a.append(q); x = x - q
        if x < 1e-12: break
        x = 1/x
    return a
import mpmath
mpmath.mp.dps = 60
def cf_mp(x, n=18):
    a = []
    for _ in range(n):
        q = int(mpmath.floor(x)); a.append(q); x = x - q; x = 1/x
    return a
T1 = dict(log3_log2=cf_mp(mpmath.log(3)/mpmath.log(2)), log5_log2=cf_mp(mpmath.log(5)/mpmath.log(2)),
          log5_log3=cf_mp(mpmath.log(5)/mpmath.log(3)), golden=cf_mp((1 + mpmath.sqrt(5))/2))
# near-closures of the (2,3) winding: |a log2 - b log3| small  (Pythagorean comma 3^12 ~ 2^19, etc.)
near = []
for b in range(1, 700):
    a = round(b*math.log(3)/math.log(2)); d = abs(a*math.log(2) - b*math.log(3))
    if not near or d < near[-1][2]: near.append((a, b, d))
T1['record_near_closures_2a_3b'] = near
# discrepancy of the winding on T^2, T^3 (box test on an m^k grid), t sampled uniformly on [0, T]
disc = {}
for P in ([2, 3], [2, 3, 5]):
    lp = np.log(np.array(P, float))/(2*np.pi); m = 8; row = []
    for T in (10, 100, 1000, 10000, 100000):
        t = rng.uniform(0, T, 200000)
        pts = np.mod(np.outer(t, lp), 1.0); idx = np.floor(pts*m).astype(int)
        flat = np.ravel_multi_index(idx.T, (m,)*len(P)); cnt = np.bincount(flat, minlength=m**len(P))/len(t)
        row.append((T, float(np.abs(cnt*m**len(P) - 1).max())))
    disc[str(P)] = row
T1['box_discrepancy_relative'] = disc
out['T1'] = T1; print('T1', T1['log3_log2'], disc, flush=True)

# ---------- T2/T3: Bagchi limit law and 'turn vs stretch' ----------
PR = primes(200000); LP = np.log(PR.astype(float))
def random_model_logabs(sigma, nsamp, a_shift=(0.0,), P=PR):
    """log|prod_p prod_j (1 - w_p p^{-(sigma+shift_j)})^{-1}| for w uniform on T (Haar)."""
    res = np.zeros(nsamp); chunk = 2000
    for i in range(0, nsamp, chunk):
        w = np.exp(2j*np.pi*rng.random((min(chunk, nsamp - i), len(P))))
        acc = np.zeros(w.shape[0])
        for sh in a_shift:
            acc -= np.log(np.abs(1 - w*(P.astype(float)**-(sigma + sh))[None, :])).sum(1)
        res[i:i + w.shape[0]] = acc
    return res
T2 = {}
for sigma in (0.6, 0.75):
    tt = rng.uniform(1e4, 3e4, 4000)
    ts = np.array([math.log(abs(Z('ZETA', complex(sigma, t)))) for t in tt])
    rm = random_model_logabs(sigma, 4000)
    qs = [0.05, 0.25, 0.5, 0.75, 0.95]
    from scipy.stats import ks_2samp
    T2[f'ZETA_{sigma}'] = dict(time_mean=float(ts.mean()), time_var=float(ts.var()), model_mean=float(rm.mean()),
        model_var=float(rm.var()), time_q=np.quantile(ts, qs).tolist(), model_q=np.quantile(rm, qs).tolist(),
        ks=float(ks_2samp(ts, rm).statistic), ks_p=float(ks_2samp(ts, rm).pvalue))
    print('T2', sigma, T2[f'ZETA_{sigma}'], flush=True)
# T3: stationarity in height: mean log|Z| on [H, 1.1H] for H = 1e2..1e5
T3 = {}
for name, sigma in (('ZETA', 0.6), ('FTH', 0.6), ('FTH', 0.85), ('DH', 0.6)):
    row = []
    for H in (1e2, 1e3, 1e4, 1e5):
        tt = rng.uniform(H, 1.1*H, 1500)
        v = np.array([math.log(max(abs(Z(name, complex(sigma, t))), 1e-300)) for t in tt])
        row.append(dict(H=H, mean=float(v.mean()), var=float(v.var()), pred_drift=float(max(0.0, 0.5 + (TH if name == 'FTH' else 0) - sigma)*math.log(1.05*H/(2*math.pi)) if name == 'FTH' else 0.0)))
    T3[f'{name}_{sigma}'] = row; print('T3', name, sigma, [(r['H'], round(r['mean'], 3), round(r['pred_drift'], 3)) for r in row], flush=True)
out['T2'] = T2; out['T3'] = T3

# ---------- T4: variance of the random model (L^2 boundary on T) ----------
T4 = {}
for name, ap2 in (('ZETA', lambda p, s: p**(-2*s)), ('FTH', lambda p, s: (p**TH + p**-TH)**2*p**(-2*s))):
    T4[name] = {str(s): [(int(P), float(np.sum(ap2(PR[PR <= P].astype(float), s)))) for P in (1e2, 1e3, 1e4, 1e5, 2e5)]
                for s in (0.55, 0.6, 0.7, 0.8, 0.9)}
# sample variance of the truncated random model for FTH at sigma = 0.6 (second factor at 0.35: diverges)
T4['FTH_0.6_sample_var_vs_P'] = [(int(P), float(random_model_logabs(0.6, 1500, (TH, -TH), PR[PR <= P]).var())) for P in (1e2, 1e3, 1e4, 1e5)]
T4['ZETA_0.6_sample_var_vs_P'] = [(int(P), float(random_model_logabs(0.6, 1500, (0.0,), PR[PR <= P]).var())) for P in (1e2, 1e3, 1e4, 1e5)]
out['T4'] = T4; print('T4', T4['FTH_0.6_sample_var_vs_P'], T4['ZETA_0.6_sample_var_vs_P'], flush=True)

# ---------- T5: strong recurrence illustration ----------
def recurrence(name, s0, r, eps, Tmax, h=0.02, nb=12):
    K = [s0 + r*complex(math.cos(2*math.pi*j/nb), math.sin(2*math.pi*j/nb)) for j in range(nb)]
    F0c = Z(name, s0); F0 = [Z(name, s) for s in K]
    taus = np.arange(1.0, Tmax, h); good = []
    for tau in taus:
        if abs(Z(name, s0 + 1j*tau) - F0c) >= eps: continue
        if max(abs(Z(name, s + 1j*tau) - f) for s, f in zip(K, F0)) < eps: good.append(float(tau))
    g = np.array(good); row = []
    for T in (Tmax/8, Tmax/4, Tmax/2, Tmax):
        row.append((T, float(np.sum(g <= T)*h/T)))
    return row, good[:50]
T5 = {}
Tm = float(sys.argv[1]) if len(sys.argv) > 1 else 4000
for name, s0 in (('ZETA', complex(0.75, 30)), ('ZETA', complex(0.6, 30)), ('FTH', complex(0.6, 30)), ('FTH', complex(0.9, 30)), ('DH', complex(0.75, 30)), ('DH', complex(0.808, 85.70)), ('ZETA', complex(0.75, 0))):
    t0 = time.time(); row, g = recurrence(name, s0 + 0j, 0.05, 0.25 if s0.imag else 0.5, Tm)
    T5[f'{name}_{s0.real}_{s0.imag}'] = dict(density_by_T=row, first_taus=g[:10], secs=time.time() - t0)
    print('T5', name, s0, row, round(time.time() - t0), flush=True)
    out['T5'] = T5; json.dump(out, open('r31_torus.json', 'w'), indent=1, default=float)
json.dump(out, open('r31_torus.json', 'w'), indent=1, default=float)
