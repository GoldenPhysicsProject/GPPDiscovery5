"""New Bot: high-precision (mpmath dps 80, S87 exact-polynomial builder hp_forms.py, unchanged, imported read-only)
test of the F_theta onset law  lambda_min^even(Q_zeta)(L) = theta^2 * kappa(L),  kappa = Q_zeta^odd(x v)/||v||^2,
v = lowest even eigvec of Q_zeta.  Derived from the PROVED identity
  Q_F(f) = Q_zeta(f cosh(theta u)) - Q_zeta(f sinh(theta u))  (builder normalisation, kernel cosh(theta r) W_zeta(r)).
HIGH-PRECISION NUMERICAL (not interval).  RH OPEN.  Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, time
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 's87'))
import hp_forms as H, mpmath as mp
DPS = 80; H.setup(DPS); NB = int(sys.argv[1]) if len(sys.argv) > 1 else 24
def forms(form, L, N, th=0):
    F = H.build(form, L, N, theta=th)
    Qe = H.add(F['Ae'], F['Pe'], 1, -1); Qo = H.add(F['Ao'], F['Po'], 1, -1)
    for p in F['pe']: Qe = [[Qe[i][j] + p[i] * p[j] for j in range(N)] for i in range(N)]
    for p in F['po']: Qo = [[Qo[i][j] - p[i] * p[j] for j in range(N)] for i in range(N)]
    return Qe, Qo
def lmin(M): return min(mp.eigsy(mp.matrix(M), eigvals_only=True))
def zeta_data(L):
    N = NB + 1; Qe, Qo = forms('ZETA', L, N)
    E, V = mp.eigsy(mp.matrix([r[:NB] for r in Qe[:NB]]))
    k = min(range(NB), key=lambda i: E[i]); v = [V[i, k] for i in range(NB)] + [mp.mpf(0)]
    Lm = mp.mpf(L); xv = [mp.mpf(0)] * N
    for i in range(NB):   # x * phi_i, x = L u, orthonormal shifted Legendre, exact three-term recurrence
        xv[i] += Lm * v[i] / 2
        xv[i + 1] += Lm * v[i] * (i + 1) / (2 * mp.sqrt((2 * i + 1) * (2 * i + 3)))
        if i > 0: xv[i - 1] += Lm * v[i] * i / (2 * mp.sqrt((2 * i - 1) * (2 * i + 1)))
    kap = mp.fsum(xv[i] * Qo[i][j] * xv[j] for i in range(N) for j in range(N))
    lo = lmin([r[:NB] for r in Qo[:NB]])
    return E[k], kap, lo
out = {"kind": "high-precision numerical (mpmath dps %d, not interval)" % DPS, "N": NB, "zeta": [], "onsets": []}
Ls = [0.40 + 0.05 * i for i in range(19)]
for L in Ls:
    t = time.time(); le, ka, lo = zeta_data(str(round(L, 2)))
    out["zeta"].append({"L": round(L, 2), "lam_even": mp.nstr(le, 12), "kappa": mp.nstr(ka, 12), "lam_odd": mp.nstr(lo, 12)})
    print(out["zeta"][-1], round(time.time() - t, 1), flush=True)
json.dump(out, open('hp_onset_N%d.json' % NB, 'w'), indent=1)
# predicted onset by log-linear interpolation of log(lam_even/kappa) = 2 log theta
g = [(z["L"], mp.log(mp.mpf(z["lam_even"]) / mp.mpf(z["kappa"]))) for z in out["zeta"]]
def predict(th):
    tgt = 2 * mp.log(th)
    for (a, ga), (b, gb) in zip(g, g[1:]):
        if ga >= tgt >= gb: return a + (b - a) * (ga - tgt) / (ga - gb)
    return None
def fe_min(L, th):
    Qe, Qo = forms('FT', mp.mpf(L), NB, th); return lmin(Qe), lmin(Qo)
for th in ['0.25', '0.1', '0.01', '1e-3', '1e-5', '1e-8', '1e-12']:
    thm = mp.mpf(th); Lp = predict(thm)
    if Lp is None: print('no prediction', th); continue
    lo, hi = max(0.35, float(Lp) - 0.05), min(1.35, float(Lp) + 0.05)
    elo, _ = fe_min(lo, thm); ehi, _ = fe_min(hi, thm)
    if not (elo > 0 and ehi < 0): print('bracket fail', th, lo, hi, elo, ehi); out["onsets"].append({"theta": th, "pred": float(Lp), "bracket_fail": [lo, hi, mp.nstr(elo, 5), mp.nstr(ehi, 5)]}); continue
    for _ in range(11):
        m = (lo + hi) / 2; em, _ = fe_min(m, thm)
        if em > 0: lo = m
        else: hi = m
    m = (lo + hi) / 2; ee, eo = fe_min(m, thm)
    r = {"theta": th, "onset_even_direct": round(m, 5), "onset_even_predicted": round(float(Lp), 5), "F_odd_min_at_onset": mp.nstr(eo, 6)}
    out["onsets"].append(r); print(r, flush=True)
    json.dump(out, open('hp_onset_N%d.json' % NB, 'w'), indent=1)
