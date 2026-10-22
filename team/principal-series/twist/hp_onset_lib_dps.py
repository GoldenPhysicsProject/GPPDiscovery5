"""New Bot: high-precision (mpmath dps 80, S87 exact-polynomial builder hp_forms.py, unchanged, imported read-only)
test of the F_theta onset law  lambda_min^even(Q_zeta)(L) = theta^2 * kappa(L),  kappa = Q_zeta^odd(x v)/||v||^2,
v = lowest even eigvec of Q_zeta.  Derived from the PROVED identity
  Q_F(f) = Q_zeta(f cosh(theta u)) - Q_zeta(f sinh(theta u))  (builder normalisation, kernel cosh(theta r) W_zeta(r)).
HIGH-PRECISION NUMERICAL (not interval).  RH OPEN.  Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
import sys, json, time
sys.path.insert(0, '/workspace/GPPDiscovery5/s87')
import hp_forms as H, mpmath as mp
DPS = int(__import__("os").environ.get("DPS", "80")); H.setup(DPS); NB = int(sys.argv[1]) if len(sys.argv) > 1 else 24
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
