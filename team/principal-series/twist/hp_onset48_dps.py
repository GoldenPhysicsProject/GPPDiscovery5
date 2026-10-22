"""New Bot S2: onset law at converged N (>= 48), both sectors, L = half-width 0.5..1.0.
Even: F_theta even form first fails where lam_min^even(Q_zeta) = theta^2 kappa_e, kappa_e = Q_zeta^odd(x v_e)/|v_e|^2.
Odd (new): where lam_min^odd(Q_zeta) = theta^2 kappa_o, kappa_o = Q_zeta^even(x v_o)/|v_o|^2.
(From the PROVED identity Q_F(f) = Q_zeta(f cosh th u) - Q_zeta(f sinh th u), first order in theta; f sinh(th u) ~ th x f flips parity.)
Direct check: lam_min of F_theta (builder 'FT') at theta = 0.9 and 1.1 x predicted theta_e / theta_o.
mpmath dps 80 (S87 builder hp_forms.py, unchanged). HIGH-PRECISION NUMERICAL. RH OPEN."""
import sys, json, os, time
sys.argv = [sys.argv[0], sys.argv[1]] + sys.argv[2:]
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hp_onset_lib_dps as H, mpmath as mp
NB = H.NB
def xmul(v, L, N):
    xv = [mp.mpf(0)] * N; Lm = mp.mpf(L)
    for i in range(len(v)):
        xv[i] += Lm * v[i] / 2
        xv[i + 1] += Lm * v[i] * (i + 1) / (2 * mp.sqrt((2 * i + 1) * (2 * i + 3)))
        if i > 0: xv[i - 1] += Lm * v[i] * i / (2 * mp.sqrt((2 * i - 1) * (2 * i + 1)))
    return xv
def low(M):
    E, V = mp.eigsy(mp.matrix(M)); k = min(range(len(M)), key=lambda i: E[i]); return E[k], [V[i, k] for i in range(len(M))]
Ls = sys.argv[2].split(',') if len(sys.argv) > 2 else ['0.5', '0.6', '0.7', '0.8', '0.9', '1.0']
for Ls_ in Ls:
    t0 = time.time(); N = NB + 1; Qe, Qo = H.forms('ZETA', mp.mpf(Ls_), N)
    le, ve = low([r[:NB] for r in Qe[:NB]]); lo, vo = low([r[:NB] for r in Qo[:NB]])
    xe = xmul(ve, Ls_, N); xo = xmul(vo, Ls_, N)
    ke = mp.fsum(xe[i] * Qo[i][j] * xe[j] for i in range(N) for j in range(N))
    ko = mp.fsum(xo[i] * Qe[i][j] * xo[j] for i in range(N) for j in range(N))
    the = mp.sqrt(le / ke) if le > 0 and ke > 0 else None; tho = mp.sqrt(lo / ko) if lo > 0 and ko > 0 else None
    r = dict(N=NB, L=Ls_, lam_even=mp.nstr(le, 8), kappa_e=mp.nstr(ke, 8), theta_e_pred=mp.nstr(the, 8) if the else None,
             lam_odd=mp.nstr(lo, 8), kappa_o=mp.nstr(ko, 8), theta_o_pred=mp.nstr(tho, 8) if tho else None)
    for tag, th, sec in (('e', the, 0), ('o', tho, 1)):
        if th is None: continue
        for f in ('0.9', '1.1'):
            QFe, QFo = H.forms('FT', mp.mpf(Ls_), NB, th * mp.mpf(f))
            r['F%s_min_at_%sxpred' % (tag, f)] = mp.nstr(H.lmin(QFe if sec == 0 else QFo), 6)
    r['sec'] = round(time.time() - t0); print(json.dumps(r), flush=True)
