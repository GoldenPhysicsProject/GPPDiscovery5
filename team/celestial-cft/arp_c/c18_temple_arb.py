"""C.18: Arb certification of the Temple-Lehmann bound (Principal's Lemma P.17.A + P.17.B) for the ODD sector, applied to the
comparison form Q'_T (C.15/C.16), at 2L = log 3 and at 2L* = log 3 + 0.03; also writes the shared-entries JSON files.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. Temple idea: Golden Eye; Lemmas P.17.A/B: Principal (P.17);
kernel/tail formalism: Prime (N.17) and S&M. RH is OPEN.

Chain (ARP C.18): lambda_1^odd(Q_L) >= lambda_1(Q'_T)   [Q_L >= Q'_T as forms: q >= q_T on |t| >= T, C.15.A step 1]
                                  >= rho - Rbar^2/(nu2 - rho)   [Temple for the bounded s.a. operator of Q'_T on odd H]
  v = sum_{j<n} c_j f_j (dyadic c_j, exact), rho = c^T A c/|c|^2 (A = Gram of Q'_T),
  |(A'-rho)v|^2 = sum_{k<Mr} r_k^2 + |P_{W_Mr} A' v|^2,  r_k = (A_{k,:} c - rho c_k)/|c|,  |P_W A'v|^2 <= (c^T M c/|c|^2) eps^abs_Mr,
  nu2 <= lambda_2(Q'_T) via: ball LDL^T of A_N - (nu2+delta) I has every pivot sign-determined and <= 1 negative pivot
  (=> lambda_2(A_N) > nu2+delta), and c_tail - s eps^abs_N/delta >= nu2 with s I - M_N > 0 (interval Cholesky),
  c_tail = q_T - eps^neg_N (min-max on {xi perp phi_1(A_N)}; ARP C.18.B).
"""
import json, time, sys
import numpy as np
from c18_entries import *
from c16_odd_arb_run import chol_pd

def bstr(x):
    s = x.str(30, radius=True)
    if '+/-' not in s: return [s.strip('[]'), '0']
    m, r = s.strip('[]').split('+/-'); return [m.strip(), r.strip()]
def fl(X): return np.array([[float(x.mid()) for x in r] for r in X])

def run(tag, twoL, twoL_str, T=150, n=64, Mr=800, nu2=None, delta=None, write=True):
    t0 = time.time(); S = [(2, arb(2).log())] + ([(3, arb(3).log())] if twoL > arb(3).log() else [])
    E = Engine(twoL, S, 'odd', T, Mr)
    sq = range(n); rect = range(n, Mr)
    NOm = E.gram('Om', range(Mr), sq); B = {p: E.gram(p, range(Mr), sq) for p, _, _ in E.S}
    A, M = E.assemble(NOm, B)            # Mr x n blocks of the Gram matrices (rows k < Mr, cols j < n)
    eps = {K: E.eps(K) for K in (n, 128, 256, Mr)}
    print("== %s: 2L = %s, odd, primes %s, T = %d, n = N = %d, Mr = %d (entries %.0fs)" % (tag, twoL_str, [p for p, _, _ in E.S], T, n, Mr, time.time() - t0))
    print("q_T =", E.QT.str(20, radius=True), "; Omega(T) =", E.OMT.str(20, radius=True))
    for K in eps: print("eps^abs_%d = %s   eps^neg_%d = %s" % (K, eps[K][0].str(12, radius=True), K, eps[K][1].str(12, radius=True)))
    print("max entry radius (A, rows < %d): %.2e" % (Mr, max(float(x.rad()) for r in A for x in r)))
    An = [r[:] for r in A[:n]]; Mn = [r[:] for r in M[:n]]
    Af = fl(An); ev, V = np.linalg.eigh(Af)
    print("float spectrum A_n: lambda_1 = %.6e  lambda_2 = %.6f ; float lambda_max(M_n) = %.4f" % (ev[0], ev[1], np.linalg.eigvalsh(fl(Mn))[-1]))
    # -- trial vector with dyadic coefficients (exact) --
    cint = [int(round(x*2**48)) for x in V[:, 0]]; c = [arb(ci)/arb(2)**48 for ci in cint]
    cc = sum((x*x for x in c), arb(0))
    Ac = [sum((A[k][j]*c[j] for j in sq), arb(0)) for k in range(Mr)]
    rho = sum((c[j]*Ac[j] for j in sq), arb(0))/cc
    sv = sum((c[i]*sum((Mn[i][j]*c[j] for j in sq), arb(0)) for i in sq), arb(0))/cc
    r2 = sum(((Ac[k] - (rho*c[k] if k < n else 0))**2 for k in range(Mr)), arb(0))/cc
    tau2 = sv*eps[Mr][0]; R2 = r2 + tau2
    print("rho = %s ; s_v = c^TMc/|c|^2 = %s" % (rho.str(15, radius=True), sv.str(8, radius=True)))
    print("sum_{k<Mr} r_k^2 = %s ; tau^2 <= s_v eps^abs_Mr = %s ; Rbar = %s" % (r2.str(8, radius=True), tau2.str(8, radius=True), R2.sqrt().str(8, radius=True)))
    # -- second eigenvalue bound --
    sN = arb(str(float(np.linalg.eigvalsh(fl(Mn))[-1])*1.0001))
    okM, pM = chol_pd([[(sN if i == j else 0) - Mn[i][j] for j in sq] for i in sq])
    ctail = E.QT - eps[n][1]; beta2 = sN*eps[n][0]
    if nu2 is None:   # float optimisation of delta, then round nu2 down
        lam2 = ev[1] - 1e-6; cf = float(ctail.mid()); bf = float(beta2.mid())
        ds = np.linspace(1e-4, lam2, 20000); val = np.minimum(lam2 - ds, cf - bf/ds); i = int(np.argmax(val))
        delta = float('%.4f' % ds[i]); nu2 = np.floor((min(lam2 - delta, cf - bf/delta) - 2e-4)*1000)/1000
    nu2 = arb(str(nu2)); delta = arb(str(delta)); sig = nu2 + delta
    okI, nneg, D = ldlt_inertia([[An[i][j] - (sig if i == j else 0) for j in sq] for i in sq])
    tailok = bool(ctail - beta2/delta - nu2 >= 0)
    print("s_V = %s: interval Cholesky s I - M_N > 0: %s" % (sN.str(8), okM))
    print("c_tail = q_T - eps^neg_N = %s ; beta^2 = s eps^abs_N = %s" % (ctail.str(10, radius=True), beta2.str(8, radius=True)))
    print("nu2 = %s, delta = %s: LDL^T(A_N - (nu2+delta) I) pivots sign-determined: %s, #negative = %s ; c_tail - beta^2/delta - nu2 = %s >= 0: %s"
          % (nu2.str(6), delta.str(6), okI, nneg, (ctail - beta2/delta - nu2).str(8, radius=True), tailok))
    gap_ok = bool(rho < nu2)
    mu = rho - R2/(nu2 - rho)
    cert = okM and okI and nneg is not None and nneg <= 1 and tailok and gap_ok and bool(mu > 0)
    print("rho < nu2: %s ; Temple: lambda_1^odd >= mu = %s ; lower endpoint %.6e ; CERTIFIED: %s  (%.0fs)"
          % (gap_ok, mu.str(12, radius=True), float(mu.lower()), cert, time.time() - t0))
    out = {'about': 'C.18 shared entries (Arb balls: [mid, rad] strings, enclosure |x - mid| <= rad). See c18_shared_entries_spec.md. '
                    'Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.',
           'twoL': twoL_str, 'L_ball': bstr(E.L), 'parity': 'odd', 'basis': 'f_k = sin(w_k x)/(sqrt(L) n_k), w_k=(k+1)pi/L, n_k=sqrt(w_k^2+1/4), k>=0',
           'T': T, 'n_square': n, 'M_rows': Mr, 'quadrature': {'panel_h': 1, 'gauss_nodes': 40, 'ellipse_b': '0.4', 'prec_bits': ctx.prec},
           'primes': [{'n': p, 'c_n': bstr(cn)} for p, _, cn in E.S], 'Omega_T': bstr(E.OMT), 'q_T': bstr(E.QT),
           'NOm_rows': [[bstr(x) for x in r] for r in NOm], 'B_rows': {str(p): [[bstr(x) for x in r] for r in B[p]] for p in B},
           'A_square': [[bstr(x) for x in r] for r in An], 'M_square': [[bstr(x) for x in r] for r in Mn],
           'eps': {str(K): {'abs': bstr(eps[K][0]), 'neg': bstr(eps[K][1])} for K in eps},
           'temple': {'c_dyadic_int': cint, 'c_scale': '2^-48', 'rho': bstr(rho), 's_v': bstr(sv), 'r2_head': bstr(r2), 'tau2': bstr(tau2),
                      'nu2': nu2.str(6), 'delta': delta.str(6), 's_V': sN.str(10), 'mu': bstr(mu), 'certified': cert}}
    if write: json.dump(out, open('c18_shared_entries_%s.json' % tag, 'w'))
    return cert, mu

if __name__ == '__main__':
    run('odd_log3', arb(3).log(), 'log 3')
    run('odd_log3p03', arb(3).log() + arb('0.03'), 'log 3 + 0.03')
    run('odd_log3p06', arb(3).log() + arb('0.06'), 'log 3 + 0.06 (extension, no JSON)', write=False)
