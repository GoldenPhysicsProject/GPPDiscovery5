"""Grok S61 H61a(4): BUDGET TABLE of Codex's energy-transport identity for the true zeta Weil-form minimizer.
Ritz space: Legendre P_j(x/L), j < M, EDGEK = 0 (f in L^2[-L,L], as S57/S59 knockouts), split by parity.
Independent assembly in Codex's form (all Arb, prec 1024):
   W = Pole - (2S_L - omega0) N + E_gamma + sum_n w_n (2N - 2 T_n),   T_n = Re g(log n), N = Gram,
   Pole_jk = l+_j l-_k + l-_j l+_k,  E_gamma = 2 int_0^inf k(u) (N - Re g(u)) du, k(u) = e^{-u/2}/(1-e^{-2u}).
Cross-check: W against the S14/S53 Gauss-digamma assembly build_all (different arch representation).
For the minimizer v (v'Nv = 1): |E|^2, |O|^2, 2S_L - omega0, E_gamma, E_prime (split by n), and the signed prime terms
-2 w_n T_n. NUMERICAL (Arb assembly, mpmath eig at dps 150). RH OPEN."""
import os, sys, json, time
os.environ['EDGEK'] = '0'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's56')); sys.path.insert(1, os.path.join(HERE, '..', 's53'))
from s14_exact import G as Gpoly, legendre
from flint import arb, acb, ctx, fmpq
import mpmath as mp
ctx.prec = 1024; mp.mp.dps = 150
def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else None
def toarb(c): return arb(c.p)/c.q
def tomp(x): return mp.mpf(x.mid().str(170, radius=False))
M = int(sys.argv[1]); Ls = sys.argv[2].split(','); XCHECK = sys.argv[3].split(',') if len(sys.argv) > 3 else []
t0 = time.time()
Gs = {}
for j in range(M):
    for k in range(j, M):
        if (j+k) % 2 == 0: g = (Gpoly(j, k)+Gpoly(k, j))/2; Gs[j, k] = Gs[k, j] = [toarb(c) for c in g.coeffs()]
D = max(len(c) for c in Gs.values())
Leg = [[toarb(c) for c in legendre(j).coeffs()] for j in range(M)]
print(json.dumps(dict(stage='G polys', M=M, maxdeg=D-1, wall_s=round(time.time()-t0))), flush=True)
def one(Lstr):
    t1 = time.time(); Lm = mp.mpf(eval(Lstr, {'log': mp.log})); h = arb(mp.nstr(Lm, 60)); hc = acb(h)
    ns = [n for n in range(2, int(mp.exp(2*Lm))+2) if vm(n) and mp.log(n) < 2*Lm]
    w = {n: arb(vm(n)).log()/arb(n).sqrt() for n in ns}; S = sum(w.values(), arb(0))
    om0 = acb(fmpq(1, 4)).digamma().real - arb.pi().log()
    # moments mu_m = int_0^2 v^(m-1) e^{hv/2}/(2 sinc(i h v)) dv  (= int_0^2 v^m k(hv) dv * h)
    mu = [None]+[acb.integral(lambda z, an, m=m: z**(m-1)*(hc*z/2).exp()/(2*(acb(0, 1)*hc*z).sinc()), 0, 2).real for m in range(1, D)]
    Kt = sum(((-(2*jj+arb(1)/2)*2*h).exp()/(2*jj+arb(1)/2) for jj in range(400)), arb(0))   # tail to e^{-800h}
    ep = [acb.integral(lambda z, an, m=m: z**m*(hc*z/2).exp(), -1, 1).real for m in range(M)]
    em = [acb.integral(lambda z, an, m=m: z**m*(-hc*z/2).exp(), -1, 1).real for m in range(M)]
    lp = [h*sum((Leg[j][m]*ep[m] for m in range(len(Leg[j]))), arb(0)) for j in range(M)]
    lm = [h*sum((Leg[j][m]*em[m] for m in range(len(Leg[j]))), arb(0)) for j in range(M)]
    res = dict(L=Lstr, Lnum=mp.nstr(Lm, 10), M=M, prime_powers=ns, S_L=mp.nstr(tomp(S), 12), omega0=mp.nstr(tomp(om0), 15),
               rev_const=mp.nstr(tomp(2*S-om0), 12))
    for sec, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        n_ = len(I)
        N = mp.matrix(n_, n_); P = mp.matrix(n_, n_); Eg = mp.matrix(n_, n_); T = {n: mp.matrix(n_, n_) for n in ns}
        for a, j in enumerate(I):
            for b, k in enumerate(I):
                co = Gs[j, k]; Nj = h*co[0]
                N[a, b] = tomp(Nj); P[a, b] = tomp(lp[j]*lm[k]+lm[j]*lp[k])
                Eg[a, b] = tomp(-2*h*sum((co[m]*mu[m] for m in range(1, len(co))), arb(0)) + 2*Nj*Kt)
                for n in ns:
                    v = arb(n).log()/h; T[n][a, b] = tomp(h*sum((co[m]*v**m for m in range(len(co))), arb(0)))
        rc = tomp(2*S-om0); wm = {n: tomp(w[n]) for n in ns}
        W = P - rc*N + Eg
        for n in ns: W += wm[n]*(2*N - 2*T[n])
        d = [1/mp.sqrt(N[i, i]) for i in range(n_)]       # Legendre: N diagonal
        offd = max(abs(N[i, j]) for i in range(n_) for j in range(n_) if i != j)
        Sm = mp.matrix(n_, n_)
        for i in range(n_):
            for j in range(n_): Sm[i, j] = d[i]*W[i, j]*d[j]
        ev, Q = mp.eigsy(Sm); order = sorted(range(n_), key=lambda i: ev[i]); i0 = order[0]
        x = mp.matrix([d[i]*Q[i, i0] for i in range(n_)])
        q = lambda X: (x.T*X*x)[0]
        lpv = mp.matrix([tomp(lp[j]) for j in I]); lmv = mp.matrix([tomp(lm[j]) for j in I])
        Ev = ((lpv+lmv).T*x)[0]/mp.sqrt(2); Ov = ((lpv-lmv).T*x)[0]/mp.sqrt(2)
        Eprime = {n: wm[n]*q(2*N-2*T[n]) for n in ns}; EP = sum(Eprime.values())
        signed = {n: -2*wm[n]*q(T[n]) for n in ns}
        egam = q(Eg); E2, O2 = Ev**2, Ov**2
        lam_rec = E2 - O2 - rc*q(N) + egam + EP
        fw = E2 + egam + EP; rv = O2 + rc
        r = dict(lam_min=mp.nstr(ev[i0], 10), lam_2=mp.nstr(ev[order[1]], 6), lam_from_budget_minus_lam=mp.nstr(lam_rec-ev[i0], 3),
                 norm=mp.nstr(q(N), 12), Gram_offdiag_max=mp.nstr(offd, 3),
                 absE2=mp.nstr(E2, 10), absO2=mp.nstr(O2, 10), E_gamma=mp.nstr(egam, 10), E_prime=mp.nstr(EP, 10),
                 E_prime_by_n={n: mp.nstr(v, 8) for n, v in Eprime.items()},
                 signed_prime_by_n={n: mp.nstr(v, 8) for n, v in signed.items()},
                 forward=mp.nstr(fw, 12), reverse=mp.nstr(rv, 12),
                 share_forward={"E2": mp.nstr(E2/fw, 5), "E_gamma": mp.nstr(egam/fw, 5), "E_prime": mp.nstr(EP/fw, 5),
                                **{"E_%d" % n: mp.nstr(Eprime[n]/fw, 5) for n in ns}},
                 Eprime_over_2S=mp.nstr(EP/(2*tomp(S)), 8), Eprime_over_4S=mp.nstr(EP/(4*tomp(S)), 8),
                 arch_total=mp.nstr(tomp(om0)+egam, 10), prime_total_signed=mp.nstr(sum(signed.values()), 10), pole=mp.nstr(E2-O2, 10),
                 coef_mass_first4=mp.nstr(sum(Q[i, i0]**2 for i in range(min(4, n_))), 5))
        if Lstr in XCHECK:
            os.chdir(os.path.join(HERE, '..', 's56'))
            _av = sys.argv; sys.argv = ['x', '0', '0', '150']
            from s53_odd_floor import build_all, weight
            sys.argv = _av; mp.mp.dps = 150
            A, P2, T2, ns2 = build_all(M, mp.nstr(Lm, 60), 1024)
            W2 = mp.matrix(n_, n_)
            for a, j in enumerate(I):
                for b, k in enumerate(I):
                    W2[a, b] = tomp(A[j, k]) + tomp(P2[j, k]) - sum(weight(n)*tomp(T2[n][j, k]) for n in ns2)
            os.chdir(HERE)
            r["xcheck_S14_maxabs_W_diff"] = mp.nstr(max(abs(W[i, j]-W2[i, j]) for i in range(n_) for j in range(n_)), 3)
            r["xcheck_S14_maxabs_W"] = mp.nstr(max(abs(W2[i, j]) for i in range(n_) for j in range(n_)), 3)
        res[sec] = r
    res['wall_s'] = round(time.time()-t1)
    return res
for Lstr in Ls:
    print(json.dumps(one(Lstr)), flush=True)
