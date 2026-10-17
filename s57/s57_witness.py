"""Grok S57 H57b(i): Arb certificate of a knockout witness.
f(x) = sum_j c_j P_j(x/L)(1-(x/L)^2)^EDGEK, j odd (odd sector) / even, supported in [-L,L], L rational.
Eval 1 (this file, Arb ball arithmetic): the S14/S53 exact assembly (acb.integral moments, exact fmpq autocorrelation polys)
gives arb matrices Arch, Pole, T_n. With RATIONAL c (rounded Ritz minimizer of the knocked-out form), we compute the arb
balls Q_ko(c) = c'(Arch+Pole - sum_{n not a power of p} w_n T_n)c, Q_zeta(c), w_p-term and Gram ||f||^2 = c'Gc
(G_jk = h*G_jk(0) exactly). Negative upper end of Q_ko => CERTIFIED witness (modulo the S14 formula, Lean-open).
usage: EDGEK=1 python s57_witness.py SECTOR p M L_rational [out.json]"""
import sys, os, json, time
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 's56')); sys.path.insert(1, os.path.join(HERE, '..', 's53'))
os.chdir(os.path.join(HERE, '..', 's56'))
from s53_odd_floor import build_all, tomp, weight, vonmangoldt
from flint import arb, arb_mat, fmpq, ctx
import mpmath as mp
sec, p, M, Ls = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
out = sys.argv[5] if len(sys.argv) > 5 else None
Lq = Fraction(Ls); mp.mp.dps = 130; t0 = time.time()
I = list(range(0 if sec == 'even' else 1, M, 2)); n = len(I)
A, P, T, ns = build_all(M, '%d/%d' % (Lq.numerator, Lq.denominator), 1024)
ctx.prec = 1024
def warb(nn):
    q = vonmangoldt(nn); return 2 * arb(q).log() / arb(nn).sqrt()
def sub(Ma): return arb_mat([[Ma[i, j] for j in I] for i in I])
Ar, Pr = sub(A), sub(P); Tr = {m: sub(T[m]) for m in ns}
ko = [m for m in ns if vonmangoldt(m) == p]
Wz = Ar + Pr
for m in ns: Wz = Wz - warb(m) * Tr[m]
Wko = Wz
for m in ko: Wko = Wko + warb(m) * Tr[m]
# Ritz minimizer (midpoint, mpmath) of the knocked-out form relative to the Gram matrix
from s14_exact import G as Gpoly
Gr = arb_mat(n, n)
for a, i in enumerate(I):
    for b, j in enumerate(I):
        g = (Gpoly(i, j) + Gpoly(j, i)) / 2; c0 = g.coeffs()[0]
        Gr[a, b] = arb(Lq.numerator) / Lq.denominator * arb(c0.p) / c0.q
def mid(X): return mp.matrix([[mp.mpf(X[i, j].mid().str(150, radius=False)) for j in range(n)] for i in range(n)])
Wm, Gm = mid(Wko), mid(Gr)
Lc = mp.cholesky(Gm); Li = mp.inverse(Lc)
S = Li * Wm * Li.T; S = (S + S.T) / 2
ev, Q = mp.eigsy(S); i0 = min(range(n), key=lambda i: ev[i])
v = Li.T * Q[:, i0]; v = v / mp.norm(v, mp.inf)
c = [Fraction(mp.nstr(v[i], 60)).limit_denominator(10**55) for i in range(n)]
ca = arb_mat([[arb(x.numerator) / x.denominator] for x in c])
def qf(X): return (ca.transpose() * X * ca)[0, 0]
Qko, Qz, Gn = qf(Wko), qf(Wz), qf(Gr)
pterm = sum((warb(m) * qf(Tr[m]) for m in ko), arb(0))
res = dict(sector=sec, p=p, M=M, EDGEK=os.environ.get('EDGEK', '4'), L=Ls, delta=float(Lq) - 0.5 * float(mp.log(p)),
           ritz_lam_ko=mp.nstr(ev[i0], 12), ritz_lam2_ko=mp.nstr(sorted(ev)[1], 12),
           Q_ko=Qko.str(15), Q_ko_upper=float((Qko.mid() + Qko.rad())), Q_ko_over_norm=mp.nstr(mp.mpf(Qko.mid().str(40, radius=False)) / mp.mpf(Gn.mid().str(40, radius=False)), 12),
           Q_zeta=Qz.str(15), p_terms_w_g=pterm.str(15), norm2=Gn.str(15),
           certified_negative=bool(Qko < 0), knocked=ko, prime_powers=ns, wall_s=round(time.time() - t0, 1),
           coeffs={str(i): str(x) for i, x in zip(I, c)})
s = json.dumps(res)
if out: open(os.path.join(HERE, out), 'w').write(s)
print(json.dumps({k: v for k, v in res.items() if k != 'coeffs'}), flush=True)
# Grok S57 (2026-10-17). RH OPEN. GRH OPEN.
