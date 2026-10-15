"""S19b (chi7): independent Arb audit of Astra's cubic mod-7 first-window signed ground-state certificate.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH OPEN, GRH OPEN.  Fixed window h <= H=(log 3)/2, REAL g only.
Audits GPP-bridge/work/2026-10-15_astra_chi7_signed_groundstate_block.md and
GPPDiscovery3/experiments/haar_parity_shift/verify_chi7_signed_block.py (Fraction intervals).
Normalization: S17/Codex  Q_eta = E_eps + (log f - c_eps) N - 2 sum Lambda(n) Re eta(n) n^{-1/2} F(log n)   (real g).
Sections: [1] characters (flint dirichlet_char + explicit generator table)   [2] constants   [3] 2x2 certificate (+ best q)
          [4] independent identity checks (T by quadrature, c0 by digamma, Q(e) two ways)
          [5] S17-basis check: Q(g) >= LB(g) for many real g (tests the reduction itself, not just arithmetic)
          [6] complex-g caveat quantified."""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "s17")); sys.path.insert(0, os.path.join(HERE, "..", "s14"))
from flint import arb, acb, ctx, fmpq, fmpq_poly, arb_mat, dirichlet_char
from s17_build import Tarb, consts, build, vM
from s14_exact import p as basis_poly
ctx.prec = 256
ok_all = True
def check(name, cond):
    global ok_all
    ok_all &= bool(cond); print(f"  [{'PASS' if cond else 'FAIL'}] {name}")
print("== [1] characters mod 7 ==")
cubic = []
for l in range(1, 7):
    c = dirichlet_char(7, l)
    if c.is_primitive():
        print(f"  conrey {l}: order {c.order()} parity {c.parity()} chi(2)={c(2)} chi(3)={c(3)}")
        if c.order() == 3: cubic.append((l, c))
check("exactly two primitive cubic characters, labels 2 and 4", sorted(l for l, _ in cubic) == [2, 4])
for l, c in cubic:
    check(f"conrey {l}: conductor 7, even (parity 0)", c.conductor() == 7 and c.parity() == 0)
    check(f"conrey {l}: Re chi(2) contains -1/2, |chi(2)|=1, chi(2)^3=1", c(2).real.contains(arb(-1)/2) and (abs(c(2))-1).contains(0) and (c(2)**3-1).contains(0))
# explicit table: 3 generates (Z/7)^*, chi(3)=w (w=e^{2pi i/3}); 2=3^2 mod 7 so chi(2)=w^2, Re=-1/2; -1=3^3 so chi(-1)=w^3=1 (even)
check("3 generates (Z/7)^*: 3^k mod 7 =", sorted(pow(3, k, 7) for k in range(6)) == [1, 2, 3, 4, 5, 6])
check("2 = 3^2 mod 7, 6 = 3^3 mod 7 (=> chi(2)=w^{+-2}, chi(-1)=1)", pow(3, 2, 7) == 2 and pow(3, 3, 7) == 6)
H = arb(3).log()/2; L = 2*H; t = arb(2).log()
pp = [n for n in range(2, 40) if vM(n) and arb(n).log() < L]
check(f"prime powers n with log n < 2H: {pp} (n=3 at boundary log3=2H: F(log3)=0 on open support)", pp == [2])
check("H < log2 < 2H (strips disjoint, n=2 active)", (H < t) and (t < L))
print("== [2] constants (Arb, prec 256) ==")
z = arb(3)**(arb(-1)/4)
TH = Tarb(H, 0)                                         # int_H^inf nu_0, nu_0=2e^{-t/2}/(1-e^{-2t})
J = lambda s: (-s/2).exp()/(1-(-2*s).exp())               # J = nu_0/2
alpha = L*J(L); alpha_cf = 9*H/(4*arb(3).sqrt())
C0, _ = consts(); kappa = t/arb(2).sqrt(); b = kappa/2
d = TH + arb(7).log() - C0
r = (L-t)/H; beta = b*r; v = b*b*r - beta*beta
for nm, x in [("H", H), ("z=3^(-1/4)", z), ("T(H)", TH), ("2(atanh z+atan z)", 2*(z.atanh()+z.atan())), ("alpha=L J(2H)", alpha), ("9H/(4 sqrt3)", alpha_cf),
              ("c0=log(8pi)+gamma+pi/2", C0), ("log7-c0", arb(7).log()-C0), ("d=T(H)+log7-c0", d), ("kappa=log2/sqrt2", kappa), ("b=kappa/2", b),
              ("r=(2H-t)/H", r), ("beta=b r", beta), ("v=b^2 r-beta^2", v)]:
    print(f"  {nm:26s} {x.str(20)}")
check("alpha closed form equals L J(2H)", (alpha-alpha_cf).contains(0))
check("T(H) closed form equals 2(atanh z+atan z)", (TH-2*(z.atanh()+z.atan())).contains(0))
print("== [3] certificate ==")
q = arb(1)/200
A = d+beta-q; D = d+alpha-b-q; det = A*D-v
for nm, x in [("A=d+beta-q", A), ("D=d+alpha-b-q", D), ("det=A D - v", det)]: print(f"  {nm:26s} {x.str(20)}")
check("A>0, D>0, det>0  =>  M > (1/200) I", A > 0 and D > 0 and det > 0)
A0, D0 = d+beta, d+alpha-b
lam = (A0+D0)/2 - (((A0-D0)/2)**2+v).sqrt()
print(f"  best block constant lambda_min(M) = {lam.str(20)}   (Q >= lambda_min ||g||^2 on real g)")
check("lambda_min(M) > 1/200", lam > q)
# Astra's conservative d uses gamma_upper=H_1000-log1000
gU = sum((arb(1)/j for j in range(1, 1001)), arb(0)) - arb(1000).log()
dA = TH + arb(7).log() - ((8*arb.pi()).log()+gU+arb.pi()/2)
AA, DA = dA+beta-q, dA+alpha-b-q
print(f"  Astra-conservative d (gamma<=H_1000-log1000) = {dA.str(17)} ; A={AA.str(17)} D={DA.str(17)} det={(AA*DA-v).str(17)}")
astra = {"d": -0.13558822271871462, "alpha": 0.7135696131722594, "b": 0.2450645358671368, "beta": 0.18089205728676117,
         "v": 0.011608291671594762, "A": 0.04030383456804654, "D": 0.32791685458640807, "det": 0.0016080149877300007}
mine = {"d": dA, "alpha": alpha, "b": b, "beta": beta, "v": v, "A": AA, "D": DA, "det": AA*DA-v}
for k in astra:
    dev = abs(float(mine[k].mid()) - astra[k]); print(f"  vs Astra {k:6s}: |diff| = {dev:.2e}")
check("Astra's printed table reproduced to < 1e-12", all(abs(float(mine[k].mid())-astra[k]) < 1e-12 for k in astra))
print("== [4] independent identity checks ==")
nu0 = lambda s: 2*(-s/2).exp()/(1-(-2*s).exp())
Tq = acb.integral(lambda s, an: 2*(-s/2).exp()/(1-(-2*s).exp()), acb(H), acb(80)).real + arb(0, 4*arb(-40).exp().mid())
check(f"T(H) by quadrature on [H,80] (+tail ball) = {Tq.str(15)} overlaps closed form", Tq.overlaps(TH))
c0dig = arb.pi().log() - acb(arb(1)/4).digamma().real
check(f"c0 = log pi - psi(1/4) = {c0dig.str(15)} overlaps log(8pi)+gamma+pi/2", c0dig.overlaps(C0))
check("T convex (T''=-nu0'>0) and J decreasing: nu0'(s)<0 on grid s in (0,6]", all((nu0(arb(k)/50+arb(10)**-8) - nu0(arb(k)/50)) < 0 for k in range(1, 301)))
# Q(e) two ways: Astra split E0(e)=(1/L) int_0^L T(s) ds ; Codex E=int nu0 (1-F_e) with F_e(s)=(1-s/L)_+
xnu = lambda s: (s/2).exp()/((acb(0, 1)*s).sinc())       # s nu_0(s), smooth at 0
E_codex = acb.integral(lambda s, an: xnu(s), 0, acb(L)).real/L + Tarb(L, 0)
E_astra = acb.integral(lambda s, an: 2*((-s/2).exp().atanh()+(-s/2).exp().atan()), arb(10)**-30, acb(L)).real/L
Qe = E_codex + arb(7).log() - C0 + beta
print(f"  E0(e): Codex form {E_codex.str(15)} | Astra split {E_astra.str(15)} ; Q(e) = {Qe.str(15)}  (>= d+beta = {(d+beta).str(12)})")
check("E0(e) Codex == Astra split (log-singular endpoint of T cut at 1e-30, contributes < 1e-28)", abs(E_codex-E_astra) < arb(10)**-20)
print("== [5] S17 basis test of the reduction: Q(g) >= LB(g) for real g ==")
M = 16
B = build(M, H)
N = B["N"]; T0 = Tarb(L, 0)
Q = B["K0"] + N*(T0 + arb(7).log() - C0) - B["T"][2]*(2*t*(arb(-1)/2)/arb(2).sqrt())
mom = []
for j in range(M):
    P = basis_poly(j).integral(); mom.append(H*(arb(P(fmpq(1)).p)/P(fmpq(1)).q - arb(P(fmpq(-1)).p)/P(fmpq(-1)).q))
sqv = v.sqrt(); sL = L.sqrt()
def evalg(c):
    qq = sum((c[i]*c[k]*Q[i, k] for i in range(M) for k in range(M)), arb(0))
    nn = sum((c[i]*c[k]*N[i, k] for i in range(M) for k in range(M)), arb(0))
    a = sum((c[i]*mom[i] for i in range(M)), arb(0))/sL
    uu = (nn - a*a); uu = uu if uu > 0 else arb(0)
    lb = (d+beta)*a*a - 2*sqv*abs(a)*uu.sqrt() + (d+alpha-b)*uu
    return qq, nn, lb, a
random.seed(19)
worst = None; cnt = 0; fails = 0
vecs = []
for trial in range(400):
    c = [arb(random.gauss(0, 1)/(1+i)**random.choice([0, 1, 2])) for i in range(M)]
    if trial % 2 == 0:
        for i in range(1, M, 2): c[i] = arb(0)      # even g emphasised (where e lives)
    vecs.append(c)
# least-squares approx of the constant (even polys) -> near-ground-state test
G = arb_mat(M, M); rhs = arb_mat(M, 1)
for i in range(M):
    rhs[i, 0] = mom[i]
    for k in range(M): G[i, k] = N[i, k] + (arb(10)**-30 if i == k else 0)
sol = G.solve(rhs)
# Ritz minimisers of Q (even / odd sectors) via S18 helper: the most stringent test vectors
sys.path.insert(0, os.path.join(HERE, "..", "s18"))
from s18_ritz_common import geig, SECT
import mpmath as mp
for sec, Iset in SECT(M):
    ev, V, Li = geig(Q, N, Iset); k = min(range(len(ev)), key=lambda i: ev[i]); w = Li.T*V[:, k]
    c = [arb(0)]*M
    for a_, i in enumerate(Iset): c[i] = arb(fmpq(int(mp.nint(w[a_]*mp.mpf(10)**30)), 10**30))
    vecs.append(c); print(f"  Ritz lambda_min ({sec} sector, M={M}) = {mp.nstr(ev[k], 12)}")
vecs.append([sol[i, 0] for i in range(M)])
for c in vecs:
    qq, nn, lb, a = evalg(c); cnt += 1
    gap = (qq - lb)/nn
    if not gap > 0: fails += 1
    if worst is None or gap.mid() < worst[0].mid(): worst = (gap, qq/nn, lb/nn, a*a/nn)
qq, nn, lb, a = evalg(vecs[-1])
print(f"  LSQ-constant g: Q/N={(qq/nn).str(10)} LB/N={(lb/nn).str(10)} a^2/N={(a*a/nn).str(10)}  (Q(e)={Qe.str(10)})")
print(f"  {cnt} real test vectors (M={M}); min (Q-LB)/N = {worst[0].str(10)} at Q/N={worst[1].str(8)}, LB/N={worst[2].str(8)}, a^2/N={worst[3].str(6)}")
check("Q(g) - LB(g) > 0 for every tested real g (certified Arb, each vector)", fails == 0)
print("== [6] complex-g caveat (single nonreal sector) ==")
# complex g: prime term -kappa(chi(2) S + conj chi(2) S*) = b(zeta S + conj zeta S*), zeta=-2chi(2), |zeta|=2 -> ||B||=2b, beta same, ||Be||^2=4b^2 r
vc = 4*b*b*r - beta*beta; Dc = d+alpha-2*b-q; Ac = A
print(f"  complex g: A={Ac.str(12)} D={Dc.str(12)} det={(Ac*Dc-vc).str(12)}  -> same block argument does {'NOT ' if not (Ac*Dc-vc) > 0 else ''}close")
print(f"\n{'PASS' if ok_all else 'FAIL'}: S19b chi7 audit. Real g, fixed window h<=log(3)/2, both primitive cubic characters mod 7. RH OPEN, GRH OPEN.")
