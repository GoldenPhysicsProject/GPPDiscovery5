"""S22 [RIGOROUS, Arb 128-bit, python-flint]: odd-sector zeta Weil coercivity on (-h,h), h <= log3/2.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH OPEN.
Certificate data (rational): C >= 0, a_k, c > 0, with u_k = 2h + k D (rational h), Dd > 0.
 m(xi) = sigma_A(xi) + C wDC(xi) - sum_k a_k cos(u_k xi),
 sigma_A = Re psi(1/4 + i xi/2) - log pi - 2 kappa cos(xi log2), kappa = log2/sqrt2,
 wDC(xi) = 2 sin(A xi)/xi * sin(b xi)/(b xi), A = 2h + Dd/2, b = Dd/2  (FT of trapezoid W, W = 1 on |t| <= 2h).
Claim 1 (checked here): m(xi) >= c for all real xi.
Claim 2 (checked here): beta := (1/2pi) int |vhat|^2/m < 1/2 (cells with exact m up to XBETA, analytic tail beyond), v = sinh(x/2) 1_(-h,h);  computed as an UPPER bound.
Theorem (proof in the bridge note): for odd g in C_c((-h,h)) (real or complex), Q_zeta(g) >= (1 - 2 beta) c ||g||^2.
Proof of claim 1: as S20 (cells with exact f, f' and global |f''| bound; digamma-monotone tail; Bernstein sup|P|), plus
 |wDC^(k)| <= int |t|^k W(t) dt <= 2 (2h+Dd)^(k+1)/(k+1) and, in the tail, |wDC(xi)| <= 4/(Dd xi^2).
Claim 2: per cell, int_cell 4 S^2 / m <= |cell| * sup_cell(4 S^2) / (c + lb_cell), S(xi) = int_0^h sinh(x/2) sin(xi x) dx
 (ball evaluation on the cell); tail chunks [X, 2X]: |S| <= (cosh(h/2)+1)/xi and m >= monotone lower bound at X."""
import json, sys, time
from flint import arb, acb, acb_poly, fmpq, ctx
ctx.prec = 128
REL = '0.01'
XBETA = 3000   # beta-only refinement: split a passing cell while its Taylor slack exceeds REL*m(mid)
def verify(fn, cfrac=None):
    d = json.load(open(fn)); t0 = time.time()
    hq = fmpq(*d["h_q"]); h = arb(hq); K = d["K"]
    aq = [fmpq(n, d["den"]) for n in d["a_num"]]; A_ = [arb(x) for x in aq]; C = arb(fmpq(d["C_num"], d["den"]))
    c = fmpq(*d["c"]); D = arb(fmpq(*d["D_q"])); Dd = arb(fmpq(*d["Dd_q"]))
    assert C >= 0
    P = acb_poly([acb(x) for x in A_]); dP = P.derivative()
    t = arb(2).log(); kap = t/arb(2).sqrt(); a0 = arb(1)/4; lpi = arb.pi().log()
    L0 = 2*h; u = [L0 + D*k for k in range(K)]
    Aw = 2*h + Dd/2; b = Dd/2
    N = 2**16; mx = arb(0)
    for j in range(N//2+1):
        z = (acb(0, 2*arb.pi()*j/N)).exp(); v = abs(P(z)).upper()
        if v > mx: mx = v
    ber = 1 - arb.pi()*(K-1)/N; assert ber > 0; M = (mx/ber).upper()
    Xi = arb(d["Xi_LP"]).floor() + 1
    def mlow_tail(X):  # lower bound for m on [X, inf)
        return acb(a0, X/2).digamma().real - lpi - 2*kap - M - C*4/(Dd*X*X)
    tail = mlow_tail(Xi) - arb(c); assert tail > 0, ("tail fails", tail)
    z3 = acb(3).zeta(acb(a0)).real
    F2 = (z3/2 + 2*kap*t*t + C*2*(2*h+Dd)**3/3 + sum(abs(x)*uk*uk for x, uk in zip(A_, u))).upper()
    def fval(x):
        x = arb(x); w = acb(a0, x/2); z = acb(0, D*x).exp(); e = acb(0, L0*x).exp(); Pz = P(z)
        s1 = (Aw*x).sin()/x; s2 = (b*x).sin()/(b*x)
        s1p = (Aw*x*(Aw*x).cos() - (Aw*x).sin())/(x*x); s2p = (b*x*(b*x).cos() - (b*x).sin())/(b*x*x)
        fv = w.digamma().real - lpi - 2*kap*(x*t).cos() + C*2*s1*s2 - (e*Pz).real - arb(c)
        Sp = (acb(0, L0)*e*Pz + e*dP(z)*acb(0, D)*z).real
        fp = -acb(2).zeta(w).imag/2 + 2*kap*t*(x*t).sin() + C*2*(s1p*s2 + s1*s2p) - Sp
        return fv, fp
    def S4(mid, r):  # sup over cell of 4 S^2 via ball evaluation
        x = arb(mid, r); out = arb(0)
        for al, sg in ((arb(1)/2, 1), (-arb(1)/2, -1)):
            zz = acb(al, x); out += sg*(((zz*h).exp() - 1)/zz).imag
        S = out/2; return (4*S*S).upper()
    XB = max(Xi, arb(XBETA))   # cells (with exact m) are run up to XB >= Xi; beyond XB the analytic tail bound is used
    nx = int(float(XB.mid())/0.01) + 1; r0 = XB/(2*nx)
    stack = [(XB*(2*i+1)/(2*nx), r0) for i in range(nx)]
    ncell = 0; worst = None; nsplit = 0; beta = arb(0); nrefine = [0]
    while stack:
        mid, r = stack.pop(); fv, fp = fval(mid); lb = fv - abs(fp)*r - F2*r*r/2; ncell += 1
        slack = abs(fp)*r + F2*r*r/2
        if lb > 0 and slack > (fv + arb(c))*arb(REL) and r > arb("1e-6"):
            nrefine[0] += 1; stack += [(mid - r/2, r/2), (mid + r/2, r/2)]; continue
        if lb > 0:
            if worst is None or lb.lower() < worst[0]: worst = (float(lb.lower()), float(mid.mid()), float(r.mid()))
            beta += 2*r*S4(mid, r)/(arb(c) + arb(lb.lower()))
        else:
            nsplit += 1; assert r > arb("1e-9"), ("cannot certify near", mid)
            stack += [(mid - r/2, r/2), (mid + r/2, r/2)]
    # |S(xi)| <= s1/xi + s2/xi^2 (two integrations by parts), s1 = sinh(h/2), s2 = (cosh(h/2) + h sinh(h/2)/2)/2
    s1 = (h/2).sinh(); s2 = ((h/2).cosh() + h*(h/2).sinh()/2)/2
    X = XB
    while X < arb(10)**7:
        beta += 4*(s1 + s2/X)**2*(1/X - 1/(2*X))/mlow_tail(X); X = 2*X
    beta += 4*(s1 + s2/X)**2/(X*mlow_tail(X))
    beta = (beta/arb.pi()).upper()   # (1/2pi) int_R = (1/pi) int_0^inf
    coer = ((1 - 2*arb(beta))*arb(c)).lower()
    out = dict(file=fn, h=f"{hq.p}/{hq.q}", c=f"{c.p}/{c.q}", C=float(C.mid()), M_upper=float(M.mid()), Xi=float(Xi.mid()),
               tail_margin=float(tail.lower()), F2=float(F2.mid()), cells=ncell, splits=nsplit, beta_refinements=nrefine[0], REL=REL, worst_cell_lb=worst[0], worst_at=worst[1],
               beta_upper=float(beta), two_beta=float(2*beta), odd_coercivity_lower=float(coer), verified=bool(2*arb(beta) < 1),
               seconds=round(time.time()-t0, 1))
    print(json.dumps(out), flush=True); return out
if __name__ == "__main__":
    res = [verify(fn) for fn in sys.argv[1:]]
    json.dump(res, open("s22_verify_results.json", "w"), indent=1)
