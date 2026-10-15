"""S20 [C]: RIGOROUS (Arb ball arithmetic, python-flint) verification of the dual cosine certificates from s20_lp.py.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH OPEN, GRH OPEN.

Certificate (rational a_k, c):  f(xi) = sigma_chi(xi) - c - sum_k a_k cos(u_k xi) >= 0 for all real xi,
   u_k = log3 + k D,  sigma_chi(xi) = Re psi((1+2eps)/4 + i xi/2) + log(f/pi) - 2 kappa chi(2) cos(xi log2), kappa = log2/sqrt2.
Then for g in C_c((-h,h)), h <= log3/2 (real OR complex g): int |ghat|^2 cos(u_k xi) dxi = 2pi Re (g*g~)(u_k) = 0 (u_k >= 2h),
   so Q_chi(g) = (1/2pi) int |ghat|^2 sigma_chi >= c ||g||^2.
Proof of f>=0:
 (T) xi >= Xi: Re psi(a+iy) is increasing in |y| (series psi(a) + sum y^2/((n+a)((n+a)^2+y^2))), the chi(2) term is >= -2kappa,
     and |sum a_k cos(u_k xi)| = |Re(e^{iL xi} P(e^{iD xi}))| <= M := sup_|z|=1 |P|; need Re psi(a+iXi/2)+log(f/pi)-2kappa-c-M > 0.
     M is bounded from an N-point grid by Bernstein: M <= m/(1 - pi (K-1)/N)  (real coeffs => half circle suffices).
 (B) 0 <= xi <= Xi: cells [m-r, m+r]: f(xi) >= f(m) - |f'(m)| r - F2 r^2/2, with the GLOBAL bound
     |f''| <= F2 := zeta(3,a)/2 + 2 kappa log(2)^2 + sum |a_k| u_k^2   (|psi''(z)| <= 2 zeta(3, Re z)).  Failing cells are bisected.
f is even, so xi >= 0 suffices."""
import json, sys, time
from flint import arb, acb, acb_poly, fmpq, ctx
ctx.prec = 128
def verify(fn):
    d = json.load(open(fn)); t0 = time.time()
    f, eps, chi2, K = d["f"], d["eps"], d["chi2"], d["K"]
    aq = [fmpq(n, d["a_den"]) for n in d["a_num"]]; c = fmpq(*d["c"]); D = arb(fmpq(*d["D"]))
    A = [arb(x) for x in aq]; P = acb_poly([acb(x) for x in A]); dP = P.derivative()
    L = arb(3).log(); t = arb(2).log(); kap = t/arb(2).sqrt(); a = arb(1+2*eps)/4; lf = (arb(f)/arb.pi()).log()
    u = [L + D*k for k in range(K)]
    # sup |P| on the unit circle
    N = 2**16; m = arb(0)
    for j in range(N//2+1):
        z = (acb(0, 2*arb.pi()*j/N)).exp(); v = abs(P(z)).upper()
        if v > m: m = v
    ber = 1 - arb.pi()*(K-1)/N; assert ber > 0
    M = (m/ber).upper()
    # tail
    Xi = arb(d["Xi_LP"]).floor() + 1
    tail = acb(a, Xi/2).digamma().real + lf - 2*kap - arb(c) - M
    assert tail > 0, ("tail fails", tail)
    # second-derivative bound
    z3 = acb(3).zeta(acb(a)).real
    F2 = (z3/2 + 2*kap*t*t + sum(abs(x)*uk*uk for x, uk in zip(A, u))).upper()
    def fval(x):
        x = arb(x); w = acb(a, x/2); z = acb(0, D*x).exp(); e = acb(0, L*x).exp()
        Pz = P(z); S = (e*Pz).real
        fv = w.digamma().real + lf - 2*kap*chi2*(x*t).cos() - arb(c) - S
        dpsi = acb(2).zeta(w)  # psi'(w)
        Sp = (acb(0, L)*e*Pz + e*dP(z)*acb(0, D)*z).real
        fp = -dpsi.imag/2 + 2*kap*chi2*t*(x*t).sin() - Sp
        return fv, fp
    # cells
    nx = int(float(Xi.mid())/0.01) + 1; r0 = Xi/(2*nx)
    stack = [(Xi*(2*i+1)/(2*nx), r0) for i in range(nx)]
    ncell = 0; worst = None; nsplit = 0
    while stack:
        mid, r = stack.pop(); fv, fp = fval(mid); lb = fv - abs(fp)*r - F2*r*r/2; ncell += 1
        if lb > 0:
            if worst is None or lb.lower() < worst[0]: worst = (float(lb.lower()), float(mid.mid()), float(r.mid()))
        else:
            nsplit += 1; assert r > arb("1e-9"), ("cannot certify near", mid)
            stack += [(mid - r/2, r/2), (mid + r/2, r/2)]
    out = dict(name=d["name"], c=f"{c.p}/{c.q}", c_float=float(arb(c).mid()), M_upper=float(M.mid()), grid_max=float(m.mid()),
               Xi=float(Xi.mid()), tail_margin=float(tail.lower()), F2=float(F2.mid()), cells=ncell, splits=nsplit,
               worst_cell_lower_bound=worst[0], worst_at=worst[1], worst_halfwidth=worst[2], verified=True, seconds=round(time.time()-t0, 1))
    print(json.dumps(out), flush=True); return out
if __name__ == "__main__":
    res = [verify(fn) for fn in sys.argv[1:]]
    json.dump(res, open("s20_verify_results.json", "w"), indent=1)
