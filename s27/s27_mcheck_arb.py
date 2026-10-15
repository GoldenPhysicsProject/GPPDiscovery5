"""S27 (copy of S26 with exact h = log3/2 support via s27_common). S26 Arb verifier: claim-1 multiplier with PL exterior DC term (h=1/2 target).
m(xi) = Re psi(1/4+i xi/2) - log pi - 2 kap cos(xi log2) + What(xi) - sum_k a_k cos(u_k xi),
W even, continuous PL: W = y0 on |t|<=2h, PL through rational knots z_0=2h<...<z_M, W(z_M)=0.
What(xi) = 2 int_0^zM W cos = -(2/xi^2) sum_j ds_j cos(xi z_j)  (ds_j = slope jump at z_j); What(0)=2 int W.
For supported g in (-h,h): form of What = y0 |<1,g>|^2 (A = B + y0 11*); u_k >= 2h invisible.
Checks m >= c on [0,inf): adaptive Taylor cells on [0,Xi], Bernstein + digamma tail beyond. RH OPEN."""
import json, sys, time, os
from flint import arb, acb, acb_poly, fmpq, ctx
from s27_common import load_h, load_W
ctx.prec = 128
def run(fn):
    d = json.load(open(fn)); t0 = time.time()
    h = load_h(d); D = arb(fmpq(*d["D_q"])); den = d["den"]; c = arb(fmpq(*d["c"]))
    A_ = [arb(fmpq(n, den)) for n in d["a_num"]]; K = len(A_)
    ys, zs = load_W(d, h)
    M = len(zs)-1
    assert len(ys) == M+1
    sl = [arb(0)] + [(ys[i+1]-ys[i])/(zs[i+1]-zs[i]) for i in range(M)] + [arb(0)]
    ds = [sl[j+1]-sl[j] for j in range(M+1)]          # jump at z_j
    W0 = 2*(ys[0]*zs[0] + sum(((ys[i]+ys[i+1])*(zs[i+1]-zs[i])/2 for i in range(M)), arb(0)))
    sumds = sum((abs(x) for x in ds), arb(0))
    # |What''| <= 2 int t^2 |W|
    W2 = 2*(abs(ys[0])*zs[0]**3/3 + sum((zs[i+1]**2*(abs(ys[i])+abs(ys[i+1]))*(zs[i+1]-zs[i])/2 for i in range(M)), arb(0)))
    L0 = 2*h; u = [L0 + D*k for k in range(K)]
    P = acb_poly([acb(x) for x in A_]); dP = P.derivative()
    tlg = arb(2).log(); kap = tlg/arb(2).sqrt(); a0 = arb(1)/4; lpi = arb.pi().log()
    NN = 2**16; mx = arb(0)
    for j in range(NN//2+1):
        z = (acb(0, 2*arb.pi()*j/NN)).exp(); v = abs(P(z)).upper()
        if v > mx: mx = v
    ber = 1 - arb.pi()*(K-1)/NN; assert ber > 0; Mbern = (mx/ber).upper()
    Xi = arb(d["Xi"]).floor() + 1
    def mlow_tail(X): return acb(a0, X/2).digamma().real - lpi - 2*kap - Mbern - 2*sumds/(X*X)
    tail = mlow_tail(Xi) - c; assert tail > 0, ("tail fails", float(tail.mid()))
    z3 = acb(3).zeta(acb(a0)).real
    F2 = (z3/2 + 2*kap*tlg*tlg + W2 + sum((abs(x)*uk*uk for x, uk in zip(A_, u)), arb(0))).upper()
    def What(x):
        if x == 0: return W0, arb(0)
        s = sum((dj*(x*zj).cos() for dj, zj in zip(ds, zs)), arb(0)); sp_ = sum((dj*zj*(x*zj).sin() for dj, zj in zip(ds, zs)), arb(0))
        return -2*s/(x*x), 4*s/(x**3) + 2*sp_/(x*x)
    def fval(x):
        w = acb(a0, x/2); z = acb(0, D*x).exp(); e = acb(0, L0*x).exp(); Pz = P(z)
        Wv, Wp = What(x)
        fv = w.digamma().real - lpi - 2*kap*(x*tlg).cos() + Wv - (e*Pz).real - c
        Sp = (acb(0, L0)*e*Pz + e*dP(z)*acb(0, D)*z).real
        fp = -acb(2).zeta(w).imag/2 + 2*kap*tlg*(x*tlg).sin() + Wp - Sp
        return fv, fp
    # What' near 0 via formula is ill-conditioned; bound |What'| <= 2 int t|W| on first cells instead
    W1 = 2*(abs(ys[0])*zs[0]**2/2 + sum((zs[i+1]*(abs(ys[i])+abs(ys[i+1]))*(zs[i+1]-zs[i])/2 for i in range(M)), arb(0)))
    r0 = arb(1)/2**7; nx = int(float(Xi.mid()))*2**6
    stack = [(r0*(2*i+1), r0) for i in range(nx)][::-1]
    ncell = nsplit = 0; worst = None
    while stack:
        mid, r = stack.pop(); ncell += 1
        if float(mid.mid()) < 1.0:   # direct derivative bound near 0 (avoid cancellation)
            w = acb(a0, mid/2); z = acb(0, D*mid).exp(); e = acb(0, L0*mid).exp(); Pz = P(z)
            fv = w.digamma().real - lpi - 2*kap*(mid*tlg).cos() + 2*(ys[0]*((mid*zs[0]).sin()/mid) ) * 0 - (e*Pz).real - c
            # exact What(mid) by segment sum (no 1/x^2 cancellation issue at x>=r0/2)
            Wv, _ = What(mid); fv = fv + Wv
            Sp = (acb(0, L0)*e*Pz + e*dP(z)*acb(0, D)*z).real
            fp = -acb(2).zeta(w).imag/2 + 2*kap*tlg*(mid*tlg).sin() - Sp
            slack = (abs(fp) + W1)*r + F2*r*r/2
        else:
            fv, fp = fval(mid); slack = abs(fp)*r + F2*r*r/2
        lb = fv - slack
        if not (lb > 0):
            nsplit += 1; assert r > arb("1e-10"), ("cannot certify", float(mid.mid()), float(fv.mid()))
            stack += [(mid+r/2, r/2), (mid-r/2, r/2)]; continue
        if worst is None or lb.lower() < worst[0]: worst = (float(lb.lower()), float(mid.mid()), float(fv.mid()))
    out = dict(file=fn, h=d.get("h_spec") or str(fmpq(*d["h_q"])), c=str(fmpq(*d["c"])), y0=float(ys[0].mid()), K=K, M=M, Mbern=float(Mbern.mid()),
               Xi=float(Xi.mid()), tail_margin=float(tail.lower()), F2=float(F2.mid()), cells=ncell, splits=nsplit,
               worst_cell_lb=worst[0], worst_at=worst[1], worst_fv=worst[2], verified=True, seconds=round(time.time()-t0, 1))
    print(json.dumps(out), flush=True); return out
if __name__ == "__main__":
    res = [run(f) for f in sys.argv[1:]]
    json.dump(res, open(os.environ.get("OUT", "s27_mcheck_results.json"), "w"), indent=1)
