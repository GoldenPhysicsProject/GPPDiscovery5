"""S20 [F]: RIGOROUS Arb verification of the complex-character (cos+sin) dual certificate from s20_lp_complex.py.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH OPEN, GRH OPEN.
Same proof as s20_verify_arb.py, but sigma is not even: cells cover [-Xi, Xi], the sup|P| grid covers the full circle,
chi(2) = -1/2 + i sqrt(3)/2 exactly, S(xi) = Re(e^{iL xi} P(e^{iD xi})) with P_k = a_k - i b_k,
F2 = zeta(3,a)/2 + 2 kappa log(2)^2 + sum (|a_k|+|b_k|) u_k^2."""
import json, sys, time
from flint import arb, acb, acb_poly, fmpq, ctx
ctx.prec = 128
d = json.load(open(sys.argv[1])); t0 = time.time()
f, eps, K = d["f"], d["eps"], d["K"]
A = [arb(fmpq(n, d["a_den"])) for n in d["a_num"]]; Bc = [arb(fmpq(n, d["a_den"])) for n in d["b_num"]]
c = fmpq(*d["c"]); D = arb(fmpq(*d["D"]))
P = acb_poly([acb(x, -y) for x, y in zip(A, Bc)]); dP = P.derivative()
L = arb(3).log(); t = arb(2).log(); kap = t/arb(2).sqrt(); a = arb(1+2*eps)/4; lf = (arb(f)/arb.pi()).log()
pq = d.get("chi2_angle", [1, 3]); th2 = 2*arb.pi()*pq[0]/pq[1]; cr, ci = th2.cos(), th2.sin()   # chi(2) = exp(2 pi i p/q) exactly
u = [L + D*k for k in range(K)]
N = 2**16; m = arb(0)
for j in range(N):
    v = abs(P((acb(0, 2*arb.pi()*j/N)).exp())).upper()
    if v > m: m = v
ber = 1 - arb.pi()*(K-1)/N; M = (m/ber).upper()
Xi = arb(d["Xi_LP"]).floor() + 1
def tailm(X): return acb(a, X/2).digamma().real + lf - 2*kap - arb(c) - M
while not tailm(Xi) > 0: Xi = Xi + 1   # Re psi increasing in |y|: enlarge the body range until the tail bound closes
tail = tailm(Xi)
F2 = (acb(3).zeta(acb(a)).real/2 + 2*kap*t*t + sum((abs(x)+abs(y))*uk*uk for x, y, uk in zip(A, Bc, u))).upper()
def fval(x):
    w = acb(a, x/2); z = acb(0, D*x).exp(); e = acb(0, L*x).exp(); Pz = P(z)
    fv = w.digamma().real + lf - 2*kap*(cr*(x*t).cos() - ci*(x*t).sin()) - arb(c) - (e*Pz).real
    fp = -acb(2).zeta(w).imag/2 + 2*kap*t*(cr*(x*t).sin() + ci*(x*t).cos()) - (acb(0, L)*e*Pz + e*dP(z)*acb(0, D)*z).real
    return fv, fp
nx = 2*(int(float(Xi.mid())/0.01) + 1); r0 = Xi/nx
stack = [(-Xi + r0*(2*i+1), r0) for i in range(nx)]; ncell = nsplit = 0; worst = None
while stack:
    mid, r = stack.pop(); fv, fp = fval(mid); lb = fv - abs(fp)*r - F2*r*r/2; ncell += 1
    if lb > 0:
        if worst is None or lb.lower() < worst[0]: worst = (float(lb.lower()), float(mid.mid()))
    else:
        nsplit += 1; assert r > arb("1e-9"), ("cannot certify near", mid); stack += [(mid - r/2, r/2), (mid + r/2, r/2)]
out = dict(name=d["name"], c=f"{c.p}/{c.q}", c_float=float(arb(c).mid()), M_upper=float(M.mid()), Xi=float(Xi.mid()),
           tail_margin=float(tail.lower()), F2=float(F2.mid()), cells=ncell, splits=nsplit, worst_cell_lower_bound=worst[0],
           worst_at=worst[1], verified=True, seconds=round(time.time()-t0, 1))
print(json.dumps(out)); json.dump(out, open(f"s20_verify_{d['name']}.json", "w"), indent=1)
