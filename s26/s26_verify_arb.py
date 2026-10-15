"""S26 [RIGOROUS, Arb]: S25 verifier with the claim-1 DC term generalised to an even continuous PL W (W=y0 on |t|<=2h;
What = -(2/xi^2) sum_j ds_j cos(xi z_j)); everything else unchanged.
S25 header follows.
S25 [RIGOROUS, Arb]: zeta odd-sector Weil coercivity, Claim 1 (m >= c, unchanged from S22/S24) + Claim 2 (beta) with
three upgrades over s24_verify_arb.py (Daniel Toupin / Golden Physics Project, with AI assistance from Grok; RH OPEN):
 (1) S(xi) = int_0^H v_w(t) sin(xi t) dt is evaluated by the exact knot-sum (one integration by parts, continuous PL on [h,H],
     a jump allowed at t = h):  S = S0 + y0 cos(xi h)/xi + xi^-2 [ s_n sin(xi H) - s_0 sin(xi h) - sum_j (s_j - s_{j-1}) sin(xi z_j) ].
 (2) Per cell [mid-r, mid+r]: S = p2 + delta, p2 the quadratic interpolant at mid-r, mid, mid+r, |delta| <= M3 r^3/(9 sqrt3),
     M3 = int_0^H t^3 |v_w| >= ||S'''||_inf.  Then  int_cell S^2 in [I2 - 2e sqrt(2r I2), I2 + 2e sqrt(2r I2) + 2r e^2],  I2 = int p2^2.
 (3) Parseval tail: int_X^inf S^2/m <= [ (pi/2)||v_w||^2_{L2(0,inf)} - (rigorous lower bound of int_0^X S^2) ] / inf_{xi>=X} m,
     inf_{xi>=X} m >= mlow_tail(X) (digamma monotone, Bernstein bound on the cosine block, |wDC| <= 4/(Dd xi^2)).
 beta = (4/pi) int_0^inf S^2/m ;  Q_zeta(g) >= (1 - 2 beta) c ||g||^2 for odd g in C_c((-h,h)) (S22 DC lift + pole Schur).
 Env: XBE (exact range X, integer), REL (m-slack refinement ratio), PREC (bits), OUT (json)."""
import json, sys, time, os
from flint import arb, acb, acb_poly, fmpq, ctx
ctx.prec = int(os.environ.get("PREC", "128"))
REL = os.environ.get("REL", "0.005")
XBETA_EXACT = int(os.environ.get("XBE", "3000"))
LOG2R0 = 7  # initial half-width r0 = 2^-7 (dyadic, so all cell points are exact)

def S0_pt(h, x):
    out = arb(0)
    for al, sg in ((arb(1)/2, 1), (-arb(1)/2, -1)):
        zz = acb(al, x); out += sg*(((zz*h).exp() - 1)/zz).imag
    return out/2

def make_S(h, ys, zs):
    n = len(ys) - 2   # ys[0..n] free values at zs[0..n]; ys[n+1] = 0 at zs[n+1] = H
    sl = [(ys[i+1]-ys[i])/(zs[i+1]-zs[i]) for i in range(n+1)]
    coef = [(zs[0], -sl[0])] + [(zs[j], -(sl[j]-sl[j-1])) for j in range(1, n+1)] + [(zs[n+1], sl[n])]
    y0 = ys[0]
    def S(x):
        if x == 0: return arb(0)
        acc = arb(0)
        for z, cf in coef: acc += cf*(x*z).sin()
        return S0_pt(h, x) + y0*(x*h).cos()/x + acc/(x*x)
    return S, sl

def verify(fn):
    d = json.load(open(fn)); t0 = time.time(); ext = d["extension"]
    hq = fmpq(*d["h_q"]); h = arb(hq); Hq = fmpq(*ext["H_q"]); H = arb(Hq); N = ext["N"]
    A_ = [arb(fmpq(n, d["den"])) for n in d["a_num"]]; d["K"] = len(A_)
    c = fmpq(*d["c"]); Dsp = arb(fmpq(*d["D_q"])); d["Xi_LP"] = d["Xi"]
    yW = [arb(fmpq(n, d["den"])) for n in d["y_num"]] + [arb(0)]; zW = [arb(fmpq(*q)) for q in d["z_q"]]; MW = len(zW)-1
    assert len(yW) == MW+1 and fmpq(*d["z_q"][0]) == 2*hq
    slW = [arb(0)] + [(yW[i+1]-yW[i])/(zW[i+1]-zW[i]) for i in range(MW)] + [arb(0)]
    dsW = [slW[j+1]-slW[j] for j in range(MW+1)]; sumdsW = sum((abs(x) for x in dsW), arb(0))
    W2 = 2*(abs(yW[0])*zW[0]**3/3 + sum((zW[i+1]**2*(abs(yW[i])+abs(yW[i+1]))*(zW[i+1]-zW[i])/2 for i in range(MW)), arb(0)))
    C = yW[0]
    ys = [arb(fmpq(n, ext["alpha_den"])) for n in ext["alpha_num"]] + [arb(0)]
    if "z_num" in ext:
        zs = [arb(fmpq(n, ext["z_den"])) for n in ext["z_num"]]
        assert fmpq(ext["z_num"][0], ext["z_den"]) == hq and fmpq(ext["z_num"][-1], ext["z_den"]) == Hq
        assert all(ext["z_num"][i] < ext["z_num"][i+1] for i in range(N+1))
    else:  # S23/S24 uniform knots
        zs = [arb(hq + (Hq-hq)*fmpq(k, N+1)) for k in range(N+2)]
    assert len(zs) == len(ys) == N+2
    zs[0] = h; zs[-1] = H
    Sfun, sl = make_S(h, ys, zs)
    # M3 >= int_0^H t^3 |v_w|,  ||v_w||^2 on (0, inf) exactly
    M3 = h**3*2*((h/2).cosh()-1) + sum(((zs[i+1]**3)*(abs(ys[i])+abs(ys[i+1]))*(zs[i+1]-zs[i])/2 for i in range(N+1)), arb(0))
    vnorm2 = ((h).sinh()-h)/2 + sum(((zs[i+1]-zs[i])*(ys[i]**2+ys[i]*ys[i+1]+ys[i+1]**2)/3 for i in range(N+1)), arb(0))
    print(f"[info] N={N} M3={float(M3.mid()):.3f} ||v_w||^2={float(vnorm2.mid()):.6f} jump={float(abs((h/2).sinh()-ys[0]).mid()):.4f}", flush=True)
    P = acb_poly([acb(x) for x in A_]); dP = P.derivative()
    tlg = arb(2).log(); kap = tlg/arb(2).sqrt(); a0 = arb(1)/4; lpi = arb.pi().log()
    L0 = 2*h; u = [L0 + Dsp*k for k in range(d["K"])]
    NN = 2**16; mx = arb(0)
    for j in range(NN//2+1):
        z = (acb(0, 2*arb.pi()*j/NN)).exp(); v = abs(P(z)).upper()
        if v > mx: mx = v
    ber = 1 - arb.pi()*(d["K"]-1)/NN; assert ber > 0; Mbern = (mx/ber).upper()
    Xi = arb(d["Xi_LP"]).floor() + 1
    def mlow_tail(X):
        return acb(a0, X/2).digamma().real - lpi - 2*kap - Mbern - 2*sumdsW/(X*X)
    tail = mlow_tail(Xi) - arb(c); assert tail > 0
    z3 = acb(3).zeta(acb(a0)).real
    F2 = (z3/2 + 2*kap*tlg*tlg + W2 + sum(abs(x)*uk*uk for x, uk in zip(A_, u))).upper()
    def fval(x):
        w = acb(a0, x/2); z = acb(0, Dsp*x).exp(); e = acb(0, L0*x).exp(); Pz = P(z)
        s = sum((dj*(x*zj).cos() for dj, zj in zip(dsW, zW)), arb(0)); spp = sum((dj*zj*(x*zj).sin() for dj, zj in zip(dsW, zW)), arb(0))
        Wv = -2*s/(x*x); Wp = 4*s/(x**3) + 2*spp/(x*x)
        fv = w.digamma().real - lpi - 2*kap*(x*tlg).cos() + Wv - (e*Pz).real - arb(c)
        Sp = (acb(0, L0)*e*Pz + e*dP(z)*acb(0, Dsp)*z).real
        fp = -acb(2).zeta(w).imag/2 + 2*kap*tlg*(x*tlg).sin() + Wp - Sp
        return fv, fp
    XBi = max(int(float(d["Xi_LP"])) + 2, XBETA_EXACT); XB = arb(XBi)
    nx = XBi*2**(LOG2R0-1); r0 = arb(1)/2**LOG2R0
    stack = [(r0*(2*i+1), r0) for i in range(nx)][::-1]
    cache = {}
    def Sc(x):
        k = int(round(float(x.mid())*2**34))
        if k not in cache: cache[k] = Sfun(x)
        return cache[k]
    ncell = nsplit = nref = 0; worst = None
    bint = arb(0); ssq_low = arb(0); rel = arb(REL); s3 = arb(3).sqrt()
    while stack:
        mid, r = stack.pop()
        fv, fp = fval(mid); slack = abs(fp)*r + F2*r*r/2; lb = fv - slack; ncell += 1
        if lb > 0 and slack > (fv+arb(c))*rel and r > arb("1e-5"):
            nref += 1; stack += [(mid+r/2, r/2), (mid-r/2, r/2)]; continue
        if not (lb > 0):
            nsplit += 1; assert r > arb("1e-9"), ("cannot certify", float(mid.mid()))
            stack += [(mid+r/2, r/2), (mid-r/2, r/2)]; continue
        if worst is None or lb.lower() < worst[0]: worst = (float(lb.lower()), float(mid.mid()))
        Sl, Sm, Sr = Sc(mid-r), Sc(mid), Sc(mid+r)
        Aq = (Sl+Sr-2*Sm)/(2*r*r); Bq = (Sr-Sl)/(2*r)
        I2 = 2*r*Sm*Sm + 4*Sm*Aq*r**3/3 + 2*Bq*Bq*r**3/3 + 2*Aq*Aq*r**5/5
        e = M3*r**3/(9*s3); Iu = I2.upper(); sq = (2*r*Iu).sqrt() if Iu > 0 else arb(0)
        bint += (I2 + 2*e*sq + 2*r*e*e)/(arb(c) + arb(lb.lower()))
        ssq_low += I2 - 2*e*sq
    T = arb.pi()/2*vnorm2 - arb(ssq_low.lower())
    T = arb(T.upper()) if T.upper() > 0 else arb(0)
    mXB = mlow_tail(XB); assert mXB > 0
    tailb = T/mXB
    beta = (4*(bint + tailb)/arb.pi()).upper()
    coer = ((1-2*arb(beta))*arb(c)).lower()
    out = dict(file=fn, h=f"{hq.p}/{hq.q}", H=f"{Hq.p}/{Hq.q}", N=N, c=f"{c.p}/{c.q}", C=float(C.mid()),
               M_upper=float(Mbern.mid()), Xi=float(Xi.mid()), tail_margin=float(tail.lower()), F2=float(F2.mid()),
               M3=float(M3.mid()), vnorm2=float(vnorm2.mid()), XB=XBi, mlow_XB=float(mXB.lower()),
               parseval_tail_S2=float(T.upper()), beta_tail_part=float((4*tailb/arb.pi()).upper()),
               cells=ncell, splits=nsplit, refinements=nref, REL=REL, prec=ctx.prec,
               worst_cell_lb=worst[0], worst_at=worst[1], beta_upper=float(beta), two_beta=float(2*beta),
               odd_coercivity_lower=float(coer), verified=bool(2*arb(beta) < 1 and coer > 0), seconds=round(time.time()-t0, 1))
    print(json.dumps(out), flush=True); return out

if __name__ == "__main__":
    res = [verify(fn) for fn in sys.argv[1:]]
    json.dump(res, open(os.environ.get("OUT", "s25_verify_results.json"), "w"), indent=1)
