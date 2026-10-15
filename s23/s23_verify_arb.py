"""S23 [RIGOROUS, Arb 128-bit]: zeta odd-sector Weil coercivity at h=9/20 with PL extension.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH OPEN.
Claim 1: m >= c (S22). Claim 2: beta for piecewise-linear odd extension of sinh(x/2).
On each cell, |S| <= max(|S(mid-r)|,|S(mid)|,|S(mid+r)|) + M2*r^2/8,
where M2 = int_0^H t^2 |v_w(t)| dt >= ||S''||_infty.
"""
import json, sys, time
from flint import arb, acb, acb_poly, fmpq, ctx
ctx.prec = 128
REL = '0.02'
XBETA_EXACT = 600

def int_lin_sin_pt(a, b, lo, hi, x):
    Fhi = -(a*hi+b)*(x*hi).cos()/x + a*(x*hi).sin()/(x*x)
    Flo = -(a*lo+b)*(x*lo).cos()/x + a*(x*lo).sin()/(x*x)
    return Fhi - Flo

def S0_pt(h, x):
    out = arb(0)
    for al, sg in ((arb(1)/2, 1), (-arb(1)/2, -1)):
        zz = acb(al, x)
        out += sg * (((zz*h).exp() - 1)/zz).imag
    return out/2

def Sw_pt(h, ys, zs, x):
    if x <= 0:
        return arb(0)
    S = S0_pt(h, x)
    for i in range(len(zs)-1):
        lo, hi = zs[i], zs[i+1]
        y0, y1 = ys[i], ys[i+1]
        D = hi - lo
        a = (y1 - y0)/D
        b = (y0*hi - y1*lo)/D
        S += int_lin_sin_pt(a, b, lo, hi, x)
    return S

def M2_bound(h, ys, zs):
    """int_0^H t^2 |v_w(t)| dt."""
    a = arb(1)/2
    # int_0^h t^2 sinh(a t) dt; integrate by parts twice
    # Known: use antiderivative via symbolic pattern
    # int t^2 sinh(at) = (t^2/a) cosh(at) - 2t/a^2 sinh(at) + 2/a^3 cosh(at)
    def Fsinh2(t):
        return (t*t/a)*(a*t).cosh() - (2*t/(a*a))*(a*t).sinh() + (2/(a*a*a))*(a*t).cosh()
    I0 = Fsinh2(h) - Fsinh2(arb(0))
    Itail = arb(0)
    for i in range(len(zs)-1):
        lo, hi = zs[i], zs[i+1]
        y0, y1 = ys[i], ys[i+1]
        D = hi - lo
        L1 = (abs(y0)+abs(y1))*D/2  # int |v|
        Itail += (hi*hi) * L1        # t^2 <= hi^2
    return I0 + Itail

def verify(fn):
    d = json.load(open(fn)); t0 = time.time()
    ext = d["extension"]
    hq = fmpq(*d["h_q"]); h = arb(hq)
    Hq = fmpq(*ext["H_q"]); H = arb(Hq)
    N = ext["N"]
    aq = [fmpq(n, d["den"]) for n in d["a_num"]]; A_ = [arb(x) for x in aq]
    C = arb(fmpq(d["C_num"], d["den"]))
    c = fmpq(*d["c"]); Dsp = arb(fmpq(*d["D_q"])); Dd = arb(fmpq(*d["Dd_q"]))
    aden = ext["alpha_den"]
    ys = [arb(fmpq(n, aden)) for n in ext["alpha_num"]] + [arb(0)]
    zs = [h + (H-h)*k/(N+1) for k in range(N+2)]
    M2 = M2_bound(h, ys, zs)
    print(f"[info] M2={float(M2.mid()):.4f} N={N}", flush=True)
    P = acb_poly([acb(x) for x in A_]); dP = P.derivative()
    tlg = arb(2).log(); kap = tlg/arb(2).sqrt(); a0 = arb(1)/4; lpi = arb.pi().log()
    L0 = 2*h; u = [L0 + Dsp*k for k in range(d["K"])]
    Aw = 2*h + Dd/2; bdc = Dd/2
    NN = 2**16; mx = arb(0)
    for j in range(NN//2+1):
        z = (acb(0, 2*arb.pi()*j/NN)).exp(); v = abs(P(z)).upper()
        if v > mx: mx = v
    ber = 1 - arb.pi()*(d["K"]-1)/NN; assert ber > 0; Mbern = (mx/ber).upper()
    Xi = arb(d["Xi_LP"]).floor() + 1
    def mlow_tail(X):
        return acb(a0, X/2).digamma().real - lpi - 2*kap - Mbern - C*4/(Dd*X*X)
    tail = mlow_tail(Xi) - arb(c); assert tail > 0
    z3 = acb(3).zeta(acb(a0)).real
    F2 = (z3/2 + 2*kap*tlg*tlg + C*2*(2*h+Dd)**3/3 + sum(abs(x)*uk*uk for x,uk in zip(A_, u))).upper()
    def fval(x):
        x = arb(x); w = acb(a0, x/2); z = acb(0, Dsp*x).exp(); e = acb(0, L0*x).exp(); Pz = P(z)
        s1 = (Aw*x).sin()/x; s2 = (bdc*x).sin()/(bdc*x)
        s1p = (Aw*x*(Aw*x).cos() - (Aw*x).sin())/(x*x)
        s2p = (bdc*x*(bdc*x).cos() - (bdc*x).sin())/(bdc*x*x)
        fv = w.digamma().real - lpi - 2*kap*(x*tlg).cos() + C*2*s1*s2 - (e*Pz).real - arb(c)
        Sp = (acb(0, L0)*e*Pz + e*dP(z)*acb(0, Dsp)*z).real
        fp = -acb(2).zeta(w).imag/2 + 2*kap*tlg*(x*tlg).sin() + C*2*(s1p*s2 + s1*s2p) - Sp
        return fv, fp
    # Claim 1+2 cells up to XB
    XB = max(Xi, arb(XBETA_EXACT))
    # use finer initial grid so r is moderate
    nx = int(float(XB.mid())/0.02) + 1; r0 = XB/(2*nx)
    stack = [(XB*(2*i+1)/(2*nx), r0) for i in range(nx)]
    ncell=0; worst=None; nsplit=0; beta=arb(0); nrefine=[0]
    s1 = (h/2).sinh()
    s2c = ((h/2).cosh() + h*(h/2).sinh()/2)/2
    jump_h = abs(s1 - ys[0])
    TV = sum((abs(ys[i+1]-ys[i]) for i in range(len(ys)-1)), arb(0))
    Scoef1 = jump_h   # leading jump of full v_w at h (v(0)=0, v(H)=0)
    Scoef2 = s1 + TV + s2c  # generous O(1/xi^2)
    Xexact = arb(XBETA_EXACT)
    while stack:
        mid, r = stack.pop()
        fv, fp = fval(mid)
        lb = fv - abs(fp)*r - F2*r*r/2
        ncell += 1
        slack = abs(fp)*r + F2*r*r/2
        if lb > 0 and slack > (fv+arb(c))*arb(REL) and r > arb("1e-5"):
            nrefine[0]+=1; stack += [(mid-r/2,r/2),(mid+r/2,r/2)]; continue
        if lb > 0:
            if worst is None or lb.lower() < worst[0]:
                worst = (float(lb.lower()), float(mid.mid()), float(r.mid()))
            if mid + r <= Xexact:
                # three-point + M2 remainder; require mid>r so left point positive
                if mid <= r:
                    nrefine[0]+=1; stack += [(mid-r/2,r/2),(mid+r/2,r/2)]; continue
                Sm = abs(Sw_pt(h,ys,zs,mid))
                Sl = abs(Sw_pt(h,ys,zs,mid-r))
                Sr = abs(Sw_pt(h,ys,zs,mid+r))
                Su = max(Sm, Sl, Sr) + M2*r*r/8
                beta += 2*r*(4*Su*Su).upper()/(arb(c)+arb(lb.lower()))
            else:
                xL = mid - r
                if xL <= 0: xL = mid
                Su = (Scoef1 + Scoef2/xL)/xL
                beta += 2*r*(4*Su*Su).upper()/(arb(c)+arb(lb.lower()))
        else:
            nsplit += 1
            assert r > arb("1e-9"), ("cannot certify", float(mid.mid()))
            stack += [(mid-r/2,r/2),(mid+r/2,r/2)]
    X = XB
    while X < arb(10)**7:
        Su = (Scoef1 + Scoef2/X)/X
        beta += 4*Su*Su*(1/X - 1/(2*X))/mlow_tail(X)
        X = 2*X
    beta += 4*((Scoef1+Scoef2/X)/X)**2/(X*mlow_tail(X))
    beta = (beta/arb.pi()).upper()
    coer = ((1-2*arb(beta))*arb(c)).lower()
    out = dict(file=fn, h=f"{hq.p}/{hq.q}", H=f"{Hq.p}/{Hq.q}", N=N,
               c=f"{c.p}/{c.q}", C=float(C.mid()), M_upper=float(Mbern.mid()),
               Xi=float(Xi.mid()), tail_margin=float(tail.lower()), F2=float(F2.mid()),
               M2=float(M2.mid()), cells=ncell, splits=nsplit, beta_refinements=nrefine[0],
               REL=REL, XBETA_EXACT=XBETA_EXACT,
               worst_cell_lb=worst[0], worst_at=worst[1],
               beta_upper=float(beta), two_beta=float(2*beta),
               odd_coercivity_lower=float(coer),
               verified=bool(2*arb(beta)<1 and coer>0),
               seconds=round(time.time()-t0,1))
    print(json.dumps(out), flush=True)
    return out

if __name__ == "__main__":
    res = [verify(fn) for fn in sys.argv[1:]]
    json.dump(res, open("s23_verify_results.json","w"), indent=1)
