"""ARP section C step tests (S103, 2026-10-28). Replays only; each block tests one named C step.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok.  RH is OPEN.

(1) C.5  F_P(sigma) = prod_{p<=P} R_p(sigma) = zeta_P(1)/sqrt(zeta_P(2 sigma) zeta_P(2-2 sigma)),
         zeta_P = truncated Euler product.  Reproduce Astra's table (Oct 26) at P = 100, 1000, 10000.
(2) C.5  Fidelity susceptibility chi_P = 4 sum (log p)^2/(p(1-1/p)^2) vs 2 (log P)^2 (mean-field data only).
(3) C.6  Astra's raw doubled kernel Q_p = d(x)d(y)h_r(x-y): h_r(v) > h_r(0) for 0<|v|<log p, so h_r is
         not positive definite (|h(v)|<=h(0) fails) and the 2x2 reflected Gram has det < 0.
(4) C.2  R5 test of the BC/primon OS step against the must-pass control G_theta' and the fake F_theta:
         K(tau) = sum_n w_n n^{-1/2-tau} with w_n = Lambda(n)(n^{i a}+n^{-i a}) (G_a, a real) and
         w_n = Lambda(n)(n^{a}+n^{-a}) (F_a).  Hankel OS Gram [K(t_i+t_j)], t_i > 1/4.
         Widder: PSD for all families <=> w_n >= 0.  G_theta' must FAIL, F_theta must PASS.
"""
import mpmath as mp
mp.mp.dps = 40

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

print("(1) mirror overlap = ratio of truncated Euler products")
for P in (100, 1000, 10000):
    ps = primes_upto(P)
    for sig in (mp.mpf("0.30"), mp.mpf("0.40"), mp.mpf("0.45"), mp.mpf("0.5")):
        prodR = mp.mpf(1); zP = lambda s: mp.fprod([1/(1-mp.mpf(p)**(-s)) for p in ps])
        for p in ps:
            prodR *= mp.sqrt((1-mp.mpf(p)**(-2*sig))*(1-mp.mpf(p)**(-2*(1-sig))))/(1-mp.mpf(1)/p)
        ratio = zP(1)/mp.sqrt(zP(2*sig)*zP(2-2*sig))
        print(f"  P={P:6d} sigma={float(sig):.2f}  prodR={mp.nstr(prodR,8):>12}  zeta_P ratio={mp.nstr(ratio,8):>12}  diff={mp.nstr(abs(prodR-ratio),3)}")

print("(2) susceptibility (uses only Mertens/PNT mean)")
for P in (10**3, 10**4, 10**5):
    ps = primes_upto(P)
    chi = sum(4*mp.log(p)**2/(p*(1-mp.mpf(1)/p)**2) for p in ps)
    print(f"  P={P:7d} chi_P={mp.nstr(chi,10)}  2(log P)^2={mp.nstr(2*mp.log(P)**2,10)}")

print("(3) raw doubled kernel h_r(v) = 1/(1+r^2-2r cosh v), r=1/p")
for p in (2, 3, 5, 17):
    r = mp.mpf(1)/p; h = lambda v: 1/(1+r*r-2*r*mp.cosh(v))
    lp = mp.log(p); x = mp.mpf("0.3")*lp; y = mp.mpf("-0.1")*lp   # sigma=0.8, tau=0.4
    d = lambda z: mp.sqrt((1-r*mp.e**(2*z))*(1-r*mp.e**(-2*z)))
    det = (d(x)*d(y))**2*(h(0)**2-h(x-y)**2)
    print(f"  p={p:2d} h(0)={mp.nstr(h(0),8)} h(x-y)={mp.nstr(h(x-y),8)} (> h(0): {h(x-y)>h(0)})  2x2 det={mp.nstr(det,6)}")

print("(4) Hankel OS Gram of the BC/primon kernel: G_theta' (must-pass) vs F_theta")
def K_G(tau, a):  # sum Lambda(n) n^{-1/2-tau} 2cos(a log n) = -2 Re zeta'/zeta(1/2+tau+ia)
    s = mp.mpf(1)/2 + tau + 1j*a
    return -2*mp.re(mp.zeta(s, derivative=1)/mp.zeta(s))
def K_F(tau, a):  # sum Lambda(n) n^{-1/2-tau}(n^a+n^-a), need 1/2+tau-a>1
    s1 = mp.mpf(1)/2 + tau - a; s2 = mp.mpf(1)/2 + tau + a
    return -(mp.zeta(s1, derivative=1)/mp.zeta(s1)) - (mp.zeta(s2, derivative=1)/mp.zeta(s2))
ts = [mp.mpf(x) for x in ("0.35", "0.6", "1.0", "1.6", "2.5", "4.0", "6.0")]
for name, K, a in (("zeta (a=0)", K_G, mp.mpf(0)), ("G_theta', a=pi/log2", K_G, mp.pi/mp.log(2)),
                   ("G_theta', a=14.134725", K_G, mp.mpf("14.134725")), ("F_theta, a=0.1", K_F, mp.mpf("0.1"))):
    tt = ts if name[0] != "F" else [t + mp.mpf("0.1") for t in ts]
    M = mp.matrix(len(tt), len(tt))
    for i, ti in enumerate(tt):
        for j, tj in enumerate(tt):
            M[i, j] = K(ti + tj, a)
    ev = mp.eigsy(M)[0]
    print(f"  {name:24s} min eig = {mp.nstr(min(ev),6):>12}   diag K(2*6)={mp.nstr(K(12+(0.2 if name[0]=='F' else 0), a),6)}")
