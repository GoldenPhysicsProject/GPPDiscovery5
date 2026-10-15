"""S17 rigorous (Arb) prime-free-window lower bounds for E_eps(g)/||g||^2 and the completed eta Weil form.
 phi(t)=1-Re F(t)/F(0) in [0,2],  E_eps = int_0^inf nu_eps phi dt,  L=2h, phi=1 on [L,inf).
 (C) Codex:    phi>=0                          -> E >= T(L)
 (A) Astra:    phi >= a(t)=1-cos(pi/(ceil(L/t)+1)) (sharp numerical radius of P_h U_t P_h) -> b_N(h)
 (S) S17 bathtub+Schur: int_0^s Re F <= min(s,h)F(0) (Schur row-sum of 1_{|x-y|<s} on (-h,h))
               => Phi(s)=int_0^s phi >= (s-h)_+ ; nu decreasing => E >= T(h)   [closed form]
 (S+A) combined: Phi >= max((s-h)_+, A_N(s)); layer cake => E >= sum_{n=3}^N a_n[T(L/n)-T(L/(n-1))] + T(h)/2 + T(s*)/2,
               s* = h + 2 A_N(h),  A_N(h)=sum_{n=3}^N a_n L(1/(n-1)-1/n).
 Completed form: zeta >= E-bound - c0 - 2(sinh h - h) (pole bound; on even g the pole is >=0 and the -2(sinh h-h) drops);
 nontrivial primitive eta >= E_eps-bound + log f - c_eps (needs Re eta(n)=0 for prime powers n<e^{2h})."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_build import Tarb, consts, arb, ctx
import mpmath as mp
ctx.prec = 256
NCH = 200
def a(n): return 1-(arb.pi()/(n+1)).cos()
def E_codex(h, e): return Tarb(2*h, e)
def E_astra(h, e, N=NCH):
    L = 2*h; s = Tarb(L, e)
    for n in range(2, N+1): s += a(n)*(Tarb(L/n, e)-Tarb(L/(n-1), e))
    return s
def E_bath(h, e): return Tarb(h, e)
def E_comb(h, e, N=NCH):
    L = 2*h; A = sum((a(n)*L*(arb(1)/(n-1)-arb(1)/n) for n in range(3, N+1)), arb(0)); ss = h+2*A
    s = sum((a(n)*(Tarb(L/n, e)-Tarb(L/(n-1), e)) for n in range(3, N+1)), arb(0))
    return s+Tarb(h, e)/2+Tarb(ss, e)/2
C0, C1 = consts()
def zeta_lb(Ef, h, even=False): return Ef(h, 0)-C0-(0 if even else 2*(h.sinh()-h))
if __name__ == "__main__":
    print("h | zeta LB: Codex | Astra(N=200) | S17 bathtub | S17 bathtub+Astra   (all g; pole bound included)")
    for hq in [(1, 32), (1, 25), (1, 16), (1, 15), (1, 12), (1, 11), (1, 10)]:
        h = arb(hq[0])/hq[1]
        print(f"{hq[0]}/{hq[1]} | " + " | ".join(zeta_lb(f, h).str(10) for f in (E_codex, E_astra, E_bath, E_comb)))
    # thresholds by bisection on midpoints (then certified by evaluating the ball at the reported h)
    def root(f, target, lo=0.01, hi=0.3):
        for _ in range(60):
            m = (lo+hi)/2
            if (f(arb(m))-target).mid() > 0: lo = m
            else: hi = m
        return lo
    print("largest h (bisection, 1e-15) with zeta LB > target; ball at that h:")
    for nm, f in (("Codex", E_codex), ("Astra", E_astra), ("bathtub", E_bath), ("bathtub+Astra", E_comb)):
        for tgt in (0.01, 0.0):
            hh = root(lambda h: zeta_lb(f, h), tgt)
            print(f"  {nm:14s} target {tgt}: h = {hh:.12f}  LB(h) = {zeta_lb(f, arb(hh)).str(8)}")
    print("Even-g zeta sector (pole >= 0): bathtub+Astra threshold h =", f"{root(lambda h: zeta_lb(E_comb,h,True),0.0):.12f}")
    print("Nontrivial characters, bathtub+Astra (valid while prime-free):")
    from s17_build import CHARS, vM
    for cn, f, e, ch in CHARS:
        lb = lambda h: E_comb(h, e)+arb(f).log()-(C0 if e == 0 else C1)
        n = 2
        while not (vM(n) and ch(n) != 0): n += 1
        hp = mp.log(n)/2
        print(f"  {cn:45s} first prime power with Re eta != 0: n={n} (window h<{mp.nstr(hp,8)}); LB>0 up to h={root(lb,0.0,0.01,1.5):.10f}; LB(1/25)={lb(arb(1)/25).str(8)}")
