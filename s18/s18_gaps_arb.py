"""S18 claims 1 + 3 + 4(analytic): CERTIFIED (Arb, prec 256).
 (1) B(h)=2 sum_{n<e^{2h}} Lambda(n)/sqrt n  <=  8 h e^h  on a grid (+ h*), looseness ratio; tighter elementary 2*(2h)*(2e^h-1) (same proof, kept -1)
     and the Rosser-Schoenfeld-based psi(x)<1.03883x partial-summation bound 4*1.03883 e^h (cited, not elementary).
 (3) T0(h*),T1(h*),c0,c1,(log2)/sqrt2 and the two Codex gaps.
 (4) analytic LB per residual primitive character: LB_T = log f + T_e(h*) - c_e - (log2)/sqrt2 |Re eta(2)|
     and LB_C = same with the S17 bathtub+Astra-chain E_comb(h*) in place of T_e(h*). Valid for every 0<h<=h* (test-space inclusion)."""
from s18_common import *
C0, C1 = consts()
print("=== (1) B(h) <= 8 h e^h  (CERTIFIED: each line is an Arb comparison) ===")
def B(h):
    s = arb(0); n = 2
    while arb(n).log() < 2*h:
        if vM(n): s += 2*arb(vM(n)).log()/arb(n).sqrt()
        n += 1
    if not (arb(n).log() > 2*h): s = None   # undecided endpoint (n = e^{2h} exactly): F(log n)=0 on open support, B uses strict <
    return s, n
grid = [arb(k)/100 for k in (1, 5, 10, 20, 30, 34, 35, 40, 50, 54, 55, 60, 69, 70, 80, 100, 120, 150, 175, 200)] + [HSTAR]
allok = True
for h in grid:
    b, n = B(h)
    if b is None:
        b, _ = B(h - arb(10)**-40); note = " (endpoint n=%d excluded: strict <)" % (n)
    else: note = ""
    ub = 8*h*h.exp(); ub2 = 4*h*(2*h.exp()-1); rs = 4*arb("1.03883")*h.exp()
    ok = bool(b < ub); allok &= ok
    print(f"h={h.str(12):34s} B={b.str(12):30s} 8he^h={ub.str(10):26s} ok={ok}  B/8he^h={(b/ub).str(6)}  "
          f"4h(2e^h-1)={ub2.str(8)} [{bool(b<ub2)}]  4.1553e^h={rs.str(8)} [{bool(b<rs)}]{note}")
print("ALL B(h) < 8he^h on grid:", allok)
# B at h* only n=2
b, _ = B(HSTAR - arb(10)**-40); print("B(h*-1e-40) [= B on (0,h*] since n=3 needs log3<2h strictly] =", b.str(20), " = 2(log2)/sqrt2 =", (2*LOG2S2).str(20), " overlap:", b.overlaps(2*LOG2S2))
print("prime powers n with log n < 2h*: ", [n for n in range(2, 4) if vM(n) and bool(arb(n).log() < 2*HSTAR)],
      "; log 3 < 2h* is", bool(arb(3).log() < 2*HSTAR), "(equality: F(log3)=0 on open support)")
print("h* < log 2 :", bool(HSTAR < arb(2).log()), " (log2 - h* =", (arb(2).log()-HSTAR).str(15), ")")
print("\n=== (3) constants and Codex gaps at h* = log3/2 (CERTIFIED) ===")
T0, T1 = T_eps(HSTAR, 0), T_eps(HSTAR, 1)
for nm, v in (("h*", HSTAR), ("T0(h*)", T0), ("T1(h*)", T1), ("c0", C0), ("c1", C1), ("(log2)/sqrt2", LOG2S2),
              ("exp(c0+(log2)/sqrt2)", (C0+LOG2S2).exp())):
    print(f"{nm:22s} {v.str(25)}")
g14 = arb(14).log()+T0-C0-LOG2S2; g13 = arb(13).log()+T0-C0-LOG2S2
g8 = arb(8).log()+T1-C1-LOG2S2; g7 = arb(7).log()+T1-C1-LOG2S2
print("even f>=14 gap log14+T0-c0-(log2)/sqrt2 =", g14.str(20), " >0.0679298:", bool(g14 > arb("0.0679298")))
print("  (f=13 would give", g13.str(12), ")")
print("odd  f>=8  gap log8 +T1-c1-(log2)/sqrt2 =", g8.str(20), " >0.0508415:", bool(g8 > arb("0.0508415")))
print("  (f=7 would give", g7.str(12), ")")
g352 = arb(352).log()-C0-LOG2S2; print("Codex sec.5 f>=352 gap (E>=0 only):", g352.str(20), " f=351:", (arb(351).log()-C0-LOG2S2).str(10))
print("\n  monotonicity: dT0/dr>0, dT1/dr = 2[1/(1-r^2)-1/(1+r^2)]>0, r=e^{-h/2} decreasing => T_e(h)>=T_e(h*) on (0,h*]  (PROVED)")
print("\n=== (3b) S17 bathtub+Astra combined E_comb(h*) (CERTIFIED, NCH=200) ===")
EC0, EC1 = E_comb(HSTAR, 0), E_comb(HSTAR, 1)
print("E_comb0(h*) =", EC0.str(15), " vs T0(h*)", T0.str(10), "  gain", (EC0-T0).str(8))
print("E_comb1(h*) =", EC1.str(15), " vs T1(h*)", T1.str(10), "  gain", (EC1-T1).str(8))
for e, (TT, EC, CC) in enumerate(((T0, EC0, C0), (T1, EC1, C1))):
    fm = (CC+LOG2S2-EC).exp(); print(f"  parity {e}: worst-case (|eta(2)|=1) threshold with E_comb: f > {fm.str(10)} ; with T: f > {(CC+LOG2S2-TT).exp().str(10)}")
    fm0 = (CC-EC).exp(); print(f"  parity {e}: eta(2)=0 threshold with E_comb: f > {fm0.str(10)} ; with T: f > {(CC-TT).exp().str(10)}")
print("\n=== (4) residual primitive characters at h* (analytic LB, CERTIFIED) ===")
print("f par label ord  eta(2)                      |Re eta2|   LB_T(worst)   LB_T(actual)    LB_C(actual)   status")
RES = {0: [5, 7, 8, 9, 11, 12, 13], 1: [3, 4, 5, 7]}
for e in (0, 1):
    TT, EC, CC = (T0, EC0, C0) if e == 0 else (T1, EC1, C1)
    for f in RES[e]:
        ch = [c for c in primitive_chars(f) if c[1] == e]
        if not ch: print(f"{f:2d} {e}   -- no primitive character of this parity --"); continue
        for l, par, o, e2, _ in ch:
            re = abs(e2.real)
            base = arb(f).log()-CC
            lbw = base+TT-LOG2S2; lbt = base+TT-LOG2S2*re; lbc = base+EC-LOG2S2*re
            st = "PROVED>0 (T)" if lbt > 0 else ("PROVED>0 (E_comb)" if lbc > 0 else ("NEG/undecided" if not lbc > 0 else "?"))
            print(f"{f:2d} {e}  {l:3d}  {o:2d}  {e2.str(8):28s} {re.str(5):10s} {lbw.str(8):14s} {lbt.str(8):14s} {lbc.str(8):14s} {st}")
print("\n=== (3c) improved universal first-prime thresholds using E_comb (|eta(2)|<=1 worst case, CERTIFIED) ===")
ge12 = arb(12).log()+EC0-C0-LOG2S2; ge11 = arb(11).log()+EC0-C0-LOG2S2
go7 = arb(7).log()+EC1-C1-LOG2S2; go6 = arb(6).log()+EC1-C1-LOG2S2
print("even f>=12: log12+E_comb0(h*)-c0-(log2)/sqrt2 =", ge12.str(15), "  (f=11:", ge11.str(10), ")")
print("odd  f>=7 : log7 +E_comb1(h*)-c1-(log2)/sqrt2 =", go7.str(15), "  (f=6:", go6.str(10), ")")
print("phi check: chain coefficient a(4) = 1-cos(pi/5) =", (1-(arb.pi()/5).cos()).str(20), " = 1-phi/2 =", (1-(1+arb(5).sqrt())/4).str(20))
