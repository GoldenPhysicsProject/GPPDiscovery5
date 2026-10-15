"""S17 self-tests: (a) Gram N vs direct quadrature; (b) zeta form vs S14 exact A (independent arch decomposition)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from s17_build import *
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s14"))
import s14_exact, mpmath as mp
M = 10
for hs in ["0.04", "0.3", "0.5"]:
    B = build(M, hs); Q = forms(B, hs)["zeta"]
    A, T = s14_exact.build(M, hs, 256)          # S14: pole+arch, and F(log2)
    if arb(2).log() < 2*arb(hs): A = A - T*(2*arb(2).log()/arb(2).sqrt())
    d = max(abs((Q[i, j]-A[i, j]).mid()) for i in range(M) for j in range(M))
    # Gram check by mpmath quadrature
    h = mp.mpf(hs); gj = lambda j: (lambda x: mp.legendre(j, x/h)*(1-(x/h)**2)**4)
    dn = max(abs(mp.quad(lambda x: gj(i)(x)*gj(j)(x), [-h, h]) - mp.mpf(B["N"][i, j].mid().str(40, radius=False))) for i in range(4) for j in range(4))
    print(f"h={hs}: max|Q_zeta(S17) - A(S14)| = {mp.nstr(mp.mpf(d.str(30,radius=False)),3)} ; max|N - quad| = {mp.nstr(dn,3)}")
