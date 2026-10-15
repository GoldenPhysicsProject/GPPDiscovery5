"""S22 [NUMERICAL sanity]: for odd basis functions b_j (S17 basis, h=2/5), compare the S17 Arb matrix (Q_zeta - P)_jk with
(1/2pi) int ghat_j conj ghat_k m(xi) dxi for the certified m (includes the DC and cosine terms, which must contribute 0)."""
import sys, os, json; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s17"))
from s17_build import *
import numpy as np, scipy.special as sp
d = json.load(open("s22_cert_h0.4.json")); h = 0.4; M = 8
B = build(M, arb(2)/5); F = forms(B, arb(2)/5); A = F["zeta"] - B["P"]
xq, wq = np.polynomial.legendre.leggauss(400); x = xq*h; w = wq*h
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's14')); from s14_exact import p as bp
bas = [np.polyval([float(int(c.p))/int(c.q) for c in bp(j).coeffs()][::-1], xq) for j in range(M)]
l2 = np.log(2); kap = l2/np.sqrt(2); a = np.array(d["a_num"])/d["den"]; C = d["C_num"]/d["den"]; u = 2*h+d["D"]*np.arange(d["K"]); Dd = d["Dd"]
xi = np.arange(0, 800, 0.002)
def m(xi):
    Aw = 2*h+Dd/2
    return sp.psi(0.25+0.5j*xi).real-np.log(np.pi)-2*kap*np.cos(xi*l2)+C*2*Aw*np.sinc(Aw*xi/np.pi)*np.sinc(Dd*xi/2/np.pi)-np.cos(np.outer(xi, u))@a
mv = m(xi); sA = sp.psi(0.25+0.5j*xi).real-np.log(np.pi)-2*kap*np.cos(xi*l2)
gh = {j: (bas[j]*w) @ np.exp(-1j*np.outer(x, xi)) for j in (1, 3, 5)}
for j, k in ((1, 1), (1, 3), (3, 5), (5, 5)):
    pr = (gh[j]*np.conj(gh[k])).real
    vm = 2*np.trapezoid(pr*mv, xi)/(2*np.pi); vs = 2*np.trapezoid(pr*sA, xi)/(2*np.pi)
    ref = float(A[j, k].mid().str(20, radius=False))
    print(f"({j},{k}) S17 {ref:.9f}  symbol sigma_A {vs:.9f}  certified m {vm:.9f}")
