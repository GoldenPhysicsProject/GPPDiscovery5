"""S26 [design, float, UNTRUSTED]: S25 jump design with the PL-DC claim-1 multiplier (s26_mult_*.json).
usage: mult.json hs H N1 N2 g1 g2 out"""
import sys, os, json, numpy as np
from fractions import Fraction as Fr
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s25"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s24"))
import s25_design as D, s24_design as S, scipy.special as sp
l2 = np.log(2); kap = l2/np.sqrt(2)
def m_pl(xi, d):
    h = float(Fr(*d["h_q"])); a = np.array(d["a_num"], float)/d["den"]; u = 2*h + float(Fr(*d["D_q"]))*np.arange(len(a))
    ys = list(np.array(d["y_num"], float)/d["den"]) + [0.0]; zs = [float(Fr(*q)) for q in d["z_q"]]; M = len(zs)-1
    sl = [0]+[(ys[i+1]-ys[i])/(zs[i+1]-zs[i]) for i in range(M)]+[0]; ds = [sl[j+1]-sl[j] for j in range(M+1)]
    out = sp.psi(0.25+0.5j*xi).real - np.log(np.pi) - 2*kap*np.cos(xi*l2) - 2*sum(dj*np.cos(xi*zj) for dj, zj in zip(ds, zs))/xi**2
    for k in range(0, len(a), 64): out -= np.cos(np.outer(xi, u[k:k+64]))@a[k:k+64]
    return out
S.m_fn = m_pl; D.S.m_fn = m_pl
if __name__ == "__main__":
    src, hs, H, N1, N2, g1, g2, outf = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
    d = json.load(open(src)); hq = Fr(hs); d["h"] = str(float(hq)); assert [hq.numerator, hq.denominator] == d["h_q"]
    zf = S.knots2(float(hq), H, N1, N2, g1, g2); zden = 10**9
    zq = [Fr(round(z*zden), zden) for z in zf]; zq[0] = hq; zq[-1] = Fr(H).limit_denominator(1000)
    zs = np.array([float(z) for z in zq]); o = D.design_jump(d, zs)
    aden = 10**12; al = [int(round(x*aden)) for x in o["alpha"]]
    d["extension"] = dict(type="piecewise_linear_odd_jump_at_h", H_q=[zq[-1].numerator, zq[-1].denominator], N=len(al)-1,
        alpha_num=al, alpha_den=aden, z_num=[int(z*zden) for z in zq], z_den=zden,
        knots=dict(N1=N1, N2=N2, g1=g1, g2=g2), two_beta_design_float=o["two_beta"], mX=o["mX"], Xi_design=o["Xi"])
    json.dump(d, open(outf, "w")); print(outf, "2beta_design=", o["two_beta"], "N=", d["extension"]["N"], flush=True)
