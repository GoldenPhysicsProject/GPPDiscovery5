"""S27 [design, float, UNTRUSTED]: S25 jump design with a PL-W claim-1 multiplier at exact h* = log3/2 (s27_mult_*.json).
Extension knots stored as rational offsets from h (off_num/off_den), H = h + last offset.
usage: mult.json H N1 N2 g1 g2 out   (RH OPEN)"""
import sys, os, json, numpy as np
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "s25")); sys.path.insert(0, os.path.join(here, "..", "s24"))
import s25_design as D, s24_design as S, scipy.special as sp
l2 = np.log(2); kap = l2/np.sqrt(2); HSTAR = np.log(3)/2
def m_pl(xi, d):
    h = HSTAR; a = np.array(d["a_num"], float)/d["den"]; u = 2*h + float(Fr(*d["D_q"]))*np.arange(len(a))
    ys = list(np.array(d["y_num"], float)/d["den"]) + [0.0]; zs = [2*h + float(Fr(*q)) for q in d["w_off_q"]]; M = len(zs)-1
    sl = [0]+[(ys[i+1]-ys[i])/(zs[i+1]-zs[i]) for i in range(M)]+[0]; ds = [sl[j+1]-sl[j] for j in range(M+1)]
    out = sp.psi(0.25+0.5j*xi).real - np.log(np.pi) - 2*kap*np.cos(xi*l2) - 2*sum(dj*np.cos(xi*zj) for dj, zj in zip(ds, zs))/xi**2
    for k in range(0, len(a), 64): out -= np.cos(np.outer(xi, u[k:k+64]))@a[k:k+64]
    return out
S.m_fn = m_pl; D.S.m_fn = m_pl
if __name__ == "__main__":
    src, H, N1, N2, g1, g2, outf = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
    d = json.load(open(src)); assert d["h_spec"] == "log3/2"; d["h"] = repr(HSTAR)
    zf = S.knots2(HSTAR, H, N1, N2, g1, g2); oden = 10**9
    on = [int(round((z - HSTAR)*oden)) for z in zf]; on[0] = 0
    assert all(on[i] < on[i+1] for i in range(len(on)-1))
    zs = np.array([HSTAR + n/oden for n in on]); o = D.design_jump(d, zs)
    aden = 10**12; al = [int(round(x*aden)) for x in o["alpha"]]
    d["extension"] = dict(type="piecewise_linear_odd_jump_at_h", N=len(al)-1, alpha_num=al, alpha_den=aden, off_num=on, off_den=oden,
        H_float=float(zs[-1]), knots=dict(N1=N1, N2=N2, g1=g1, g2=g2), two_beta_design_float=o["two_beta"], mX=o["mX"], Xi_design=o["Xi"])
    json.dump(d, open(outf, "w")); print(outf, "2beta_design=", o["two_beta"], "N=", d["extension"]["N"], flush=True)
