"""S22 [design, float]: beta = (1/2pi) int |vhat|^2 / m  for v = sinh(x/2) 1_(-h,h) (no extension), m from s22_cert_h*.json."""
import numpy as np, scipy.special as sp, json, sys
from scipy.integrate import quad
l2 = np.log(2); kap = l2/np.sqrt(2)
def S(xi, h):
    out = 0
    for al, sg in ((0.5, 1), (-0.5, -1)):
        z = al+1j*xi; out = out + sg*((np.exp(z*h)-1)/z).imag
    return out/2
d = json.load(open(sys.argv[1])); h = np.log(3)/2 if d["h"] == "hstar" else float(d["h"])
a = np.array(d["a_num"])/d["den"]; C = d["C_num"]/d["den"]; u = 2*h+d["D"]*np.arange(d["K"]); Dd = d["Dd"]
def m(xi):
    A = 2*h+Dd/2
    return sp.psi(0.25+0.5j*xi).real-np.log(np.pi)-2*kap*np.cos(xi*l2)+C*2*A*np.sinc(A*xi/np.pi)*np.sinc(Dd*xi/2/np.pi)-np.cos(np.outer(np.atleast_1d(xi), u))@a
xs = np.arange(0, 3000, 0.002); mv = m(xs); print("min m on grid", mv.min(), "at", xs[np.argmin(mv)])
f = 4*S(xs, h)**2/mv
core = np.trapezoid(f, xs)/np.pi
# tail beyond 3000: |S|<= (e^{h/2}+1)/xi roughly -> bound int 4 S^2/m <= 4*(cosh(h/2)+1)^2/(m_tail*3000)
tail = 4*(np.cosh(h/2)+1)**2/(np.log(1500)-np.log(np.pi)-2*kap-3.2)/3000/np.pi
print(json.dumps(dict(h=d["h"], beta_core=core, tail_bound=tail, two_beta=2*(core+tail), c_LP=d["c_LP"])))
