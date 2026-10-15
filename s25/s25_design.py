"""S25 [design, float, UNTRUSTED]: PL odd extension of v = sinh(x/2)·1_(-h,h) with a FREE value just outside h
(jump allowed at |x| = h), graded knots (S24 knots2), honest Parseval tail as in S24 design_tail.
beta = (4/pi) int_0^inf S^2/m,  S(xi) = int_0^inf v_w(t) sin(xi t) dt."""
import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "s24"))
import s24_design as S
def design_jump(d, zs, dxi=0.01, Xi=1500.0, ridge=1e-12):
    h = float(d["h"]); xs = np.arange(dxi/2, Xi, dxi); mv = S.m_fn(xs, d)
    xt = np.linspace(Xi, 4*Xi, 20001); mX = float(S.m_fn(xt, d).min())
    W = 4/np.pi*(1/mv - 1/mX)*dxi
    n = len(zs) - 2                      # free values y_0..y_n at zs[0..n]; y_{n+1}=0 at zs[-1]=H
    base = S.S0(xs, h)
    Phi = np.column_stack([S.tent(xs, zs[0], zs[0], zs[1])] + [S.tent(xs, zs[i-1], zs[i], zs[i+1]) for i in range(1, n+1)])
    A = (Phi*W[:, None]).T@Phi; b = (Phi*W[:, None]).T@base; c0 = float(base*W@base)
    dz = np.diff(zs); Mf = np.zeros((n+2, n+2))
    for i in range(n+1):
        Mf[i, i] += dz[i]/3; Mf[i+1, i+1] += dz[i]/3; Mf[i, i+1] += dz[i]/6; Mf[i+1, i] += dz[i]/6
    k = 2/mX
    A += k*Mf[:n+1, :n+1]; c0 += k*(np.sinh(h)-h)/2
    al = -np.linalg.solve(A + ridge*np.eye(n+1), b); bc = c0 + 2*b@al + al@A@al
    return dict(beta=float(bc), two_beta=float(2*bc), alpha=al.tolist(), zs=list(map(float, zs)), mX=mX, Xi=Xi)
if __name__ == "__main__":
    d = json.load(open(sys.argv[1])); h = float(d["h"]); res = []
    for H in map(float, sys.argv[2].split(",")):
        for N1, N2 in [tuple(map(int, p.split(":"))) for p in sys.argv[3].split(",")]:
            for g1 in map(float, sys.argv[4].split(",")):
                for g2 in map(float, sys.argv[5].split(",")):
                    t0 = time.time(); zs = S.knots2(h, H, N1, N2, g1, g2); o = design_jump(d, zs)
                    print(f"H={H} N1={N1} N2={N2} g1={g1} g2={g2}: 2beta={o['two_beta']:.5f} y0={o['alpha'][0]:.4f} ({time.time()-t0:.0f}s)", flush=True)
                    res.append(dict(H=H, N1=N1, N2=N2, g1=g1, g2=g2, two_beta=o["two_beta"], y0=o["alpha"][0]))
    json.dump(res, open(sys.argv[6] if len(sys.argv) > 6 else "s25_design_sweep.json", "w"), indent=1)
