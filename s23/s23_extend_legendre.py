"""S23 [design, float]: Legendre / knot exterior extension to drive 2β below 1.
v_w = sinh on (-h,h); on [h,H] = sum_{k=0}^{N-1} a_k P_k(2(x-h)/L - 1), zero outside.
Also: continuous-at-h knot (piecewise-linear) variant.
"""
import numpy as np, scipy.special as sp, json, sys, time
from numpy.linalg import solve, lstsq

l2 = np.log(2); kap = l2/np.sqrt(2)

def S0(xi, h):
    a = 0.5
    num = a*np.cosh(a*h)*np.sin(xi*h) - xi*np.sinh(a*h)*np.cos(xi*h)
    den = a*a + xi*xi
    return np.where(den < 1e-30, 0.0, num/den)

def m_fn(xi, h, a, C, u, Dd):
    A = 2*h + Dd/2
    return (sp.psi(0.25 + 0.5j*xi).real - np.log(np.pi)
            - 2*kap*np.cos(xi*l2)
            + C*2*A*np.sinc(A*xi/np.pi)*np.sinc(Dd*xi/(2*np.pi))
            - np.cos(np.outer(np.atleast_1d(xi), u)) @ a)

def legendre_phi(xi, h, H, k, nquad=64):
    """sine transform of P_k(2(x-h)/L-1) on [h,H] via Gauss-Legendre quad."""
    L = H - h
    # nodes on [-1,1]
    nodes, wts = np.polynomial.legendre.leggauss(nquad)
    # x = h + L*(t+1)/2, t in [-1,1]
    t = nodes
    x = h + L*(t+1)/2
    Pk = sp.eval_legendre(k, t)
    # int_h^H Pk sin(xi x) dx = (L/2) int_{-1}^1 Pk(t) sin(xi x(t)) dt
    xi = np.asarray(xi, float)
    # broadcast: (nxi, nquad)
    sx = np.sin(xi[:, None] * x[None, :])
    return (L/2) * (sx * (Pk * wts)[None, :]).sum(axis=1)

def tent_phi(xi, x_left, x_mid, x_right):
    """sine transform of the tent peaking at 1 at x_mid, 0 at x_left/right (support [x_left,x_right])."""
    xi = np.asarray(xi, float)
    # piecewise linear: use explicit antiderivative of (ax+b) sin(xi x)
    def int_lin_sin(a, b, lo, hi, xi):
        # int (a x + b) sin(xi x) dx
        # = - (a x + b) cos(xi x)/xi + a sin(xi x)/xi^2
        def F(x):
            out = np.empty_like(xi)
            small = np.abs(xi) < 1e-12
            # int (a x + b) x^0 ~ for xi=0: int (ax+b) = a(x^2)/2 + b x
            out[small] = a*(hi*hi - lo*lo)/2 + b*(hi - lo)  # whole interval — handled outside
            xs = xi[~small]
            def f(x, xs):
                return -(a*x + b)*np.cos(xs*x)/xs + a*np.sin(xs*x)/(xs*xs)
            out[~small] = f(hi, xs) - f(lo, xs)
            return out
        # fix small xi separately for the interval
        out = np.empty_like(xi, dtype=float)
        small = np.abs(xi) < 1e-12
        out[small] = a*(hi*hi - lo*lo)/2 + b*(hi - lo)
        xs = xi[~small]
        def f(x, xs):
            return -(a*x + b)*np.cos(xs*x)/xs + a*np.sin(xs*x)/(xs*xs)
        out[~small] = f(hi, xs) - f(lo, xs)
        return out
    # left piece [x_left, x_mid]: value = (x - x_left)/(x_mid - x_left)
    # a = 1/(x_mid-x_left), b = -x_left/(x_mid-x_left)
    out = np.zeros_like(xi, dtype=float)
    if x_mid > x_left:
        a = 1.0/(x_mid - x_left); b = -x_left/(x_mid - x_left)
        out = out + int_lin_sin(a, b, x_left, x_mid, xi)
    if x_right > x_mid:
        # right piece [x_mid, x_right]: value = (x_right - x)/(x_right - x_mid)
        a = -1.0/(x_right - x_mid); b = x_right/(x_right - x_mid)
        out = out + int_lin_sin(a, b, x_mid, x_right, xi)
    return out

def design_legendre(cert, H=3.0, N=24, dxi=0.004, Xi=2500.0, ridge=1e-10):
    t0 = time.time()
    d = json.load(open(cert))
    h = float(d["h"]); a = np.array(d["a_num"], float)/d["den"]
    C = d["C_num"]/d["den"]; u = 2*h + d["D"]*np.arange(d["K"]); Dd = d["Dd"]
    xs = np.arange(dxi, Xi, dxi)
    mv = m_fn(xs, h, a, C, u, Dd)
    W = 4.0/(mv*np.pi)*dxi
    s0 = S0(xs, h)
    Phi = np.column_stack([legendre_phi(xs, h, H, k) for k in range(N)])
    A = (Phi*W[:,None]).T @ Phi
    b = (Phi*W[:,None]).T @ s0
    c0 = float(np.dot(s0*W, s0))
    alpha = -solve(A + ridge*np.eye(N), b)
    beta_core = c0 + 2*float(b@alpha) + float(alpha@A@alpha)
    # tail numerical
    xt = np.arange(Xi, Xi*2.5, 0.05)
    Phit = np.column_stack([legendre_phi(xt, h, H, k) for k in range(N)])
    Swt = S0(xt,h) + Phit@alpha
    beta_tail = float(np.trapezoid(4*Swt**2/(m_fn(xt,h,a,C,u,Dd)*np.pi), xt))
    Xlast = float(xt[-1]); s1=np.sinh(h/2)
    m_low = np.log(Xlast/2)-np.log(np.pi)-2*kap-4
    beta_rem = 4*(s1+1)**2/(Xlast*max(m_low,0.1))/np.pi
    beta = beta_core + beta_tail + beta_rem
    out = dict(method="legendre", h=h, H=H, N=N, beta_core=beta_core, beta_tail=beta_tail,
               beta_rem=beta_rem, beta=beta, two_beta=2*beta, alpha=alpha.tolist(),
               seconds=round(time.time()-t0,2))
    print(json.dumps({k:round(v,8) if isinstance(v,float) else v for k,v in out.items() if k!='alpha'}), flush=True)
    return out

def design_knots(cert, H=3.0, N=20, dxi=0.004, Xi=2500.0, continuous=True, ridge=1e-10):
    """N interior free knots on (h,H), plus optional fixed value at h = sinh(h/2), value at H = 0.
    If continuous=False, value at h is also free (N+1 free including edge)."""
    t0 = time.time()
    d = json.load(open(cert))
    h = float(d["h"]); aa = np.array(d["a_num"], float)/d["den"]
    C = d["C_num"]/d["den"]; u = 2*h + d["D"]*np.arange(d["K"]); Dd = d["Dd"]
    xs = np.arange(dxi, Xi, dxi)
    mv = m_fn(xs, h, aa, C, u, Dd)
    W = 4.0/(mv*np.pi)*dxi
    s0 = S0(xs, h)
    # nodes: h = z0 < z1 < ... < z_{N+1} = H  with N interior
    zs = np.linspace(h, H, N+2)
    cols = []; fixed = np.zeros_like(xs)
    # continuous: y0 = sinh(h/2) fixed via tent at z0 (half-tent from z0 to z1)
    # free y_i at z_i for i=1..N; y_{N+1}=0 fixed
    if continuous:
        y0 = np.sinh(h/2)
        # contribution of left half-tent at boundary: tent with peak at z0, but only right half exists
        # equivalent: linear from y0 at h to 0 at z1, times... actually the basis tent at z0
        # supported [z0,z1] peaking at 1: value (z1-x)/(z1-z0)
        fixed = y0 * tent_phi(xs, zs[0], zs[0], zs[1])  # degenerate left=mid=z0 → only right piece
        # Wait tent_phi with x_left=x_mid gives only right piece from mid to right with value (right-x)/...
        # At x=mid, value = 1. Good.
        for i in range(1, N+1):
            cols.append(tent_phi(xs, zs[i-1], zs[i], zs[i+1]))
    else:
        # free value at h as well: tent at z0 (right half) + tents at interiors
        cols.append(tent_phi(xs, zs[0], zs[0], zs[1]))
        for i in range(1, N+1):
            cols.append(tent_phi(xs, zs[i-1], zs[i], zs[i+1]))
        fixed = 0
    Phi = np.column_stack(cols)
    n = Phi.shape[1]
    # S = s0 + fixed + Phi @ alpha
    s_base = s0 + fixed
    A = (Phi*W[:,None]).T @ Phi
    b = (Phi*W[:,None]).T @ s_base
    c0 = float(np.dot(s_base*W, s_base))
    alpha = -solve(A + ridge*np.eye(n), b)
    beta_core = c0 + 2*float(b@alpha) + float(alpha@A@alpha)
    xt = np.arange(Xi, Xi*2.5, 0.05)
    # rebuild fixed/phi on tail
    if continuous:
        y0 = np.sinh(h/2)
        fixed_t = y0 * tent_phi(xt, zs[0], zs[0], zs[1])
        Phit = np.column_stack([tent_phi(xt, zs[i-1], zs[i], zs[i+1]) for i in range(1, N+1)])
    else:
        fixed_t = 0
        cols_t = [tent_phi(xt, zs[0], zs[0], zs[1])]
        cols_t += [tent_phi(xt, zs[i-1], zs[i], zs[i+1]) for i in range(1, N+1)]
        Phit = np.column_stack(cols_t)
    Swt = S0(xt,h) + fixed_t + Phit@alpha
    beta_tail = float(np.trapezoid(4*Swt**2/(m_fn(xt,h,aa,C,u,Dd)*np.pi), xt))
    Xlast=float(xt[-1]); s1=np.sinh(h/2)
    m_low=np.log(Xlast/2)-np.log(np.pi)-2*kap-4
    beta_rem=4*(s1+1)**2/(Xlast*max(m_low,0.1))/np.pi
    beta=beta_core+beta_tail+beta_rem
    out=dict(method="knots", continuous=continuous, h=h, H=H, N=N, n_free=n,
             beta_core=beta_core, beta_tail=beta_tail, beta_rem=beta_rem,
             beta=beta, two_beta=2*beta, alpha=alpha.tolist(), zs=zs.tolist(),
             seconds=round(time.time()-t0,2))
    print(json.dumps({k:round(v,8) if isinstance(v,float) else v for k,v in out.items() if k not in ('alpha','zs')}), flush=True)
    return out

if __name__ == "__main__":
    cert = sys.argv[1]
    mode = sys.argv[2] if len(sys.argv)>2 else "both"
    results=[]
    if mode in ("legendre","both"):
        for H in (2.0, 3.0, 4.0):
            for N in (16, 24, 32, 48):
                o=design_legendre(cert, H=H, N=N)
                results.append(dict(method="leg", H=H, N=N, two_beta=o["two_beta"]))
    if mode in ("knots","both"):
        for H in (2.0, 3.0, 4.0):
            for N in (12, 20, 30, 40):
                for cont in (True, False):
                    o=design_knots(cert, H=H, N=N, continuous=cont)
                    results.append(dict(method="knots", H=H, N=N, cont=cont, two_beta=o["two_beta"]))
    print("BEST", min(results, key=lambda r: r["two_beta"]))
    json.dump(results, open("out_sweep2.json","w"), indent=1)
