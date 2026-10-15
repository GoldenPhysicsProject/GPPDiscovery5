"""S23 [design, float]: optimize odd extension of v=sinh(x/2) outside (-h,h) to minimize beta.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH OPEN.
v_w = sinh(x/2) on |x|<h; on [h,H]:
  alpha0 * (H-x)/L + sum_{j=1}^J alpha_j * sin(j*pi*(x-h)/L),  L=H-h, mirrored oddly.
(linear ramp free at the window edge + sine corrections; v=0 for |x|>H).
beta = (1/pi) int_0^inf 4 S_w^2 / m = c0 + 2 b·α + αᵀ A α;  α = -A^{-1} b.
"""
import numpy as np, scipy.special as sp, json, sys, time
from numpy.linalg import solve

l2 = np.log(2); kap = l2/np.sqrt(2)

def S0(xi, h):
    a = 0.5
    num = a*np.cosh(a*h)*np.sin(xi*h) - xi*np.sinh(a*h)*np.cos(xi*h)
    den = a*a + xi*xi
    return np.where(den < 1e-30, 0.0, num/den)

def _isinc(w, L):
    """vectorized int_0^L sin(w t) dt = (1-cos(w L))/w, = L if w=0."""
    out = np.empty_like(w, dtype=float)
    small = np.abs(w) < 1e-12
    out[small] = L
    ws = w[~small]
    out[~small] = (1 - np.cos(ws*L))/ws
    return out

def _icosc(w, L):
    """int_0^L cos(w t) dt = sin(wL)/w, = L if w=0."""
    out = np.empty_like(w, dtype=float)
    small = np.abs(w) < 1e-12
    out[small] = L
    ws = w[~small]
    out[~small] = np.sin(ws*L)/ws
    return out

def phi_linear(xi, h, H):
    """sine transform of (H-x)/L on [h,H]. t=x-h in [0,L], f=(L-t)/L.
    int_0^L (1 - t/L) sin(xi(t+h)) dt."""
    xi = np.asarray(xi, float)
    L = H - h
    # int_0^L sin(xi(t+h)) = cos(xi h)*_isinc(xi,L) +? wait
    # sin(xi(t+h)) = sin(xit)cos(xih)+cos(xit)sin(xih)
    I_sin = _isinc(xi, L)           # int sin(xi t)
    I_cos = _icosc(xi, L)           # int cos(xi t)
    # int t sin(xi t): integrate by parts
    # int_0^L t sin(xt) = -L cos(xL)/x + int cos(xt)/x = (sin(xL)/x^2 - L cos(xL)/x) ; x=0 -> L^2/2
    def int_t_sin(x, L):
        out = np.empty_like(x)
        small = np.abs(x) < 1e-12
        out[small] = L*L/2
        xs = x[~small]
        out[~small] = np.sin(xs*L)/(xs*xs) - L*np.cos(xs*L)/xs
        return out
    def int_t_cos(x, L):
        out = np.empty_like(x)
        small = np.abs(x) < 1e-12
        out[small] = 0.0  # int_0^L t = wait int t cos(0)=L^2/2 ... actually at 0: L^2/2
        # int_0^L t cos(xt) = L sin(xL)/x - int sin(xt)/x = L sin/x + (cos(xL)-1)/x^2
        out[small] = L*L/2
        xs = x[~small]
        out[~small] = L*np.sin(xs*L)/xs + (np.cos(xs*L) - 1)/(xs*xs)
        return out
    # int (1-t/L) sin(xi(t+h)) = cos(xih)*int(1-t/L)sin(xit) + sin(xih)*int(1-t/L)cos(xit)
    Is = I_sin - int_t_sin(xi, L)/L
    Ic = I_cos - int_t_cos(xi, L)/L
    return Is*np.cos(xi*h) + Ic*np.sin(xi*h)

def phi_sine(xi, h, H, j):
    """sine transform of sin(j*pi*(x-h)/L) on [h,H], vectorized."""
    xi = np.asarray(xi, float)
    L = H - h
    om = j * np.pi / L
    # I1 = int_0^L sin(om t) sin(xi t)
    # = [sin((om-xi)L)/(2(om-xi)) - sin((om+xi)L)/(2(om+xi))] with limits
    def I_ss(om, xi, L):
        out = np.empty_like(xi)
        d = om - xi; s = om + xi
        # near om≈xi
        near = np.abs(d) < 1e-10
        out[near] = L/2 - np.sin(2*om*L)/(4*om) if abs(om) > 1e-14 else 0.0
        # for near elements individually if om varies — here om scalar
        if np.any(near):
            out[near] = 0.5*L - (np.sin(2*om*L)/(4*om) if abs(om)>1e-14 else 0.0)
        dn = ~near
        out[dn] = np.sin(d[dn]*L)/(2*d[dn]) - np.sin(s[dn]*L)/(2*s[dn])
        return out
    def I_sc(om, xi, L):
        # int sin(om t) cos(xi t) = 0.5*(isinc(om+xi)+isinc(om-xi))
        return 0.5*(_isinc(om+xi, L) + _isinc(om-xi, L))
    I1 = I_ss(om, xi, L)
    I2 = I_sc(om, xi, L)
    return I1*np.cos(xi*h) + I2*np.sin(xi*h)

def m_fn(xi, h, a, C, u, Dd):
    A = 2*h + Dd/2
    return (sp.psi(0.25 + 0.5j*xi).real - np.log(np.pi)
            - 2*kap*np.cos(xi*l2)
            + C*2*A*np.sinc(A*xi/np.pi)*np.sinc(Dd*xi/(2*np.pi))
            - np.cos(np.outer(np.atleast_1d(xi), u)) @ a)

def design(cert_path, H=1.5, J=16, dxi=0.005, Xi=2000.0, save=None, match_edge=True):
    t0 = time.time()
    d = json.load(open(cert_path))
    h = float(d["h"]) if d["h"] != "hstar" else np.log(3)/2
    a = np.array(d["a_num"], float)/d["den"]
    C = d["C_num"]/d["den"]
    u = 2*h + d["D"]*np.arange(d["K"])
    Dd = d["Dd"]
    xs = np.arange(dxi, Xi, dxi)  # skip 0
    mv = m_fn(xs, h, a, C, u, Dd)
    assert mv.min() > 0, ("m not positive", mv.min(), xs[np.argmin(mv)])
    W = (4.0 / (mv * np.pi)) * dxi
    s0 = S0(xs, h)
    cols = []
    names = []
    if match_edge:
        cols.append(phi_linear(xs, h, H)); names.append("lin")
    for j in range(1, J+1):
        cols.append(phi_sine(xs, h, H, j)); names.append(f"s{j}")
    Phi = np.column_stack(cols)
    n = Phi.shape[1]
    A = (Phi * W[:, None]).T @ Phi
    b = (Phi * W[:, None]).T @ s0
    c0 = float(np.dot(s0 * W, s0))
    # mild ridge
    alpha = -solve(A + 1e-12*np.eye(n), b)
    beta_core = c0 + 2*float(np.dot(b, alpha)) + float(alpha @ A @ alpha)
    # numerical tail
    xt = np.arange(Xi, min(Xi*3, 12000), 0.05)
    s0t = S0(xt, h)
    cols_t = []
    if match_edge:
        cols_t.append(phi_linear(xt, h, H))
    for j in range(1, J+1):
        cols_t.append(phi_sine(xt, h, H, j))
    Swt = s0t + np.column_stack(cols_t) @ alpha
    mvt = m_fn(xt, h, a, C, u, Dd)
    beta_tail_num = float(np.trapezoid(4*Swt**2/(mvt*np.pi), xt))
    Xlast = float(xt[-1])
    s1 = np.sinh(h/2)
    L = H - h
    # generous O(1/xi^2) coeff
    s2 = (np.cosh(h/2) + h*np.sinh(h/2)/2)/2 + 10*np.sum(np.abs(alpha))
    m_low = np.log(Xlast/2) - np.log(np.pi) - 2*kap - 4.0
    beta_tail_rem = 4*(s1 + s2/Xlast)**2 / (Xlast * max(m_low, 0.1)) / np.pi
    beta = beta_core + beta_tail_num + beta_tail_rem
    out = dict(
        h=h, H=H, J=J, n_coeff=n, match_edge=match_edge, dxi=dxi, Xi=Xi,
        c_LP=d["c_LP"], c=d.get("c"),
        beta_core=beta_core, beta_tail_num=beta_tail_num, beta_tail_rem=beta_tail_rem,
        beta=beta, two_beta=2*beta,
        beta0_no_ext_core=c0, two_beta0_core=2*c0,
        alpha=alpha.tolist(), names=names,
        edge_value=float(alpha[0]) if match_edge else None,
        sinh_h2=float(np.sinh(h/2)),
        m_min_grid=float(mv.min()),
        seconds=round(time.time()-t0, 2),
    )
    summary = {k: (round(v, 8) if isinstance(v, float) else v)
               for k, v in out.items() if k not in ("alpha", "names")}
    print(json.dumps(summary, indent=1), flush=True)
    print("alpha=", np.round(alpha, 6).tolist(), flush=True)
    if save:
        json.dump(out, open(save, "w"), indent=1)
        print("saved", save, flush=True)
    return out

if __name__ == "__main__":
    cert = sys.argv[1] if len(sys.argv) > 1 else "s22_cert_h0.45.json"
    # sweep mode
    if len(sys.argv) > 2 and sys.argv[2] == "sweep":
        results = []
        for H in (1.0, 1.2, 1.5, 2.0, 2.5, 3.0):
            for J in (8, 16, 24, 32):
                o = design(cert, H=H, J=J, dxi=0.005, Xi=2000.0,
                           save=f"out_extend_H{H}_J{J}.json")
                results.append(dict(H=H, J=J, two_beta=o["two_beta"], beta=o["beta"]))
        print("SWEEP", json.dumps(results, indent=1))
        json.dump(results, open("out_sweep.json", "w"), indent=1)
    else:
        H = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
        J = int(sys.argv[3]) if len(sys.argv) > 3 else 24
        Xi = float(sys.argv[4]) if len(sys.argv) > 4 else 2000.0
        design(cert, H=H, J=J, Xi=Xi, save=f"out_extend_H{H}_J{J}.json")
