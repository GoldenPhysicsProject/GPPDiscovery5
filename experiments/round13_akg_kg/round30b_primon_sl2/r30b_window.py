"""Round 30-b, part B: the Mazac-Paulos 'window' mechanism transplanted to the principal series Re Delta = 1/2.
MP (AFB II, arXiv:1811.10646, abstract and sec. 5): crossing + unitarity (a_O >= 0) + Regge boundedness force at least
one primary in every window [2D_phi+2n-2, 2D_phi+2n+2] for large n.  Transplant: the completed function
Lambda_M(1/2+it) is REAL for real t whenever Lambda_M(D) = Lambda_M(1-D) with real coefficients (shadow + reality).
Window functional, positive kernel k_T(t) = cos^2(pi (t-T)/(2H)) on |t-T| < H:
     R_M(T,H) = | int Z_M k_T | / int |Z_M| k_T ,    Z_M(t) = Lambda_M(1/2+it)/|Q^s G(s)|.
If Z_M kept one sign on the window, R = 1.  R < 1  =>  Z_M changes sign  =>  an ODD-order zero ON the line in the window
(PROVED logic; the values are NUMERICAL, mpmath).  No zeros are inputs: only L-values on the line.
Models (as in rounds 13-29): ZETA; FTH = xi(s-th)xi(s+th), th = 1/4; EPS = Epstein x^2+5y^2;
DH = Davenport-Heilbronn (chi mod 5, kappa); BEUR has no functional equation -> no real Z (structural)."""
import json, math, sys, numpy as np, mpmath as mp
mp.mp.dps = 20
TH = mp.mpf('0.25'); KAP = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
chi5 = [0, 1, 1j, -1j, -1]; chi5b = [0, 1, -1j, 1j, -1]
chim4 = [0, 1, 0, -1]; chi5r = [0, 1, -1, -1, 1]
chim20 = [chim4[n % 4]*chi5r[n % 5] for n in range(20)]
def Lfun(s, c): return mp.dirichlet(s, c)
def Z(name, t):
    s = mp.mpc(0.5, t)
    if name == 'ZETA':
        g = mp.pi**(-s/2)*mp.gamma(s/2); v = g*mp.zeta(s)
    elif name == 'FTH':
        g1 = mp.pi**(-(s - TH)/2)*mp.gamma((s - TH)/2); g2 = mp.pi**(-(s + TH)/2)*mp.gamma((s + TH)/2)
        g = g1*g2; v = g*mp.zeta(s - TH)*mp.zeta(s + TH)
    elif name == 'DH':
        g = (5/mp.pi)**(s/2)*mp.gamma((s + 1)/2)
        v = g*((1 - 1j*KAP)*Lfun(s, chi5) + (1 + 1j*KAP)*Lfun(s, chi5b))/2
    elif name == 'EPS':
        g = (mp.sqrt(20)/(2*mp.pi))**s*mp.gamma(s)
        v = g*(mp.zeta(s)*Lfun(s, chim20) + Lfun(s, chim4)*Lfun(s, chi5r))/2
    return v/abs(g)

def eps_coeff_check(N=300):
    r = np.zeros(N + 1)
    for x in range(-20, 21):
        for y in range(-10, 11):
            q = x*x + 5*y*y
            if 0 < q <= N: r[q] += 1
    a = r/2
    def conv(f, g):
        h = np.zeros(N + 1)
        for d in range(1, N + 1):
            h[d::d] += f[d]*g[1:N//d + 1]
        return h
    one = np.ones(N + 1); one[0] = 0
    c20 = np.array([0] + [chim20[n % 20] for n in range(1, N + 1)], float)
    c4 = np.array([0] + [chim4[n % 4] for n in range(1, N + 1)], float)
    c5 = np.array([0] + [chi5r[n % 5] for n in range(1, N + 1)], float)
    b = (conv(one, c20) + conv(c4, c5))/2
    return float(np.abs(a[1:] - b[1:]).max())

if __name__ == '__main__':
    out = dict(eps_decomposition_max_err=eps_coeff_check())
    T0, T1, step = 10.0, float(sys.argv[1]) if len(sys.argv) > 1 else 300.0, 0.05
    ts = np.arange(T0, T1 + step/2, step)
    names = ['ZETA', 'FTH', 'EPS', 'DH']
    vals = {}
    for nm in names:
        z = np.array([complex(Z(nm, t)) for t in ts])
        vals[nm] = z
        out.setdefault('max_imag_over_abs', {})[nm] = float(np.max(np.abs(z.imag))/np.max(np.abs(z)))
        print(nm, 'done', out['max_imag_over_abs'][nm], flush=True)
    res = {}
    for H in (0.5, 1.0, 2.0, 4.0):
        res[str(H)] = {}
        for nm in names:
            z = vals[nm].real; cert = []; Rs = []
            for T in np.arange(T0 + H, T1 - H + 1e-9, H):
                m = np.abs(ts - T) < H; k = np.cos(np.pi*(ts[m] - T)/(2*H))**2
                R = abs(np.trapezoid(z[m]*k, ts[m]))/np.trapezoid(np.abs(z[m])*k, ts[m])
                Rs.append(R); cert.append(R < 1 - 1e-6)
            Rs = np.array(Rs); cert = np.array(cert)
            # sign changes on the sampling grid (diagnostic only)
            sc = int(np.sum(np.sign(z[1:]) != np.sign(z[:-1])))
            res[str(H)][nm] = dict(windows=int(len(Rs)), certified=int(cert.sum()), frac=float(cert.mean()),
                                   R_max=float(Rs.max()), R_median=float(np.median(Rs)), grid_sign_changes=sc,
                                   first_uncertified_T=float(np.arange(T0 + H, T1 - H + 1e-9, H)[~cert][0]) if (~cert).any() else None,
                                   last_uncertified_T=float(np.arange(T0 + H, T1 - H + 1e-9, H)[~cert][-1]) if (~cert).any() else None)
            print(H, nm, res[str(H)][nm], flush=True)
    out['windows'] = res; out['T_range'] = [T0, T1]
    out['min_Z_FTH'] = float(vals['FTH'].real.min())
    np.savez('r30b_window_Z.npz', ts=ts, **{k: v.real for k, v in vals.items()})
    json.dump(out, open('r30b_window.json', 'w'), indent=1)
