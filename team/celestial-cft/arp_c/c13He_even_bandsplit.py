"""C.14.6 follow-up (NUMERICAL): failing C.14.K profiles, and the even band-split coercivity C.13.H^e at L0.
R = P + Omega(0) ||.||_H^2 (hat Galerkin on H^1_0(-L0, L0), L0 = (1/2) log 2, uniform mesh, M elements).
(i)  m0^{e/o}(L0) = min R/||.||_H^2 per sector (Galerkin upper bound for the true min; converges from above).
(ii) mu_hi*(mu_lo) = min [R - mu_lo ||.^{<T}||_H^2] / ||.^{>T}||_H^2 per sector: C.13.H^e needs mu_hi* > c2 for some mu_lo > 0.
     Galerkin values are upper bounds on the true infimum (NUMERICAL, not certified).
(iii) profiles of the extremal pair (xi_b, xi_s) for C.14.K.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
import numpy as np
from scipy.linalg import toeplitz, eigh
from c14K_gamma_test import toeplitz_P, gamma, gstar, Om, Om0, a
c2 = a/np.sqrt(2)
def H_cols(M, h, T, nGL=600):
    hcol = np.zeros(M+1); hcol[0] = 2/h + h/6; hcol[1] = -1/h + h/24
    x, wt = np.polynomial.legendre.leggauss(nGL); t = 0.5*T*(x+1); wt = 0.5*T*wt
    hat2 = (h*np.sinc(t*h/(2*np.pi))**2)**2
    low = np.array([np.sum(wt*(0.25+t*t)*hat2*np.cos(t*k*h))/np.pi for k in range(M+1)])
    return hcol, low
def sector_basis(n, sector):
    R = np.eye(n)[::-1]; F = np.eye(n) + (1 if sector == 'even' else -1)*R
    q, r = np.linalg.qr(F); keep = np.abs(np.diag(r)) > 1e-10; return q[:, keep]
def mats(M, T):
    L = a/2; h = 2*L/M; n = M-1
    col, colhi = toeplitz_P(M, h, T); hcol, hlow = H_cols(M, h, T)
    Rm = toeplitz(col[:n]) + Om0*toeplitz(hcol[:n])
    H = toeplitz(hcol[:n]); Hlo = toeplitz(hlow[:n]); Hhi = H - Hlo
    return Rm, H, Hlo, Hhi
if __name__ == "__main__":
    print("c2 =", c2)
    print("== (i) m0(L0) = min R/||.||_H^2 (Galerkin, M elements)")
    for M in (40, 80, 160, 320):
        Rm, H, Hlo, Hhi = mats(M, 14.0)
        row = []
        for s in ('even', 'odd'):
            B = sector_basis(M-1, s)
            w, V = eigh(B.T@Rm@B, B.T@H@B)
            v = B@V[:, 0]; lof = (v@Hlo@v)/(v@H@v)
            row.append(f"{s}: m0={w[0]:+.6f} (low-band H-fraction of minimiser at T=14: {lof:.4f})")
        print(f"  M={M}: " + "; ".join(row))
    print("== (ii) C.13.H^e: mu_hi*(mu_lo) per sector (need > c2 = %.4f)" % c2)
    for T in (12.875, 14.0, 16.0):
        for M in (80, 160, 320):
            Rm, H, Hlo, Hhi = mats(M, T)
            out = []
            for s in ('even', 'odd'):
                B = sector_basis(M-1, s)
                vals = []
                for mlo in (0.0, 0.1, 0.25, 0.3):
                    w, V = eigh(B.T@(Rm - mlo*Hlo)@B, B.T@Hhi@B)
                    vals.append(f"{mlo}:{w[0]:+.4f}")
                out.append(s + " " + " ".join(vals))
            print(f"  T={T} M={M}: " + " | ".join(out))
    print("== (iii) C.14.K extremal profiles (T=14, n_s=16)")
    for s in ('even', 'odd'):
        for lt in (0.0301, 0.0603, 0.1155):
            M0 = round(a*2/lt)+2; ns = 16; M = M0*ns//2
            r = gamma(M, ns, 14.0, s, retvec=True)
            xb, xs = r['xi_b'], r['xi_s']; n = M-1; hh = r['h']; L = r['L']
            x = -L + hh*np.arange(1, M)
            col, colhi = toeplitz_P(M, hh, 14.0); G = toeplitz(colhi[:n])
            w = xb + xs; ratio = (w@G@w)/(xb@G@xb + xs@G@xs)
            alpha = xb[ns-1]; jb = np.argmax(np.abs(xb)); js = np.argmax(np.abs(xs))
            # bulk profile: fraction of xi_b sup norm within one strip-width of the strip edge
            near = np.abs(x) > L - 2*r['ell']
            print(f"  {s} ell={r['ell']:.4f} gamma={r['gamma']:.4f} (gamma*={gstar(14.0):.4f}); P_hi(xi_b+xi_s)/[P_hi(xi_b)+P_hi(xi_s)]={ratio:.4f}"
                  f"; edge value alpha={alpha:+.3f}, max|xi_b|={abs(xb[jb]):.3f} at x={x[jb]:+.4f}; "
                  f"max|xi_s|={abs(xs[js]):.3f} at x={x[js]:+.4f} (strip [{-L:.4f},{-L+r['ell']:.4f}]); "
                  f"share of sum xi_b^2 within 2 ell of +-L: {np.sum(xb[near]**2)/np.sum(xb**2):.3f}")
