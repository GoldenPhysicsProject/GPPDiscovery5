# N.17: high-band H-energy fraction of the existing N.15/N.16 near-null vectors at 2L = log 3 (N = 320, P1 hats; the same
# vectors as n16_nearnull.py, re-formed, no new scan). f(T) = (1/pi) int_T^inf (1/4+t^2)|xi^|^2 dt / ||xi||_H^2. NUMERICAL.
# Test of Principal's (H3^sigma)(T, eta, mu > 0) on this vector: needs rho >= c0*(1+eta) f(T).
import numpy as np
src = open(__file__.replace('n17_hiband.py', 'n15_P2_margins.py')).read().split('c2, c3, c4 =')[0]; exec(src)
c2 = np.log(2)/np.sqrt(2); L = np.log(3)/2; Nn = 320; c0 = 0.90093
x, h, GR, GH = grams(L, Nn); n = x.size; A = GR - 2*c2*Bgram(x, h, L, np.log(2)); Rf = np.eye(n)[::-1]
def xihat(c, t): return (np.exp(1j*np.outer(t, x)) @ c)*h*np.sinc(t*h/2/np.pi)**2
for sgn, nm in ((+1, "even"), (-1, "odd")):
    P = (np.eye(n) + sgn*Rf)/2; U, s, _ = np.linalg.svd(P); B = U[:, s > 0.5]
    Lc = np.linalg.cholesky(B.T @ GH @ B); Ci = np.linalg.inv(Lc)
    e, V = np.linalg.eigh(Ci @ (B.T @ A @ B) @ Ci.T)
    cvec = B @ (Ci.T @ V[:, 0]); Hn = cvec @ GH @ cvec
    out = []
    for T, eta in [(16.0, 0.0373), (19.35, 0.248), (38.1, 1.0)]:
        dt = 0.001; t = np.arange(dt/2, T, dt)
        low = np.sum((0.25 + t**2)*np.abs(xihat(cvec, t))**2)*dt/np.pi/Hn
        f = 1 - low; need = c0*(1 + eta)*f
        out.append(f"T={T:5.2f}: f={f:.3e}, c0*(1+eta)f={need:.3e} {'<= rho ok' if need <= e[0] else '> rho: (H3) FAILS on this vector'}")
    print(f"{nm} #0  rho={e[0]:.4e}  " + " | ".join(out), flush=True)
