# N.16 (mechanism). (i) Evaluate the derived near-null model  m^sigma(L) ~ w * (1 - lambda_{2N+2+sigma}(L*gamma_{N+1}))
#     (prolate leakage after killing the first N zeros), against the existing N.15 Ritz numbers (no new scan);
# (ii) the ONE check: N = 320 Ritz near-null vectors at 2L = log 3, their xi^ at the first zeros, and the
#     explicit-formula zero sum 2*sum_{gamma>0}(1/4+gamma^2)|xi^(gamma)|^2 / ||xi||_H^2 vs the Ritz value.
# Zero ordinates: mpmath.zetazero, n <= 400 (NUMERICAL input; all lie below the verified height 3e12, Platt-Trudgian 2021).
import numpy as np, json, sys
Z = np.array(json.load(open(__file__.replace('n16_nearnull.py', 'n16_zeros400.json'))))
def lam(c, M=160):
    x, w = np.polynomial.legendre.leggauss(M); D = x[:, None] - x[None, :]
    K = np.where(D == 0, c/np.pi, np.sin(c*D)/(np.pi*np.where(D == 0, 1, D)))
    S = np.sqrt(w); return np.sort(np.linalg.eigvalsh(S[:, None]*K*S[None, :]))[::-1]
print("(i) model: margin^sigma(L) ~ min_N log(g_{N+1}/2pi) * (1 - lambda_{2N+2+sigma}(c = L g_{N+1}))")
N15 = {"log2": (0.0730, 0.5493), "log2+0.1": (0.0168, 0.2452), "log2+0.3": (4.5e-4, 0.0142), "log3": (2.5e-4, 1.34e-3)}
for lab, twoL in [("log2", np.log(2)), ("log2+0.1", np.log(2)+0.1), ("log2+0.3", np.log(2)+0.3), ("log3", np.log(3)), ("log3.5", np.log(3.5)), ("log4", np.log(4))]:
    L = twoL/2; mu = np.exp(twoL); row = []
    for sg in (0, 1):
        best = (np.inf, None)
        for N in range(0, 8):
            g = Z[N]; c = L*g; k = 2*N + 2 + sg; lk = lam(c)
            val = np.log(g/(2*np.pi))*max(1 - lk[k], 1e-16)
            if val < best[0]: best = (val, N, c, k)
        row.append(best)
    D = L*2*np.pi*mu/np.pi - ((mu)*(twoL - 1) + 7/8) - 1
    ref = N15.get(lab, ("-", "-"))
    print(f"  2L={lab:9s} mu={mu:.3f} D(t*)={D: .3f}  even ~ {row[0][0]:.2e} (kill N={row[0][1]}, c={row[0][2]:.2f}, k={row[0][3]})"
          f"   odd ~ {row[1][0]:.2e} (N={row[1][1]}, c={row[1][2]:.2f}, k={row[1][3]})   N.15 Ritz <= {ref}")
if len(sys.argv) > 1 and sys.argv[1] == "nocheck": sys.exit()
print("(ii) ONE check at 2L = log 3, N = 320 (Ritz vectors of R - 2c2 B2, the full form there)")
src = open(__file__.replace('n16_nearnull.py', 'n15_P2_margins.py')).read().split('c2, c3, c4 =')[0]; exec(src)
c2 = np.log(2)/np.sqrt(2); L = np.log(3)/2; Nn = 320
x, h, GR, GH = grams(L, Nn); n = x.size; A = GR - 2*c2*Bgram(x, h, L, np.log(2)); Rf = np.eye(n)[::-1]
def xihat(c, t): return (np.exp(1j*np.outer(t, x)) @ c)*h*np.sinc(t*h/2/np.pi)**2
for sgn, nm in ((+1, "even"), (-1, "odd")):
    P = (np.eye(n) + sgn*Rf)/2; U, s, _ = np.linalg.svd(P); B = U[:, s > 0.5]
    Lc = np.linalg.cholesky(B.T @ GH @ B); Ci = np.linalg.inv(Lc)
    e, V = np.linalg.eigh(Ci @ (B.T @ A @ B) @ Ci.T)
    for j in range(2):
        cvec = B @ (Ci.T @ V[:, j]); Hn = cvec @ GH @ cvec          # = 1
        wz = (0.25 + Z**2)*np.abs(xihat(cvec, Z))**2/Hn
        Tm = 300/h; dt = 0.02; tail = 0.0
        for t0 in np.arange(Z[-1], Tm, 1e5*dt):      # smoothed density (1/2pi)log(t/2pi) per sign: NUMERICAL heuristic
            tt = np.arange(t0, min(t0 + 1e5*dt, Tm), dt)
            tail += np.sum((1/np.pi)*np.log(tt/(2*np.pi))*(0.25 + tt**2)*np.abs(xihat(cvec, tt))**2)*dt/Hn
        tg = np.linspace(0, 40, 4001); spec = (0.25 + tg**2)*np.abs(xihat(cvec, tg))**2/Hn
        print(f"  {nm} eig#{j}: Ritz value {e[j]:.4e} | zero sum: 2*sum_(n<=400) = {2*wz.sum():.4e}, density tail above g_400 = {tail:.4e},"
              f" total = {2*wz.sum()+tail:.4e}")
        print(f"     (1/4+g^2)|xi^(g)|^2 at g1..g5 = " + ", ".join(f"{v:.2e}" for v in wz[:5]) +
              f";  max over t<=40 of (1/4+t^2)|xi^|^2 = {spec.max():.3f} at t = {tg[spec.argmax()]:.2f}")
