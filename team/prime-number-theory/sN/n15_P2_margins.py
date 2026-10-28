# N.15 NUMERICAL (Rayleigh-Ritz, P1 hats): H-norm margins of the Weil forms on windows (-L, L).
#   Q^{(<=X)}(xi) = R(xi) - sum_{n<=X} 2 c_n B_n(xi),  c_n = Lambda(n)/sqrt(n),
#   B_n(xi) = int xi'(x) xi'(x+log n) + 1/4 xi(x) xi(x+log n) dx,   R = (1/2pi) int (1/4+t^2) Omega |xi^|^2,
#   norm ||xi||_H^2 = int (xi'^2 + xi^2/4).   Outputs inf Q/||.||_H^2 per parity (Ritz = UPPER bounds of the inf,
#   valid as true upper bounds only up to quadrature error; NUMERICAL, not certified).
import numpy as np, sys
def psi(z, N=12):
    z = np.asarray(z, complex); acc = np.zeros_like(z)
    for k in range(N): acc -= 1/(z + k)
    w = z + N; w2 = 1/(w*w)
    return acc + np.log(w) - 0.5/w - w2*(1/12 - w2*(1/120 - w2*(1/252 - w2*(1/240 - w2/132))))
Om = lambda t: np.real(psi(0.25 + 0.5j*t)) - np.log(np.pi)
lg = np.log
def grams(L, N):
    h = 2*L/N; x = -L + h*np.arange(1, N); n = x.size
    Tmax = 300.0/h; dt = min(0.02, 0.6/(2*L + 1/h*0 + 2*L)); rowR = np.zeros(n); dd = np.arange(n)*h
    for t0 in np.arange(0, Tmax, 2e5*dt):        # chunked to bound memory
        t = np.arange(t0 + dt/2, min(t0 + 2e5*dt, Tmax), dt)
        wR = (0.25 + t**2)*Om(t)*h*h*np.sinc(t*h/2/np.pi)**4*dt/np.pi
        rowR += np.cos(np.outer(dd, t)) @ wR
    tl = (16/h**2)*(np.log(Tmax/(2*np.pi)) + 1)/Tmax/np.pi
    for k, m in [(0, 3/8), (1, -1/4), (2, 1/16)]: rowR[k] += tl*m
    I = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    GR = rowR[I]
    GH = np.where(I == 0, 2/h + h/6, np.where(I == 1, -1/h + h/24, 0.0))
    return x, h, GR, GH
def hatval(x, h, y):   # values and slopes of all hats at points y
    r = (y[None, :] - x[:, None])/h
    v = np.clip(1 - np.abs(r), 0, None)
    s = np.where((r > -1) & (r < 0), 1/h, np.where((r >= 0) & (r < 1), -1/h, 0.0))
    return v, s
def Bgram(x, h, L, sft, M=40):
    # exact for P1: integrate on each sub-cell between merged breakpoints with M-point Gauss
    if sft >= 2*L: return np.zeros((x.size, x.size))
    bp = np.unique(np.concatenate([np.r_[-L, x, L], np.r_[-L, x, L] - sft]))
    bp = bp[(bp >= -L) & (bp <= L - sft)]
    g, w = np.polynomial.legendre.leggauss(4)
    ys, ws = [], []
    for lo, hi in zip(bp[:-1], bp[1:]):
        ys.append((lo + hi)/2 + (hi - lo)/2*g); ws.append((hi - lo)/2*w)
    y = np.concatenate(ys); wt = np.concatenate(ws)
    v0, s0 = hatval(x, h, y); v1, s1 = hatval(x, h, y + sft)
    M_ = (s0*wt) @ s1.T + 0.25*(v0*wt) @ v1.T
    return (M_ + M_.T)/2
def margins(L, N, primes):
    x, h, GR, GH = grams(L, N); n = x.size; A0 = GR.copy()
    for (p, c) in primes: A0 -= 2*c*Bgram(x, h, L, lg(p))
    Rf = np.eye(n)[::-1]; out = []
    for sgn in (+1, -1):
        P = (np.eye(n) + sgn*Rf)/2; U, s, _ = np.linalg.svd(P); B = U[:, s > 0.5]
        Ci = np.linalg.inv(np.linalg.cholesky(B.T @ GH @ B))
        out.append(np.linalg.eigvalsh(Ci @ (B.T @ A0 @ B) @ Ci.T)[0])
    return out
c2, c3, c4 = lg(2)/np.sqrt(2), lg(3)/np.sqrt(3), lg(2)/2
print("check: Bgram unit test, c2 =", c2, " c3 =", c3)
print("[a] strip-supported margin r(ell) = inf R/||.||_H^2 on H^1_0 of an interval of length ell (min over parity)")
for ell in [0.03, 0.06, 0.119, 0.2, 0.3, lg(1.5), 0.55, lg(2)]:
    r = margins(ell/2, 40, []); print(f"  ell={ell:.4f}  even {r[0]: .4f}  odd {r[1]: .4f}   (compare c2={c2:.4f}, c2+c3={c2+c3:.4f})")
print("[b] P2 form Q^(<=2) = R - 2c2 B2 (Ritz N=80 / N=160): parity margins in H^1_1/4")
for lab, twoL in [("log2", lg(2)), ("log2+0.028", lg(2)+0.028), ("log2+0.1", lg(2)+0.1), ("log2+0.3", lg(2)+0.3), ("log3", lg(3)), ("log3.5", lg(3.5)), ("log4", lg(4))]:
    r1 = margins(twoL/2, 80, [(2, c2)]); r2 = margins(twoL/2, 160, [(2, c2)])
    print(f"  2L={lab:11s} even {r1[0]: .5f} / {r2[0]: .5f}   odd {r1[1]: .5f} / {r2[1]: .5f}")
print("[c] P3 form R - 2c2 B2 - 2c3 B3 (+ 2c4 B4 at 2L > log 4) (N=80 / N=160)")
for lab, twoL in [("log3", lg(3)), ("log3+0.028", lg(3)+0.028), ("log3+0.1", lg(3)+0.1), ("log3.5", lg(3.5)), ("log4", lg(4)), ("log4+0.1", lg(4)+0.1)]:
    pr = [(2, c2), (3, c3)] + ([(4, c4)] if twoL > lg(4) else [])
    r1 = margins(twoL/2, 80, pr); r2 = margins(twoL/2, 160, pr)
    print(f"  2L={lab:11s} even {r1[0]: .5f} / {r2[0]: .5f}   odd {r1[1]: .5f} / {r2[1]: .5f}")
