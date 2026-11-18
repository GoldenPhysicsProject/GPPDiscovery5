"""Locate zeros of the Epstein zeta Z_Q, Q = x^2 + 5y^2 (disc -20, h = 2), for test data only.
Z_Q(s) = zeta_K(s) + L(s,chi_-4) L(s,chi_5)  with zeta_K = zeta(s) L(s,chi_-20)  (forms x^2+5y^2 and 2x^2+2xy+3y^2; w=2).
Checked against direct representation numbers r_Q(n). On-line zeros: sign changes of the real function
Lambda_Q(1/2+it) e^{...} (FE: (sqrt20/2pi)^s Gamma(s) Z_Q(s) symmetric). Off-line: |Z_Q| minima on a grid, refined by findroot.
"""
import mpmath as mp, numpy as np, json
mp.mp.dps = 15
def kron(D, n):
    return int(mp.re(mp.mpf(0))) if False else None
def chi_table(D, q):
    # Kronecker symbol (D/n) for n mod q, D fundamental discriminant
    from sympy import jacobi_symbol
    from sympy.ntheory import factorint
    def k(n):
        if n == 0: return 0
        r = 1; m = n
        for p, e in factorint(n).items():
            if p == 2:
                if D % 2 == 0: return 0
                r *= (1 if D % 8 in (1, 7) else -1)**e
            else:
                r *= jacobi_symbol(D % p, p)**e
        return r
    return [k(n) for n in range(q)]
c20 = chi_table(-20, 20); c4 = chi_table(-4, 4); c5 = chi_table(5, 5)
Z = lambda s: mp.zeta(s)*mp.dirichlet(s, c20) + mp.dirichlet(s, c4)*mp.dirichlet(s, c5)
# coefficient check
N = 200
r = np.zeros(N+1, int)
for x in range(-15, 16):
    for y in range(-7, 8):
        v = x*x+5*y*y
        if 0 < v <= N: r[v] += 1
chk = abs(Z(mp.mpf(3)) - mp.fsum(r[n]*mp.mpf(n)**-3 for n in range(1, N+1)))
A = mp.sqrt(20)/(2*mp.pi)
Lam = lambda s: A**s*mp.gamma(s)*Z(s)
out = dict(coef_check_at_s3=float(chk), FE_check=float(abs(Lam(mp.mpf('0.3')+5j) - Lam(mp.mpf('0.7')-5j))))
# on-line zeros: Lambda(1/2+it) is real
ts = np.arange(0.2, 60.0, 0.05); vals = [float(mp.re(Lam(0.5+1j*t))) for t in ts]
online = []
for i in range(len(ts)-1):
    if vals[i] == 0 or vals[i]*vals[i+1] < 0:
        online.append(float(mp.findroot(lambda t: mp.re(Lam(0.5+1j*t)), (ts[i], ts[i+1]), solver='anderson')))
# off-line search
off = []
for sig in np.arange(0.55, 1.5, 0.05):
    row = [abs(Z(sig+1j*t)) for t in np.arange(0.5, 60, 0.25)]
    for j in range(1, len(row)-1):
        if row[j] < row[j-1] and row[j] < row[j+1] and row[j] < 0.15:
            try:
                z = mp.findroot(Z, mp.mpc(sig, 0.5+0.25*j))
                if abs(Z(z)) < 1e-10 and z.real > 0.5+1e-6 and 0 < z.imag < 60 and all(abs(z - mp.mpc(*o)) > 1e-6 for o in off):
                    off.append((float(z.real), float(z.imag)))
            except Exception:
                pass
out.update(online=online, n_online=len(online), offline=sorted(off, key=lambda z: z[1]))
json.dump(out, open('epstein_zeros.json', 'w'), indent=1); print(json.dumps(out)[:2000])
