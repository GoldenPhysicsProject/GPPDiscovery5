# N.14 NUMERICAL test of Celestial Conjecture C.14.K (team/celestial-cft a7f943a).
# Piecewise-linear xi in H^1_0(-L,L); split V = V_b (linear on each strip) (+) V_s (strip interiors).
# P_{>T}(u,v) = (1/pi) int_T^inf (1/4+t^2)(Omega(t)-Omega(0)) Re(u^ conj v^) dt,  u^ = -(sum_k s_k e^{-itx_k})/t^2
# (s_k = slope jumps). Best gamma = sup |P(u,v)|/sqrt(P(u)P(v)) = ||Pbb^{-1/2} Pbs Pss^{-1/2}||, which is the
# sharp constant in |2P(u,v)| <= gamma[P(u)+P(v)]. Compare with gamma_crit = Omega(T)/(Omega(T)+|Omega(0)|).
import numpy as np, sys
def psi(z, N=12):
    z = np.asarray(z, complex); acc = np.zeros_like(z)
    for k in range(N): acc -= 1/(z + k)
    w = z + N; w2 = 1/(w*w)
    return acc + np.log(w) - 0.5/w - w2*(1/12 - w2*(1/120 - w2*(1/252 - w2*(1/240 - w2/132))))
def sqrtm(A):
    e, U = np.linalg.eigh((A + A.T)/2); return (U*np.sqrt(e)) @ U.T
a = np.log(2.0)
def Om(t): return np.real(psi(0.25 + 0.5j*np.asarray(t))) - np.log(np.pi)
Om0 = float(Om(0.0))

def basis(L, ell, ns, hb):
    # nodes: strip S- = [-L,-L+ell], bulk, strip S+ ; returns lists of (jump positions, jump values)
    sm = np.linspace(-L, -L+ell, ns+1); sp = np.linspace(L-ell, L, ns+1)
    nb = max(2, int(np.ceil((2*L-2*ell)/hb))); bk = np.linspace(-L+ell, L-ell, nb+1)
    Vb, Vs = [], []
    def hat(xl, x0, xr): return ([xl, x0, xr], [1/(x0-xl), -1/(x0-xl)-1/(xr-x0), 1/(xr-x0)])
    for i in range(1, ns): Vs.append(hat(sm[i-1], sm[i], sm[i+1])); Vs.append(hat(sp[i-1], sp[i], sp[i+1]))
    # V_b: bulk hats, endpoints of bulk carry ramps over the whole strip
    for j in range(0, nb+1):
        xl = bk[j-1] if j > 0 else -L
        xr = bk[j+1] if j < nb else L
        Vb.append(hat(xl, bk[j], xr))
    return Vb, Vs

TMAX = 4000.0
def gram(funcs, T, Tmax=None, dt=0.02):
    Tmax = Tmax or TMAX
    n = len(funcs); G = np.zeros((n, n)); t0 = T
    while t0 < Tmax:
        t = np.arange(t0, min(t0+400.0, Tmax), dt) + dt/2
        w = (0.25 + t**2)*(Om(t) - Om0)/t**4*dt/np.pi
        S = np.zeros((n, t.size), complex)
        for i, (xs, ss) in enumerate(funcs):
            for x, s in zip(xs, ss): S[i] += s*np.exp(-1j*t*x)
        G += np.real((S*w) @ S.conj().T); t0 += 400.0
    # tail t > Tmax: oscillatory cross terms average out; keep coincident-node terms (log t/t^2 weight)
    tail = (np.log(Tmax) + 1 - np.log(np.pi) + 0.0)/Tmax/np.pi  # int_Tmax^inf (log(t/2pi)... ) ~ upper-ish, small
    D = {}
    for i, (xs, ss) in enumerate(funcs):
        for x, s in zip(xs, ss): D.setdefault(round(x, 12), {})[i] = D.get(round(x, 12), {}).get(i, 0) + s
    for x, d in D.items():
        for i, si in d.items():
            for j, sj in d.items(): G[i, j] += tail*si*sj
    return G

def gamma(L, ell, ns, hb, T):
    Vb, Vs = basis(L, ell, ns, hb); nb = len(Vb)
    G = gram(Vb + Vs, T)
    A, Bm, C = G[:nb, :nb], G[:nb, nb:], G[nb:, nb:]
    Ai = np.linalg.inv(np.real(sqrtm(A))); Ci = np.linalg.inv(np.real(sqrtm(C)))
    M = Ai @ Bm @ Ci
    # parity split: x -> -x maps the basis to itself (symmetric grids); report full-space sup
    return np.linalg.svd(M, compute_uv=False)[0], np.linalg.eigvalsh(G)[0]

if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "refine":
    TMAX = 24000.0; T = 12.8; gc = Om(T)/(Om(T) + abs(Om0))
    print("sweep 3: refinement at T = 12.8, gamma_crit =", round(float(gc), 4), " Tmax =", TMAX)
    for ell in [0.01, 0.02, 0.025]:
        L = (a + ell)/2
        for ns, hb in [(12, 0.025), (24, 0.0125), (48, 0.00625)]:
            g, _ = gamma(L, ell, ns, hb, T)
            print(f"{ell:5.3f} {T:6.1f} {ns:4d} {hb:7.5f}  {g:8.4f}  {'FAILS' if g >= gc else 'holds'}"); sys.stdout.flush()
    sys.exit()
if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "fine":
    print("Omega(0) =", Om0, " (sweep 2: small ell, refinement, quadrature check)")
    print(" ell     T    ns   hb    Tmax   gamma_best  gamma_crit  verdict")
    L = (a + 0.02)/2
    for Tm in [4000.0, 16000.0]:
        TMAX = Tm; g, _ = gamma(L, 0.02, 12, 0.025, 12.8)
        print(f"0.020   12.8   12  0.025 {Tm:7.0f}  {g:10.4f}  (quadrature check)")
    TMAX = 12000.0
    for ell in [0.005, 0.01, 0.02, 0.03, 0.04]:
        L = (a + ell)/2
        for T in [12.8, 16.0, 20.0]:
            gc = Om(T)/(Om(T) + abs(Om0))
            for ns, hb in [(12, 0.025), (24, 0.0125)]:
                g, _ = gamma(L, ell, ns, hb, T)
                print(f"{ell:5.3f} {T:6.1f} {ns:4d} {hb:6.4f} {TMAX:7.0f}  {g:10.4f}  {gc:10.4f}  {'FAILS' if g >= gc else 'holds'}  ell*T={ell*T:.3f}")
            sys.stdout.flush()
    sys.exit()
print("Omega(0) =", Om0)
print(" ell     T    ns   hb    gamma_best  gamma_crit  minEig(G)  verdict")
for ell in [0.02, 0.05, 0.1, 0.119]:
    L = (a + ell)/2
    for T in [12.8, 20.0, 40.0, 100.0]:
        gc = Om(T)/(Om(T) + abs(Om0))
        for ns, hb in [(6, 0.05), (12, 0.025)]:
            g, me = gamma(L, ell, ns, hb, T)
            print(f"{ell:5.3f} {T:6.1f} {ns:4d} {hb:6.3f}  {g:10.4f}  {gc:10.4f}  {me:9.2e}  {'FAILS' if g >= gc else 'holds'}")
        sys.stdout.flush()
# sweep 3 (run as: python3 n14_c14K_gamma.py refine): one more refinement level at T = 12.8
