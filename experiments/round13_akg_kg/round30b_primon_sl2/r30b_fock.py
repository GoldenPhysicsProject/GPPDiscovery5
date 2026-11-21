"""Round 30-b, part A: the primon Fock space, its positive Gram, the sl(2) no-go for L0 = log N, the per-prime
metaplectic (weight 1/2, 'spinor') oscillator, the Planck <-> Matsubara (KMS) duality that supplies the Weyl element
(shadow Delta <-> 1-Delta = functional equation), and the norm of the shadow intertwiner chi(Delta).
No zeta zeros are used anywhere. Labels: see ROUND30B.md."""
import json, math, numpy as np, mpmath as mp
mp.mp.dps = 30
out = {}

# A1. Fock Gram <v_s|v_s'> = sum_n n^{-conj(s)-s'} = zeta(conj(s)+s') on Re s > 1/2: PSD (it is a Gram matrix).
rng = np.random.default_rng(30)
def fock_gram(points, coeff=None):
    k = len(points); G = mp.matrix(k, k)
    for i, a in enumerate(points):
        for j, b in enumerate(points):
            w = mp.conj(a) + b
            G[i, j] = mp.zeta(w) if coeff is None else coeff(w)
    return G
rows = []
for sig_lo in (0.51, 0.6, 0.8):
    pts = [mp.mpc(sig_lo + 0.4*rng.random(), 30*(rng.random() - 0.5)) for _ in range(12)]
    G = fock_gram(pts)
    ev = mp.eighe(G)[0] if hasattr(mp, 'eighe') else mp.eigh(G)[0]
    rows.append(dict(sigma_min=sig_lo, n=12, min_eig=float(min(ev)), max_eig=float(max(ev))))
out['A1_fock_gram_PSD'] = rows

# A1b. boundary Toeplitz kernel zeta(1+2e+i(t'-t)) on a t-grid: PSD for every e > 0 (Bochner: FT of the positive
# measure sum_n n^{-1-2e} delta_{log n}); the diagonal blows up like 1/(2e) -- the pole = Plancherel delta.
tb = []
for e in (0.25, 0.1, 0.03, 0.01):
    ts = np.linspace(0, 40, 41); K = np.array([[complex(mp.zeta(1 + 2*e + 1j*(b - a))) for b in ts] for a in ts])
    w = np.linalg.eigvalsh(K); tb.append(dict(eps=e, min_eig=float(w[0]), diag=float(K[0, 0].real)))
out['A1b_boundary_toeplitz'] = tb

# A2. no-go: [L0, L1] = -L1 with L0|n> = log n |n>  =>  <m|L1|n>(log m - log n + 1) = 0  =>  L1 = 0 unless m/n = e^{-1}.
# e is irrational, so L1 = L-1 = 0, then [L1, L-1] = 2 L0 forces L0 = 0: contradiction on |n>, n >= 2.  (PROVED)
# With L0 = log N / log q the ladder needs m/n = q: possible only for q a positive rational; for an sl(2)
# lowest-weight ladder the step must be a FIXED multiplicative factor, and the free-boson structure then lives on
# single-prime chains n, p n, p^2 n, ...  Illustration: closest |log(m/n) - 1| over 1 <= n < m <= 2000.
best = (9, None)
for n in range(1, 2001):
    m = round(n*math.e)
    for mm in (m - 1, m, m + 1):
        if mm > n and mm <= 6000:
            d = abs(math.log(mm/n) - 1)
            if d < best[0]: best = (d, (mm, n))
out['A2_closest_ratio_to_e'] = dict(min_gap=best[0], pair=best[1], note='never 0: e irrational (PROVED)')

# A3. per-prime metaplectic oscillator: K0 = (N + 1/2)/2, K+ = (a^dag)^2/2, K- = a^2/2 on the truncated Fock space
# of one prime (occupation N = v_p(n)).  Check sl(2) relations and the Casimir K0^2 - K0 - K+K- = k(k-1) on the even
# (k = 1/4) and odd (k = 3/4) sectors.  k = 1/4 and 3/4 are SHADOWS (k + k' = 1) with the SAME Casimir -3/16.
M = 60; a = np.diag(np.sqrt(np.arange(1, M)), 1); ad = a.T; Nn = ad @ a
K0 = (Nn + 0.5*np.eye(M))/2; Kp = ad @ ad/2; Km = a @ a/2
c1 = np.abs((K0 @ Kp - Kp @ K0) - Kp)[:M-4, :M-4].max()
c2 = np.abs((Kp @ Km - Km @ Kp) + 2*K0)[:M-4, :M-4].max()
Cas = K0 @ K0 - K0 - Kp @ Km
cas_even = np.diag(Cas)[0:M-4:2]; cas_odd = np.diag(Cas)[1:M-4:2]
out['A3_metaplectic'] = dict(comm_K0_Kp_err=float(c1), comm_Kp_Km_err=float(c2),
    casimir_even=[float(cas_even.min()), float(cas_even.max())], casimir_odd=[float(cas_odd.min()), float(cas_odd.max())],
    lowest_weights=[float(K0[0, 0]), float(K0[1, 1])], shadow_sum=float(K0[0, 0] + K0[1, 1]))
# Euler factor = trace of the oscillator semigroup on the prime's Fock space, split into the two shadow sectors:
# (1-x)^{-1} = (1-x^2)^{-1} [even, k=1/4: squares] + x(1-x^2)^{-1} [odd, k=3/4]; globally zeta(b) = zeta(2b) * sum |mu(n)| n^{-b}.
N = 200000; mu2 = np.ones(N + 1, bool); mu2[0] = False
for p in range(2, int(N**0.5) + 1): mu2[p*p::p*p] = False
sf = float(np.sum(np.arange(1, N + 1, dtype=float)[mu2[1:]]**-2.3))
out['A3_square_split'] = dict(beta=2.3, zeta=float(mp.zeta(2.3)), zeta2b_times_squarefree=float(mp.zeta(4.6))*sf)

# A4. Planck <-> Matsubara (KMS): 1/(e^w-1) = 1/w - 1/2 + sum_{k>=1} 2w/(w^2+(2 pi k)^2).
# Mellin: for -1 < Re s < 0, int_0^inf w^{s-1}[1/(e^w-1) - 1/w + 1/2] dw = Gamma(s) zeta(s)   (thermal side, n >= 1)
#                              = pi (2 pi)^{s-1} zeta(1-s) / cos(pi s/2)                         (Matsubara side, k >= 1)
# i.e. the functional equation (Weyl element / shadow) = the KMS duality between occupation n and Matsubara k.
a4 = []
for s in (mp.mpc(-0.5, 0), mp.mpc(-0.3, 4), mp.mpc(-0.7, 11)):
    def br(w):   # 1/(e^w-1) - 1/w + 1/2, Bernoulli series near 0 to avoid cancellation
        if w < 0.2: return sum(mp.bernoulli(2*k)*w**(2*k - 1)/mp.factorial(2*k) for k in range(1, 14))
        return 1/mp.expm1(w) - 1/w + mp.mpf(1)/2
    f = lambda w: w**(s - 1)*br(w)
    # [0,1] numerically; on [1,inf) the elementary pieces are integrated in closed form (Re s < 0)
    I0 = sum(mp.bernoulli(2*k)/mp.factorial(2*k)*mp.mpf('0.2')**(s + 2*k - 1)/(s + 2*k - 1) for k in range(1, 14))
    I = I0 + mp.quad(f, mp.linspace(0.2, 1, 9)) + mp.quad(lambda w: w**(s - 1)/mp.expm1(w), [1, 10, 60]) + 1/(s - 1) - 1/(2*s)
    a4.append(dict(s=str(s), mellin=complex(I), thermal=complex(mp.gamma(s)*mp.zeta(s)),
                   matsubara=complex(mp.pi*(2*mp.pi)**(s - 1)*mp.zeta(1 - s)/mp.cos(mp.pi*s/2))))
out['A4_planck_matsubara'] = a4

# A5. shadow intertwiner zeta(D) = chi(D) zeta(1-D): |chi(1/2+it)| = 1 (unitary on the principal series),
# |chi(s+it)| ~ (t/2pi)^{1/2-s} off it; for F_theta, chi_th(D) = chi(D-th) chi(D+th) is ALSO unitary on the line.
chi = lambda s: 2**s*mp.pi**(s - 1)*mp.sin(mp.pi*s/2)*mp.gamma(1 - s)
a5 = []
for t in (14, 50, 200):
    for sg in (0.5, 0.6, 0.75):
        s = mp.mpc(sg, t); a5.append(dict(t=t, sigma=sg, abs_chi=float(abs(chi(s))),
              asym=float((t/(2*mp.pi))**(0.5 - sg)), abs_chi_FTH=float(abs(chi(s - 0.25)*chi(s + 0.25)))))
out['A5_intertwiner'] = a5
json.dump(out, open('r30b_fock.json', 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
