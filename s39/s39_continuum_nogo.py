#!/usr/bin/env python3
"""H39a.4: Numeric witness for continuum-only no-go.
L > log2/2, f0>=0 smooth with g0(log2)>0; f_T=e^{i T x} f0;
find T where active log-phases nearly align → P_L stays O(1)>0 while H_L,|Pole| shrink.
Label NUMERICAL (Kronecker step not proved here).
"""
import json, math
import numpy as np
from numpy.linalg import norm

L = 1.0  # > log2/2 ≈ 0.34657; support 2L=2 > log2, so n=2 active
assert 2*L > math.log(2)

# Smooth nonnegative bump: f0(x) = exp(-1/(1-x^2)) on (-1,1), else 0  (C^∞_c)
# Use slightly smaller support inside (-L,L): scale to (-0.95,0.95)
def make_f0_grid(L, N=2048):
    xs = np.linspace(-L, L, N, endpoint=False)
    dx = xs[1]-xs[0]
    r = xs / (0.95*L)
    f = np.zeros_like(xs)
    mask = np.abs(r) < 1
    f[mask] = np.exp(-1.0 / (1.0 - r[mask]**2))
    # normalize L2 ~ 1
    f /= np.sqrt(np.sum(f**2)*dx)
    return xs, f, dx

xs, f0, dx = make_f0_grid(L, N=4096)

def autocorr_g(f, dx, u):
    # g(u)=∫ f(x) conj(f(x+u)) dx  via shift on grid
    shift = int(round(u/dx))
    if shift >= len(f) or shift <= -len(f):
        return 0j
    if shift >= 0:
        return np.dot(f[:len(f)-shift], np.conj(f[shift:])) * dx
    else:
        s = -shift
        return np.dot(f[s:], np.conj(f[:len(f)-s])) * dx

# Active n with log n < 2L
ns = [n for n in range(2, int(math.exp(2*L))+2) if math.log(n) < 2*L]
# von Mangoldt
def Lambda(n):
    # n = p^k
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            if m == 1:
                return math.log(p)
            return 0.0
    return 0.0

weights = [(n, Lambda(n)/math.sqrt(n), math.log(n)) for n in ns if Lambda(n)>0]
print('active', [(n,w,u) for n,w,u in weights])

g0_log2 = autocorr_g(f0, dx, math.log(2))
print('g0(log2)=', g0_log2)

def P_L(f, dx):
    s = 0.0
    for n,w,u in weights:
        g = autocorr_g(f, dx, u)
        s += 2*w*np.real(g)
    return float(s)

def ell_pm(f, xs, dx, sign):
    return np.sum(f * np.exp(sign*xs/2)) * dx

def Pole(f, xs, dx):
    ep = ell_pm(f, xs, dx, +1)
    em = ell_pm(f, xs, dx, -1)
    return float(2*np.real(ep*np.conj(em)))

def H_L_quad(f, xs, dx):
    # ⟨f, H f⟩ = ∬ conj(f(x)) e^{-|x-y|/2} f(y) dx dy
    # = ∫ |û(t)|^2 / (t²+1/4) dt / 2π with û = FT of f extended by 0
    # Discrete: use double sum on grid (O(N²)) — N=4096 is 16M, ok
    # Better: FFT of kernel
    N = len(f)
    # kernel on differences: for lag k, K = exp(-|k|dx/2)
    # Toeplitz quadratic form via embed circulant FFT
    lags = np.arange(-N+1, N)
    K = np.exp(-np.abs(lags)*dx/2)
    # build first row of circulant embedding length 2N-1 → pad to 2N
    M = 2*N
    circ = np.zeros(M)
    circ[0] = K[N-1]  # lag 0 at index N-1 in lags
    for k in range(1, N):
        circ[k] = K[N-1+k]       # lag +k
        circ[M-k] = K[N-1-k]     # lag -k
    Kf = np.fft.fft(circ)
    fp = np.zeros(M, dtype=complex)
    fp[:N] = f
    Hf = np.fft.ifft(Kf * np.fft.fft(fp))[:N]
    return float(np.real(np.vdot(f, Hf) * dx * dx / dx))  # vdot with dx: careful
    # Actually discrete: sum_i sum_j conj(f_i) K_ij f_j * dx*dx
    # Toeplitz via FFT: (dx)* Hf_i ≈ sum_j K_{i-j} f_j dx, then sum conj(f) (Hf) dx
    # Fix:

def H_L_quad_fixed(f, xs, dx):
    N = len(f)
    M = 2*N
    circ = np.zeros(M)
    circ[0] = 1.0  # lag 0
    for k in range(1, N):
        val = math.exp(-k*dx/2)
        circ[k] = val
        circ[M-k] = val
    Kf = np.fft.fft(circ)
    fp = np.zeros(M, dtype=complex); fp[:N] = f
    conv = np.fft.ifft(Kf * np.fft.fft(fp))[:N]  # sum_j K_{i-j} f_j
    return float(np.real(np.vdot(f, conv)) * dx * dx)

H0 = H_L_quad_fixed(f0, xs, dx)
P0 = P_L(f0, dx)
Pole0 = Pole(f0, xs, dx)
print(f'f0: P={P0:.6f} H={H0:.6f} Pole={Pole0:.6f} g(log2)={g0_log2}')

# Diophantine: find T minimizing sum_p w_p |e^{-i T log p} - 1|^2 ≈ maximizing Re phases
# Use continued search: sample T or use LLL-ish grid on torus
logs = np.array([u for n,w,u in weights])
ws = np.array([w for n,w,u in weights])

def phase_score(T):
    # want e^{-i T u} ≈ 1 for all u; score = sum w (1-cos(T u))  small is good
    return float(np.sum(ws * (1 - np.cos(T * logs))))

# Coarse then refine: Kronecker — search T in [0, Tmax]
rng = np.random.default_rng(42)
best = []
Tmax = 5000.0
# random + grid
candidates = list(np.linspace(0, Tmax, 2000)) + list(rng.uniform(0, Tmax, 5000))
for T in candidates:
    sc = phase_score(T)
    best.append((sc, T))
best.sort()
# refine top 20 with local golden search
def local_min(T0, span=2.0, steps=80):
    Ts = np.linspace(T0-span, T0+span, steps)
    return min((phase_score(T), T) for T in Ts)

refined = []
for sc,T in best[:30]:
    refined.append(local_min(T, span=5.0, steps=200))
    refined.append(local_min(T, span=0.5, steps=200))
refined.sort()

witnesses = []
for sc, T in refined[:15]:
    fT = f0 * np.exp(1j * T * xs)
    PT = P_L(fT, dx)
    HT = H_L_quad_fixed(fT, xs, dx)
    PoleT = Pole(fT, xs, dx)
    # C_L = Pole - H  (Astra identity Pointwise ⇒ quadratic)
    CT = PoleT - HT  # should equal continuum K quadratic; for modulated →0
    ratio = PT / max(HT, 1e-30)
    witnesses.append({
        'T': float(T), 'phase_score': float(sc),
        'P_L': PT, 'H_L': HT, 'Pole': PoleT, 'C_L_via_Pole_minus_H': CT,
        'P_over_H': ratio,
        'P0': P0, 'H0': H0, 'Pole0': Pole0,
    })
    print(f'T={T:.4f} sc={sc:.4e} P={PT:.5f} H={HT:.5e} Pole={PoleT:.5e} P/H={ratio:.3e}')

# Also try larger Tmax for better alignment
print('--- larger T search ---')
Tmax2 = 50000.0
candidates2 = list(rng.uniform(0, Tmax2, 20000))
best2 = sorted((phase_score(T), T) for T in candidates2)[:20]
for sc,T in best2:
    sc2,T2 = local_min(T, span=2.0, steps=100)
    fT = f0 * np.exp(1j * T2 * xs)
    PT = P_L(fT, dx); HT = H_L_quad_fixed(fT, xs, dx); PoleT = Pole(fT, xs, dx)
    ratio = PT / max(HT, 1e-30)
    witnesses.append({'T': float(T2), 'phase_score': float(sc2), 'P_L': PT, 'H_L': HT,
                      'Pole': PoleT, 'P_over_H': ratio, 'search': 'Tmax50k'})
    print(f'T={T2:.4f} sc={sc2:.4e} P={PT:.5f} H={HT:.5e} Pole={PoleT:.5e} P/H={ratio:.3e}')

# Best witness: max P/H with P>0.5*P0
pos = [w for w in witnesses if w['P_L'] > 0.3*P0]
pos.sort(key=lambda w: -w['P_over_H'])
best_w = pos[0] if pos else max(witnesses, key=lambda w: w['P_over_H'])

out = {
    'L': L, 'n_active': len(weights), 'weights': [(n,w,u) for n,w,u in weights],
    'g0_log2_real': float(np.real(g0_log2)),
    'f0_P_H_Pole': {'P': P0, 'H': H0, 'Pole': Pole0},
    'best_witness': best_w,
    'n_witnesses': len(witnesses),
    'verdict': 'NUMERICAL',
    'note': 'Kronecker recurrence not proved here; phases found by random+local search on torus',
    'all_window_domination': 'OPEN',
    'RH': 'OPEN',
}
json.dump(out, open('/workspace/GPPDiscovery5/s39/h39a_continuum_nogo.json','w'), indent=2)
print('BEST', json.dumps(best_w, indent=2))
print('wrote h39a_continuum_nogo.json')
