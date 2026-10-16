"""Grok S30 H30a secondary: truncated Toeplitz / trig-moment dual toward A*=λ_max(P_L).
Visible freqs Ω_vis = {±log n : log n < 2L}. Maximize Σ c_ω Re μ_ω over trig moment
sequences with Toeplitz PSD of size M (Fejér–Riesz: nonnegative trig poly ↔ Toeplitz PSD).
No cvxpy/MOSEK: use projected gradient on Cholesky factor of Toeplitz (or scipy SLSQP
with eig-min≥0). Compare dual value to λ_max. NUMERICAL. RH OPEN."""
from __future__ import annotations
import numpy as np, math, json, sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's28'))
from s28_pmax import mangoldt_terms, Pmat

def visible_coeffs(L):
    """Return list of (ω, c) for cos-terms: B(t)=Σ c cos(ω t), ω=log n, c=2Λ(n)/√n (n>1)."""
    terms = []
    for n, lp in mangoldt_terms(L):
        w = math.log(n)
        if w <= 0 or w >= 2 * L:
            continue
        terms.append((w, 2.0 * lp / math.sqrt(n)))
    return terms

def toeplitz_from_moments(r):
    """Hermitian Toeplitz from moments r[0..M-1] (r[k]=E e^{ikθ} on circle model).
    Here we use a 1D circle proxy (not full Bohr), so this is a RELAXATION / lower bound probe."""
    M = len(r)
    T = np.zeros((M, M), dtype=complex)
    for i in range(M):
        for j in range(M):
            T[i, j] = r[abs(i - j)]
            if i < j:
                T[i, j] = np.conj(T[i, j])
    # enforce Hermitian
    T = 0.5 * (T + T.conj().T)
    return T

def dual_circle_proxy(L, M=48, steps=400, lr=0.05, seed=0):
    """Crude circle Fejér–Riesz proxy: treat frequencies as commensurate on a fine grid.
    Maximize Σ_vis c_k · Re(ρ(ω_k)) where ρ is the density of a nonnegative trig poly
    of degree <M (Toeplitz PSD moments), ∫ρ=1.
    This is NOT the true Bohr dual — only a diagnostic. If even this stays near λ_max,
    it supports strong duality; if it overshoots raw A_L it is meaningless."""
    rs = np.random.default_rng(seed)
    vis = visible_coeffs(L)
    # Represent nonnegative trig poly p(θ)=|q(e^{iθ})|^2 via q coeffs (Fejér–Riesz).
    # q degree m-1; moments r_k = Σ_{j} q_j conj(q_{j+k}).
    m = M
    # parameterize q as complex vector; normalize ‖q‖² = r_0 = 1/(2π) for probability on [0,2π]
    # Actually for torus measure normalized ∫dμ=1, r_0=1.
    q = rs.normal(size=m) + 1j * rs.normal(size=m)
    q /= np.linalg.norm(q)  # => sum |q|^2 = 1 = r_0

    def moments(q):
        # r_k = Σ_j q_j conj(q_{j+k}) for k=0..m-1 (rest 0)
        r = np.zeros(m, dtype=complex)
        for k in range(m):
            r[k] = np.dot(q[:m - k], np.conj(q[k:]))
        return r

    def objective(q):
        # On the CIRCLE we need to map visible ω to angles — use a large period Tper
        # embedding ω ↦ ω * α mod 2π with α chosen so frequencies are distinct.
        # Better diagnostic: evaluate the trig poly density on a fine θ-grid and
        # compute ∫ B(θ) ρ(θ) dθ where ρ = |q(e^{iθ})|^2 / 2π, and B uses
        # cos(ω t) with t identified to θ somehow — this is the weak point.
        # Instead: use discrete Bohr sample of length Ngrid approximating the dual
        # value under Fejér–Riesz densities of order m on ONE random 1-D winding.
        return None

    # --- Practical diagnostic used here ---
    # Sample density ρ on [0,2π), set t = R·θ for large random R (irrational winding),
    # compute E_ρ[B(t)] = Σ c_k Re( mean e^{i ω_k R θ} ρ ).
    # Maximize over q. This probes whether bandlimited Fejér–Riesz densities can
    # approach λ_max along a 1-parameter subgroup (they cannot reach full Bohr dual,
    # so value ≤ A*_Bohr; useful as a lower-bound probe only if ≥ λ_max).
    Ngrid = 4096
    theta = np.linspace(0, 2 * math.pi, Ngrid, endpoint=False)
    # Fixed winding scale: map θ → t so that the smallest ω has many oscillations
    R = 1.0  # t=θ; frequencies ω enter as cos(ω θ) — only valid as a toy if we pretend
    # commensurate. Better: optimize the complex Fourier coeffs at exact visible
    # frequencies subject to being a positive measure's Fourier transform via
    # truncated Toeplitz on a refined frequency lattice.

    # Refined approach: embed all visible ±ω and many invisible frequencies into a
    # frequency lattice with spacing Δ, build Toeplitz of size M, maximize
    # Σ c_j m_{k(ω_j)} with m = moments of PSD Toeplitz, m_0=1.
    # Index map: take Ω_max = 4L, Δ = Ω_max/(M-1), bin each ω to nearest index.
    Omega_max = 4.0 * L
    Delta = Omega_max / (M - 1)
    def bin_w(w):
        return int(round(abs(w) / Delta))

    idx_vis = [(bin_w(w), c) for (w, c) in vis]
    # Fejér–Riesz param: q real for cosine (even measure) — use real q
    q = rs.normal(size=M)
    q /= np.linalg.norm(q)

    def value_and_grad(q):
        # r_k = Σ_j q_j q_{j+k} (real autocorrelation)
        r = np.correlate(q, q, mode='full')[M - 1:]  # r[0..M-1], r[0]=‖q‖²
        # normalize r[0]=1
        nrm = r[0]
        r = r / nrm
        val = 0.0
        # d val / d r_k
        dr = np.zeros(M)
        for k, c in idx_vis:
            if 0 < k < M:
                val += c * r[k]   # cos moment ≈ r_k for even measure (Re μ_k)
                dr[k] += c
            elif k == 0:
                val += c * r[0]
                dr[0] += c
        # chain rule through normalization and correlate
        # r_raw = correlate; r = r_raw/r_raw[0]
        # ∂val/∂q via finite diff for robustness (M small)
        return val, r

    hist = []
    best = -1e9
    for step in range(steps):
        val, r = value_and_grad(q)
        # finite-diff grad
        eps = 1e-5
        g = np.zeros(M)
        for i in range(M):
            q2 = q.copy(); q2[i] += eps
            v2, _ = value_and_grad(q2)
            g[i] = (v2 - val) / eps
        q = q + lr * g
        q /= (np.linalg.norm(q) + 1e-15)
        if val > best:
            best = val
        if step % 50 == 0 or step == steps - 1:
            hist.append((step, val))
            print(f'  step {step}: dual_proxy={val:.6f} best={best:.6f}', flush=True)
    # λ_max for comparison
    P = Pmat(L, 80)
    lam = float(np.linalg.eigvalsh(P)[-1])
    raw_A = sum(c for (_, c) in vis)
    return dict(
        L=L, M=M, dual_proxy=best, lam_max=lam, raw_A=raw_A,
        hist=hist,
        note=('NUMERICAL truncated-Toeplitz circle/lattice PROXY only — not true Bohr dual. '
              'If dual_proxy ≈ λ_max, consistent with strong duality on this relaxation; '
              'if dual_proxy ≪ λ_max, lattice too coarse or method weak. '
              'cvxpy/MOSEK unavailable; skipped full SDP.'),
        status='DIAGNOSTIC',
    )

if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS'] = '1'
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    out = {}
    for L in (1.0, 1.19):
        print(f'\n=== Toeplitz proxy L={L} ===', flush=True)
        out[f'L{L}'] = dual_circle_proxy(L, M=64, steps=300, lr=0.08)
        print(f"RESULT L={L}: proxy={out[f'L{L}']['dual_proxy']:.5f} λ={out[f'L{L}']['lam_max']:.5f} "
              f"rawA={out[f'L{L}']['raw_A']:.5f}", flush=True)
    json.dump(out, open('h30a_toeplitz_proxy.json', 'w'), indent=1)
    print('wrote h30a_toeplitz_proxy.json')
