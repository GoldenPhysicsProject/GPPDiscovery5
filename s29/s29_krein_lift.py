"""Grok S29 H29a: near-optimal invisible lift from discrete P_L top eigenfunction (Krein/extension dual).
1) Galerkin top eigenpair (lam, v) of P_L (reuse s28_pmax construction).
2) Synthesize f, FT F(t); dual weight w=|F|^2.
3) Invisible characters: per-prime high harmonics + cross-prime |k|_1 <= D with |k·log p|>=2L.
4) Chebyshev LP on dense t-grid (optionally densified at w-peaks): minimise A s.t. A >= comb - m.
5) Also: weighted LS fit of (comb - lam) by invisible m on w-measure, then measure true max.
NUMERICAL. Compare A vs lam_max(P_L) and vs per-prime Sigma kappa_p. RH OPEN."""
import numpy as np, math, json, sys, os
from numpy.polynomial.legendre import legval, leggauss
from scipy.optimize import linprog, lsq_linear

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's28'))
from s28_pmax import mangoldt_terms, phi, Pmat

def primes_below(X):
    ps = []
    for n in range(2, int(X)+1):
        if all(n % p for p in ps if p*p <= n):
            ps.append(n)
    return ps

def comb_terms(L):
    """All von Mangoldt cosine coefficients for log n < +inf (for evaluating comb on R).
    For the lift we only need the full arithmetic comb; envelope uses all n."""
    # Use n up to a generous cutoff for evaluating comb(t); the one-stroke uses Psi ~ log
    # Here for A_lift we need the finite trigonometric part matching s27: n with log n contributing.
    # Per s27/s28: comb contribution on the window is sum_{log n < 2L} Lambda(n)n^{-1/2}*2 cos(t log n)
    # PLUS the invisible lift cancels peaks. The quantity being maximised is
    # sum_{ALL primes powers} c_n * 2? Wait — look at s27 carefully.
    #
    # Actually the one-stroke envelope is about Psi(x)-x etc. The finite A_L in s28 is
    # A_L = sum_{log n < 2L} 2 Lambda(n)/sqrt(n) = max of the visible trigonometric comb
    # (all cosines align at t=0).
    #
    # For lift: we replace visible+invisible per prime. Here for cross-prime LP we take
    # B(t) = sum_{log n < 2L} Lambda(n)n^{-1/2} * 2 cos(t log n)   [visible arithmetic]
    # and m(t) = sum_{invisible omega} a_omega cos(omega t), minimise max(B-m).
    out = []
    for n, lp in mangoldt_terms(L):  # log n < 2L
        out.append((math.log(n), lp / math.sqrt(n)))
    return out  # (freq, coeff) with B = sum 2*c*cos? 
    # mangoldt_terms gives (n, log p); c = lp/sqrt(n).
    # At t=0, sum_n 2*c = A_L in s28_pmax: sum 2*lp/sqrt(n). Yes B(t)=sum 2*c*cos(t*log n).

def invisible_freqs(L, D=3, kmax_mult=8):
    """Invisible frequencies: per-prime k log p >= 2L, and cross |k|_1<=D with |k·log|>=2L."""
    ps = [p for p in primes_below(math.exp(2*L)+1) if math.log(p) < 2*L]
    freqs = []  # list of omega > 0
    # per-prime
    for p in ps:
        lp = math.log(p)
        k0 = int(math.ceil(2*L / lp - 1e-12))
        for k in range(k0, k0 * kmax_mult + 1):
            freqs.append(k * lp)
    # cross-prime multi-indices |k|_1 <= D, not pure axis already covered
    # enumerate small integer vectors
    n = len(ps)
    logs = [math.log(p) for p in ps]
    # D=2 pairs
    if D >= 2:
        for i in range(n):
            for j in range(i+1, n):
                for a in range(-D, D+1):
                    for b in range(-D, D+1):
                        if abs(a)+abs(b) == 0 or abs(a)+abs(b) > D: continue
                        if abs(a)+abs(b) < 2: continue  # pure axis handled above if one is 0
                        if a == 0 or b == 0: continue
                        om = a*logs[i] + b*logs[j]
                        if abs(om) >= 2*L - 1e-12:
                            freqs.append(abs(om))
    if D >= 3 and n >= 3:
        # limited triples: |k|_1=3 with support 3
        for i in range(min(n, 6)):
            for j in range(i+1, min(n, 6)):
                for k in range(j+1, min(n, 6)):
                    for a in (-1,1):
                        for b in (-1,1):
                            for c in (-1,1):
                                om = a*logs[i]+b*logs[j]+c*logs[k]
                                if abs(om) >= 2*L - 1e-12:
                                    freqs.append(abs(om))
    # unique sorted, drop near-duplicates
    freqs = sorted(set(round(f, 10) for f in freqs))
    freqs = [f for f in freqs if f >= 2*L - 1e-9]
    return freqs

def synthesize_ft(L, v, nt=8000, t_max=None):
    """v = Legendre coeffs of top eigenfunction; return t-grid and |F(t)|^2."""
    N = len(v)
    # spatial grid
    xg, wg = leggauss(max(4*N, 200))
    x = L * xg; w = L * wg
    f = np.zeros_like(x)
    for i, ci in enumerate(v):
        f += ci * phi(i, x, L)
    # FT F(t) = int f(x) e^{-itx} dx
    if t_max is None:
        t_max = 80.0
    t = np.linspace(0, t_max, nt)
    # F(t) = sum w_j f_j exp(-i t x_j)
    # vectorized
    phase = np.exp(-1j * np.outer(t, x))
    F = phase @ (f * w)
    wgt = np.abs(F)**2
    return t, wgt, f, x, w

def eval_B(t, terms):
    B = np.zeros_like(t, dtype=float)
    for freq, c in terms:
        B += 2 * c * np.cos(t * freq)
    return B

def chebyshev_lp(terms, omegas, t_grid, w_weight=None, densify_peaks=True, n_peak=2000):
    """min A s.t. A >= B(t) - sum a_j cos(omega_j t) on t_grid (and peak densification)."""
    B = eval_B(t_grid, terms)
    if densify_peaks and w_weight is not None:
        # add extra points at top weight percentiles
        idx = np.argsort(w_weight)[-n_peak:]
        t_extra = t_grid[idx]
        # merge
        t_all = np.unique(np.concatenate([t_grid, t_extra]))
        B = eval_B(t_all, terms)
    else:
        t_all = t_grid
    nA = len(omegas)
    # variables: a[0..nA), A
    # A - B + sum a cos >= 0  =>  -sum a cos - A <= -B
    Cmat = np.cos(np.outer(t_all, omegas))  # (npts, nA)
    Aub = np.hstack([-Cmat, -np.ones((len(t_all), 1))])
    bub = -B
    cobj = np.r_[np.zeros(nA), 1.0]
    # mild bound on coefficients to avoid wild oscillation
    bnds = [(-50, 50)] * nA + [(None, None)]
    r = linprog(cobj, A_ub=Aub, b_ub=bub, bounds=bnds, method='highs')
    if not r.success:
        return None, None, dict(success=False, message=r.message)
    a = r.x[:nA]; A = float(r.x[-1])
    # true max on fine grid
    m = Cmat @ a if len(t_all)==len(B) else np.cos(np.outer(t_grid, omegas)) @ a
    # recompute true max on original t_grid densely
    t_fine = np.linspace(0, t_grid[-1], 100000)
    Bf = eval_B(t_fine, terms)
    mf = np.cos(np.outer(t_fine, omegas)) @ a
    true_max = float(np.max(Bf - mf))
    return a, A, dict(success=True, lp_A=A, true_max=true_max, n_constraints=len(t_all), n_omega=nA)

def weighted_ls_lift(terms, omegas, t_grid, w_weight, lam):
    """Fit m ≈ B - lam in weighted LS on w; return coeffs and true max(B-m)."""
    B = eval_B(t_grid, terms)
    target = B - lam
    Cmat = np.cos(np.outer(t_grid, omegas))
    sw = np.sqrt(np.maximum(w_weight, 0))
    # solve min ||sw * (C a - target)||
    A_ls = Cmat * sw[:, None]
    b_ls = target * sw
    res = lsq_linear(A_ls, b_ls, bounds=(-50, 50))
    a = res.x
    t_fine = np.linspace(0, t_grid[-1], 100000)
    Bf = eval_B(t_fine, terms)
    mf = np.cos(np.outer(t_fine, omegas)) @ a
    return a, float(np.max(Bf - mf)), float(np.sum(w_weight * (B - Cmat @ a)) / np.sum(w_weight))

def run(L, N=100, D=3, t_max=60.0, nt=12000):
    print(f'=== H29a L={L} N={N} D={D} ===', flush=True)
    P = Pmat(L, N)
    ev, V = np.linalg.eigh(P)
    lam = float(ev[-1]); v = V[:, -1]
    print(f'lambda_max(P_L)={lam:.6f}', flush=True)
    terms = comb_terms(L)
    A_raw = sum(2*c for _,c in terms)
    print(f'A_raw (visible t=0)={A_raw:.6f}', flush=True)
    omegas = invisible_freqs(L, D=D)
    print(f'# invisible freqs D={D}: {len(omegas)}', flush=True)
    t, wgt, f, x, w = synthesize_ft(L, v, nt=nt, t_max=t_max)
    # Rayleigh check via FT: int comb |F|^2 / int |F|^2
    B = eval_B(t, terms)
    # trapezoid
    dt = t[1]-t[0]
    # F was computed for t>=0; full line ~ 2*int_0^inf for even |F|^2 of real f... 
    # Actually F(-t)=conj(F(t)) for real f, |F| even. comb even. So int_{-inf}^{inf} = 2 int_0 - endpoint.
    num = 2*np.sum(B * wgt)*dt  # but B=sum 2c cos already; careful with (1/2pi)
    # <f,P f> = (1/2pi) int comb |F|^2. Our Galerkin lam should match.
    # Skip exact FT normalisation; Galerkin lam is authoritative.
    a_lp, A_lp, info_lp = chebyshev_lp(terms, omegas, t, w_weight=wgt)
    print(f'LP: {info_lp}', flush=True)
    a_ls, A_ls, dual_mean = weighted_ls_lift(terms, omegas, t, wgt, lam)
    print(f'WLS: A_true_max={A_ls:.6f}  dual_mean(B-m)={dual_mean:.6f} (target lam={lam:.6f})', flush=True)
    # per-prime reference
    per_prime = {1.0: 3.047608, 1.19: 3.919408, 1.25: 4.582996}.get(round(L,2))
    torus = {1.0: 2.95, 1.19: 3.602}.get(round(L,2))
    out = dict(
        L=L, N=N, D=D, lam_max=lam, A_raw=A_raw,
        n_omega=len(omegas),
        lp=info_lp, lp_A=info_lp.get('true_max') if info_lp.get('success') else None,
        wls_A=A_ls, wls_dual_mean=dual_mean,
        per_prime_A=per_prime, torus_D3_A=torus,
        T1_star=2*math.pi*math.exp(lam),
        T1_lp=2*math.pi*math.exp(info_lp['true_max']) if info_lp.get('success') else None,
        T1_per_prime=2*math.pi*math.exp(per_prime) if per_prime else None,
        gap_lp_to_floor=(info_lp['true_max']-lam) if info_lp.get('success') else None,
        gap_pp_to_floor=(per_prime-lam) if per_prime else None,
        # save coeffs for possible cert
        omegas=omegas[:500],  # truncate for JSON size
        a_lp=list(map(float, a_lp[:500])) if a_lp is not None else None,
        a_ls=list(map(float, a_ls[:500])),
    )
    return out

if __name__ == '__main__':
    results = {}
    for L in (1.0, 1.19):
        for D in (2, 3):
            key = f'L{L}_D{D}'
            try:
                results[key] = run(L, N=80 if L<=1.0 else 100, D=D, t_max=50.0, nt=10000)
            except Exception as e:
                results[key] = dict(error=str(e))
                print('ERROR', key, e)
    json.dump(results, open('h29a_krein_lift.json','w'), indent=2)
    print('WROTE h29a_krein_lift.json')
    for k,v in results.items():
        if 'error' in v: print(k, 'ERR', v['error']); continue
        print(f"{k}: lam={v['lam_max']:.4f} LP_A={v.get('lp_A')} WLS_A={v['wls_A']:.4f} "
              f"per_prime={v['per_prime_A']} torus={v['torus_D3_A']} gap_lp={v.get('gap_lp_to_floor')}")
