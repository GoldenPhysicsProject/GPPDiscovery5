"""Grok S30 H30a: Fejér–Riesz / autocorrelation construction toward A ≈ λ_max(P_L).
Reconstruct top eigenfunction f* of P_L, form dual weight |F*(t)|², fit invisible
trigonometric m by (1) |F*|²-weighted WLS to B−λ and (2) torus LP cutting planes
seeded by F*-peak phases (not dense D=5 planes). Verify almost-periodic max via
≥1e5 Bohr samples + local ascent. Compare A_true vs λ_max vs per-prime vs S28 torus.
NUMERICAL. Success: A ≤ 1.05·λ_max with verified max, OR sharp DEAD END. RH OPEN."""
from __future__ import annotations
import numpy as np, math, json, sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 's28'))
from s28_pmax import Pmat, mangoldt_terms, phi
from s28_torus_lp import setup, chars, Bval, Fval, Fgrad, local_max
from scipy.optimize import linprog, nnls
from numpy.polynomial.legendre import leggauss

REF = {
    1.0:  dict(lam=1.9305635972086237, per_prime=3.047608, s28_torus=2.95,  s29=2.936),
    1.19: dict(lam=2.6612525390517066, per_prime=3.919408, s28_torus=3.602, s29=3.602),
}

def top_eigen(L, N=120):
    P = Pmat(L, N)
    ev, V = np.linalg.eigh(P)
    return float(ev[-1]), V[:, -1].copy(), P

def reconstruct_f(v, L, nq=None):
    N = len(v)
    nq = nq or max(4 * N, 256)
    xg, wg = leggauss(nq)
    x = L * xg; w = L * wg
    f = np.zeros_like(x)
    for i, ci in enumerate(v):
        f += ci * phi(i, x, L)
    return x, w, f

def autocorrelation(x, w, f, L, nu=400):
    """g(u) = ∫ f(x) f(x+u) dx for u in (-2L, 2L) via GL on overlap."""
    us = np.linspace(-2 * L + 1e-9, 2 * L - 1e-9, nu)
    g = np.zeros(nu)
    # denser GL for overlap integrals
    xg, wg = leggauss(max(len(x), 200))
    for j, u in enumerate(us):
        a, b = max(-L, -L - u), min(L, L - u)
        if b <= a:
            continue
        xx = 0.5 * (b - a) * xg + 0.5 * (b + a)
        ww = 0.5 * (b - a) * wg
        # interpolate f at xx and xx+u (f supported on [-L,L], samples on x)
        f1 = np.interp(xx, x, f, left=0.0, right=0.0)
        f2 = np.interp(xx + u, x, f, left=0.0, right=0.0)
        g[j] = float(np.dot(ww, f1 * f2))
    return us, g

def FT_weight(x, w, f, t_max=120.0, nt=24000):
    t = np.linspace(0.0, t_max, nt)
    # F(t) = ∫ f(x) e^{-itx} dx
    F = np.exp(-1j * np.outer(t, x)) @ (f * w)
    wgt = np.abs(F) ** 2
    return t, F, wgt

def phase_seeds(t, F, wgt, lp, n_seeds=8000):
    """Map F*-peak times to Bohr torus; also keep complex phase for Fejér–Riesz cue."""
    idx = np.argsort(wgt)[-n_seeds:]
    t_star = t[idx]
    th = np.mod(np.outer(t_star, lp) + math.pi, 2 * math.pi) - math.pi
    phases = np.angle(F[idx])
    return th, phases, t_star, wgt[idx]

def wls_fit_invisible(th_w, w_w, vis, K, lam, ridge=1e-8):
    """Weighted LS: m(θ) ≈ B(θ)−λ on weight samples; m = Σ a_j cos(K_j·θ)."""
    B = Bval(th_w, vis)
    target = B - lam
    C = np.cos(th_w @ K.T)  # (npts, nK)
    # weighted normal equations: (C^T W C + ridge I) a = C^T W target
    sw = np.sqrt(np.maximum(w_w, 0.0))
    Cw = C * sw[:, None]
    tw = target * sw
    G = Cw.T @ Cw + ridge * np.eye(C.shape[1])
    rhs = Cw.T @ tw
    try:
        a = np.linalg.solve(G, rhs)
    except np.linalg.LinAlgError:
        a = np.linalg.lstsq(G, rhs, rcond=None)[0]
    resid = target - C @ a
    A_on_support = float(np.max(B - C @ a))  # max on the weight sample only
    dual_mean = float(np.average(B - C @ a, weights=np.maximum(w_w, 0)))
    return a, A_on_support, dual_mean, float(np.sqrt(np.average(resid ** 2, weights=np.maximum(w_w, 0))))

def torus_true_max(vis, K, a, n_rand=200000, n_local=8000, seed=0):
    """Almost-periodic max estimate: random Bohr samples + local ascent."""
    rs = np.random.default_rng(seed)
    n = len(vis) and (max(i for (i, _, _) in vis) + 1)
    # better: dim from K
    n = K.shape[1]
    pts = rs.uniform(-math.pi, math.pi, (n_rand, n))
    Fc = Fval(pts, vis, K, a)
    top_idx = np.argsort(-Fc)[:n_local]
    climbed = local_max(pts[top_idx], vis, K, a, iters=250, lr=0.03)
    Fc2 = Fval(climbed, vis, K, a)
    return float(max(Fc.max(), Fc2.max())), climbed[np.argmax(Fc2)]

def seeded_lp(L, D, lam, th_seed, vis, K, npts0=12000, rounds=14, amax=25.0,
              seed=2, gap=5e-4, a_init=None):
    """Cutting-plane LP for min A s.t. B−m ≤ A, seeded by F* phases + optional WLS warmstart."""
    rs = np.random.default_rng(seed)
    n = K.shape[1]
    pts = np.vstack([th_seed, rs.uniform(-math.pi, math.pi, (npts0, n))])
    a = np.zeros(len(K)) if a_init is None else np.asarray(a_init, dtype=float).copy()
    if len(a) != len(K):
        a = np.zeros(len(K))
    hist = []
    s = smax = 1e9
    for r in range(rounds):
        C = np.cos(pts @ K.T)
        Bp = Bval(pts, vis)
        Aub = np.hstack([-C, -np.ones((len(pts), 1))])
        res = linprog(np.r_[np.zeros(len(K)), 1.0], A_ub=Aub, b_ub=-Bp,
                      bounds=[(-amax, amax)] * len(K) + [(None, None)], method='highs')
        if not res.success:
            print(f'  LP fail round {r}: {res.message}', flush=True)
            break
        a = res.x[:-1]; s = float(res.x[-1])
        cand = rs.uniform(-math.pi, math.pi, (60000, n))
        cand = np.vstack([cand, th_seed + rs.normal(0, 0.08, th_seed.shape),
                          pts[np.argsort(-(Bp - C @ a))[:2500]]])
        cand = local_max(cand, vis, K, a)
        Fc = Fval(cand, vis, K, a); smax = float(Fc.max())
        hist.append((s, smax))
        print(f'  round {r}: LP s={s:.6f} true>={smax:.6f} pts={len(pts)} chars={len(K)}', flush=True)
        if smax - s < gap:
            break
        pts = np.vstack([pts, cand[np.argsort(-Fc)[:3500]]])
    return a, s, smax, hist

def run_one(L, D, N=120, t_max=150.0, nt=30000, n_seeds=8000, budget_s=1200):
    t0 = time.time()
    print(f'\n=== H30a Fejér–Riesz L={L} D={D} N={N} ===', flush=True)
    lam, v, _ = top_eigen(L, N)
    # refine λ with N already at 120; also report N=80 for continuity with S28/S29
    ref = REF.get(round(L, 2), {})
    print(f'λ_max(P_L)={lam:.8f} (N={N}); ref S28 N120={ref.get("lam")}', flush=True)
    x, w, f = reconstruct_f(v, L)
    us, g = autocorrelation(x, w, f, L)
    g0 = float(np.interp(0.0, us, g))
    print(f'‖f‖²≈g(0)={g0:.6f}  (should be ~1 for orthonormal basis combo)', flush=True)
    t, F, wgt = FT_weight(x, w, f, t_max=t_max, nt=nt)
    # Parseval check: (1/2π)∫|F|² ≈ ‖f‖² — discrete proxy
    dt = t[1] - t[0]
    parseval_proxy = float(np.sum(wgt) * dt / (2 * math.pi))
    print(f'Parseval proxy ∫|F|² dt/2π ≈ {parseval_proxy:.4f} (partial t∈[0,{t_max}])', flush=True)

    ps, lp, vis = setup(L)
    th_seed, phases, t_star, w_seed = phase_seeds(t, F, wgt, lp, n_seeds=n_seeds)
    kmax_axis = 6 * (max(k for (_, k, _) in vis) + 1)
    K = chars(lp, L, D, kmax_axis)
    print(f'primes={ps} vis={len(vis)} invisible chars={len(K)} seeds={len(th_seed)}', flush=True)

    # (1) WLS Fejér–Riesz fit on |F*|²-weighted torus samples
    # reweight seeds by wgt
    a_wls, A_wls_supp, dual_mean, rms = wls_fit_invisible(th_seed, w_seed, vis, K, lam)
    A_wls_true, _ = torus_true_max(vis, K, a_wls, n_rand=150000, n_local=6000, seed=11)
    print(f'WLS: dual_mean={dual_mean:.6f} (≈λ?)  A_on_seed={A_wls_supp:.6f}  A_true≈{A_wls_true:.6f}  rms={rms:.4f}', flush=True)

    # (2) LP cutting planes seeded by F* phases, warm-started from WLS
    remain = budget_s - (time.time() - t0)
    rounds = 12 if remain > 600 else (8 if remain > 300 else 5)
    a_lp, s_lp, smax_lp, hist = seeded_lp(
        L, D, lam, th_seed, vis, K, rounds=rounds, a_init=a_wls, gap=5e-4)
    # independent true-max audit ≥1e5 Bohr + ascent
    A_true, _ = torus_true_max(vis, K, a_lp, n_rand=250000, n_local=10000, seed=42)
    A_final = max(smax_lp, A_true)

    per_prime = ref.get('per_prime') or float(sum(c for (_, _, c) in vis))
    s28 = ref.get('s28_torus')
    ratio = A_final / lam
    status = 'NEW' if ratio <= 1.05 else ('PROGRESS' if (s28 and A_final < s28 - 0.02) else 'NO_PROGRESS')
    out = dict(
        L=L, D=D, N=N, lam_max=lam, g0=g0, parseval_proxy=parseval_proxy,
        n_primes=len(ps), primes=ps, n_vis=len(vis), n_chars=len(K),
        wls=dict(dual_mean=dual_mean, A_on_seed=A_wls_supp, A_true=A_wls_true, rms=rms),
        lp=dict(A_lp=s_lp, A_max_est=smax_lp, hist=hist),
        A_true=A_final, A_true_indep=A_true,
        per_prime=per_prime, s28_torus=s28, s29_ref=ref.get('s29'),
        ratio_to_floor=ratio, gap_to_floor=A_final - lam,
        status=status,
        T1_star=2 * math.pi * math.exp(lam),
        T1_A=2 * math.pi * math.exp(A_final),
        T1_per_prime=2 * math.pi * math.exp(per_prime),
        elapsed_s=time.time() - t0,
        note='NUMERICAL Fejér–Riesz / |F*|²-seeded torus; not B&B certified.',
    )
    # keep coefficients if useful for later cert
    out['a_lp'] = a_lp.tolist()
    out['K'] = K.tolist()
    out['vis'] = [(int(i), int(k), float(c)) for (i, k, c) in vis]
    print(f'RESULT L={L} D={D}: A_true≈{A_final:.5f} λ={lam:.5f} ratio={ratio:.3f} '
          f'per_prime={per_prime:.4f} S28={s28} status={status} [{out["elapsed_s"]:.0f}s]', flush=True)
    return out

if __name__ == '__main__':
    os.environ['OMP_NUM_THREADS'] = '1'
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    results = {}
    # Primary order: L=1 D=2,3,4 then L=1.19 D=2,3,4 — stop branch if clearly no better than S28
    plan = [(1.0, 2), (1.0, 3), (1.0, 4), (1.19, 2), (1.19, 3), (1.19, 4)]
    for L, D in plan:
        key = f'L{L}_D{D}'
        try:
            results[key] = run_one(L, D, N=120 if L <= 1.0 else 120, budget_s=900)
        except Exception as e:
            results[key] = dict(error=str(e), status='ERROR')
            print(f'ERROR {key}: {e}', flush=True)
        # early-exit a D-branch if still stuck at S28 level after D=3
        json.dump(results, open('h30a_fejer_riesz.json', 'w'), indent=1)
        if L == 1.0 and D == 3:
            r = results[key]
            if r.get('status') == 'NO_PROGRESS' and r.get('ratio_to_floor', 9) > 1.4:
                print('L=1 D=3 no progress toward floor; still run D=4 once then move to L=1.19', flush=True)
    # summary verdict
    best = min((v.get('A_true', 9e9), k, v) for k, v in results.items() if 'A_true' in v)
    summary = dict(
        best_key=best[1], best_A=best[0],
        best_ratio=best[2].get('ratio_to_floor'),
        all_status={k: v.get('status') for k, v in results.items()},
        conclusion='',
    )
    # Decide DEAD END vs progress
    any_new = any(v.get('status') == 'NEW' for v in results.values())
    any_prog = any(v.get('status') == 'PROGRESS' for v in results.values())
    if any_new:
        summary['conclusion'] = 'NEW: A within 5% of λ_max (strong-duality near-hit).'
    elif any_prog:
        summary['conclusion'] = 'PROGRESS: beat S28 torus slightly; floor still open.'
    else:
        summary['conclusion'] = (
            'DEAD END / no progress: Fejér–Riesz |F*|² seeding + WLS+LP at D≤4 does not '
            'close gap to λ_max beyond S28/S29 torus (~2.9 at L=1, ~3.6 at 1.19). '
            'Strong duality H28 remains HYPOTHESIS. Do not circle.'
        )
        summary['label'] = 'DEAD END'
    results['_summary'] = summary
    json.dump(results, open('h30a_fejer_riesz.json', 'w'), indent=1)
    print('\nSUMMARY:', summary['conclusion'], flush=True)
