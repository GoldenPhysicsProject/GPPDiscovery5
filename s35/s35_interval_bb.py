"""Grok S35 H35a: INDEPENDENT libm-free re-verification of torus upper envelopes
   A >= sup_{theta in T^n} F(theta),  F = sum_vis c_{p,k} cos(k theta_p) - sum_j a_j cos(K_j . theta).
Independent of s28_torus_bb (different code path, different arithmetic model):
 * Box centres lie on the exact grid theta = 2*pi*j/N (j integer, N = 2^M). Phases K.theta = 2*pi*(K.j mod N)/N
   are computed in EXACT integer arithmetic; cos/sin come from a table built in Arb (python-flint, 128 bits),
   each entry the nearest double to an Arb ball of radius < 2^-100. NO libm cos/sin is used anywhere.
 * Visible coefficients c = 2 log p / p^(k/2) are recomputed in Arb from (p,k); invisible a_j are the decimal
   strings a_q stored in the certificate (parsed exactly by Arb, then rounded to nearest double).
 * Float error of the table-lookup sums is bounded a priori by standard IEEE-754 analysis valid for ANY
   summation order / FMA (Higham Thm 3.1 style): |F~ - F| <= 1.01 (n+4) u sum|c|, |g~_i - g_i| <= 1.01 (n+5) u sum|c||k_i|,
   u = 2^-53 (round-to-nearest only; no libm, no directed rounding needed).
 * Box bound (Taylor, Lagrange remainder): F(x) <= F(c) + sum_i |dF_i(c)| d_i + 1/2 sum_il H_il d_i d_l,
   H_il = sum |coef||k_i||k_l| (since |d2/dt2 cos| <= 1 along the segment);  d_i = w_i * 2*pi_up/N >= true half-width.
   Plus an absolute slack 1e-9 covering float evaluation of the <= O(n^2) nonnegative bound terms.
 * Symmetry F(-x)=F(x): theta_1 in [0,pi], others in [-pi,pi]; start from ONE box (different partition from s28).
 * Invisibility of every K_j re-checked independently with mpmath (60 digits) via log of the exact integer ratio.
Usage: python s35_interval_bb.py CERT.json TARGET [M=20]   RH OPEN."""
import sys, json, math, time
import numpy as np
from flint import arb, ctx
import mpmath as mp

def load_cert(fn):
    d = json.load(open(fn)); ps = d['primes']; n = len(ps)
    rows, cs = [], []
    for (i, k, _c) in d['vis']:
        kv = [0] * n; kv[i] = int(k); rows.append(kv); cs.append(2 * arb(ps[i]).log() / arb(ps[i]) ** (arb(int(k)) / 2))
    nvis = len(rows)
    for kv, x in zip(d['K'], d['a_q']):
        rows.append([int(t) for t in kv]); cs.append(-arb(x))
    return d, ps, np.array(rows, np.int64), cs, nvis

def check_invisible(ps, K, nvis, L):
    mp.mp.dps = 60; twoL = 2 * mp.mpf(L); gmin = None
    for r, kv in enumerate(K):
        P = 1; Q = 1
        for p, k in zip(ps, kv):
            if k > 0: P *= p ** int(k)
            elif k < 0: Q *= p ** int(-k)
        u = abs(mp.log(mp.mpf(P) / Q)); g = u - twoL
        if r < nvis: assert g < -mp.mpf('1e-30'), ('visible char not visible', kv)
        else:
            assert g > mp.mpf('1e-30'), ('INVISIBLE CHECK FAILED', kv, g)
            gmin = g if gmin is None or g < gmin else gmin
    return float(gmin)

def table(N):
    ctx.prec = 128
    tc = np.empty(N); ts = np.empty(N); twopi = 2 * arb.pi(); rmax = 0.0
    for m in range(N):
        x = twopi * m / N; c = x.cos(); s = x.sin()
        tc[m] = float(c.mid()); ts[m] = float(s.mid()); rmax = max(rmax, float(c.rad()), float(s.rad()))
    assert rmax < 2.0 ** -100, rmax
    return tc, ts

def run(fn, target, M=20, L=None, chunk=20000, log=print):
    t0 = time.time()
    d, ps, K, cs, nvis = load_cert(fn); n = len(ps); nt = len(cs)
    L = L if L is not None else d['L']
    gmin = check_invisible(ps, K, nvis, L)
    log(f'{fn}: n={n} terms={nt} (vis {nvis}); invisibility (mpmath, independent) OK at L={L}, min gap {gmin:.4e}')
    coef = np.array([float(c.mid()) for c in cs]); absc = np.abs(coef) * (1 + 2.0 ** -50)
    u = 2.0 ** -53; Kabs = np.abs(K).astype(float)
    EF = 1.01 * (nt + 4) * u * absc.sum()
    Eg = 1.01 * (nt + 5) * u * (absc @ Kabs)
    H = (Kabs.T * absc) @ Kabs
    N = 2 ** M; tc, ts = table(N); log(f'  Arb table N=2^{M} built [{time.time()-t0:.0f}s]')
    hstep = 2 * 3.1415926535897936 / N          # 2*pi_up/N >= 2*pi/N (pi_up = next double above pi)
    Kf = K.astype(float)
    # one start box: theta_1 centre N/4 half-width N/4 ; others centre 0 half-width N/2 (grid units)
    J = np.zeros((1, n), np.int32); W = np.full((1, n), N // 2, np.int32); J[0, 0] = N // 4; W[0, 0] = N // 4
    best_lo = -1e300; nproc = 0; it = 0; maxub_last = None; minw = N
    while len(J):
        it += 1; nxtJ = []; nxtW = []; kept = 0; maxub = -1e300
        for s in range(0, len(J), chunk):
            j = J[s:s+chunk]; w = W[s:s+chunk]
            idx = np.mod(j.astype(np.int64) @ K.T, N)
            F = tc[idx] @ coef
            g = -(ts[idx] * coef) @ Kf
            dd = w * hstep
            ub = F + EF + ((np.abs(g) + Eg) * dd).sum(1) + 0.5 * np.einsum('bi,il,bl->b', dd, H, dd) + 1e-9
            best_lo = max(best_lo, float((F - EF).max())); maxub = max(maxub, float(ub.max()))
            keep = ub > target
            if keep.any():
                jk = j[keep]; wk = w[keep]; ddk = dd[keep]; gk = np.abs(g[keep])
                sc = ddk * (gk + ddk @ H)
                sc[wk < 2] = -1                       # cannot split below grid resolution
                ax = np.argmax(sc, 1)
                if (wk[np.arange(len(wk)), ax] < 2).any(): raise RuntimeError('grid too coarse: raise M')
                w2 = wk.copy(); w2[np.arange(len(wk)), ax] //= 2
                off = np.zeros_like(wk); off[np.arange(len(wk)), ax] = w2[np.arange(len(wk)), ax]
                nxtJ += [jk - off, jk + off]; nxtW += [w2, w2]; kept += int(keep.sum())
                minw = min(minw, int(w2.min()))
        nproc += len(J)
        log(f'  it {it}: boxes {len(J)} kept {kept} best_lo {best_lo:.6f} max ub {maxub:.6f} min w {minw}')
        if best_lo > target: return dict(ok=False, reason='sample exceeds target', best_lo=best_lo)
        if kept == 0: break
        J = np.vstack(nxtJ); W = np.vstack(nxtW); nxtJ = nxtW = None
    return dict(ok=True, cert=fn, target=target, L=str(L), best_lower=best_lo, boxes=nproc, iters=it, M=M,
                min_halfwidth_rad=minw * hstep, EF=EF, Eg_max=float(Eg.max()), min_invisible_gap=gmin,
                elapsed_s=time.time() - t0,
                method='s35_interval_bb: exact-integer phases + Arb cos/sin table (no libm), a-priori IEEE sum bounds, 2nd-order Taylor B&B')

if __name__ == '__main__':
    fn = sys.argv[1]; target = float(sys.argv[2]); M = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    L = sys.argv[4] if len(sys.argv) > 4 else None
    r = run(fn, target, M, L)
    print(json.dumps(r, indent=1))
