"""Grok S82 H80c step 2: speed levers for the Ozaki-sliced cell products (support 2.84 even, L=71/50, N=2376, prec 512).
Same Arb inputs as s81/sliced_cell.py (S71 mori_cert code), one cell, nch consecutive quadrature chunks of 64 nodes.
Levers: (1) batching chunks so the inner dim is k = 64*m (fixed-point scale common to the whole batch);
(2) Cp2 = F diag(w) F^T = G G^T with G = F diag(sqrt w) (Gauss-Legendre weights 2hh*w_k > 0, asserted in Arb),
    via dsyrk (diagonal slice pairs) + dsyr2k (pairs i<j on each anti-diagonal), upper triangle only;
(3) amortised carry: every slice product (< 2^53, exact) is added into int64 base-2^b limb arrays
    (headroom 2^10 products per limb), one carry per accumulated batch-group instead of per chunk.
Exactness: every BLAS call has output entries bounded by (inner dim)*2^(2b) <= 2^53 with integer operands -> exact.
Radius (entrywise, rigorous, float64 x (1+1e-9)):
  dgemm  A W^T : (rA+uA)(|W|+rW)^T + |A|(rW+uW)^T
  syrk   G G^T : |G|E^T + E|G|^T + E E^T,  E = rG + uG
Validation vs Arb (same inputs): |mid_arb - sliced| <= rad_arb + rad_sliced for every entry (exact dyadic diff)."""
import sys, os, time, json, math, argparse
import numpy as np
from scipy.linalg import blas
from threadpoolctl import threadpool_limits
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's71'))
import mori_cert as mc
from flint import arb, arb_mat, ctx

def pick_b(P, k):
    for b in range(26, 1, -1):
        S = math.ceil((P + 1) / b)
        if 2 * b + math.ceil(math.log2(S * k)) <= 53: return b, S

def fixed(M, P):
    """arb_mat -> (python int mids at scale 2^(P-e), e, rad float array, u)"""
    n, k = M.nrows(), M.ncols()
    mx = 0.0; mids = [[None] * k for _ in range(n)]; rad = np.zeros((n, k))
    for i in range(n):
        for j in range(k):
            x = M[i, j]; m = x.mid(); mids[i][j] = m.man_exp(); rad[i, j] = float(x.rad())
            mx = max(mx, abs(float(m)))
    e = math.ceil(math.log2(mx)); sh = P - e
    Z = [[0] * k for _ in range(n)]
    for i in range(n):
        for j in range(k):
            m, ex = mids[i][j]; m = int(m); t = int(ex) + sh
            Z[i][j] = (m << t) if t >= 0 else ((m + (1 << (-t - 1))) >> (-t))
    return Z, e, rad * (1 + 1e-12), 2.0 ** (-(sh + 1))

def slice_words(Z, P, b, S):
    """exact slices: |x| << pad as 512-bit little-endian words, then b-bit windows; returns list of S float arrays
    (slice 0 = most significant), lo exponent (= -pad), |Z| as float, sign"""
    n, k = len(Z), len(Z[0]); pad = S * b - (P + 1); assert 0 <= pad and S * b <= 512
    buf = bytearray(); sg = np.empty((n, k))
    for i in range(n):
        for j in range(k):
            x = Z[i][j]; sg[i, j] = -1.0 if x < 0 else 1.0
            buf += (abs(x) << pad).to_bytes(64, 'little')
    Wd = np.frombuffer(bytes(buf), dtype=np.uint64).reshape(n, k, 8)
    out = []; mask = np.uint64((1 << b) - 1)
    for s in range(S):
        lo = (S - 1 - s) * b; wi, off = lo // 64, lo % 64
        v = Wd[:, :, wi] >> np.uint64(off)
        if off + b > 64: v = v | (Wd[:, :, wi + 1] << np.uint64(64 - off))
        out.append(sg * (v & mask).astype(np.float64))
    absf = np.array([[float(abs(x)) for x in r] for r in Z])
    return out, -pad, absf

def acc_gemm(D, As, Bs, S, k):
    """D[2S-2-t] += sum_{i+j=t} As[i] Bs[j]^T  (one dgemm per anti-diagonal, concatenated)"""
    AC = np.concatenate(As, axis=1)
    for t in range(2 * S - 1):
        i0, i1 = max(0, t - S + 1), min(S - 1, t)
        Bb = np.concatenate([Bs[t - i] for i in range(i0, i1 + 1)], axis=1)
        D[2 * S - 2 - t] += (AC[:, i0 * k:(i1 + 1) * k] @ Bb.T).astype(np.int64)

def acc_syrk(D, Gs, S, k):
    """upper triangle of D[2S-2-t] += sum_{i+j=t} Gs[i] Gs[j]^T via dsyr2k (i<j) + dsyrk (i=j)"""
    GF = [np.asfortranarray(g) for g in Gs]
    for t in range(2 * S - 1):
        i0, i1 = max(0, t - S + 1), min(S - 1, t)
        prs = [(i, t - i) for i in range(i0, i1 + 1) if i < t - i]
        if prs:
            A = np.asfortranarray(np.concatenate([Gs[i] for i, _ in prs], axis=1))
            B = np.asfortranarray(np.concatenate([Gs[j] for _, j in prs], axis=1))
            D[2 * S - 2 - t] += blas.dsyr2k(1.0, A, B).astype(np.int64)
        if t % 2 == 0:
            D[2 * S - 2 - t] += blas.dsyrk(1.0, GF[t // 2]).astype(np.int64)

def carry(D, b):
    m = (1 << b) - 1
    for d in range(D.shape[0] - 1):
        cy = D[d] >> b; D[d] &= m; D[d + 1] += cy

def limbs_to_int(D, b, rows, cols):
    """exact integers from carried limbs (all but top limb in [0,2^b), top signed)"""
    nl = D.shape[0]; nb = (nl - 1) * b; nw = (nb + 63) // 64 + 1
    Wd = np.zeros(D.shape[1:] + (nw,), dtype=np.uint64)
    for d in range(nl - 1):
        p = d * b; wi, off = p // 64, p % 64
        v = D[d].astype(np.uint64)
        Wd[..., wi] |= v << np.uint64(off)
        if off + b > 64: Wd[..., wi + 1] |= v >> np.uint64(64 - off)
    top = D[nl - 1]
    out = {}
    for i, j in zip(rows, cols):
        out[(i, j)] = int.from_bytes(Wd[i, j].tobytes(), 'little') + (int(top[i, j]) << nb)
    return out

def validate(D, b, lo, eA, eB, P, radS, Ref, sym=False, ctl=None):
    n = D.shape[1]
    if sym:
        iu = np.triu_indices(n); rows, cols = iu
    else:
        rows, cols = np.indices((n, n)).reshape(2, -1)
    vals = limbs_to_int(D, b, rows, cols)
    se = lo - (P - eA) - (P - eB)
    worst = 0.0; md = 0.0; mra = 0.0; mrs = 0.0; bad = 0
    for (i, j), g in vals.items():
        for (ii, jj) in ((i, j), (j, i)) if (sym and i != j) else ((i, j),):
            x = Ref[ii, jj]; m, ex = x.mid().man_exp(); m = int(m); ex = int(ex)
            em = min(se, ex); di = (g << (se - em)) - (m << (ex - em))
            dd = abs(math.ldexp(float(di), em)) if di else 0.0
            ra = float(x.rad()); rs = radS[i, j]
            q = dd / (ra + rs); worst = max(worst, q); md = max(md, dd); mra = max(mra, ra); mrs = max(mrs, rs)
            if q > 1: bad += 1
    return dict(maxdiff=md, max_rad_arb=mra, max_rad_sliced=mrs, worst_ratio=worst, n_fail=bad, ok=bad == 0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cell', type=int, default=6); ap.add_argument('--c0', type=int, default=8)
    ap.add_argument('--nch', type=int, default=16); ap.add_argument('--nm', type=int, default=1188)
    ap.add_argument('--P', type=int, default=480); ap.add_argument('--ks', default='64,256,1024')
    ap.add_argument('--out', default=''); ap.add_argument('--control', action='store_true')
    ap.add_argument('--threads', default='1,8'); ap.add_argument('--novalidate', action='store_true')
    a = ap.parse_args()
    prec = 512; ctx.prec = prec
    L = mc.parseL('71/50'); N = 2376; p = 0
    mc.G.update(prec=prec, L=L, N=N, p=p)
    mc.G['nrm'] = [((2 * d + 1) / (2 * L)).sqrt() for d in range(N)]
    e2L = (2 * L).exp(); sh = []; n = 2
    while arb(n) < e2L:
        lam = mc.mangoldt(n)
        if lam is not None:
            c = arb(n).log() / arb(n).sqrt() if lam == math.log(n) else arb([q for q in range(2, n + 1) if n % q == 0][0]).log() / arb(n).sqrt()
            sh += [(arb(n).log(), c), (-arb(n).log(), c)]
        n += 1
    cand = []
    for (ell, _) in sh:
        if float(ell.mid()) > 0:
            z = ell - L if float((ell - L).mid()) > 0 else L - ell
            if 0 < float(z.mid()) < float(L.mid()) and all(abs(float((z - y).mid())) > 1e-30 for y in cand): cand.append(z)
    bpa = [arb(0)] + sorted(cand, key=lambda z: float(z.mid())) + [L]
    ic = a.cell; x0, x1 = bpa[ic], bpa[ic + 1]; degs = list(range(p, N, 2))[:a.nm]; Np = len(degs)
    hh = (x1 - x0) / 2; xc = (x0 + x1) / 2
    ks = list(range(a.c0 * 64, (a.c0 + a.nch) * 64)); K = len(ks)
    t0 = time.time()
    V = arb_mat(Np, K); F = arb_mat(Np, K); FW = arb_mat(K, Np); G = arb_mat(Np, K)
    for c, kk in enumerate(ks):
        r, w = arb.legendre_p_root(N, kk, weight=True); x = xc + hh * r; wk = 2 * hh * w
        assert wk > 0; sw = wk.sqrt()
        v = mc.legv(x / L, N); f = [arb(0)] * N
        for (aa, ca) in sh:
            y = x - aa
            if (y > -L) and (y < L):
                va = mc.legv(y / L, N); f = [f[d] + ca * va[d] for d in range(N)]
        for i, d in enumerate(degs):
            V[i, c] = v[d]; F[i, c] = f[d]; FW[c, i] = f[d] * wk; G[i, c] = f[d] * sw
    t_gen = time.time() - t0; print('gen', round(t_gen, 1), 'K', K, 'Np', Np, flush=True)
    ths = [int(x) for x in a.threads.split(',')]
    arbt = {}
    for th in ths:
        ctx.threads = th; t0 = time.time(); Cp = V * FW; Cp2 = F * FW; arbt[th] = time.time() - t0
        print('arb', th, round(arbt[th], 2), flush=True)
    ctx.threads = 1
    t0 = time.time()
    Vz, eV, rV, uV = fixed(V, a.P); Wz, eW, rW, uW = fixed(FW.transpose(), a.P); Gz, eG, rG, uG = fixed(G, a.P)
    del G, V, F, FW
    import gc; gc.collect()
    if a.control: Vz[0][0] += 1 << 20; Gz[0][0] += 1 << 20   # negative controls for Cp and Cp2
    t_fix = time.time() - t0; print('fixed', round(t_fix, 1), flush=True)
    for kb in map(int, a.ks.split(',')):
        b, S = pick_b(a.P, kb); m = K // kb
        Vs = Ws = Gs = None; gc.collect()
        t0 = time.time()
        Vs, lo, Va = slice_words(Vz, a.P, b, S); Ws, _, Wa = slice_words(Wz, a.P, b, S); Gs, _, Ga = slice_words(Gz, a.P, b, S)
        t_slice = time.time() - t0
        sV, sW, sG = 2.0 ** -(a.P - eV), 2.0 ** -(a.P - eW), 2.0 ** -(a.P - eG)
        for th in ths:
            res = {}
            with threadpool_limits(limits=th):
                for name in ('Cp', 'Cp2_gemm', 'Cp2_syrk'):
                    D = np.zeros((2 * S, Np, Np), dtype=np.int64)
                    t0 = time.time()
                    for q in range(m):
                        sl = slice(q * kb, (q + 1) * kb)
                        if name == 'Cp': acc_gemm(D, [x[:, sl] for x in Vs], [x[:, sl] for x in Ws], S, kb)
                        elif name == 'Cp2_gemm':
                            # F itself is not stored in slices; F W^T == G G^T exactly in truth; gemm variant uses G slices too
                            acc_gemm(D, [x[:, sl] for x in Gs], [x[:, sl] for x in Gs], S, kb)
                        else: acc_syrk(D, [x[:, sl] for x in Gs], S, kb)
                    t_mm = time.time() - t0; t0 = time.time(); carry(D, b); t_cy = time.time() - t0
                    # carry per chunk-batch (no amortisation) would cost m * t_cy
                    r = dict(name=name, k=kb, batches=m, K=K, b=b, S=S, P=a.P, threads=th, Np=Np, t_arb_both=round(arbt[th], 3),
                             t_mm=round(t_mm, 3), t_carry_once=round(t_cy, 3), t_carry_if_per_batch=round(m * t_cy, 3),
                             t_carry_if_per_chunk64=round(K // 64 * t_cy, 3))
                    if not a.novalidate and th == ths[0]:
                        t0 = time.time()
                        if name == 'Cp':
                            radS = ((rV + uV) @ (Wa * sW + rW).T + (Va * sV) @ (rW + uW).T) * (1 + 1e-9)
                            r.update(validate(D, b, 2 * lo, eV, eW, a.P, radS, Cp))
                        else:
                            E = rG + uG; Gab = Ga * sG
                            radS = (Gab @ E.T + E @ Gab.T + E @ E.T) * (1 + 1e-9)
                            r.update(validate(D, b, 2 * lo, eG, eG, a.P, radS, Cp2, sym=(name == 'Cp2_syrk')))
                        r['t_validate'] = round(time.time() - t0, 1)
                    r.update(t_gen=round(t_gen, 1), t_fix=round(t_fix, 1), t_slice_py=round(t_slice, 2), control=a.control, cell=ic, c0=a.c0)
                    print(json.dumps(r), flush=True)
                    if a.out:
                        with open(a.out, 'a') as fh: fh.write(json.dumps(r) + '\n')
                    del D

if __name__ == '__main__':
    main()
