"""Grok S81 H80c step 1: Ozaki-sliced (error-free float64 dgemm) replacement for the x-space cell products
Cp_cell = V * FW and Cp2_cell = F * FW of s71/mori_cert.py work_cell, validated against the Arb products on the SAME
Arb-generated V, F, FW (support 2.84 even: L=71/50, N=2376, prec 512).
Method: X_mid -> signed integer at scale 2^(P-e) (e = ceil log2 max|X|), split |X| into S slices of b bits
(2b + ceil log2 k <= 53 so every slice dgemm is exact in float64), accumulate per diagonal t=i+j in int64 (exact),
recombine in Python big ints. Optional dropping of pairs i+j > S-1+keep, with a rigorous bound.
Rigorous sliced radius (entrywise): sum_c (rV+uV)(|FW|+rFW) + |Vr|(rFW+uFW) + dropped, u = 2^-(P-e+1) rounding,
r = Arb input radii; evaluated in float64 times (1+1e-9) margin.
Check: |mid_arb - sliced| <= rad_arb + rad_sliced entrywise (both balls contain the true product)."""
import sys, os, time, json, math, argparse
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's71'))
import mori_cert as mc
from flint import arb, arb_mat, ctx

def to_fixed(M, P):
    """arb_mat -> (signed int object array of mids at scale 2^(P-e), e, radius float array, rounding u)"""
    n, k = M.nrows(), M.ncols()
    mids = [[M[i, j].mid() for j in range(k)] for i in range(n)]
    mx = max(abs(float(x)) for r in mids for x in r)
    e = math.ceil(math.log2(mx)) if mx > 0 else 0
    sh = P - e
    out = np.empty((n, k), dtype=object); rad = np.zeros((n, k))
    for i in range(n):
        for j in range(k):
            m, ex = mids[i][j].man_exp(); m = int(m); ex = int(ex); t = ex + sh
            if t >= 0: v = m << t
            else:   # round to nearest
                v = (m + (1 << (-t - 1))) >> (-t)
            out[i, j] = v; rad[i, j] = float(M[i, j].rad())
    return out, e, rad, 2.0 ** (-(sh + 1))

def slices(Z, P, b, S):
    sgn = np.where(np.array([z < 0 for z in Z.ravel()]).reshape(Z.shape), -1.0, 1.0)
    A = np.vectorize(abs, otypes=[object])(Z)
    out = []
    for s in range(S):
        sh = P - (s + 1) * b + 1   # top slice includes bit P (values may reach 2^P after rounding)
        mask = (1 << b) - 1
        if sh >= 0: f = np.vectorize(lambda x: float((x >> sh) & mask), otypes=[float])(A)
        else: f = np.vectorize(lambda x: float((x << (-sh)) & mask), otypes=[float])(A)
        out.append(sgn * f)
    return out, P + 1 - S * b   # lowest-bit exponent offset (<=0)

TIM = {}
CAT = bool(os.environ.get('H80C_CAT'))
def sliced_prod(Va, Wa, P, b, S, keep=None):
    """exact Va @ Wa.T-like product (Va n x k, Wa k x n given as n x k transposed already: Va @ Wb) via slices."""
    k = Va.shape[1]; assert 2 * b + math.ceil(math.log2(k)) <= 53
    t0 = time.time(); Vs, lo = slices(Va, P, b, S); Ws, _ = slices(Wa, P, b, S); TIM['slice'] = TIM.get('slice', 0) + time.time() - t0
    t0 = time.time()
    n1, n2 = Va.shape[0], Wa.shape[1] if False else Wa.shape[0]
    acc = {}; dropped_terms = 0
    tmax = 2 * S - 2 if keep is None else min(2 * S - 2, S - 1 + keep)
    if CAT:   # one dgemm per anti-diagonal t: [Vs_i]_{i} concatenated along k (inner dim <= S k; needs 2b+log2(Sk) <= 53)
        assert 2 * b + math.ceil(math.log2(S * k)) <= 53
        VC = np.concatenate(Vs, axis=1); WC = np.concatenate(Ws, axis=1)
        for t in range(tmax + 1):
            i0, i1 = max(0, t - S + 1), min(S - 1, t)
            Vblk = VC[:, i0 * k:(i1 + 1) * k]
            Wblk = np.concatenate([Ws[t - i] for i in range(i0, i1 + 1)], axis=1)
            acc[t] = (Vblk @ Wblk.T).astype(np.int64)
        dropped_terms = sum(1 for i in range(S) for j in range(S) if i + j > tmax)
    for i in (range(S) if not CAT else []):
        for j in range(S):
            t = i + j
            if t > tmax: dropped_terms += 1; continue
            C = Vs[i] @ Ws[j].T
            Ci = C.astype(np.int64)
            acc[t] = acc[t] + Ci if t in acc else Ci
    TIM['dgemm'] = TIM.get('dgemm', 0) + time.time() - t0
    # production recombination: base-2^b limbs, vectorized carry propagation (int64 only)
    t0 = time.time(); D = np.zeros((2 * S,) + acc[0].shape, dtype=np.int64)
    for t, Ct in acc.items(): D[2 * S - 2 - t] += Ct
    for d in range(2 * S - 1):
        cy = D[d] >> b; D[d] &= (1 << b) - 1; D[d + 1] += cy
    TIM['carry'] = TIM.get('carry', 0) + time.time() - t0
    # value = sum_t acc[t] * 2^{(S-1-i)b + (S-1-j)b} * 2^{2 lo}  with i+j=t
    G = np.zeros(acc[0].shape, dtype=object)
    for t, Ct in acc.items():
        G = G + (Ct.astype(object) << ((2 * S - 2 - t) * b))
    # scale: entry = G * 2^{2 lo} in units of 2^-P each factor -> true = G * 2^{2lo} * 2^{-(P-eV)} * 2^{-(P-eW)}
    for jj in range(G.shape[1]):   # limb path == bigint path (row 0 and last row)
        for ii in (0, G.shape[0] - 1):
            assert sum(int(D[d, ii, jj]) << (d * b) for d in range(2 * S)) == G[ii, jj]
    # dropped bound: each dropped slice pair contributes < k * 2^{2b} * 2^{(2S-2-t)b} units; sum bounded below
    drop_units = 0
    for i in range(S):
        for j in range(S):
            if i + j > tmax: drop_units += k * (1 << (2 * b)) * (1 << ((2 * S - 2 - i - j) * b))
    return G, 2 * lo, drop_units

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cells', default='0,5,11'); ap.add_argument('--chunks', default='0,19')
    ap.add_argument('--nm', type=int, default=300); ap.add_argument('--P', type=int, default=480)
    ap.add_argument('--keep', type=int, default=-1); ap.add_argument('--out', default='')
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
    mc.G['shifts'] = sh
    cand = []
    for (ell, _) in sh:
        if float(ell.mid()) > 0:
            z = ell - L if float((ell - L).mid()) > 0 else L - ell
            if 0 < float(z.mid()) < float(L.mid()) and all(abs(float((z - y).mid())) > 1e-30 for y in cand): cand.append(z)
    bpa = [arb(0)] + sorted(cand, key=lambda z: float(z.mid())) + [L]; mc.G['bpa'] = bpa
    print('cells', len(bpa) - 1, flush=True)
    chunk = 64; gl = {}
    nm = a.nm; keep = None if a.keep < 0 else a.keep
    b = (53 - math.ceil(math.log2(chunk))) // 2; S = math.ceil((a.P + 1) / b)
    if CAT:
        b = 21
        while 2 * b + math.ceil(math.log2(math.ceil((a.P + 1) / b) * chunk)) > 53: b -= 1
        S = math.ceil((a.P + 1) / b)
    res = []
    for ic in map(int, a.cells.split(',')):
        for c0 in map(int, a.chunks.split(',')):
            ks = list(range(c0 * chunk, min(N, (c0 + 1) * chunk)))
            for kk in ks:
                if kk not in gl: gl[kk] = arb.legendre_p_root(N, kk, weight=True)
            mc.G['gl'] = gl
            t0 = time.time()
            x0, x1 = bpa[ic], bpa[ic + 1]; degs = list(range(p, N, 2))[:nm]; Np = len(degs)
            hh = (x1 - x0) / 2; xc = (x0 + x1) / 2
            V = arb_mat(Np, len(ks)); F = arb_mat(Np, len(ks)); FW = arb_mat(len(ks), Np)
            for c, kk in enumerate(ks):
                r, w = gl[kk]; x = xc + hh * r; wk = 2 * hh * w
                v = mc.legv(x / L, N); f = [arb(0)] * N
                for (aa, ca) in sh:
                    y = x - aa
                    if (y > -L) and (y < L):
                        va = mc.legv(y / L, N); f = [f[d] + ca * va[d] for d in range(N)]
                for i, d in enumerate(degs):
                    V[i, c] = v[d]; F[i, c] = f[d]; FW[c, i] = f[d] * wk
            t_gen = time.time() - t0
            t0 = time.time(); Cp = V * FW; Cp2 = F * FW; t_arb = time.time() - t0
            t0 = time.time()
            Vz, eV, rV, uV = to_fixed(V, a.P); Fz, eF, rF, uF = to_fixed(F, a.P); Wz, eW, rW, uW = to_fixed(FW.transpose(), a.P)
            t_conv = time.time() - t0
            if os.environ.get('H80C_CONTROL'):   # negative control: perturb one V entry by 2^20 ulps (fixed-point units)
                Vz[0, 0] += 1 << 20
            t0 = time.time()
            G1, lo1, d1 = sliced_prod(Vz, Wz, a.P, b, S, keep); G2, lo2, d2 = sliced_prod(Fz, Wz, a.P, b, S, keep)
            t_sl = time.time() - t0; tim = dict(TIM); TIM.clear()
            out = {}
            for name, Gm, lo, dr, Az, eA, rA, uA, Ref in (('Cp', G1, lo1, d1, Vz, eV, rV, uV, Cp), ('Cp2', G2, lo2, d2, Fz, eF, rF, uF, Cp2)):
                sA = 2.0 ** -(a.P - eA); sW = 2.0 ** -(a.P - eW)
                Aabs = np.abs(Az.astype(float)) * sA; Wabs = np.abs(Wz.astype(float)) * sW
                radS = (rA + uA) @ (Wabs + rW).T + Aabs @ (rW + uW).T
                radS = (radS + dr * 2.0 ** lo * sA * sW) * (1 + 1e-9)
                ctx.prec = 2048
                scale_exp = lo - (a.P - eA) - (a.P - eW)
                worst = 0.0; maxdiff = 0.0; maxradA = 0.0; ok = True
                for i in range(Np):
                    for j in range(Np):
                        val = arb(int(Gm[i, j])) * arb(2) ** scale_exp
                        dd = abs(float((val - Ref[i, j].mid()).mid()))
                        ra = float(Ref[i, j].rad())
                        maxdiff = max(maxdiff, dd); maxradA = max(maxradA, ra)
                        q = dd / (ra + radS[i, j]); worst = max(worst, q)
                        if q > 1: ok = False
                ctx.prec = prec
                out[name] = dict(maxdiff=maxdiff, max_rad_arb=maxradA, max_rad_sliced=float(radS.max()), worst_ratio=worst, ok=ok, dropped_pairs=dr > 0)
            r = dict(cell=ic, chunk=c0, Np=Np, k=len(ks), P=a.P, b=b, S=S, keep=a.keep, cat=CAT, threads=os.environ.get('OPENBLAS_NUM_THREADS', 'default'), t_gen=round(t_gen, 2), t_arb=round(t_arb, 3), t_conv=round(t_conv, 3), t_sliced_total_incl_bigint_check=round(t_sl, 3), t_slice=round(tim['slice'], 3), t_dgemm=round(tim['dgemm'], 3), t_carry=round(tim['carry'], 3), speedup_core=round(t_arb / (tim['dgemm'] + tim['carry']), 2), **out)
            print(json.dumps(r), flush=True); res.append(r)
            if a.out:
                with open(a.out, 'a') as fh: fh.write(json.dumps(r) + '\n')

if __name__ == '__main__':
    main()
