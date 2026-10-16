"""S45 H45b: p=3 (k=1) single-prime coupling window spot — reuse S44 gen-eig, T3 not T2.
c_true = 2 log3 / √3. Hold all else at ζ. Outer gaps at L∈{0.55,0.70,1.0,1.3}, M=40.
CERTIFIED edges via rational Ritz vector (same as s44_cert). Rigidity ≠ RH evidence (S44 DEAD END).
"""
import sys, json, time
sys.path.insert(0, '.')
from s14_exact import *
import mpmath as mp
from flint import arb, fmpq, ctx, arb_mat, acb

def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return (p if m == 1 else None)
    return None

def build_all(M, h, prec):
    ctx.prec = prec
    ha = arb(h)
    Gs = {}
    for j in range(M):
        for k in range(j, M):
            if (j + k) % 2 == 0:
                g = (G(j, k) + G(k, j)) / 2
                Gs[j, k] = Gs[k, j] = g
    D = max(g.degree() for g in Gs.values()) + 1
    def mom(phi):
        return [acb.integral(lambda z, an, m=m: z ** m * phi(z), 0, 2).real for m in range(D + 1)]
    hc = acb(ha)
    m1 = mom(lambda z: (-hc * z).exp() / (acb(0, 1) * hc * z).sinc())
    m2 = mom(lambda z: (-hc * z / 4).exp() * (acb(0, 1) * hc * z * 3 / 4).sinc() * fmpq(3, 4) / (acb(0, 1) * hc * z).sinc())
    m3 = mom(lambda z: (hc * z / 2).cosh())
    c0 = -arb.const_euler() - arb.pi().log() - (1 - (-4 * ha).exp()).log()
    ns = [n for n in range(2, int(mp.exp(2 * mp.mpf(h))) + 2)
          if vonmangoldt(n) and mp.log(n) < 2 * mp.mpf(h)]
    A = arb_mat(M, M)
    T = {n: arb_mat(M, M) for n in ns}
    for (j, k), g in Gs.items():
        co = [arb(c.p) / c.q for c in g.coeffs()]
        G0 = co[0]
        q = [-c for c in co[1:]]
        pole = 4 * ha * ha * sum(co[m] * m3[m] for m in range(len(co)))
        arch = (ha * c0 * G0 + ha * sum(q[m] * m1[m] for m in range(len(q)))
                - 2 * ha * ha * sum(co[m] * m2[m] for m in range(len(co))))
        A[j, k] = pole + arch
        for n in ns:
            v = arb(n).log() / ha
            T[n][j, k] = ha * sum(co[m] * v ** m for m in range(len(co)))
    return A, T, ns

def tomp(Ma, I):
    return mp.matrix([[mp.mpf(Ma[i, j].mid().str(mp.mp.dps + 20, radius=False)) for j in I] for i in I])

def weight(n):
    p = vonmangoldt(n)
    return 2 * mp.log(p) / mp.sqrt(n)

def run_cert(M, h, prec, n0=3):
    t0 = time.time()
    A, T, ns = build_all(M, h, prec)
    if n0 not in ns:
        return {"M": M, "h": h, "n0": n0, "note": f"n={n0} not in support (need 2h>log {n0})", "sec": 0}
    W = arb_mat(M, M)
    for i in range(M):
        for j in range(M):
            W[i, j] = A[i, j] - sum((2 * arb(vonmangoldt(n)).log() / arb(n).sqrt()) * T[n][i, j] for n in ns)
    c_true = weight(n0)
    row = {"M": M, "h": h, "n0": n0, "c_true": mp.nstr(c_true, 20), "prime_powers": ns}
    for nm, I in (("even", list(range(0, M, 2))), ("odd", list(range(1, M, 2)))):
        w = tomp(W, I)
        t = tomp(T[n0], I)
        try:
            Lc = mp.cholesky(w)
        except Exception:
            row[nm] = {"note": "W_zeta NOT PD on Ritz space"}
            continue
        Li = mp.inverse(Lc)
        S = Li * t * Li.T
        S = (S + S.T) / 2
        mu, V = mp.eigsy(S)
        out = {"mu_max": mp.nstr(max(mu), 12), "mu_min": mp.nstr(min(mu), 12),
               "gap_up_ritz": mp.nstr(1 / max(mu), 12) if max(mu) > 0 else "inf",
               "gap_down_ritz": mp.nstr(-1 / min(mu), 12) if min(mu) < 0 else "inf"}
        for side, idx in (("up", max(range(len(I)), key=lambda i: mu[i])),
                          ("down", min(range(len(I)), key=lambda i: mu[i]))):
            y = Li.T * V[:, idx]
            sc = max(abs(v) for v in y)
            y = [v / sc for v in y]
            xs = [fmpq(int(mp.nint(v * mp.mpf(10) ** 70)), 10 ** 70) for v in y]
            xa = [arb(q.p) / q.q for q in xs]
            wv = sum(xa[a] * xa[b] * W[I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
            tv = sum(xa[a] * xa[b] * T[n0][I[a], I[b]] for a in range(len(I)) for b in range(len(I)))
            r = wv / tv
            out[side] = {"Tnx_sign": "+" if tv > 0 else ("-" if tv < 0 else "?"),
                         "r_cert": r.str(14), "ritz_gap": mp.nstr(1 / mu[idx], 14)}
        row[nm] = out
    row["sec"] = round(time.time() - t0, 1)
    return row

if __name__ == "__main__":
    M = int(sys.argv[1])
    hs = [float(x) for x in sys.argv[2].split(',')]
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 110
    prec = int(sys.argv[4]) if len(sys.argv) > 4 else 768
    mp.mp.dps = dps
    for h in hs:
        print(json.dumps(run_cert(M, h, prec)), flush=True)
