"""Grok S90: H89a -- where do prime-power integrality and Lambda>=0 enter the in-band cancellation?
For the zeta parity window forms Q_e = A - P + p p^T, Q_o = A - P - q q^T (S87 builder, sign conventions CHECKED in s89/calib.py):
 (1) anatomy of the near-null vector v: v^T Q v = arch(in band) + arch(out of band) + pole - sum_n Lambda(n) v^T P_n v (per prime power);
 (2) atom-energy profile g(x) = v^T P_x v / c(x) for a unit atom at continuous x in (1, e^{2L}) -- where v 'wants' primes;
 (3) exact perturbations of the prime side: remove / sign-flip one prime (all powers), drop prime powers k>=2,
     Beurling-type jitter x -> x(1+delta) of all prime powers, global scale, insert composite 6 as a 'prime',
     and the smooth (PNT-main-term) density dpsi = dx replacing the primes.
HIGH-PRECISION NUMERICAL (mpmath dps 90, N=24), not interval. RH OPEN."""
import json, sys, os, mpmath as mp
from math import comb, factorial
from multiprocessing import Pool
sys.path.insert(0, "../s87"); import hp_forms as H

def mono_atom(x, c, L, N, fac, acc):
    ln = mp.log(x); sg = ln / L; Mh, Mh2, Mt = acc
    if sg < 2:
        if sg <= 1:
            for k in range(N):
                for l in range(N): Mh[k][l] += c * fac[k] * fac[l] / fac[k + l + 1] * sg ** (k + l + 1)
        else:
            t = 2 - sg
            for k in range(N):
                for l in range(N): Mh2[k][l] += c * fac[k] * fac[l] / fac[k + l + 1] * t ** (k + l + 1)
    if sg < 1:
        for k in range(N):
            for l in range(N):
                Mt[k][l] += 2 * c * mp.fsum(comb(k, a) / mp.mpf(a + l + 1) * sg ** (k - a) * (1 - sg) ** (a + l + 1) for a in range(k + 1))

def Pmats(atoms, L, N):  # atoms: list of (x, weight) ; coefficient = weight/sqrt(x)
    fac = [mp.mpf(factorial(n)) for n in range(2 * N + 2)]
    acc = [[[mp.mpf(0)] * N for _ in range(N)] for _ in range(3)]
    for x, w in atoms:
        if w != 0: mono_atom(mp.mpf(x), w / mp.sqrt(x), L, N, fac, acc)
    C = H.legC(N); D = [(-1) ** i for i in range(N)]
    P1 = H.CMCt(C, acc[0]); P2 = H.CMCt(C, acc[1]); Pt = H.CMCt(C, acc[2])
    Bp = H.sym([[P1[i][j] + D[i] * D[j] * P2[i][j] for j in range(N)] for i in range(N)]); Tp = H.sym(Pt)
    return mp.matrix(H.add(Tp, Bp)), mp.matrix(H.add(Tp, Bp, 1, -1))

def primepowers(X):
    out = []
    for n in range(2, int(X) + 1):
        m = H.mang(n)
        if m != 0:
            p = next(q for q in range(2, n + 1) if n % q == 0); out.append((n, p, m))
    return out

def qmin(M):
    E, V = mp.eigsy(M); i = min(range(M.rows), key=lambda j: E[j]); return E[i], V[:, i]

def task(L):
    H.setup(90); N = int(os.environ.get("NB", "24")); L = mp.mpf(L)
    d = H.build('ZETA', L, N); X = mp.exp(2 * L)
    pp = primepowers(X); primes = sorted({p for _, p, _ in pp})
    res = {"L": float(L), "N": N, "prime_powers": [n for n, _, _ in pp]}
    Pn = {n: Pmats([(n, 1)], L, N) for n, _, _ in pp}           # unit-Lambda atom matrices (e, o)
    comp = [n for n in (6, 10, 12) if n <= X]
    Pc = {n: Pmats([(n, 1)], L, N) for n in comp}
    def Qform(par, Pe_o):
        A = mp.matrix(d["A" + par]); sgn = 1 if par == 'e' else -1; Q = A - Pe_o
        for p in d["p" + par]: pv = mp.matrix(p); Q += sgn * pv * pv.T
        return Q
    for pi_, par in enumerate("eo"):
        A = mp.matrix(d["A" + par]); sgn = 1 if par == 'e' else -1
        P0 = sum((m * Pn[n][pi_] for n, _, m in pp), mp.zeros(N, N))
        Q0 = Qform(par, P0); q0, v = qmin(Q0)
        EA, VA = mp.eigsy(A)
        c = [(v.T * VA[:, i])[0] for i in range(N)]
        arch_in = mp.fsum(EA[i] * c[i] ** 2 for i in range(N) if EA[i] < 0)
        arch_out = mp.fsum(EA[i] * c[i] ** 2 for i in range(N) if EA[i] >= 0)
        pole = sgn * mp.fsum((mp.matrix(p).T * v)[0] ** 2 for p in d["p" + par])
        per_n = {n: -m * (v.T * Pn[n][pi_] * v)[0] for n, _, m in pp}
        per_comp = {n: -(v.T * Pc[n][pi_] * v)[0] for n in comp}
        r = {"qmin": mp.nstr(q0, 8), "arch_in": mp.nstr(arch_in, 8), "arch_out": mp.nstr(arch_out, 8), "pole": mp.nstr(pole, 8),
             "per_n": {str(n): mp.nstr(x, 8) for n, x in per_n.items()},
             "per_unit_comp": {str(n): mp.nstr(x, 6) for n, x in per_comp.items()},
             "check_sum": mp.nstr(arch_in + arch_out + pole + mp.fsum(per_n.values()) - q0, 5)}
        # profile g(x) = -v^T P_x v per unit weight (contribution to v^T Q v of a unit-Lambda atom at x)
        xs = [mp.e ** (2 * L * (k + 0.5) / 40) for k in range(40)]
        r["profile"] = [(float(x), float(-(v.T * Pmats([(x, 1)], L, N)[pi_] * v)[0])) for x in xs]
        V = {}
        for p in primes:
            V[f"remove_{p}"] = Qform(par, P0 - sum((m * Pn[n][pi_] for n, q, m in pp if q == p), mp.zeros(N, N)))
            V[f"flip_{p}"] = Qform(par, P0 - 2 * sum((m * Pn[n][pi_] for n, q, m in pp if q == p), mp.zeros(N, N)))
        V["primes_only"] = Qform(par, P0 - sum((m * Pn[n][pi_] for n, q, m in pp if n != q), mp.zeros(N, N)))
        for eps in ("0.01", "-0.01"):
            V[f"scale_{eps}"] = Qform(par, (1 + mp.mpf(eps)) * P0)
        for dl in ("0.001", "-0.001", "0.01", "-0.01", "0.05", "-0.05"):
            V[f"jitter_{dl}"] = Qform(par, Pmats([(n * (1 + mp.mpf(dl)), m) for n, _, m in pp], L, N)[pi_])
        for n in comp:
            V[f"add_comp_{n}"] = Qform(par, P0 + mp.log(n) * Pc[n][pi_])
        # smooth PNT density: dpsi = dx on (1, e^{2L}) -> atoms at x=e^u with weight e^u du (Gauss-Legendre on [0,L],[L,2L])
        GL = H.gnodes(7); at = []
        for a0 in (0, 1):
            for k in range(8):
                for x, w in GL:
                    u = L * (a0 + (k + x) / 8); at.append((mp.e ** u, L / 8 * w * mp.e ** u))
        Psm = Pmats(at, L, N)[pi_]
        V["smooth_density"] = Qform(par, Psm)
        r["P_minus_smooth_norm"] = mp.nstr(mp.mnorm(P0 - Psm, 2) if False else mp.norm(P0 - Psm), 6)
        r["variants"] = {}
        for k, M in V.items():
            e, w = qmin(M)
            r["variants"][k] = {"qmin": mp.nstr(e, 6), "overlap_v": float(abs((w.T * v)[0]))}
        res[par] = r
        print(float(L), par, "qmin", r["qmin"], "arch_in", r["arch_in"], "arch_out", r["arch_out"], "pole", r["pole"], "chk", r["check_sum"], flush=True)
        print("   per_n", r["per_n"], "comp", r["per_unit_comp"], flush=True)
        print("   ", {k: x["qmin"] for k, x in r["variants"].items()}, flush=True)
    return res

if __name__ == "__main__":
    Ls = [float(x) for x in sys.argv[1:]] or [0.6, 0.8, 1.0, 1.2]
    with Pool(len(Ls)) as pool:
        out = pool.map(task, Ls)
    json.dump(out, open("prime_anatomy.json", "w"), indent=1)
