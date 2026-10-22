"""Grok S89: explicit band split of the Weil window form, Q = [[Q_GG, Q_G+],[Q_+G, Q_++]] in the eigenbasis of the
arch-only form A (G = arch-negative 'ghost' span, + = arch-positive complement).
  Q >= 0  <=>  (i) out-of-band  Q_++ > 0   and   (ii) in-band  S_G = Q_GG - Q_G+ Q_++^{-1} Q_+G >= 0   (Schur/Haynsworth).
Also the pole-free out-of-band ratio rho_out = lambda_max(P_++, A_++) (primes vs Gamma-positive tail; rho_out < 1 <=> R_++ > 0).
Forms: ZETA, DH, FT (theta), and HT(t) = backward heat flow H_t: Q_t = Q_zeta + Delta_t (zero-side swap of the first M zeros,
zero lists from s64/hneg). Sign conventions CHECKED against the zero side (calib.py): Q_e = R_e + sum p p^T, Q_o = R_o - sum q q^T.
HIGH-PRECISION NUMERICAL (mpmath dps 90), not interval. RH OPEN."""
import json, sys, os, mpmath as mp
from multiprocessing import Pool
sys.path.insert(0, "../s87"); import hp_forms as H
sys.path.insert(0, "../s88"); from ghost_overlap import nodes
ZZ = None
def zeta_zeros(M):
    return [mp.mpf(l.split()[1]) for _, l in zip(range(M), open('/workspace/GPPVerify/discovery/cutkosky_weil/zeta_zeros_100k.txt'))]
def ht_zeros(t):
    d = '/workspace/GPPDiscovery5/s64/hneg/'
    re_ = [mp.mpf(repr(x)) / 2 for x in json.load(open(d + f'scan_t{t}.json'))['real_zeros_x']]
    cp = [mp.mpc(repr(a), repr(b)) for a, b in json.load(open(d + f'cplx_t{t}.json'))['gamma_form']]
    zs = sorted([(z, 1) for z in re_] + [(z, 2) for z in cp], key=lambda q: mp.re(q[0]))
    return zs
def I(z, i, L):  # int_0^1 P_i(2s-1) e^{i z L s} ds = e^{i a} i^i j_i(a), a = zL/2
    a = z * L / 2
    j = mp.sqrt(mp.pi / 2) * mp.besselj(i + mp.mpf(1) / 2, a) / mp.sqrt(a)  # consistent principal branches: = a^i * entire
    return mp.exp(1j * a) * (1j) ** i * j
def uvec(z, L, N, par):
    out = []
    for i in range(N):
        c = 2 * mp.sqrt(L) * mp.sqrt(2 * i + 1)
        if par == 'e': out.append(c * (I(z, i, L) + I(-z, i, L)) / 2)
        else: out.append(c * (I(z, i, L) - I(-z, i, L)) / 2j)
    return out
def delta_ht(t, L, N, par, M=int(os.environ.get("MZ", "50"))):
    zs = ht_zeros(t); sel, cnt = [], 0
    for z, m in zs:
        if cnt + m > M: break
        sel.append((z, m)); cnt += m
    Mz = cnt; D = mp.zeros(N, N)
    for z, m in sel:
        u = uvec(z, L, N, par)
        for a in range(N):
            for b in range(N):
                D[a, b] += mp.re(u[a] * u[b]) * (1 if m == 1 else 2)
    for g in zeta_zeros(Mz):
        u = [mp.re(x) for x in uvec(g, L, N, par)]
        for a in range(N):
            for b in range(N): D[a, b] -= u[a] * u[b]
    return D, Mz
def analyse(Q, A, P, poles, sgn, N):
    EA, VA = mp.eigsy(A)
    gi = [i for i in range(N) if EA[i] < 0]; pi_ = [i for i in range(N) if EA[i] >= 0]
    VG = mp.matrix(N, len(gi)); VP = mp.matrix(N, len(pi_))
    for c, i in enumerate(gi):
        for r in range(N): VG[r, c] = VA[r, i]
    for c, i in enumerate(pi_):
        for r in range(N): VP[r, c] = VA[r, i]
    Qpp = VP.T * Q * VP; Qgg = VG.T * Q * VG; Qgp = VG.T * Q * VP
    Epp, Vpp = mp.eigsy(Qpp); out_min = min(Epp)
    App = VP.T * A * VP; Ppp = VP.T * P * VP
    # rho_out = lambda_max(App^{-1/2} Ppp App^{-1/2})
    Ea, Va = mp.eigsy(App); Ais = Va * mp.diag([1 / mp.sqrt(x) for x in Ea]) * Va.T
    rho = max(mp.eigsy(Ais * Ppp * Ais, eigvals_only=True))
    S = Qgg - Qgp * mp.inverse(Qpp) * Qgp.T if len(gi) else None
    in_min = min(mp.eigsy(S, eigvals_only=True)) if len(gi) else None
    EQ, VQ = mp.eigsy(Q); iq = min(range(N), key=lambda i: EQ[i]); v = VQ[:, iq]
    span = sum(float((v.T * VA[:, i])[0] ** 2) for i in gi)
    # failing / tightest vector inside the + block
    jp = min(range(len(pi_)), key=lambda i: Epp[i]); wv = VP * Vpp[:, jp]
    # pole weight out of band
    pole_out = [float(mp.norm(VP.T * mp.matrix(p)) ** 2 / mp.norm(mp.matrix(p)) ** 2) for p in poles]
    return {"nG": len(gi), "qmin": mp.nstr(EQ[iq], 8), "nodes_v": nodes([v[i] for i in range(N)]), "span_v": span,
            "out_min": mp.nstr(out_min, 8), "nodes_out": nodes([wv[i] for i in range(N)]), "rho_out": mp.nstr(rho, 10),
            "in_schur_min": mp.nstr(in_min, 8) if in_min is not None else None, "pole_out_frac": pole_out, "sgn": sgn}
def task(args):
    FORM, th, L, N = args
    H.setup(90); base = 'ZETA' if FORM == 'HT' else FORM
    d = H.build(base, L, N, th if FORM == 'FT' else 0)
    out = {"form": FORM, "theta": th, "L": L, "N": N}
    for par in "eo":
        A = mp.matrix(d["A" + par]); P = mp.matrix(d["P" + par]); Q = A - P
        sgn = 1 if par == 'e' else -1
        for p in d["p" + par]:
            pv = mp.matrix(p); Q += sgn * pv * pv.T
        if FORM == 'HT':
            D, Mz = delta_ht(th, mp.mpf(L), N, par); Q = Q + D; out["M_zeros"] = Mz
        out[par] = analyse(Q, A, P, d["p" + par], sgn, N)
    return out
if __name__ == "__main__":
    Ls = [round(0.4 + 0.1 * i, 1) for i in range(11)]
    which = sys.argv[1] if len(sys.argv) > 1 else 'main'
    if which == 'main':
        tasks = [(F, th, L, 24) for F, th in (("ZETA", 0), ("DH", 0), ("FT", 0.1)) for L in Ls]
    else:
        tasks = [("HT", t, L, int(which.split(':')[1])) for t in ("-1", "-2", "-4") for L in Ls]
    res = []
    with Pool(int(os.environ.get('NP', '6'))) as pool:
        for r in pool.imap_unordered(task, tasks):
            res.append(r)
            print(r["form"], r["theta"], r["L"], "| E", {k: r["e"][k] for k in ("qmin", "out_min", "in_schur_min", "rho_out", "nG", "nodes_v", "nodes_out")},
                  "span %.3f" % r["e"]["span_v"], "| O", {k: r["o"][k] for k in ("qmin", "out_min", "in_schur_min", "rho_out", "nG", "nodes_v", "nodes_out")}, "span %.3f" % r["o"]["span_v"], flush=True)
    res.sort(key=lambda r: (r["form"], str(r["theta"]), r["L"]))
    json.dump(res, open(f"band_{which.replace(':', '_')}.json", "w"), indent=1)
