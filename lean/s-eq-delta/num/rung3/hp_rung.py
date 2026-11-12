"""High-precision NUMERICAL (python-flint Arb at `prec` bits for the matrix build, mpmath for
eigen-solves) Rayleigh-Ritz of the TRUE zeta Weil form Q on the window [-H, H], split as
V_old = L^2[-h, h]  (+)  annulus h < |x| <= H.  No zeros are used; no R_T lower form either:
Q is assembled directly from the in-tree time-domain definitions
  Q(f) = poleTerm + archimedeanTerm - primeSide  on  g = autocorr f,
  poleTerm = 2 A+ A-,  A+- = int f e^{+-x/2};
  archimedeanTerm = C_inf g(0) + 2 int_0^inf e^{-u/2} (g(0) - g(u)) / (1 - e^{-2u}) du;
  primeSide = sum_{log n < 2H} Lambda(n)/sqrt(n) (g(log n) + g(-log n)).
Basis: orthonormal Legendre modes on each piece.  The symmetrised cross-correlation of two basis
functions is a piecewise polynomial in u with breakpoints at differences of piece endpoints; it is
evaluated exactly (Gauss in y with enough nodes) at Gauss nodes in u on panels between breakpoints,
where the archimedean integrand is analytic.  Quadrature errors are not bounded rigorously, hence
NUMERICAL (high precision), not CERTIFIED."""
import sys, json, time, math
from flint import arb, arb_mat, ctx
import mpmath as mp

def gauss_nodes(m, dps):
    with mp.workdps(dps + 20):
        xs, ws = [], []
        for i in range(1, m + 1):
            x = mp.cos(mp.pi * (i - mp.mpf(1) / 4) / (m + mp.mpf(1) / 2))
            for _ in range(100):
                p0, p1 = mp.mpf(1), x
                for k in range(2, m + 1):
                    p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
                dp = m * (x * p1 - p0) / (x * x - 1)
                dx = p1 / dp; x -= dx
                if abs(dx) < mp.mpf(10) ** (-(dps + 10)): break
            p0, p1 = mp.mpf(1), x
            for k in range(2, m + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = m * (x * p1 - p0) / (x * x - 1)
            xs.append(x); ws.append(2 / ((1 - x * x) * dp * dp))
        return [arb(mp.nstr(x, dps + 15)) for x in xs], [arb(mp.nstr(w, dps + 15)) for w in ws]

def mangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return (p if m == 1 else None)

class Setup:
    def __init__(self, h, H, K0, K1, prec=192, mu=64, panel=0.12):
        ctx.prec = prec; self.dps = int(prec * 0.30)
        self.h, self.H = h, H
        cA = (h + H) / 2; wA = (H - h) / 2
        self.pieces = [(arb(0), h, K0), (cA, wA, K1), (-cA, wA, K1)]
        self.modes = [(pi, k) for pi, (_, _, K) in enumerate(self.pieces) for k in range(K)]
        self.N = len(self.modes); self.K0 = K0
        self.Kmax = max(K0, K1)
        self.gy = gauss_nodes(self.Kmax + 2, self.dps)
        self.gu = gauss_nodes(mu, self.dps)
        self.panel = panel
        self.ends = [-H, -h, h, H]
        n = 2; self.terms = []
        while True:
            ln = arb(n).log()
            if not (ln < 2 * H): break
            p = mangoldt(n)
            if p: self.terms.append((n, arb(p).log() / arb(n).sqrt(), ln))
            n += 1
    def piece_of(self, x):
        return None
    def eval_all(self, ys):
        """E[i][j] = e_i(y_j) for interior points (each y strictly inside one piece or outside)."""
        N = self.N; E = arb_mat(N, len(ys))
        for j, y in enumerate(ys):
            for pi, (c, w, K) in enumerate(self.pieces):
                z = (y - c) / w
                if z > -1 and z < 1:
                    base = sum(self.pieces[q][2] for q in range(pi))
                    p0, p1 = arb(1), z
                    for k in range(K):
                        val = p0 if k == 0 else p1 if k == 1 else None
                        if k >= 2:
                            p0, p1 = p1, ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; val = p1
                        E[base + k, j] = ((2 * k + 1) / (2 * w)).sqrt() * val
                    break
        return E
    def gsym(self, u):
        """symmetrised cross-correlation matrix G_ij(u) = (g_ij(u)+g_ij(-u))/2, u >= 0."""
        H = self.H
        pts = sorted(set([float(e) for e in self.ends] + [float(e + u) for e in self.ends]))
        pts = [p for p in pts if p >= float(-H + u) - 1e-30 and p <= float(H) + 1e-30]
        # rebuild exact arb endpoints from the float order (all are endpoints or endpoints+u)
        cand = list(self.ends) + [e + u for e in self.ends]
        cand = [c for c in cand if not (c < -H + u) and not (c > H)]
        cand.sort(key=lambda a: float(a))
        uniq = []
        for c in cand:
            if not uniq or float(c) - float(uniq[-1]) > 1e-25: uniq.append(c)
        ys, ws = [], []
        for a, b in zip(uniq[:-1], uniq[1:]):
            half = (b - a) / 2; mid = (a + b) / 2
            for x, w in zip(*self.gy):
                ys.append(mid + half * x); ws.append(half * w)
        if not ys:
            return arb_mat(self.N, self.N)
        E1 = self.eval_all(ys); E2 = self.eval_all([y - u for y in ys])
        W = arb_mat(len(ys), len(ys))
        for j, w in enumerate(ws): W[j, j] = w
        G = E1 * W * E2.transpose()
        return (G + G.transpose()) * arb(0.5)
    def build(self, verbose=True):
        t0 = time.time(); H, h = self.H, self.h; N = self.N
        G0 = self.gsym(arb(0))
        # pole moments
        Ap, Am = arb_mat(N, 1), arb_mat(N, 1)
        for pi, (c, w, K) in enumerate(self.pieces):
            xs = [c + w * x for x in self.gu[0]]; E = self.eval_all(xs)
            for i in range(N):
                sp = arb(0); sm = arb(0)
                for j, (x, wt) in enumerate(zip(xs, self.gu[1])):
                    if E[i, j] != 0:
                        sp += w * wt * E[i, j] * (x / 2).exp(); sm += w * wt * E[i, j] * (-x / 2).exp()
                Ap[i, 0] += sp; Am[i, 0] += sm
        pole = Ap * Am.transpose() + Am * Ap.transpose()
        self.cvec = (Ap + Am) * arb(0.5); self.svec = (Ap - Am) * arb(0.5)
        # archimedean
        Cinf = -arb.const_euler() - arb.pi().log() - arb.pi() / 2 - 3 * arb(2).log()
        bps = sorted(set([0.0, float(H - h), float(2 * h), float(H + h), float(2 * H)]))
        exact = {0.0: arb(0), float(H - h): H - h, float(2 * h): 2 * h, float(H + h): H + h, float(2 * H): 2 * H}
        I = arb_mat(N, N); nodes = 0
        for a, b in zip(bps[:-1], bps[1:]):
            A_, B_ = exact[a], exact[b]
            npan = max(1, int(math.ceil((b - a) / self.panel)))
            for q in range(npan):
                pa = A_ + (B_ - A_) * q / npan; pb = A_ + (B_ - A_) * (q + 1) / npan
                half = (pb - pa) / 2; mid = (pa + pb) / 2
                for x, w in zip(*self.gu):
                    u = mid + half * x
                    ker = (-u / 2).exp() / (1 - (-2 * u).exp())
                    I += (G0 - self.gsym(u)) * (half * w * ker); nodes += 1
            if verbose: print(f"  arch panel [{a:.4f},{b:.4f}] done, nodes={nodes}, {time.time()-t0:.0f}s", flush=True)
        tail = arb(0); m = 0
        while True:
            c = arb(1) / 2 + 2 * m; term = (-2 * H * c).exp() / c; tail += term; m += 1
            if term < arb(10) ** (-self.dps - 5): break
        arch = G0 * Cinf + (I + G0 * tail) * 2
        primes = {}
        for n, coef, ln in self.terms:
            primes[n] = self.gsym(ln) * (2 * coef)
        Q = pole + arch
        for n in primes: Q -= primes[n]
        self.parts = dict(pole=pole, arch=arch, **{f"prime_{n}": M for n, M in primes.items()})
        self.Q = Q; self.G0 = G0
        if verbose: print(f"  built N={N} in {time.time()-t0:.0f}s; comb={[n for n,_,_ in self.terms]}", flush=True)
        return Q

def tomp(M):
    return mp.matrix([[mp.mpf(M[i, j].mid().str(60, radius=False)) for j in range(M.ncols())] for i in range(M.nrows())])

def analyse(S, label):
    mp.mp.dps = S.dps
    Q = tomp(S.Q); G0 = tomp(S.G0); N = S.N; n0 = S.K0
    # orthonormality check
    orth = max(abs(G0[i, j] - (1 if i == j else 0)) for i in range(N) for j in range(N))
    A = Q[0:n0, 0:n0]; B = Q[0:n0, n0:N]; D = Q[n0:N, n0:N]
    la, Va = mp.eigsy(A); i0 = min(range(n0), key=lambda i: la[i]); lam_old = la[i0]; phi = Va[:, i0]
    lq, Vq = mp.eigsy(Q); j0 = min(range(N), key=lambda i: lq[i]); lam_new = lq[j0]
    ld, Vd = mp.eigsy(D); Dmin = min(ld)
    # D = L L^T ; c_full = 1 - lam_min(L^-1 S L^-T)
    L = mp.cholesky(D); Li = mp.inverse(L)
    S_ = D - B.T * (mp.inverse(A) * B)
    M = Li * S_ * Li.T; M = (M + M.T) / 2
    lm, Vm = mp.eigsy(M); deficit = min(lm)
    bv = B.T * phi; y = mp.lu_solve(D, bv); r_phi = (bv.T * y)[0] / lam_old
    psi = y / mp.norm(y); b = (phi.T * B * psi)[0]; d = (psi.T * D * psi)[0]
    split = {}
    for k, P in S.parts.items():
        Pm = tomp(P); split[k] = float((phi.T * Pm[0:n0, n0:N] * psi)[0])
    sgn = {"pole": 1, "arch": 1}
    contrib = {k: (v if k in sgn else -v) for k, v in split.items()}
    out = dict(label=label, h=float(S.h), H=float(S.H), K0=S.K0, N=N, comb=[n for n, _, _ in S.terms],
               orth_err=mp.nstr(orth, 3), lam_old=mp.nstr(lam_old, 12), lam_new=mp.nstr(lam_new, 12),
               D_min=mp.nstr(Dmin, 8), one_minus_c_full=mp.nstr(deficit, 12), c_full=mp.nstr(1 - deficit, 15),
               r_phi=mp.nstr(r_phi, 12), b=mp.nstr(b, 10), d=mp.nstr(d, 10),
               b_contrib={k: f"{v:.6e}" for k, v in contrib.items()},
               new_kick_over_b=(f"{contrib[f'prime_{S.terms[-1][0]}']/float(b):.3f}" if S.terms else None),
               ratio_new_old=mp.nstr(lam_new / lam_old, 8))
    # channel shares s_j = (B^T phi_j)^T D^{-1} (B^T phi_j) / lam_j for the lowest modes of A
    Di = mp.inverse(D); order = sorted(range(n0), key=lambda i: la[i])
    shares = []
    for i in order[:8]:
        bj = B.T * Va[:, i]; shares.append([mp.nstr(la[i], 8), mp.nstr((bj.T * Di * bj)[0] / la[i], 8)])
    out["channel_shares_lam_s"] = shares
    # saturating annulus direction w* (top of D^{-1/2} B^T A^{-1} B D^{-1/2}) and the phi_0 share there
    kmin = min(range(len(lm)), key=lambda i: lm[i]); wst = Li.T * Vm[:, kmin]
    Ai = mp.inverse(A); tot = (wst.T * B.T * Ai * B * wst)[0]; p0 = (phi.T * B * wst)[0] ** 2 / lam_old
    out["phi0_share_at_saturation"] = mp.nstr(p0 / tot, 10)
    Phi = Vq[:, j0]; u = Phi[0:n0]; w = Phi[n0:N]
    out["new_ground_state"] = dict(annulus_weight=mp.nstr(mp.norm(w) ** 2, 8),
                                   overlap2_with_phi0=mp.nstr((phi.T * u)[0] ** 2 / mp.norm(u) ** 2, 10))
    # Toeplitz + rank-2 structure: Q = T + 2 c c^T - 2 s s^T, T = arch - primes (truncated convolution)
    T = tomp(S.parts["arch"])
    for k_, P_ in S.parts.items():
        if k_.startswith("prime_"): T -= tomp(P_)
    cv = tomp(S.cvec); sv = tomp(S.svec)
    # parity operator in this basis
    Pm = mp.zeros(N, N); K0_, K1_ = S.pieces[0][2], S.pieces[1][2]
    for k_ in range(K0_): Pm[k_, k_] = (-1) ** k_
    for k_ in range(K1_):
        i1 = K0_ + k_; i2 = K0_ + K1_ + k_; Pm[i1, i2] = (-1) ** k_; Pm[i2, i1] = (-1) ** k_
    lp, Vp = mp.eigsy(Pm)
    ev_idx = [i for i in range(N) if lp[i] > 0]; od_idx = [i for i in range(N) if lp[i] < 0]
    def restrict(M, idx):
        V_ = mp.matrix(N, len(idx))
        for a_, i_ in enumerate(idx):
            for r_ in range(N): V_[r_, a_] = Vp[r_, i_]
        return V_.T * M * V_, V_
    diag = {}
    for name, idx in (("even", ev_idx), ("odd", od_idx)):
        Tr, V_ = restrict(T, idx); Qr = V_.T * Q * V_
        lt, Vt = mp.eigsy(Tr); lt = sorted([lt[i] for i in range(len(idx))])
        lq_ = mp.eigsy(Qr, eigvals_only=True); lq_ = sorted([lq_[i] for i in range(len(idx))])
        diag[name] = dict(T_negatives=sum(1 for x in lt if x < 0), T_lowest=[mp.nstr(x, 8) for x in lt[:3]],
                          Q_lowest=[mp.nstr(x, 10) for x in lq_[:3]])
    # Q without the newest prime power on the full window
    if S.terms:
        nn = S.terms[-1][0]; Qwo = Q + tomp(S.parts[f"prime_{nn}"])
        lw = mp.eigsy(Qwo, eigvals_only=True); diag["Q_without_newest_prime_min"] = mp.nstr(min(lw), 8)
    # pole+arch only (no primes at all) on the full window
    lpa = mp.eigsy(tomp(S.parts["pole"]) + tomp(S.parts["arch"]), eigvals_only=True)
    diag["pole_plus_arch_min"] = mp.nstr(min(lpa), 8)
    out["structure"] = diag
    out["Q_matrix_dps40"] = [[mp.nstr(Q[i, j], 40) for j in range(N)] for i in range(N)]
    print(json.dumps({k: v for k, v in out.items() if k != "Q_matrix_dps40"}), flush=True)
    print(json.dumps(out), flush=True)
    return out

if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--h"); ap.add_argument("--H"); ap.add_argument("--K0", type=int, default=36)
    ap.add_argument("--K1", type=int, default=14); ap.add_argument("--prec", type=int, default=192)
    ap.add_argument("--mu", type=int, default=64); ap.add_argument("--label", default="")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    ctx.prec = a.prec
    def parse(s):
        if s.startswith("log"):  # logN/2
            n = int(s[3:].split("/")[0]); return arb(n).log() / 2
        if "/" in s:
            p, q = s.split("/"); return arb(int(p)) / int(q)
        return arb(s)
    S = Setup(parse(a.h), parse(a.H), a.K0, a.K1, a.prec, a.mu)
    S.build()
    import pickle
    dump = dict(h=str(S.h), H=str(S.H), K0=S.K0, pieces=[(str(c), str(w), K) for c, w, K in S.pieces],
                terms=[n for n, _, _ in S.terms],
                mats={k: [[M[i, j].mid().str(80, radius=False) for j in range(M.ncols())] for i in range(M.nrows())]
                      for k, M in dict(Q=S.Q, G0=S.G0, cvec=S.cvec, svec=S.svec, **S.parts).items()})
    if a.out: pickle.dump(dump, open(a.out.replace('.json', '.pkl'), 'wb'))
    o = analyse(S, a.label)
    if a.out: json.dump(o, open(a.out, "w"), indent=1)
