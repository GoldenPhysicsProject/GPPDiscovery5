"""Queen S1c, step (c): Ritz eigenvalues of R = Q - pole (= A - P) for zeta, even and odd sectors, on the
S87 high-precision window forms (s87/hp_forms.py; orthonormal shifted Legendre on [0,L] in s = log x, parity-extended;
half-width L, primes n <= e^{2L}). Quantities: lam1(R_e) (<0), lam2(R_e), lam1(R_o); also lam_min(Q_e), lam_min(Q_o).
Ritz values in nested polynomial spaces are non-increasing in N (Cauchy interlacing), so each value is an upper bound
for the operator eigenvalue (min-max). HIGH-PRECISION NUMERICAL (mpmath, not interval). RH OPEN."""
import sys, json, time, mpmath as mp
sys.path.insert(0, "/workspace/queen/GPPDiscovery5/s87"); import hp_forms as H
def task(a):
    L, N, dps = a
    t0 = time.time(); H.setup(dps); H._GN.clear()
    d = H.build('ZETA', mp.mpf(L), N)
    Re = mp.matrix(H.add(d["Ae"], d["Pe"], 1, -1)); Ro = mp.matrix(H.add(d["Ao"], d["Po"], 1, -1))
    pe = mp.matrix(d["pe"][0]); po = mp.matrix(d["po"][0])
    Qe = Re + pe * pe.T; Qo = Ro - po * po.T
    Ee, Ve = mp.eigsy(Re); Eo = mp.eigsy(Ro, eigvals_only=True)
    eQe = mp.eigsy(Qe, eigvals_only=True); eQo = mp.eigsy(Qo, eigvals_only=True)
    ie = sorted(range(N), key=lambda i: Ee[i])
    se = [Ee[i] for i in ie]; so = sorted(Eo[i] for i in range(N))
    v2 = [mp.nstr(Ve[j, ie[1]], 40) for j in range(N)]
    nrm = max(abs(x) for x in se)
    return dict(L=L, N=N, dps=dps, lam1e=mp.nstr(se[0], 30), lam2e=mp.nstr(se[1], 30), lam3e=mp.nstr(se[2], 30),
                lam1o=mp.nstr(so[0], 30), lam2o=mp.nstr(so[1], 30), qmin_e=mp.nstr(min(eQe), 30), qmin_o=mp.nstr(min(eQo), 30),
                normR=mp.nstr(nrm, 10), v2e=v2, secs=round(time.time() - t0, 1))
if __name__ == "__main__":
    from multiprocessing import Pool
    Ls = [float(x) for x in sys.argv[1].split(",")]; Ns = [int(x) for x in sys.argv[2].split(",")]; dps = int(sys.argv[3]); out = sys.argv[4]
    tasks = sorted([(L, N, dps) for L in Ls for N in Ns], key=lambda t: -t[1])
    res = []
    with Pool(int(sys.argv[5]) if len(sys.argv) > 5 else 7) as pool:
        for r in pool.imap_unordered(task, tasks):
            res.append(r); print({k: v for k, v in r.items() if k != "v2e"}, flush=True)
            json.dump(res, open(out, "w"))
