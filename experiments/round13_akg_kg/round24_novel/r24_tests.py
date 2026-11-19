"""Round 24 numerical tests (integers + gamma factors only; zeros appear ONLY in rows marked 'illustration').
Usage: python3 r24_tests.py PART   with PART in {audit, exponent, cells, levy, cramer}."""
import sys, json, math, time, numpy as np, mpmath as mp
import r22_models as RM
from r22_fejer import rho
from r22_models import lam_table

def DM(name, L, lam):
    """D_M(L) = (pi L/2) rho_M(e^L, 0): the gamma = 0 slice, explicit formula (R22/R23)."""
    if L == 0: return 0.0
    return float(np.pi*L/2*rho(name, math.exp(L), np.array([0.0]), lam)[0])

def set_theta(th):
    RM.TH = th
    RM.MODELS['FTH'] = dict(logQ=-math.log(math.pi), arch=[(0.5, 0.25 - th/2, 0.5), (0.5, 0.25 + th/2, 0.5)], poles=[1 + th, 1 - th])

out = {}
part = sys.argv[1]
if part == 'audit':
    # paper's claim "local terms always >= 0": split zeta's D(L) = pole + arch - prime
    mp.mp.dps = 30
    c1 = float((mp.digamma(0.25) - mp.log(mp.pi))/2); C2 = float(mp.psi(1, 0.25)/4)
    lam = lam_table('ZETA', 10**6 + 1)
    rows = []
    for L in [0.3, 0.69, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 13.8]:
        k = np.arange(0, 400); T = float(np.sum(np.exp(-(2*k + 0.5)*L)/(2*k + 0.5)**2))
        pole = 4*(math.exp(L/2) - 2 + math.exp(-L/2)); arch = c1*L + C2 - T
        X = math.exp(L); n = np.arange(2, int(X) + 1); lm = lam[2:int(X) + 1]; nz = lm != 0
        S = float(np.sum(lm[nz]*n[nz]**-0.5*(L - np.log(n[nz]))))
        rows.append(dict(L=L, pole=pole, arch=arch, prime=-S, total=pole + arch - S))
        print(rows[-1], flush=True)
    out['split'] = rows
elif part == 'exponent':
    # Theorem 4 sharpness on F_theta: negative envelope of D_FTH grows like e^{theta L}
    res = {}
    for th in (0.0, 0.125, 0.25):
        name = 'ZETA' if th == 0 else 'FTH'
        if th: set_theta(th)
        lam = lam_table(name, int(math.exp(11.6)) + 2)
        Ls = np.arange(3.0, 11.5, 0.005); D = np.array([DM(name, L, lam) for L in Ls])
        cen, amp, mins = [], [], []
        for a in np.arange(3.0, 11.0, 1.0):
            w = (Ls >= a) & (Ls < a + 1)
            cen.append(a + 0.5); amp.append(math.log((D[w].max() - D[w].min())/2)); mins.append(float(D[w].min()))
        slope = float(np.polyfit(cen, amp, 1)[0])
        res[str(th)] = dict(model=name, windows=cen, log_half_range=amp, window_min=mins, fitted_exponent=slope)
        print(th, res[str(th)], flush=True)
    out['FTH_exponent'] = res
elif part == 'cells':
    # N x N additive-cell Weil Gram in multiplicative coordinates
    Nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    res = {}
    for name in ('ZETA', 'FTH', 'EPS', 'DH'):
        if name == 'FTH': set_theta(0.25)
        t0 = time.time(); lam = lam_table(name, Nmax + 3)
        cache = {}
        def Dr(p, q):          # D(|log(p/q)|)
            if p == q: return 0.0
            if p < q: p, q = q, p
            g = math.gcd(p, q); key = (p//g, q//g)
            if key not in cache: cache[key] = DM(name, math.log(key[0]/key[1]), lam)
            return cache[key]
        n = np.arange(1, Nmax + 1)
        M = np.zeros((Nmax, Nmax))
        for i in range(Nmax):
            a = i + 1
            for j in range(i, Nmax):
                b = j + 1
                v = Dr(b + 1, a) + Dr(b, a + 1) - Dr(b, a) - Dr(b + 1, a + 1)
                M[i, j] = M[j, i] = v
        d = np.sqrt(np.diag(M)); C = M/np.outer(d, d)
        rows = {}
        first = None
        for N in [10, 20, 30, 50, 75, 100, 150, 200, 300, 400, 600, 800]:
            if N > Nmax: break
            w = np.linalg.eigvalsh(C[:N, :N]); rows[N] = float(w[0])
        lo, hi = 2, Nmax
        if np.linalg.eigvalsh(C)[0] < -1e-9:
            while lo < hi:
                mid = (lo + hi)//2
                if np.linalg.eigvalsh(C[:mid, :mid])[0] < -1e-9: hi = mid
                else: lo = mid + 1
            first = lo
        res[name] = dict(min_eig_corr=rows, first_indefinite_N=first, diag_min=float(d.min()**2), n_D_evals=len(cache), secs=time.time() - t0)
        print(name, res[name], flush=True)
    out['cells'] = res
elif part == 'levy':
    # e^{-tD} Toeplitz (Schoenberg) and Schur powers G^{o k} (tensor powers) on the R23 grid
    Lmax = math.log(1e5); h = 0.05
    res = {}
    for name in ('ZETA', 'FTH', 'EPS', 'DH'):
        if name == 'FTH': set_theta(0.25)
        lam = lam_table(name, 10**5 + 2)
        allL = np.arange(0, int(Lmax/h) + 1)*h
        Dv = np.array([DM(name, L, lam) for L in allL])
        n = len(allL) - 1; I = np.arange(1, n + 1)
        G = Dv[I][:, None] + Dv[I][None, :] - Dv[np.abs(I[:, None] - I[None, :])]
        r = {}
        for k in (1, 2, 3, 4):
            Gk = G**k; w = np.linalg.eigvalsh(Gk); r[f'schur_k{k}_min_over_max'] = float(w[0]/w[-1])
        J = np.arange(0, n + 1)
        for t in (0.3, 1.0, 3.0, 10.0, 30.0):
            E = np.exp(-t*Dv[np.abs(J[:, None] - J[None, :])]); w = np.linalg.eigvalsh(E)
            r[f'expD_t{t}_min_eig'] = float(w[0])
        res[name] = r; print(name, r, flush=True)
    out['levy'] = res
elif part == 'cramer':
    # candidate Omega from the integers: E(x) = e^{-x/2}(psi(e^x) - main(e^x)); Bohr mean of E^2 and autocorrelation
    res = {}
    for name, N in (('ZETA', 10**8), ('FTH', 10**6), ('EPS', 10**6), ('DH', 10**6)):
        if name == 'FTH': set_theta(0.25)
        t0 = time.time()
        if name == 'ZETA':
            s = np.ones(N + 1, bool); s[:2] = False
            for p in range(2, int(N**0.5) + 1):
                if s[p]: s[p*p::p] = False
            ps = np.nonzero(s)[0]; pts = [ps]; lv = [np.log(ps)]
            for p in ps[ps <= int(N**0.5)]:
                q = int(p)**2
                while q <= N: pts.append(np.array([q])); lv.append(np.array([math.log(int(p))])); q *= int(p)
            pos = np.concatenate(pts); val = np.concatenate(lv); o = np.argsort(pos); pos = pos[o]; val = val[o]
        else:
            lam = lam_table(name, N + 1); pos = np.nonzero(lam)[0]; pos = pos[pos >= 2]; val = lam[pos]
        cum = np.cumsum(val)
        xmax = math.log(N); dx = 1e-4; xs = np.arange(0, xmax, dx); ts = np.exp(xs)
        idx = np.searchsorted(pos, ts, side='right'); psi = np.where(idx > 0, cum[np.maximum(idx - 1, 0)], 0.0)
        if name in ('ZETA', 'EPS'): main = ts
        elif name == 'FTH': main = ts**1.25/1.25 + ts**0.75/0.75
        else: main = 0*ts
        E = np.exp(-xs/2)*(psi - main)
        means = {}
        for X in [5, 7, 9, 11, 13, 15, 17, 18.4]:
            if X <= xmax: w = xs < X; means[str(X)] = float(np.mean(E[w]**2))
        # autocorrelation over the last part of the range (Bohr mean estimate), y in [0, 30]
        x0 = 4.0; w = xs >= x0; Ew = E[w]
        ys = np.arange(0, 30.0, 0.05); ac = []
        for y in ys:
            k = int(round(y/dx)); ac.append(float(np.mean(Ew[:len(Ew) - k]*Ew[k:])) if k < len(Ew) - 1000 else None)
        res[name] = dict(N=N, B2_mean_by_X=means, autocorr_y=ys[::20].tolist(), autocorr=ac[::20], secs=time.time() - t0)
        if name == 'ZETA':
            gam = np.array([float(mp.zetazero(k).imag) for k in range(1, 401)])
            pred = [float(np.sum(2*np.cos(gam*y)/(0.25 + gam**2))) for y in ys[::20]]
            res[name]['zero_side_illustration_400'] = pred
            res[name]['sum_1_over_abs_rho2_exact'] = float(2 + mp.euler - mp.log(4*mp.pi))
        print(name, {k: v for k, v in res[name].items() if k in ('B2_mean_by_X', 'secs', 'sum_1_over_abs_rho2_exact')}, flush=True)
    out['cramer'] = res
json.dump(out, open(f'r24_{part}.json', 'w'), indent=1)
