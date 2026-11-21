"""Round 32: discrete Beurling (DMV-type) control, exactly the Round 13e harness construction (primes >= 20 tracking
dPsi_C = (1 - sum_q x^{rho_q - 1}) dx, quartet rho_q = 1/2 +- 0.1 +- 28i), to NMAX (default 1e8).
Outputs: Moebius-floor accumulators (via r32_beur_dfs.c), Lambda_B accumulators, |Pi_B - Pi_C| check. No zeta zeros."""
import sys, math, json, subprocess, numpy as np
NMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 1e8
H = 0.005; x0 = 20.0
q = [0.6+28j, 0.6-28j, 0.4+28j, 0.4-28j]
lx = np.linspace(math.log(x0), math.log(NMAX)+0.01, 2_000_000); x = np.exp(lx)
dens = 1 - sum((x**complex(r-1)).real for r in q); assert dens.min() > 0
dPi = dens/lx*x
Pi = np.concatenate([[0], np.cumsum(0.5*(dPi[1:]+dPi[:-1])*np.diff(lx))])
pb = Pi.copy()
for it in range(8):
    corr = np.zeros_like(Pi); m = 2
    while lx[-1]/m >= lx[0]:
        corr += np.floor(np.interp(lx/m, lx, pb, left=0.0)+0.5)/m; m += 1
    pb = Pi - corr
cnt = np.floor(pb+0.5); cnt = np.maximum.accumulate(np.maximum(cnt, 0))
jumps = np.nonzero(np.diff(cnt) > 0)[0]+1
primes = []
for j in jumps:
    for lev in np.arange(cnt[j-1]+1, cnt[j]+1):
        a, b = pb[j-1], pb[j]; th = 0.5 if b == a else min(max((lev-0.5-a)/(b-a), 0.0), 1.0)
        primes.append(math.exp(lx[j-1] + th*(lx[j]-lx[j-1])))
P = np.sort(np.array(primes)); P = P[P <= NMAX]
# Pi_B (prime powers weighted 1/m) vs Pi_C on the grid
PiB = np.zeros_like(lx); m = 1
while True:
    pm = P**m; pm = pm[pm <= NMAX]
    if not len(pm): break
    PiB += np.searchsorted(np.log(pm), lx, side='right')/m; m += 1
dev = (PiB - Pi)[lx <= math.log(NMAX)]
P.astype(np.float64).tofile('beur_primes.bin')
out = subprocess.run(['./r32_beur_dfs', 'beur_primes.bin', repr(NMAX), repr(H), 'beur_bins.bin'], capture_output=True, text=True).stdout
print(out, flush=True)
A = np.fromfile('beur_bins.bin').reshape(3, -1)
nb = A.shape[1]
# Lambda_B binned
lam_n, lam_c = [], []
for p in P:
    mm = p
    while mm <= NMAX: lam_n.append(mm); lam_c.append(math.log(p)); mm *= p
lam_n = np.array(lam_n); lam_c = np.array(lam_c)
j = np.ceil(np.log(lam_n)/H - 1e-12).astype(int)
L0 = np.bincount(j, weights=lam_c*lam_n**-0.5, minlength=nb); L1 = np.bincount(j, weights=lam_c*lam_n**-0.5*np.log(lam_n), minlength=nb)
res = dict(NMAX=NMAX, H=H, nprimes=int(len(P)), first_primes=P[:8].tolist(), dfs=out.strip(),
           PiB_minus_PiC=dict(max_abs=float(np.abs(dev).max()), mean=float(dev.mean()), at_end=float(dev[-1])),
           L=(np.arange(nb)*H).tolist(), A0=np.cumsum(A[0]).tolist(), A1=np.cumsum(A[1]).tolist(), A2=np.cumsum(A[2]).tolist(),
           P0=np.cumsum(L0[:nb]).tolist(), P1=np.cumsum(L1[:nb]).tolist())
json.dump(res, open('r32_beurling.json', 'w'))
print('nprimes', len(P), 'max|PiB-PiC|', float(np.abs(dev).max()), flush=True)
