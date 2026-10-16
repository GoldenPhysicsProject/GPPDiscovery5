# Grok S51 H51b(2): finite-N Hankel detection depth for a planted off-line quartet (NUMERICAL).
# H_N[i,j] = sum_nodes x^{i+j+2} + 2 Re w^{i+j+2},  x=gamma^-2 (on-line), w = 1/(-(delta+i*gamma)^2) (planted pair).
# Planted quartet REPLACES on-line zero k. Detection depth N* = first N with a negative LDL pivot
# (= first N at which H_N is not PD; #neg eigenvalues of H_N = #neg pivots d_1..d_N, Jacobi/Sylvester).
# Tail beyond the K listed nodes is omitted (it is PSD, so omission can only make detection EASIER:
# reported N* are lower bounds for the true-tail problem).
import sys, mpmath as mp
dps = int(sys.argv[1]) if len(sys.argv) > 1 else 300
mode = sys.argv[2] if len(sys.argv) > 2 else 'zeta'
K = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
NMAX = int(sys.argv[4]) if len(sys.argv) > 4 else 70
mp.mp.dps = dps
if mode == 'zeta':
    gam = [mp.mpf(l.split()[0]) for l in open('zeros_arb_3000.txt')][:K]
else:  # lattice model gamma_j = c*j
    c = mp.mpf(mode.split(':')[1]); gam = [c*j for j in range(1, K+1)]
s = gam[0]**2   # congruence scaling y = x*s (signature unchanged)
ys = [s/g**2 for g in gam]
P = 2*NMAX + 2
pw = [mp.fsum(y**p for y in ys) for p in range(P+1)]
def pivots(entries, N):
    A = [[entries[i+j+2] for j in range(N)] for i in range(N)]
    d = []
    for k in range(N):
        piv = A[k][k]; d.append(piv)
        if piv == 0: break
        for i in range(k+1, N):
            f = A[i][k]/piv
            for j in range(k+1, N): A[i][j] -= f*A[k][j]
    return d
def depth(k, delta):
    g = gam[k-1]
    w = s/(-(mp.mpc(delta, g))**2)
    e = [pw[p] - ys[k-1]**p + 2*mp.re(w**p) for p in range(P+1)]
    d = pivots(e, NMAX)
    for n, piv in enumerate(d, 1):
        if piv < 0: return n, piv, g, w
    return None, None, g, w
def spacing(g): return 2*mp.pi/mp.log(g/(2*mp.pi))
targets = [int(t) for t in sys.argv[5].split(',')] if len(sys.argv) > 5 else ([4, 10, 21, 37] if mode == 'zeta' else [8, 15, 25, 40])
DL = sys.argv[6].split(',') if len(sys.argv) > 6 else ['0.05sp','0.2sp','0.45sp','abs0.1','abs0.4']
print('mode', mode, 'dps', dps, 'K', K, 'NMAX', NMAX)
print('k  gamma   N_below+2  spacing | delta : N*  (pivot sign at N*, scaled)')
for k in targets:
    g = gam[k-1]; sp = spacing(g) if mode == 'zeta' else gam[1]-gam[0]
    row = []
    for lab in DL:
        dl = mp.mpf(lab[3:]) if lab.startswith('abs') else mp.mpf(lab[:-2])*sp
        n, piv, _, _ = depth(k, dl)
        flag = '' if dl < 0.5 else '(off-strip model)'
        row.append('%s=%s:%s%s' % (lab, mp.nstr(dl, 3), n, flag))
    print(k, mp.nstr(g, 6), k-1+2, mp.nstr(sp, 4), '|', '  '.join(row), flush=True)
