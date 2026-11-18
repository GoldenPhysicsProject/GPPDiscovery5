"""Round 15 part (b): full Hasse Gram  G[(a,i),(b,j)] = W(tau_{s_a} phi_i, tau_{s_b} phi_j) = B_{s_b - s_a}[i,j]
on a shift set S (Z-analogue of span{1, F, F^2, ...}), base window l.  Families:
  P2/P3/P5 : S = {k log p, k = 0..K}   (one Frobenius p: block-Toeplitz, literal Z[F_p])
  INT      : S = {log n, n = 1..K}      (identity + all integer dilations: the multiplicative monoid N)
  PP       : S = {0} U {log p^k <= log K} (identity + prime powers)
Generalised eigenproblem G v = lam N v, N[(a,i),(b,j)] = <tau_{s_a}phi_i, tau_{s_b}phi_j>; N-null directions (eig < 1e-9 max) dropped.
Ghost = eigenvalue < -1e-7 max|G| (14B convention). Span lies in width L = l + max S, so ghosts <= 14B ghosts at L (Cauchy interlacing).
Usage: python3 hasse_gram.py MODEL FAMILY Kmax [l]"""
import sys, json, math, numpy as np, hasse_core as hc
name, fam, Kmax = sys.argv[1], sys.argv[2], int(sys.argv[3]); l = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
fmax = 125.0 if name == 'DH' else 45.0
mo = hc.MODELS[name](); M = hc.Mof(l, fmax); cache = {}
def B(s):
    k = round(s, 12)
    if k not in cache:
        if k < 0: cache[k] = dict(B=B(-s)['B'].T, N=B(-s)['N'].T)
        else: D = hc.Bs(mo, l, M, s); cache[k] = dict(B=D['B'], N=D['N'])
    return cache[k]
def S_of(K):
    if fam.startswith('P') and fam[1:].isdigit(): p = int(fam[1:]); return [k*math.log(p) for k in range(K+1)]
    if fam == 'INT': return [math.log(n) for n in range(1, K+1)]
    if fam == 'PP':
        out = [0.0]
        for n in range(2, K+1):
            m = n; q = next(q for q in range(2, n+1) if n % q == 0)
            while m % q == 0: m //= q
            if m == 1: out.append(math.log(n))
        return out
rows = []
Smax = math.log(mo['n'].max()) - l - 0.05
for K in range(1, Kmax+1):
    S = S_of(K)
    if S[-1] > Smax: break
    if rows and len(S) == rows[-1]['dim']//M: continue
    n = len(S); G = np.zeros((n*M, n*M)); N = np.zeros_like(G)
    for a in range(n):
        for b in range(n):
            D = B(S[b]-S[a]); G[a*M:(a+1)*M, b*M:(b+1)*M] = D['B']; N[a*M:(a+1)*M, b*M:(b+1)*M] = D['N']
    G = 0.5*(G+G.T); N = 0.5*(N+N.T)
    en, U = np.linalg.eigh(N); keep = en > 1e-9*en.max(); U = U[:, keep]/np.sqrt(en[keep])
    ev = np.linalg.eigvalsh(U.T@G@U); sc = np.abs(G).max()
    gh = int(np.sum(ev < -1e-7*sc))
    rows.append(dict(K=K, dim=n*M, rank=int(keep.sum()), width=float(l+S[-1]), lam_min=float(ev[0]), n_ghosts=gh, scale=float(sc)))
    print(name, fam, K, n*M, rows[-1]['rank'], round(l+S[-1], 3), '%.3e' % ev[0], gh, flush=True)
    json.dump(dict(model=name, family=fam, l=l, M=M, fmax=fmax, rows=rows), open(f'gram_{name}_{fam}_l{l}.json', 'w'), indent=1)
json.dump(dict(model=name, family=fam, l=l, M=M, fmax=fmax, rows=rows), open(f'gram_{name}_{fam}_l{l}.json', 'w'), indent=1)
