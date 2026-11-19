"""Round 19, N2: Hamburger rigidity, numerically (mpmath).
Unknowns a_0..a_N (real). Class H_1 (conductor 1, Gamma_R(s), even): theta(x)=a_0+2 sum a_n e^{-pi n^2 x},
constraint theta(1/x) = sqrt(x) theta(x) on x in [xmin,1]. Hamburger: solution space = R * (1,1,1,...).
Class H_5^odd (conductor 5, Gamma_R(s+1)): theta(x)=sum a_n n e^{-pi n^2 x/5}, constraint theta(1/x)=x^{3/2} theta(x).
Kaczorowski-Perelli: fixed-conductor degree-1 extended Selberg class is finite-dimensional (combinations of L(s,chi)). We report singular values of the column-normalized collocation matrix
and the null vector(s). Truncation N is matched so the neglected tail is < 10^-(dps/2).
NUMERICAL. No zeros used."""
import mpmath as mp, json
mp.mp.dps = 120
def colloc(q, N, xs, odd):
    rows = []
    for x in xs:
        r = []
        for n in range(0 if not odd else 1, N+1):
            if odd:
                r.append(n*(mp.e**(-mp.pi*n*n/(q*x)) - x**1.5*mp.e**(-mp.pi*n*n*x/q)))
            else:
                w = 1 if n == 0 else 2
                r.append(w*(mp.e**(-mp.pi*n*n/x) - mp.sqrt(x)*mp.e**(-mp.pi*n*n*x)))
        rows.append(r)
    A = mp.matrix(rows)
    sc = [mp.sqrt(sum(A[i, j]**2 for i in range(A.rows))) for j in range(A.cols)]
    for j in range(A.cols):
        for i in range(A.rows): A[i, j] /= sc[j]
    return A, sc
def run(q, N, odd, xmin, M=60):
    xs = [mp.mpf(xmin) + (1-mp.mpf(xmin))*(k+mp.mpf(1)/2)/M for k in range(M)]
    A, sc = colloc(q, N, xs, odd)
    U, S, V = mp.svd_r(A)
    S = sorted([S[i] for i in range(len(S))], reverse=True)
    # null vector = last row of V (rows of V are right singular vectors), undo scaling
    idx = min(range(len(S)), key=lambda i: S[i])
    U2, S2, V2 = mp.svd_r(A)
    k = min(range(len(S2)), key=lambda i: S2[i])
    v = [V2[k, j]/sc[j] for j in range(A.cols)]
    v = [vi/v[0] for vi in v]
    return [float(mp.log10(s)) for s in S], [float(vi) for vi in v]
out = {}
xi = (mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
dh = {0: 0, 1: 1, 2: xi, 3: -xi, 4: -1}
# H_1 even, conductor 1
for N in (6, 8, 10, 12):
    xmin = 0.35
    lg, v = run(1, N, False, xmin)
    out[f'H1_N{N}'] = dict(log10_singular_values=lg, null_vector_normalized_a0=v)
    print('H1', N, [round(x, 1) for x in lg], [round(x, 12) for x in v[:6]], flush=True)
# H_5 odd, conductor 5 (matched truncation N5 ~ sqrt5 * N1)
import math
for N in (13, 18, 22, 27, 34):
    lg, v = run(5, N, True, 0.35, M=90)
    dhv = [float(dh[n % 5]) for n in range(1, N+1)]
    err = max(abs(a-b) for a, b in zip(v[:min(N, 12)], dhv[:min(N, 12)]))
    out[f'H5odd_N{N}'] = dict(log10_singular_values=lg, null_vector_first12=v[:12], max_dev_from_DH_first12=err)
    print('  null vector vs DH, max dev on n<=12:', err)
    print('H5', N, [round(x, 1) for x in lg], flush=True)
    out[f'H5odd_N{N}']['count_below_1e-40'] = sum(1 for x in lg if x < -40)
# control: odd class at conductor 1 (no odd character mod 1): expect NO isolated null direction
for N in (8, 12):
    lg, v = run(1, N, True, 0.35)
    out[f'H1odd_N{N}'] = dict(log10_singular_values=lg); print('H1odd', N, [round(x, 1) for x in lg], flush=True)
# DH coefficients check in H_5^odd
def thDH(x): return mp.fsum(dh[n % 5]*n*mp.e**(-mp.pi*n*n*x/5) for n in range(1, 400))
res = [float(abs(thDH(1/x) - x**1.5*thDH(x))/abs(thDH(1/x))) for x in (mp.mpf('0.4'), mp.mpf('0.7'), mp.mpf('0.9'))]
out['DH_in_H5odd_relative_residual'] = res; print('DH residual', res)
json.dump(out, open('r19_hamburger.json', 'w'), indent=1)
