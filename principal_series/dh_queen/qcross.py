exec(open('parity_split.py').read().split('out = []')[0])
for L in [1.85, 1.9, 1.95]:
    T = Q(L, lambda u: S(u, L), c, Om, nu); print(L, T['Q_over_norm'], flush=True)
