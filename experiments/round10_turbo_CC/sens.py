import numpy as np, json
exec(open('assembly_identity.py').read().split("res = {}")[0])
lam = 5.24628; res = {}
def ok(m):
    M = A - m*CC/lam
    return ydPD(M) > 0 and schurE(M) > 1.3e-8
for k in [1, 2, 3, 4, 5, 6, 8]:
    M = A - k*CC/lam; res[f'CC_x{k}'] = dict(schur=schurE(M), YdPD=ydPD(M))
lo, hi = 1.0, 8.0
for _ in range(40):
    m = (lo+hi)/2
    (lo, hi) = (m, hi) if ok(m) else (lo, m)
res['max_multiplier_keeping_beta2'] = lo; res['lamZ_breakeven'] = lam/lo
print(json.dumps(res, indent=1)); json.dump(res, open('sens_results.json', 'w'), indent=1)
