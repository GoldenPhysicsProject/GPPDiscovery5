"""Round 7 float study (NUMERICAL): bordered Schur onto E_12 with Y = E_{12+K} ⊕ {delta_l}; extra loss from modes beyond 12+K."""
import numpy as np, json, os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
src = open('p2_boundary_float.py').read().split("out = dict")[0].replace("../round5_turbo","../round5_turbo")
exec(src)
N = 12
def gram(J, with_delta):
    F = [S(j) for j in J]; s = [sig(j) for j in J]
    if with_delta: F.append(Fd); s.append(sd)
    F = np.array(F); s = np.array(s)
    return (F*meas) @ F.T - 2*bI2**2*np.outer(s, s)
def schur(M):
    M11, M12, M22 = M[:N,:N], M[:N,N:], M[N:,N:]
    return M11 - M12 @ np.linalg.solve(M22, M12.T), np.linalg.eigvalsh(M22)[0]
res = {}
for Kh in [40, 60, 100, 200, 400]:
    J = list(range(1, N+Kh+1))
    for wd in [False, True]:
        M = gram(J, wd); Sc, l22 = schur(M)
        res[f"K{Kh}_delta{int(wd)}"] = dict(lam_min_schur=float(np.linalg.eigvalsh(Sc)[0]), lam_min_D=float(l22))
        print(Kh, wd, res[f"K{Kh}_delta{int(wd)}"], flush=True)
json.dump(res, open('r7_float_results.json','w'), indent=1)
