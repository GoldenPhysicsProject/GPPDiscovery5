"""Round 14C: matched-local-data pair test + Golden Eye falsifiers for prolate/archimedean/semilocal T.
Candidates (all from data at a set of places S, never zeros):
  T_inf : ||T_inf g||^2 := (arch + pole)(g)   [Connes-Consani 2021: archimedean Weil form, prolate/Sonin proof for
          support in [2^-1/2, 2^1/2], i.e. width L <= log 2]; exists as an operator iff arch+pole >= 0.
  T_S   : ||T_S g||^2 := (arch + pole - prime_S)(g), prime_S = atoms p^k with p in S  [semilocal trace formula form].
Identity test: W(g) - ||T g||^2 = exact defect (computed), on the SAME g for both members of a pair.
Pairs: (ZETA, EPS), (CHI5, DH), (ZETA, BEUR2); controls F_theta (theta=0.2)."""
import sys, os, math, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'round14b_euler')); sys.path.insert(0, os.path.join(HERE, '..'))
import weil_galerkin as wg, models14b as md
out = {}
def chi5():
    Z = md.zeta(); D = md.dh(); chi = {0: 0, 1: 1, 2: 1j, 3: -1j, 4: -1}
    re = np.array([chi[int(n) % 5].real if not isinstance(chi[int(n) % 5], int) else chi[int(n) % 5] for n in Z['n']], float)
    m = np.abs(re) > 0
    return dict(name='CHI5', n=Z['n'][m], c=(Z['c']*re)[m], C=D['C'], w=D['w'], poles=[])
def pp_mask(model):  # atoms that are prime powers, and their prime
    n = model['n'].astype(np.int64); P = np.zeros(len(n), np.int64)
    for i, v in enumerate(n):
        for p in range(2, int(v**0.5)+2):
            if v % p == 0:
                q = v
                while q % p == 0: q //= p
                P[i] = p if q == 1 else 0; break
        else: P[i] = v
    return P
def sub(model, mask): return dict(model, n=model['n'][mask], c=model['c'][mask])

Z, E, B, D, F, X = md.zeta(), md.eps(), md.beur2(), md.dh(), md.fth(0.2), chi5()

# 1. T_inf on zeta across L (Golden Eye: break just past log 2?)
rows = []
for L in [0.3, 0.5, 0.6, 0.65, math.log(2)-0.01, math.log(2)+0.05, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0]:
    M = max(24, int(math.ceil(45*L/math.pi))+8)
    a, p, q = wg.operators(Z, L, M); Ainf = a+p; W = Ainf - q
    ev = np.linalg.eigvalsh(Ainf); evW = np.linalg.eigvalsh(W)
    rows.append(dict(L=L, M=M, lam_min_arch_pole=float(ev[0]), lam_min_W=float(evW[0]), defect_norm=float(np.abs(q).max()),
                     defect_is_prime_term=bool(np.allclose(W-Ainf, -q)), atoms=[int(n) for n in Z['n'][np.log(Z['n']) < L]]))
    print('Tinf', rows[-1], flush=True)
out['T_inf_zeta'] = rows

# 2. pair tests on the fake's ghost vectors (same g for both members)
def pair(m1, m2, L, fmax, tag):
    M = int(math.ceil(fmax*L/math.pi))+8
    a1, p1, q1 = wg.operators(m1, L, M); a2, p2, q2 = wg.operators(m2, L, M)
    W1, W2 = a1+p1-q1, a2+p2-q2
    ev, V = np.linalg.eigh(W2); neg = np.nonzero(ev < -1e-7*np.abs(W2).max())[0]
    res = dict(L=L, M=M, T_inf_identical=bool(np.allclose(a1+p1, a2+p2, atol=1e-12)),
               lam_min_Tinf=[float(np.linalg.eigvalsh(a1+p1)[0]), float(np.linalg.eigvalsh(a2+p2)[0])], ghosts=[])
    for k in neg[:6]:
        v = V[:, k]
        res['ghosts'].append(dict(W_member1=float(v@W1@v), W_member2=float(v@W2@v), Tinf1=float(v@(a1+p1)@v), Tinf2=float(v@(a2+p2)@v),
                                  d_arch=float(v@(a2-a1)@v), d_pole=float(v@(p2-p1)@v), d_prime=float(-v@(q2-q1)@v)))
    out.setdefault('pairs', {}).setdefault(tag, []).append(res); print(tag, res, flush=True)
    return res, V[:, neg], (a1, p1, q1, a2, p2, q2)
for L in [4.0, 6.0]: pair(Z, E, L, 45.0, 'ZETA_vs_EPS')
for L in [4.5, 6.0]: pair(Z, B, L, 45.0, 'ZETA_vs_BEUR2')
pair(Z, F, 2.0, 45.0, 'ZETA_vs_FTH0.2')
# CHI5 vs DH, with the DH prime-side difference split by atom type
for L in [4.0, 6.0]:
    res, Vg, (a1, p1, q1, a2, p2, q2) = pair(X, D, L, 125.0, 'CHI5_vs_DH')
    M = Vg.shape[0]
    Pd = pp_mask(D); Px = pp_mask(X)
    groups = {'DH_atoms_5|n': D['n'].astype(np.int64) % 5 == 0, 'DH_atoms_prime_power_p!=5': (Pd > 0) & (Pd != 5),
              'DH_atoms_composite_not_pp_5!|n': (Pd == 0) & (D['n'].astype(np.int64) % 5 != 0)}
    split = {}
    for gname, msk in groups.items():
        qq = wg.operators(dict(sub(D, msk), C=None, poles=[]), L, M)[2]
        split[gname] = [float(-v@qq@v) for v in Vg.T[:4]]
    split['CHI5_atoms(all prime powers)'] = [float(-v@q1@v) for v in Vg.T[:4]]
    out['pairs']['CHI5_vs_DH'][-1]['prime_term_split_on_ghosts'] = split; print(split, flush=True)

# 3. Golden Eye: unchanged construction on EPS L=4 and F_theta L=2, M ~ 60-100: identity W = ||T g||^2 must fail
ge = {}
for mo, L, M in [(E, 4.0, 66), (E, 4.0, 100), (F, 2.0, 60), (F, 2.0, 100)]:
    a, p, q = wg.operators(mo, L, M); evA = np.linalg.eigvalsh(a+p); evW = np.linalg.eigvalsh(a+p-q)
    ge[f"{mo['name']}_L{L}_M{M}"] = dict(lam_min_Tinf=float(evA[0]), lam_min_W=float(evW[0]), n_neg_W=int(np.sum(evW < -1e-7*np.abs(a+p-q).max())))
    print(ge, flush=True)
out['golden_eye'] = ge

# 4. semilocal T_S on zeta: S = {inf} u {p <= P}; T_S exists iff W_S >= 0; defect = omitted atoms
sl = []
for L in [2.0, 3.0, 4.0]:
    M = int(math.ceil(45*L/math.pi))+8; Pz = pp_mask(Z)
    for Pmax in [2, 3, 5, 7, 11, 23, 53]:
        msk = (Pz > 0) & (Pz <= Pmax)
        a, p, qS = wg.operators(sub(Z, msk), L, M)
        ev = np.linalg.eigvalsh(a+p-qS)
        sl.append(dict(L=L, Pmax=Pmax, lam_min_W_S=float(ev[0]), n_neg=int(np.sum(ev < -1e-7*np.abs(a+p-qS).max())),
                       S_complete=bool(all(p <= Pmax for p in range(2, int(math.exp(L))+1) if all(p % r for r in range(2, int(p**0.5)+1))))))
        print('semilocal', sl[-1], flush=True)
out['semilocal_zeta'] = sl
json.dump(out, open(os.path.join(HERE, 'pairtest.json'), 'w'), indent=1)
