"""Round 21, Task 3: 'the zero is the connections passing through the boundary'.
(a) Boundary-flux split (Riemann 1859):  M(s) = A(s) + A(1-s),  A(s) = int_1^oo x^{s-1} F(x) dx  (flux of the x > 1 bulk),
    A(1-s) = int_0^1 x^{s-1} F(x) dx (its J-partner, the x < 1 bulk).  Zero <=> A(rho) = -A(1-rho)  (destructive interference
    at the boundary).  On the line 1-s = conj(s), so the partner amplitude is automatically conj A(s) (TFD pairing) and
    T(t) := M(1/2+it) = 2 Re A(1/2+it) is REAL: an on-line zero is a sign flip of the transmitted amplitude.
(b) Turing-type count: number of sign flips of T on (0,T] versus the argument-principle count N(T) of xi_M = s(s-1)M(s)
    (DH: entire, xi = M) in [-1,2] x (0,T].  RH for M up to T  <=>  flips = N(T).
Zeros are NOT inputs: flips and N(T) are computed from M alone.  Two zeros (one zeta, one EPS off-line) are located by Newton
from printed labels ONLY to illustrate (a)."""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from r21_models import *
import mpmath as mp
mp.mp.dps = 20
t0 = time.time(); out = {}
# ---------- (a) boundary flux ----------
def A_zeta(s):
    return mp.quad(lambda x: x**(s - 1) * mp_thplus(x * x), [1, 2, 4, mp.inf]) - 1 / (2 * (1 - s))
def A_eps(s):
    r5 = mp.sqrt(5)
    bulk = lambda t: (1 + 2 * mp_thplus(t / r5)) * (1 + 2 * mp_thplus(r5 * t)) - 1
    return mp.quad(lambda t: t**(s - 1) * bulk(t), [1, 3, 10, 30, mp.inf]) - 1 / (1 - s)
A = {'ZETA': A_zeta, 'EPS': A_eps}; M = {'ZETA': M_zeta, 'EPS': M_eps}
chk = {}
for nm in A:
    rows = []
    for s in (mp.mpc(0.5, 3), mp.mpc(0.3, 10), mp.mpc(0.8, -7)):
        lhs = A[nm](s) + A[nm](1 - s); rows.append(dict(s=mp.nstr(s, 5), err=mp.nstr(abs(lhs - M[nm](s)), 3)))
    chk[nm] = rows
out['flux_split_check'] = chk; print(chk, flush=True)
illus = {}
r1 = mp.zetazero(1)                                       # illustration only
r2 = mp.findroot(ZQ, mp.mpc(0.933, 15.67))                # EPS off-line zero, Newton from Round 17 label (illustration)
for nm, rho in (('ZETA', r1), ('EPS', r2)):
    a, b = A[nm](rho), A[nm](1 - rho)
    illus[nm] = dict(rho=mp.nstr(rho, 12), A_right=mp.nstr(a, 10), A_left=mp.nstr(b, 10),
                     cancellation_abs_sum=mp.nstr(abs(a + b), 3), scale_abs_A=mp.nstr(abs(a), 6),
                     partner_is_conjugate_abs_diff=mp.nstr(abs(b - mp.conj(a)), 6),
                     arg_A_over_pi=mp.nstr(mp.arg(a) / mp.pi, 10))
    print(nm, illus[nm], flush=True)
out['zero_as_boundary_cancellation'] = illus
# ---------- (b) sign flips vs argument principle ----------
def L_dh(s):
    return sum(DH_A[r] * mp.zeta(s, mp.mpf(r) / 5) for r in range(1, 5)) * mp.mpf(5)**(-s)
def M_dh(s): return (mp.mpf(5) / mp.pi)**((s + 1) / 2) * mp.gamma((s + 1) / 2) * L_dh(s) / 2
XI = {'ZETA': lambda s: s * (s - 1) * M_zeta(s), 'EPS': lambda s: s * (s - 1) * M_eps(s), 'DH': M_dh}
TMAX = {'ZETA': 100.0, 'EPS': 50.0, 'DH': 95.0}
def flips(f, T, h):
    k = 0; prev = None; t = h; pts = []
    while t <= T:
        v = mp.re(f(mp.mpc(0.5, t)))
        if prev is not None and (v > 0) != (prev > 0): k += 1; pts.append(round(float(t), 3))
        prev = v; t += h
    return k, pts
def arg_count(f, T):
    path = [mp.mpc(2, 1e-3), mp.mpc(2, T), mp.mpc(-1, T), mp.mpc(-1, 1e-3), mp.mpc(2, 1e-3)]
    tot = mp.mpf(0)
    for z0, z1 in zip(path[:-1], path[1:]):
        L = abs(z1 - z0); m = int(L / 0.05) + 1; zs = [z0 + (z1 - z0) * j / m for j in range(m + 1)]
        vals = [f(z) for z in zs]; j = 0
        while j < len(zs) - 1:
            d = mp.arg(vals[j + 1] / vals[j])
            if abs(d) > 0.4 and abs(zs[j + 1] - zs[j]) > 1e-6:
                zm = (zs[j] + zs[j + 1]) / 2; zs.insert(j + 1, zm); vals.insert(j + 1, f(zm)); continue
            tot += d; j += 1
    return tot / (2 * mp.pi)
res = {}
for nm, f in XI.items():
    T = TMAX[nm]
    # nudge T so the top edge is not near an on-line zero
    while abs(mp.re(f(mp.mpc(0.5, T)))) < 1e-3 * abs(mp.re(f(mp.mpc(0.5, T - 0.5)))): T += 0.137
    fl, pts = flips(f, T, 0.01)
    Nc = arg_count(f, T)
    # real axis: sign of xi on [-1,2] (bottom edge must avoid zeros)
    rs = [mp.re(f(mp.mpf(-1) + 3 * (mp.mpf(j) + 0.5) / 300)) for j in range(300)]
    res[nm] = dict(T=T, sign_flips_on_line=fl, argument_principle_N=mp.nstr(Nc, 8), deficit=round(float(Nc)) - fl,
                   real_axis_sign_changes=sum(1 for a, b in zip(rs, rs[1:]) if (a > 0) != (b > 0)), first_flips=pts[:12])
    print(nm, res[nm], flush=True)
out['turing_count'] = res
# right-of-line count for EPS/DH, as a cross-check of the deficit (off-line pairs rho, 1-conj(rho))
def rect_count(f, s0, s1, T):
    path = [mp.mpc(s1, 0.01), mp.mpc(s1, T), mp.mpc(s0, T), mp.mpc(s0, 0.01), mp.mpc(s1, 0.01)]
    tot = mp.mpf(0)
    for z0, z1 in zip(path[:-1], path[1:]):
        m = int(abs(z1 - z0) / 0.05) + 1; zs = [z0 + (z1 - z0) * j / m for j in range(m + 1)]; vals = [f(z) for z in zs]; j = 0
        while j < len(zs) - 1:
            d = mp.arg(vals[j + 1] / vals[j])
            if abs(d) > 0.4 and abs(zs[j + 1] - zs[j]) > 1e-6:
                zm = (zs[j] + zs[j + 1]) / 2; zs.insert(j + 1, zm); vals.insert(j + 1, f(zm)); continue
            tot += d; j += 1
    return tot / (2 * mp.pi)
out['right_of_line_count'] = {nm: mp.nstr(rect_count(XI[nm], 0.52, 2, res[nm]['T']), 6) for nm in ('ZETA', 'EPS', 'DH')}
print(out['right_of_line_count'])
out['meta'] = dict(dps=mp.mp.dps, seconds=round(time.time() - t0, 1))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'r21_transmission.json'), 'w'), indent=1)
