"""Grok S61 H61a(2,3): numerical audit of Codex's all-L energy-transport identity (PR #11 comment 6231393772)
   Q_L(f) = |E|^2-|O|^2 - (2S_L-omega0)||f||^2 + E_gamma(f) + E_prime(f)
against (a) the independent Fourier-side Weil explicit formula (digamma integral over r, pole via fhat(+-i/2),
primes via direct quadrature of g(log n)), and (b) the zero side sum_rho |fhat(gamma)|^2 over 3000 Arb zeros + tail.
Test f(x) = q(x/L)(1-(x/L)^2)^K on [-L,L] (K=6 so |fhat|^2 ~ t^-14), q rational, possibly complex.
Conventions: fhat(t)=int f e^{itx}, g(u)=int f(x+u)conj f(x)dx, U_u f(x)=f(x+u), w_n=Lambda(n)/sqrt n. NUMERICAL."""
import sys, json, time
import mpmath as mp
from fractions import Fraction as Fr
mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 45
K = 6
def vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else None
def pmul(a, b):
    c = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] += x*y
    return c
def ppow(a, k):
    r = [1]
    for _ in range(k): r = pmul(r, a)
    return r
def peval(c, s): return mp.polyval(list(reversed([mp.mpmathify(x) if not isinstance(x, complex) else mp.mpc(x) for x in c])), s)
EDGE = ppow([1, 0, -1], K)
def mkf(qre, qim):   # returns list of complex Fraction pairs (re, im) coefficients in s
    a = pmul(qre, EDGE); b = pmul(qim, EDGE) if qim else [0]
    n = max(len(a), len(b)); a += [0]*(n-len(a)); b += [0]*(n-len(b))
    return [(Fr(x), Fr(y)) for x, y in zip(a, b)]
def cval(c):
    return [mp.mpc(mp.mpf(x.numerator)/x.denominator, mp.mpf(y.numerator)/y.denominator) for x, y in c]
def binom(n, k): return mp.binomial(n, k)
def autoc_re(coef):
    """Exact (Fraction) polynomial R(v) = Re int_{-1}^{1-v} F(s+v) conj F(s) ds, v in [0,2]. Returns coeff list in v."""
    from math import comb
    deg = len(coef)-1; out = [Fr(0)]*(2*deg+2)
    for a in range(deg+1):
        ar, ai = coef[a]
        for b in range(deg+1):
            br, bi = coef[b]
            cr = ar*br + ai*bi            # Re(c_a conj c_b)
            if cr == 0: continue
            for i in range(a+1):           # (s+v)^a = sum C(a,i) s^i v^(a-i)
                e = i+b+1; co = cr*comb(a, i)/e
                # int_{-1}^{1-v} s^(i+b) ds = [(1-v)^e - (-1)^e]/e
                for m in range(e+1):       # (1-v)^e = sum C(e,m)(-v)^m
                    out[a-i+m] += co*comb(e, m)*(-1)**m
                out[a-i] -= co*(-1)**e
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out
def fhat(c, L, t):
    """fhat(t) = L int_{-1}^1 F(s) e^{i t L s} ds: Taylor/moment series for |tL| < 12, else closed form by repeated
    integration by parts; both at dps+40 internally (the boundary form cancels badly for small |tL|)."""
    with mp.workdps(mp.mp.dps + 40):
        t = mp.mpmathify(t); w = t*L
        if abs(w) < 12:
            tot = mp.mpc(0); cabs = 2*sum(abs(x) for x in c)
            for m in range(0, 600):
                mom = sum(c[k]*mp.mpf(1-(-1)**(k+m+1))/(k+m+1) for k in range(len(c)))   # (S61 fix: was int/int float division)
                term = (1j*w)**m/mp.factorial(m)*mom; tot += term
                if m > 40 and abs(w)**m/mp.factorial(m)*cabs < mp.mpf(10)**(-mp.mp.dps-5): break   # S61 fix: parity zeros made term==0 and broke early
            res = L*tot
        else:
            d = list(c); tot = mp.mpc(0); iw = 1j*w; k = 0
            while any(x != 0 for x in d):
                val = lambda s: mp.polyval(list(reversed(d)), s)
                tot += (-1)**k*(val(1)*mp.exp(iw) - val(-1)*mp.exp(-iw))/iw**(k+1)
                d = [i*d[i] for i in range(1, len(d))]; k += 1
            res = L*tot
    return +res
def run(name, qre, qim, Ls, zeros):
    coef = mkf(qre, qim); c = cval(coef); L = mp.mpf(Fr(Ls).numerator)/Fr(Ls).denominator
    F = lambda s: mp.polyval(list(reversed(c)), s)
    R = autoc_re(coef); Rm = [mp.mpf(x.numerator)/x.denominator for x in R]
    N = L*Rm[0]                                         # ||f||^2
    ns = [n for n in range(2, int(mp.exp(2*L))+2) if vm(n) and mp.log(n) < 2*L]
    w = {n: mp.log(vm(n))/mp.sqrt(n) for n in ns}; S = sum(w.values())
    om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)
    # --- Codex real-space form ---
    lp = L*mp.quad(lambda s: F(s)*mp.exp(L*s/2), [-1, 0, 1]); lm = L*mp.quad(lambda s: F(s)*mp.exp(-L*s/2), [-1, 0, 1])
    E = (lp+lm)/mp.sqrt(2); O = (lp-lm)/mp.sqrt(2)
    pole_codex = abs(E)**2 - abs(O)**2; pole_lr = 2*mp.re(lp*mp.conj(lm))
    qv = [-x for x in Rm[1:]]                            # (R(0)-R(v))/v
    kup = lambda u: mp.exp(-u/2)*u/(1-mp.exp(-2*u)) if u != 0 else mp.mpf(1)/2   # k(u)*u
    Eg = 2*L*mp.quad(lambda v: kup(L*v)*mp.polyval(list(reversed(qv)), v), [0, 0.5, 1, 1.5, 2])  # 2 int k(u)(N-Re g)du, u=Lv
    Kt = mp.nsum(lambda j: mp.exp(-(2*j+mp.mpf(1)/2)*2*L)/(2*j+mp.mpf(1)/2), [0, mp.inf])
    Eg += 2*N*Kt
    Rg = lambda u: L*mp.polyval(list(reversed(Rm)), u/L)
    Ep = {n: w[n]*(2*N-2*Rg(mp.log(n))) for n in ns}; EP = sum(Ep.values())
    Q_codex = pole_codex - (2*S-om0)*N + Eg + EP
    # --- independent Fourier-side explicit formula ---
    fh = lambda t: fhat(c, L, t)
    pole_ef = 2*mp.re(fh(1j/2)*mp.conj(fh(-1j/2)))      # ghat(i/2)+ghat(-i/2), ghat(t)=fhat(t)conj(fhat(conj t))
    gdir = lambda u: mp.quad(lambda x: F((x+u)/L)*mp.conj(F(x/L)), [-L, (L-u)/2 - L/2 + 0, L-u]) if u < 2*L else 0
    prim = -2*sum(w[n]*mp.re(gdir(mp.log(n))) for n in ns)
    Rmax = 1500; br = [mp.mpf(x) for x in range(0, Rmax+1, 3)]
    integrand = lambda r: (abs(fh(r))**2 + abs(fh(-r))**2)*mp.re(mp.digamma(mp.mpf(1)/4 + 1j*r/2))
    archI = mp.quad(integrand, br) + mp.quad(integrand, [Rmax, mp.inf])
    N_planch = None
    arch_ef = -mp.log(mp.pi)*N + archI/(2*mp.pi)
    Q_ef = pole_ef + arch_ef + prim
    # --- zero side ---
    Zs = mp.fsum(abs(fh(g))**2 + abs(fh(-g))**2 for g in zeros)
    T = zeros[-1]
    tail = mp.quad(lambda t: (abs(fh(t))**2 + abs(fh(-t))**2)*mp.log(t/(2*mp.pi))/(2*mp.pi), [T, 2*T, mp.inf])
    Qz = Zs + tail
    # arch check of codex arch vs Fourier arch: arch_codex = om0*N + Eg
    out = dict(name=name, L=Ls, K=K, q_re=[str(x) for x in qre], q_im=[str(x) for x in qim] if qim else None,
               norm2=mp.nstr(N, 25), S_L=mp.nstr(S, 20), omega0=mp.nstr(om0, 20), prime_powers=ns,
               ellp=mp.nstr(lp, 20), ellm=mp.nstr(lm, 20), absE2=mp.nstr(abs(E)**2, 20), absO2=mp.nstr(abs(O)**2, 20),
               pole_codex=mp.nstr(pole_codex, 25), pole_minus_2Re_lp_conj_lm=mp.nstr(pole_codex-pole_lr, 3),
               pole_ef_minus_codex=mp.nstr(pole_ef-pole_codex, 3),
               E_gamma=mp.nstr(Eg, 25), E_prime=mp.nstr(EP, 25), E_prime_by_n={n: mp.nstr(v, 15) for n, v in Ep.items()},
               arch_codex_minus_ef=mp.nstr(om0*N+Eg-arch_ef, 3), prime_codex_minus_direct=mp.nstr(EP-2*S*N-prim, 3),
               Q_codex=mp.nstr(Q_codex, 30), Q_explicit_formula=mp.nstr(Q_ef, 30), Q_zeros=mp.nstr(Qz, 30),
               zero_tail=mp.nstr(tail, 4), n_zeros=len(zeros),
               diff_codex_ef=mp.nstr(Q_codex-Q_ef, 3), diff_codex_zeros=mp.nstr(Q_codex-Qz, 3),
               rel_diff_codex_zeros=mp.nstr((Q_codex-Qz)/Qz, 3))
    return out
if __name__ == "__main__":
    t0 = time.time()
    zl = [mp.mpf(line.split()[0]) for line in open('../s51/zeros_arb_3000.txt')]
    nz = int(sys.argv[2]) if len(sys.argv) > 2 else 3000; zl = zl[:nz]
    # S39 cross-check (no edge factor): f=x and f=1 on [-1,1]
    for nm, cc in (("S39_odd_f=x", [0, 1]), ("S39_even_f=1", [1])):
        Fv = lambda s, cc=cc: mp.polyval(list(reversed(cc)), s)
        lp = mp.quad(lambda s: Fv(s)*mp.exp(s/2), [-1, 1]); lm = mp.quad(lambda s: Fv(s)*mp.exp(-s/2), [-1, 1])
        E2 = (lp+lm)**2/2; O2 = (lp-lm)**2/2
        print(json.dumps(dict(name=nm, L=1, pole=mp.nstr(2*lp*lm, 20), absE2_minus_absO2=mp.nstr(E2-O2, 20))), flush=True)
    tests = [("even", [1, 0, Fr(3, 10)], None), ("odd", [0, 1, 0, Fr(-1, 2)], None),
             ("complex", [1, 2, 0, Fr(1, 3)], [Fr(1, 2), Fr(-1, 3), 1])]
    for Ls in sys.argv[3].split(',') if len(sys.argv) > 3 else ["1/2", "1", "27/20"]:
        for nm, qr, qi in tests:
            r = run(nm, qr, qi, Ls, zl); r["wall_s"] = round(time.time()-t0); print(json.dumps(r), flush=True)
