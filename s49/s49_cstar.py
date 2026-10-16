# refine two-mode critical c*: smallest c>1 with min_{u>=0} -m'_c(u) < 0, F_c = cosh2z + c cosh z
import mpmath as mp, json
mp.mp.dps = 40
def mm(c, U):
    s = mp.sqrt(U); C1, S1, C2, S2 = mp.cosh(s), mp.sinh(s), mp.cosh(2*s), mp.sinh(2*s)
    # f(u)=C2 + c C1 ; derivatives in u via s: d/du = (1/(2s)) d/ds
    f = C2 + c*C1
    fs = 2*S2 + c*S1; fss = 4*C2 + c*C1
    f1 = fs/(2*s); f2 = (fss/(2*s) - fs/(2*s**2))/(2*s)
    return (f1**2 - f*f2)/f**2
def minval(c):
    g = lambda U: mm(c, U)
    # coarse then golden refine
    Us = [mp.mpf(x)/4 for x in range(4, 800)]
    v = [g(U) for U in Us]; i = min(range(len(v)), key=lambda j: v[j])
    a, b = Us[max(i-1, 0)], Us[min(i+1, len(Us)-1)]
    for _ in range(80):
        m1, m2 = a + (b-a)*0.382, a + (b-a)*0.618
        if g(m1) < g(m2): b = m2
        else: a = m1
    return g((a+b)/2), (a+b)/2
lo, hi = mp.mpf(300), mp.mpf(500)
for _ in range(60):
    mid = (lo+hi)/2
    if minval(mid)[0] < 0: hi = mid
    else: lo = mid
v, U = minval(hi)
print('c* =', mp.nstr(hi, 20), ' e^6 =', mp.nstr(mp.e**6, 20), ' log c* =', mp.nstr(mp.log(hi), 15), ' touching u* =', mp.nstr(U, 15))
json.dump(dict(cstar=mp.nstr(hi, 20), log_cstar=mp.nstr(mp.log(hi), 15), u_touch=mp.nstr(U, 15), e6=mp.nstr(mp.e**6, 20),
               label='NUMERICAL (mpmath 40 digits, bisection+golden section)'), open('s49_cstar.json', 'w'), indent=1)
