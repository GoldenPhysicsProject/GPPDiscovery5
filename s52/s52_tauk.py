# Grok S52 H52b(1): unconditional complete-monotonicity window for the zero-free heat trace
# K(tau)=sum_{Im rho>0} exp((rho-1/2)^2 tau).  PROVED bound (constants evaluated in mpmath, dps 50):
#  (-1)^k K^(k)(tau) >= g1^{2k} e^{-g1^2 tau} - T_k,   T_k <= g(H)[H log H + (log H+1)/a + 1/(H a^2)],
#  g(t)=e^{tau/4}(t^2+1/4)^k e^{-t^2 tau},  a=2H(tau - k/H^2) > 0,
# using: zeros with 0<gamma<=H on the line (Platt-Trudgian 2021, H=3e12 < 3.0000175e12),
#        N(t) <= t log t for t>=H, |w|<=gamma^2+1/4, Re w <= 1/4-gamma^2, g decreasing on [H,inf) when H^2 tau>k.
import mpmath as mp
mp.mp.dps = 50
H = mp.mpf('3e12'); g1 = mp.mpf('14.134725141734693')   # lower bound of gamma_1
def lhs_minus_T(k, tau):
    a = 2*H*(tau - k/H**2)
    if a <= 0: return None
    lg = tau/4 + k*mp.log(H**2+mp.mpf(1)/4) - H**2*tau
    lT = lg + mp.log(H*mp.log(H) + (mp.log(H)+1)/a + 1/(H*a**2))
    lL = 2*k*mp.log(g1) - g1**2*tau
    return lL - lT     # >0  => CM inequality of order k holds at tau
def tauk(k):
    lo, hi = mp.mpf(k)/H**2*(1+mp.mpf(10)**-30), mp.mpf(1)
    for _ in range(300):
        mid = (lo+hi)/2
        v = lhs_minus_T(k, mid)
        if v is not None and v > 0: hi = mid
        else: lo = mid
    return hi
print('H =', H)
for k in [0, 1, 2, 5, 10, 100, 1000, 10**4, 10**6]:
    t = tauk(k); print('k=%-8d tau_k=%s   tau_k*H^2=%s' % (k, mp.nstr(t, 8), mp.nstr(t*H**2, 8)))
slope = 2*mp.log(H/g1)
print('asymptotic slope tau_k*H^2/k -> log((H^2+1/4)/g1^2) + 1 ~', mp.nstr(slope, 8), '(+1 from k/H^2 shift)')
# monotonicity in tau: v(k,tau) increasing in tau for tau>=tau_k? check a few multiples
for k in [0, 10, 1000]:
    t = tauk(k); print(' k=%d checks at tau_k*{1.01,2,10,1e6}:' % k, [mp.nstr(lhs_minus_T(k, t*m), 4) for m in [1.01, 2, 10, 10**6]])
