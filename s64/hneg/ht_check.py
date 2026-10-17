import sys; sys.argv=['x','0','1','1','/dev/null']
exec(open('ht_scan.py').read().split('t0 = time.time()')[0])
import mpmath as mp
for x in [1, 10, 28.2694, 100, 250]:
    xi = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    print(x, mp.nstr(mp.mpf(Ht(repr(x)).mid().str(60, radius=False)),20), mp.nstr(xi(mp.mpf(1)/2+1j*mp.mpf(x)/2).real/8,20))
