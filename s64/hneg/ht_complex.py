"""Grok S64 H64c stage 2: locate the non-real zeros of H_t (t<0) near the scan candidates by complex Newton (Arb, NUMERICAL)."""
import sys, json
sys.argv = ['x', sys.argv[1], '300', '0.02', '/dev/null'] + sys.argv[2:]
src = open('ht_scan.py').read().split('t0 = time.time()')[0]; exec(src)
from flint import acb
def F(z): return sum((c*(z*u).cos() for u, c in zip(U, C)), acb(0))
def dF(z): return sum((-c*u*(z*u).sin() for u, c in zip(U, C)), acb(0))
d = json.load(open(f'scan_t{sys.argv[1]}.json')); found = []
for x0 in d['absmin_nosign_candidates']:
    xa = arb(repr(x0)); h0 = Ht(xa); h2 = sum((-c*u*u*(xa*u).cos() for u, c in zip(U, C)), arb(0))
    r = float((2*h0/h2).mid()); yest = r**0.5 if r > 0 else 0.3
    for y0 in (yest, 0.7*yest, 1.4*yest, 0.3, 0.8):
        z = acb(x0, y0)
        for _ in range(200):
            dz = F(z)/dF(z); z = (z - dz).mid()
            if abs(complex(dz.mid())) < 1e-40: break
        zc = complex(z.mid())
        if abs(complex(dz.mid())) < 1e-30 and abs(zc.imag) > 1e-6 and 0 < zc.real < 320:
            zc = complex(zc.real, abs(zc.imag))
            if all(abs(zc - f) > 1e-6 for f in found): found.append(zc)
            break
found.sort(key=lambda z: z.real)
json.dump({"t": sys.argv[1], "complex_zeros_x_upper": [[z.real, z.imag] for z in found],
           "gamma_form": [[z.real/2, z.imag/2] for z in found]}, open(f'cplx_t{sys.argv[1]}.json', 'w'))
print(sys.argv[1], len(found), [(round(z.real/2, 4), round(z.imag/2, 5)) for z in found])
