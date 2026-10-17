"""Grok S64 H64c: complete the non-real zero list of H_t (x<=300) by complex Newton from a 2D start grid; count check vs scan deficit."""
import sys, json
tt = sys.argv[1]; sys.argv = ['x', tt, '300', '0.02', '/dev/null']
exec(open('ht_scan.py').read().split('t0 = time.time()')[0])
from flint import acb
def F(z): return sum((c*(z*u).cos() for u, c in zip(U, C)), acb(0))
def dF(z): return sum((-c*u*(z*u).sin() for u, c in zip(U, C)), acb(0))
found = [complex(a, b) for a, b in json.load(open(f'cplx_t{tt}.json'))['complex_zeros_x_upper']]
for x0 in range(2, 300, 3):
    for y0 in (0.4, 1.0, 1.8, 2.8):
        z = acb(x0, y0); ok = False
        for _ in range(40):
            dz = F(z)/dF(z); z = (z - dz).mid()
            if abs(complex(dz.mid())) < 1e-35: ok = True; break
            if abs(complex(z.mid()).imag) > 8: break
        zc = complex(z.mid())
        if ok and abs(zc.imag) > 1e-6 and 0 < zc.real <= 300:
            zc = complex(zc.real, abs(zc.imag))
            if all(abs(zc - f) > 1e-6 for f in found): found.append(zc)
found.sort(key=lambda z: z.real)
d = json.load(open(f'scan_t{tt}.json')); deficit = (52 - d['nreal'])//2
json.dump({"t": tt, "complex_zeros_x_upper": [[z.real, z.imag] for z in found], "gamma_form": [[z.real/2, z.imag/2] for z in found],
           "expected_pairs_from_count": deficit, "method": "grid Newton"}, open(f'cplx_t{tt}.json', 'w'))
print(tt, 'found', len(found), 'expected', deficit, [(round(z.real/2, 3), round(z.imag/2, 4)) for z in found])
