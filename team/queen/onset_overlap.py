"""Queen S1c: overlap of Goddess of Prime Number Theory's DH even onset witness (L = 3421/2000 half-width, K = 266,
cosine modes b_m(u) = L^{-1/2} cos((2m+1) pi u/(2L)) on [-L,L], mode weights |v_m|^2 from her runs/cos/dh_witness.json;
signs not deposited) with our quartet vector f0 = e^{-i g0 u} sinh(b0 u) on [-L,L].  <u,v> = int_{-L}^{L} u conj(v) du.
Bound: |<v,f0>|/||f0|| <= sum_{m<8} sqrt(w_m)|<b_m,f0>|/||f0|| + sqrt(1 - sum w_m)  (Cauchy-Schwarz, sign-free). NUMERICAL."""
import mpmath as mp, json
mp.mp.dps = 30
L = mp.mpf(3421)/2000; g0 = mp.mpf("85.699348"); b0 = mp.mpf("0.308517")
w = {0: 0.7987, 1: 0.1969, 2: 0.0027, 3: 0.0017, 4: 0.0001, 5: 0.0, 6: 0.0, 7: 0.0}
f0 = lambda u: mp.e**(-1j*g0*u)*mp.sinh(b0*u)
nf = mp.sqrt(mp.quad(lambda u: mp.sinh(b0*u)**2, [-L, L]))
pts = mp.linspace(-L, L, 200)
c = {}
for m in range(8):
    bm = lambda u, m=m: mp.cos((2*m+1)*mp.pi*u/(2*L))/mp.sqrt(L)
    c[m] = abs(mp.quad(lambda u: f0(u)*bm(u), pts))/nf
bound = sum(mp.sqrt(w[m])*c[m] for m in range(8)) + mp.sqrt(max(0, 1 - sum(w.values())))
res = dict(L=str(L), norm_f0=mp.nstr(nf, 10), cos_overlap_per_mode={m: mp.nstr(c[m], 6) for m in c}, overlap_upper_bound=mp.nstr(bound, 6),
           overlap_bound_top8_only=mp.nstr(sum(mp.sqrt(w[m])*c[m] for m in range(8)), 6))
print(json.dumps(res, indent=1)); json.dump(res, open("onset_overlap.json", "w"), indent=1)
