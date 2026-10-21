"""Grok S67: RIGOROUS (Arb) Gauss-panel quadrature error bound for the band parts of J in mori_route/mori_cert
(analogue of Mori VII Lemma 5.1 / Mori VIII Lemma 4.1, our own proof route).
Panels as in mori_route: edges 0,1/2,1,2,4,..,pw then width pw up to Omega; ng Gauss nodes per panel.
For a panel [t0,t1] (centre tc, half-length h) and rho>1 put E = image of the Bernstein ellipse E_rho.
Integrands, for any real unit x in one parity (b(z) = (b_d(z))_d, entire):
   F_x(z) = (1/pi) S(z) (x.b(z))^2,  S(z) = (psi(1/4+iz/2)+psi(1/4-iz/2))/2 - log pi - gamma      (gamma-weighted band)
   G_x(z) = (1/pi) (x.b(z)) b_0(z)                                                                (bb0 = B_band v_p)
Bessel: b_d(z) = +-int_I v_d(x) e^{izx} dx, so sum_d |b_d(z)|^2 <= int_{-L}^{L} e^{-2 Im z x} dx <= 2L e^{2L|Im z|}.
Trefethen ATAP Thm 19.3: |int_{-1}^1 f - Q_n f| <= (64/15) M rho^{-2n} / (rho^2-1), M = sup_{E_rho} |f|.
Since the error matrix is real symmetric, ||E_J||_2 = sup_x |x^T E_J x| <= sum_panels h * (64/15) rho^{-2n}/(rho^2-1)
   * (1/pi) * 2L e^{2L hb} * sup_{boundary} |S|,  hb = h(rho-1/rho)/2 (max |Im z| on E);  ||E_bb||_2 likewise without |S|.
S is holomorphic on E iff no pole z = +-i(2k+1/2) lies in the closed ellipse; checked rigorously per panel.
sup|S| on E = sup on its boundary (max modulus), enclosed by covering the boundary with K acb balls.
Output JSON: eps_J (band gamma-weighted part), eps_bb (vector), per-panel data.  Status: CERTIFIED (Arb), our proof."""
import sys, json, argparse
from flint import arb, acb, ctx
ap = argparse.ArgumentParser()
ap.add_argument('--L', default='71/50'); ap.add_argument('--Omega', default='536'); ap.add_argument('--s', default='763/200')
ap.add_argument('--m', default='1/2'); ap.add_argument('--pw', default='16'); ap.add_argument('--ng', type=int, default=96)
ap.add_argument('--K', type=int, default=4096); ap.add_argument('--out', default='')
a = ap.parse_args(); ctx.prec = 256
def q(s):
    p, _, r = s.partition('/'); return arb(int(p)) / (int(r) if r else 1)
L = q(a.L); Om = q(a.Omega); g = q(a.s) + q(a.m); pw = q(a.pw); n = a.ng; pi = arb.pi(); lp = pi.log()
edges = [arb(0), arb(1) / 2, arb(1)]
while edges[-1] < pw: edges.append(edges[-1] * 2)
while edges[-1] < Om: edges.append(edges[-1] + pw)
edges[-1] = Om
def S(z): return (acb(arb(1) / 4) + acb(0, 1) * z / 2).digamma() / 2 + (acb(arb(1) / 4) - acb(0, 1) * z / 2).digamma() / 2 - lp - g
def poles_out(tc, A_, B_, marg):
    # poles of S: z = +-i y, y = 2k+1/2 (k >= 0); z = i y in closed ellipse iff tc^2/A^2 + y^2/B^2 <= 1
    k = 0
    while True:
        y = 2 * arb(k) + arb(1) / 2
        if bool(y > B_ * 2) and bool(y > 1): return True     # y^2/B^2 > 4 > marg for all further k
        if not bool(tc * tc / (A_ * A_) + y * y / (B_ * B_) > marg): return False
        k += 1
def Ssup_arc(th0, th1, tc, A_, B_, depth=0):
    thm = (th0 + th1) / 2; zc = acb(tc + A_ * thm.cos(), B_ * thm.sin()); r = A_ * (th1 - th0) / 2
    zb = acb(arb(zc.real.mid(), r + zc.real.rad()), arb(zc.imag.mid(), r + zc.imag.rad()))
    v = abs(S(zb))
    if v.is_finite() and float(v.rad()) < 0.25 * (float(v.mid()) + 1): return arb(v.upper())
    assert depth < 25, 'cover failed'
    return max(Ssup_arc(th0, thm, tc, A_, B_, depth + 1), Ssup_arc(thm, th1, tc, A_, B_, depth + 1))
tot_J = arb(0); tot_b = arb(0); panels = []
for i in range(len(edges) - 1):
    t0, t1 = edges[i], edges[i + 1]; tc = (t0 + t1) / 2; h = (t1 - t0) / 2
    # rho: optimum of 2L h (rho-1/rho)/2 - 2n log rho, capped so the ellipse avoids the poles
    best = None
    for k in range(1, 101):   # rho <= 6 (cap: keeps the ellipses moderate; the bound is ample)
        rho = 1 + arb(k) / 20
        A_ = h * (rho + 1 / rho) / 2; B_ = h * (rho - 1 / rho) / 2
        # pole check: z = +-i y_k, y_k = 2k+1/2: inside closed ellipse iff tc^2/A^2 + y^2/B^2 <= 1
        if not poles_out(tc, A_, B_, arb(13) / 10): break
        lg = 2 * L * B_ - 2 * n * rho.log() - (rho * rho - 1).log()
        if best is None or float(lg.mid()) < best[0]: best = (float(lg.mid()), rho)
    rho = best[1]; A_ = h * (rho + 1 / rho) / 2; B_ = h * (rho - 1 / rho) / 2
    # rigorous pole exclusion (only y = 1/2 can matter: all others have y >= 5/2 > B_ or Re z > 0 region)
    assert poles_out(tc, A_, B_, arb(1))     # RIGOROUS: no pole of S in the closed ellipse
    # sup |S| on boundary: z(th) = tc + A cos th + i B sin th, arcs of length <= A*dth (|z'| <= A)
    Smax = arb(0); K = a.K; dth = 2 * pi / K
    for j in range(K):
        Smax = max(Smax, Ssup_arc(arb(j) * dth, arb(j + 1) * dth, tc, A_, B_))
    base = h * arb(64) / 15 * rho ** (-2 * n) / (rho * rho - 1) / pi * 2 * L * (2 * L * B_).exp()
    eJ = base * Smax; eb = base
    tot_J += eJ; tot_b += eb
    print(i, rho.str(3), eJ.str(3), file=sys.stderr, flush=True)
    panels.append(dict(t0=t0.str(6), t1=t1.str(6), rho=rho.str(5), semi_minor=B_.str(5), Smax=Smax.str(5), epsJ=eJ.str(5), epsbb=eb.str(5)))
res = dict(L=a.L, Omega=a.Omega, s=a.s, m=a.m, pw=a.pw, ng=n, K_cover=a.K, n_panels=len(panels),
           eps_J=arb(tot_J.upper()).str(6), eps_bb=arb(tot_b.upper()).str(6), eps_J_upper=float(tot_J.upper()), eps_bb_upper=float(tot_b.upper()), panels=panels)
print(json.dumps(res, indent=1))
if a.out: json.dump(res, open(a.out, 'w'), indent=1)
