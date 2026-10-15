"""S26 [design, UNTRUSTED float LP]: enriched claim-1 for zeta odd sector.
m(xi) = sigma_A(xi) + sum_j C_j wDC(xi; Dd_j) + sum_p v_p wPL_p(xi) - sum_k a_k cos(u_k xi) >= c
- Multi-DC: several trapezoids with different Dd (C_j >= 0 by default).
- PL exterior: W=1 on |t|<=2h, free heights at exterior knots (signed), 0 at end.
- Cosine lattice u_k = 2h + k D, k=0..K-1, with Bernstein cap on the trig poly.
For ODD g on (-h,h): int |ghat|^2 wDC = int |ghat|^2 wPL = 0 and cos(u>=2h) invisible.
UNTRUSTED float; Arb verify separately. RH OPEN.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok."""
from __future__ import annotations
import numpy as np, scipy.special as sp, json, sys, time, argparse
from scipy.optimize import linprog

l2 = np.log(2)
kap = l2 / np.sqrt(2)


def sigA(xi):
    return sp.psi(0.25 + 0.5j * xi).real - np.log(np.pi) - 2 * kap * np.cos(xi * l2)


def wDC(xi, h, Dd):
    """FT of even trapezoid: 1 on |t|<=2h, linear to 0 at 2h+Dd."""
    a = 2 * h + Dd / 2
    xi = np.asarray(xi, float)
    return 2 * a * np.sinc(a * xi / np.pi) * np.sinc(Dd * xi / 2 / np.pi)


def _seg_ft(a, b, ya, yb, xi):
    """2 * int_a^b linear(ya,yb) cos(xi t) dt  (even-function convention uses factor 2 outside).
    Returns the contribution to int_0^inf for one segment; caller doubles for even FT,
    OR we return 2*int directly here."""
    xi = np.asarray(xi, float)
    out = np.zeros_like(xi, dtype=float)
    dz = b - a
    # W(t) = ya + (yb-ya)/dz * (t-a)
    # I = int_a^b W cos = ya * S0 + slope * S1
    # S0 = int_a^b cos = (sin(b xi)-sin(a xi))/xi
    # int_a^b (t-a) cos(xi t) dt
    small = np.abs(xi) < 1e-12
    xs = xi[~small]
    if xs.size:
        sa, sb = np.sin(a * xs), np.sin(b * xs)
        ca, cb = np.cos(a * xs), np.cos(b * xs)
        S0 = (sb - sa) / xs
        # int (t-a) cos(xt) dt = [(t-a)sin(xt)/x + cos(xt)/x^2]_a^b - 0
        # = (b-a)sin(b x)/x + (cos(b x)-cos(a x))/x^2
        S1 = dz * sb / xs + (cb - ca) / xs**2
        slope = (yb - ya) / dz
        out[~small] = 2.0 * (ya * S0 + slope * S1)  # factor 2 for even FT
    # xi~0: 2 * int_a^b W = 2 * dz * (ya+yb)/2 = dz*(ya+yb)
    out[small] = dz * (ya + yb)
    return out


def wPL_basis(xi, h, knots, which):
    """FT contribution of unit height at exterior knot index `which` (1..M-1),
    with y[0]=1 at knots[0]=2h fixed, y[M]=0 at knots[M] fixed.
    Returns the FT of the PL function that is 1 on [0,2h] and has y[which]=1, other free y=0.
    Actually we want columns = d(FT)/d y_p for p=1..M-1 of the FULL W (including the fixed y0=1 piece as separate).

    Cleaner: return FT of W with y0=1, yM=0, and y_p = e_p (unit vector at free knots).
    For the 'box+exterior' we split:
      W0: y0=1, all free=0, yM=0  -> base column (like a single tent-to-zero from 2h if M=1)
      Wp: y0=0, y_p=1, others 0, yM=0 -> free columns
    But W must be 1 on [0,2h]. So W0 includes the box [0,2h] plus linear from 1 at 2h toward y1
    (or to 0 if M=1).

    Implementation: build W with given y_ext (length M-1), y0=1, yM=0, plus box.
    """
    raise NotImplementedError  # use wPL_columns below


def wPL_columns(xi, h, knots):
    """knots: array z[0]=2h < z[1] < ... < z[M], M>=1.
    Returns (W0, [W1..W_{M-1}]) where
      W0 = FT of (1 on [0,2h], linear 1->y1_overlap with free=0 so linear 1->0 on first segment if we set path via zeros... )
    Better concrete:
    Full W(y) with y = (y1,...,y_{M-1}):
      on [0,2h]: 1
      on [z_i, z_{i+1}]: linear from y_i to y_{i+1}, with y0=1, yM=0.
    FT(W) = FT(box[0,2h]) + sum segments
           = box_ft + seg(z0,z1,1,y1) + seg(z1,z2,y1,y2) + ... + seg(z_{M-1},zM,y_{M-1},0)
    Affine in y: FT = A0 + sum_p y_p A_p
    Returns A0 and list A_p (p=1..M-1).
    """
    z = np.asarray(knots, float)
    assert z[0] == 2 * h or abs(z[0] - 2 * h) < 1e-14
    M = len(z) - 1  # number of segments
    xi = np.asarray(xi, float)
    # box [0, 2h]
    small = np.abs(xi) < 1e-12
    box = np.empty_like(xi, dtype=float)
    box[small] = 2 * (2 * h)  # 2*int_0^{2h} 1
    box[~small] = 2 * np.sin(2 * h * xi[~small]) / xi[~small]

    A0 = box.copy()
    # segment 0: y0=1 to y1 (y1 free or 0 if M=1). Contribution to A0 from y0=1 side; y1 goes to A_1.
    # seg(a,b,ya,yb) = ya * S_a + yb * S_b  (linear response)
    def seg_ya_yb(a, b, xi):
        """Return (coeff of ya, coeff of yb) for 2*int_a^b linear cos."""
        out_a = np.zeros_like(xi, dtype=float)
        out_b = np.zeros_like(xi, dtype=float)
        dz = b - a
        small = np.abs(xi) < 1e-12
        xs = xi[~small]
        if xs.size:
            sa, sb = np.sin(a * xs), np.sin(b * xs)
            ca, cb = np.cos(a * xs), np.cos(b * xs)
            S0 = (sb - sa) / xs
            S1 = dz * sb / xs + (cb - ca) / xs**2  # int (t-a) cos
            # W = ya*(1 - u/dz) + yb*(u/dz) = ya + (yb-ya)/dz * u
            # I = ya*S0 + ((yb-ya)/dz)*S1 = ya*(S0 - S1/dz) + yb*(S1/dz)
            out_a[~small] = 2.0 * (S0 - S1 / dz)
            out_b[~small] = 2.0 * (S1 / dz)
        out_a[small] = dz  # 2 * int ya*(1-u/dz) = 2*ya*dz/2 = ya*dz; here return coeff
        out_b[small] = dz
        return out_a, out_b

    Ap = [np.zeros_like(xi, dtype=float) for _ in range(M - 1)]  # free y1..y_{M-1}
    for i in range(M):
        a, b = z[i], z[i + 1]
        ca, cb = seg_ya_yb(a, b, xi)
        # y_i and y_{i+1}; y0=1 fixed, yM=0 fixed
        if i == 0:
            A0 += ca  # y0=1
            if M == 1:
                pass  # y1 = yM = 0
            else:
                Ap[0] += cb  # y1
        elif i == M - 1:
            # y_{M-1} free (unless M==1 handled), yM=0
            if M - 1 >= 1:
                Ap[M - 2] += ca
            # yM=0 contributes nothing
        else:
            Ap[i - 1] += ca  # y_i
            Ap[i] += cb      # y_{i+1}
    return A0, Ap


def choose_Xi(B, nph=8, margin=0.05, xs_max=20000.0, xs_step=0.02):
    """Smallest Xi where envelope of sigA exceeds B/cos(pi/nph)+margin."""
    thr = B / np.cos(np.pi / nph) + margin
    xs = np.arange(0.0, xs_max, xs_step)
    # lower envelope of sigA: Re psi - log pi - 2 kap
    low = sp.psi(0.25 + 0.5j * xs).real - np.log(np.pi) - 2 * kap
    ok = low > thr
    if not np.any(ok):
        return float(xs_max)
    return float(xs[np.argmax(ok)])


def lp_enriched(
    h,
    U=12.0,
    D=0.05,
    Dd_list=(0.5,),
    B=3.0,
    useDC=True,
    C_signed=False,
    pl_knots=None,
    pl_signed=True,
    xi_step=0.005,
    nph=8,
    Xi=None,
    method="highs",
    bernstein_samples=1500,
):
    """Solve max c s.t. m >= c on [0,Xi] and Bernstein on a[].
    Dd_list: multi-DC widths (used if useDC).
    pl_knots: optional exterior knot array starting at 2h; if set, PL columns included.
    """
    K = int(round((U - 2 * h) / D)) + 1
    u = 2 * h + D * np.arange(K)
    if Xi is None:
        Xi = choose_Xi(B, nph=nph)
    xi = np.arange(0.0, Xi, xi_step)
    s = sigA(xi)

    cols = [np.ones((xi.size, 1))]  # c
    nd = 0
    npl = 0
    # DC columns: constraint is c - sum C_j wDC_j + sum a cos <= sigA
    # i.e. A_ub row: [1, -wDC_1, ..., -wDC_nd, (-A0 if pl base absorbed?), ...]
    # Wait: m = sigA + sum C wDC + sum v A_pl - sum a cos >= c
    # => c - sum C wDC - sum v A_pl + sum a cos <= sigA
    if useDC:
        for Dd in Dd_list:
            cols.append(-wDC(xi, h, Dd)[:, None])
        nd = len(Dd_list)
    pl_A0 = None
    if pl_knots is not None:
        A0, Ap = wPL_columns(xi, h, pl_knots)
        # Include A0 with free coefficient C0>=0 (or signed): m += C0 * A0 + sum v_p Ap
        cols.append(-A0[:, None])
        for Ap_p in Ap:
            cols.append(-Ap_p[:, None])
        npl = 1 + len(Ap)
        pl_A0 = True
    cols.append(np.cos(np.outer(xi, u)))
    A1 = np.hstack(cols)

    # Bernstein: sum a_k cos(u_k theta + ph) <= B  (no c, C, v in these rows)
    th = np.linspace(0, 2 * np.pi, bernstein_samples, endpoint=False)
    n_left = 1 + nd + npl
    A2 = np.hstack(
        [
            np.zeros((bernstein_samples * nph, n_left)),
            np.vstack(
                [
                    np.cos(np.outer(th, np.arange(K)) + ph)
                    for ph in 2 * np.pi * np.arange(nph) / nph
                ]
            ),
        ]
    )
    # bounds: c free; C_j >=0 or free; v_p free (or >=0); a free
    if C_signed:
        bC = [(None, None)] * nd
    else:
        bC = [(0, None)] * nd
    if pl_knots is not None:
        if pl_signed:
            bPL = [(None, None)] * npl
        else:
            bPL = [(0, None)] * npl
    else:
        bPL = []
    bnds = [(None, None)] + bC + bPL + [(None, None)] * K
    cobj = np.r_[-1.0, np.zeros(nd + npl + K)]
    r = linprog(
        cobj,
        A_ub=np.vstack([A1, A2]),
        b_ub=np.r_[s, np.full(bernstein_samples * nph, B)],
        bounds=bnds,
        method=method,
    )
    meta = dict(h=h, U=U, D=D, Dd_list=list(Dd_list), B=B, useDC=useDC, C_signed=C_signed,
                pl_knots=None if pl_knots is None else list(map(float, pl_knots)),
                pl_signed=pl_signed, xi_step=xi_step, Xi=Xi, K=K, nd=nd, npl=npl,
                n_rows=A1.shape[0] + A2.shape[0], n_cols=A1.shape[1])
    return r, u, meta


def pack_result(r, u, meta):
    out = dict(meta)
    out["status"] = int(r.status)
    out["message"] = str(r.message)
    if r.status == 0:
        x = r.x
        nd, npl, K = meta["nd"], meta["npl"], meta["K"]
        out["c_LP"] = float(-r.fun)
        out["C"] = x[1 : 1 + nd].tolist()
        out["v_pl"] = x[1 + nd : 1 + nd + npl].tolist()
        out["a"] = x[1 + nd + npl :].tolist()
        out["u0"] = float(u[0])
    else:
        out["c_LP"] = None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("h", help="h value or hstar")
    ap.add_argument("--U", type=float, default=12.0)
    ap.add_argument("--D", type=float, default=0.05)
    ap.add_argument("--Dd", type=str, default="0.5", help="comma list of Dd for multi-DC")
    ap.add_argument("--B", type=float, default=3.0)
    ap.add_argument("--xi-step", type=float, default=0.005)
    ap.add_argument("--no-DC", action="store_true")
    ap.add_argument("--C-signed", action="store_true")
    ap.add_argument("--pl", type=str, default="", help="PL exterior: 'M,L' => M free knots over length L past 2h")
    ap.add_argument("--pl-knots", type=str, default="", help="explicit knots comma list starting at 2h")
    ap.add_argument("--pl-nonneg", action="store_true")
    ap.add_argument("--method", default="highs")
    ap.add_argument("--out", default="")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    h = np.log(3) / 2 if args.h == "hstar" else float(args.h)
    Dd_list = tuple(float(x) for x in args.Dd.split(",") if x.strip())
    pl_knots = None
    if args.pl_knots:
        pl_knots = np.array([float(x) for x in args.pl_knots.split(",")])
    elif args.pl:
        M, L = args.pl.split(",")
        M, L = int(M), float(L)
        # M = number of segments; knots: 2h = z0 < z1 < ... < zM = 2h+L; free heights M-1
        pl_knots = 2 * h + np.linspace(0, L, M + 1)
    t0 = time.time()
    r, u, meta = lp_enriched(
        h,
        U=args.U,
        D=args.D,
        Dd_list=Dd_list,
        B=args.B,
        useDC=not args.no_DC,
        C_signed=args.C_signed,
        pl_knots=pl_knots,
        pl_signed=not args.pl_nonneg,
        xi_step=args.xi_step,
        method=args.method,
    )
    out = pack_result(r, u, meta)
    out["time_s"] = time.time() - t0
    out["tag"] = args.tag
    line = (
        f"h={args.h} tag={args.tag} DC={Dd_list} pl={args.pl or args.pl_knots or '-'} "
        f"U={args.U} D={args.D} B={args.B} xi={args.xi_step}: "
        f"status={out['status']} c_LP={out['c_LP']} Xi={out['Xi']:.1f} "
        f"K={out['K']} nd={out['nd']} npl={out['npl']} ({out['time_s']:.0f}s)"
    )
    print(line, flush=True)
    path = args.out or f"lp_h{args.h}_{args.tag or 'run'}.json"
    json.dump(out, open(path, "w"), indent=1)
    print("wrote", path, flush=True)


if __name__ == "__main__":
    main()
