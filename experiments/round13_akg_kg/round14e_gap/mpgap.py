"""Round 14E2 (Golden Eye AUDIT14 6.8 falsifier): the 14E quotient question for zeta at 80+ digits.
S = window-truncated range of E_F, built exactly as in gap.py::radical_space (same family layout c, w, k, phi; same pole constraints;
same concentration cut rho >= 1-eps; same Parseval route through F on the critical line), except:
  * every step is in 90-digit arithmetic (python-flint arb for the t-integrals, mpmath for all eigenproblems);
  * the bump profile is replaced by a Gaussian of s.d. SIG*w in v = log y (centre c, half-width w as in gap.py), because the bump's
    Mellin transform decays only like exp(-c sqrt t), which caps any t-truncated Parseval integral near 1e-12; the Gaussian's
    transform is closed form and decays like exp(-(SIG w)^2 t^2/2), so the T = 260 trapezoid (h = 0.1, entire integrand) is exact
    far below 1e-80 (checked against h = 0.2).
  * whitening cut TAU on the family Gram (gap.py: 1e-13 relative) is a parameter.
Outputs, for each cut eps: dim S, bilinear residual ||W S||/||W||, lambda_min of the RESTRICTION Q^T W Q (Q = L2-complement of S,
the 14E pipeline) and of the TRUE QUOTIENT seminorm q([x]) = inf_{s in S} W(x+s) on S^perp (Schur complement W/W_SS), the
interlacing bracket lambda_1(W) <= both <= lambda_{dimS+1}(W), and principal angles between S and the bottom-dimS eigenvectors of W.
Integer data only (Dirichlet coefficients through F(1/2+it) values, gamma factor, poles); no zeros.
Usage: python3 mpgap.py MODEL L F TAU   (needs mp80/W_<MODEL>_<L>_<M>.pkl and mp80/line_<MODEL>.pkl)."""
import sys, os, json, math, pickle, time, mpmath as mp, flint
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpW
DPS = 90; mp.mp.dps = DPS; flint.ctx.dps = DPS + 10
SIG = mp.mpf(1)/4
EPSS = [mp.mpf('1e-10'), mp.mpf('1e-8'), mp.mpf('1e-6'), mp.mpf('1e-4')]
def M_of(F, L): return int(math.ceil(F*L/math.pi)) + 8
def family(L, F, KF=1.1, widths=(1.5, 1.0)):           # gap.py::family, deduplicated
    out = []
    for wmax in widths:
        w = min(wmax, L/2*0.98); cs = []
        if L/2 - w > -L/2 + w:
            x = -L/2 + w
            while x <= L/2 - w + 1e-9: cs.append(round(x, 12)); x += w/2
        else: cs = [0.0]
        dk = math.pi/(2*w); ks = []; k = 0.0
        while k <= KF*F + 1e-9: ks.append(round(k, 2)); k += dk
        for c in cs:
            for k in ks:
                for ph in ((0.0,) if k == 0 else (0.0, math.pi/2)):
                    f = (c, w, k, ph)
                    if f not in out: out.append(f)
    return out
def A2F(x): return flint.arb(mp.nstr(x, DPS+5))
def F2M(x): return mp.mpf(x.mid().str(DPS+5, radius=False))
def Hhat_acb(f, s):          # Mellin of h(v) = exp(-(v-c)^2/(2 sg^2)) cos(k(v-c)+phi): int h e^{s v} dv, s complex (acb)
    c, w, k, ph = f; sg = flint.arb(SIG_STR)*flint.arb(repr(w)); k = flint.arb(repr(k)); c = flint.arb(repr(c))
    I = flint.acb(0, 1); e = flint.acb(0, ph).exp() if ph else flint.acb(1)
    pref = (s*c).exp()*flint.arb(2*flint.arb.pi()).sqrt()*sg/2
    return pref*(e*(sg*sg*(s+I*k)**2/2).exp() + (sg*sg*(s-I*k)**2/2).exp()/e)
SIG_STR = '0.25'
def build(model, L, F, h_step=1):
    fam = family(L, F); M = M_of(F, L)
    d = pickle.load(open(os.path.join(HERE, 'mp80', f'line_{model}.pkl'), 'rb'))
    idx = list(range(0, d['NT'], h_step)); h = flint.arb('0.1')*h_step
    ts = [flint.arb(k)*flint.arb('0.1') for k in idx]
    Fl = [flint.acb(flint.arb(d['F'][k][0]), flint.arb(d['F'][k][1])) for k in idx]
    if os.environ.get('CONTROL'):   # same-band NON-radical control (as 14E bilinear.py): F' = F - 0.1*2^{-s}, s = 1/2+it
        Fl = [Fv - flint.arb('0.1')*(-(flint.acb(flint.arb(1)/2, t))*flint.arb(2).log()).exp() for Fv, t in zip(Fl, ts)]
    wt = [h]*len(ts); wt[0] = h/2; wt[-1] = h/2
    sw = [x.sqrt() for x in wt]
    Lf = flint.arb(repr(L)); nrm = (2/Lf).sqrt(); I = flint.acb(0, 1)
    # G[i][t] = Hhat_i(1/2+it) F(1/2+it) * sqrt(wt);  P[j][t] = phihat_j(t) * sqrt(wt)
    Gr, Gi = [], []
    for f in fam:
        row = [Hhat_acb(f, flint.acb(flint.arb(1)/2, t))*Fv*s_ for t, Fv, s_ in zip(ts, Fl, sw)]
        Gr.append([z.real for z in row]); Gi.append([z.imag for z in row])
    Pr, Pi = [], []
    for j in range(1, M+1):
        kj = flint.arb(j)*flint.arb.pi()/Lf; sgn = -1 if j % 2 else 1
        row = [(-I*t*Lf/2).exp()*nrm*kj*(1 - sgn*(I*t*Lf).exp())/(kj*kj - t*t)*s_ for t, s_ in zip(ts, sw)]
        Pr.append([z.real for z in row]); Pi.append([z.imag for z in row])
    Gr, Gi, Pr, Pi = map(flint.arb_mat, (Gr, Gi, Pr, Pi))
    pi = flint.arb.pi()
    Gtot = (Gr*Gr.transpose() + Gi*Gi.transpose())*(1/pi)       # ||E H||^2 Gram (real parts; H real)
    A = (Pr*Gr.transpose() + Pi*Gi.transpose())*(1/pi)          # a_ij = <phi_i, E_F H_j>
    poles = [mp.mpf(1)/2] if model == 'ZETA' else [mp.mpf(1)/2 + mp.mpf('0.2'), mp.mpf(1)/2 - mp.mpf('0.2')]
    Cn = mp.matrix([[F2M(Hhat_acb(f, flint.acb(A2F(mp.mpf(1)/2 + a))).real) for f in fam] for a in poles])
    tomp = lambda X: mp.matrix([[F2M(X[i, j]) for j in range(X.ncols())] for i in range(X.nrows())])
    return fam, M, tomp(Gtot), tomp(A), Cn, max(float(abs(x.rad())) for x in Gtot.entries())
def eigsy_sorted(X):
    E, U = mp.eigsy(X); o = sorted(range(len(E)), key=lambda i: E[i])
    return [E[i] for i in o], U.extract(list(range(U.rows)), o) if hasattr(U, 'extract') else mp.matrix([[U[r, i] for i in o] for r in range(U.rows)])
def cols(U, js): return mp.matrix([[U[r, j] for j in js] for r in range(U.rows)])
def nullspace_rows(C):       # orthonormal basis of {x : C x = 0}
    n = C.cols; r = C.rows
    Q, R = mp.qr(C.T, mode='full'); return cols(Q, list(range(r, n)))
def complement(S):
    Q, R = mp.qr(S, mode='full'); return cols(Q, list(range(S.cols, S.rows)))
def run(model, L, F, tau):
    t0 = time.time()
    fam, M, Gtot, A, Cn, rad = build(model, L, F)
    W = mpW.load(model, L, M); mp.mp.dps = DPS
    evW, UW = eigsy_sorted(W); nW = max(abs(evW[0]), abs(evW[-1]))
    B = nullspace_rows(Cn); Gt = B.T*Gtot*B; At = A*B
    eg, Ug = eigsy_sorted(Gt); gmax = eg[-1]; keep = [i for i in range(len(eg)) if eg[i] > tau*gmax]
    Z = mp.matrix([[Ug[r, i]/mp.sqrt(eg[i]) for i in keep] for r in range(Ug.rows)])
    AZ = At*Z; rho, V = eigsy_sorted(AZ.T*AZ)
    base = dict(model=model, L=L, F=F, M=M, tau=float(tau), sig_over_w=float(SIG), nfam=len(fam), rank_family=len(keep),
                gram_ball_radius_max=rad, rho_top=[mp.nstr(1-x, 6) for x in rho[-12:][::-1]], normW=float(nW),
                low_W=[mp.nstr(x, 8) for x in evW[:10]], n_neg_W=sum(1 for x in evW if x < -1e-7*nW))
    rows = []
    for eps in EPSS:
        sel = [i for i in range(len(rho)) if rho[i] >= 1 - eps]
        r = dict(base, eps=float(eps), nconc=len(sel))
        if not sel:
            r.update(dimS=0, lam_min_restr=mp.nstr(evW[0], 8), lam_min_quot=mp.nstr(evW[0], 8), interlace_upper=mp.nstr(evW[0], 8)); rows.append(r); continue
        vecs = AZ*cols(V, sel)
        es, Us = eigsy_sorted(vecs.T*vecs); smax = es[-1]
        good = [i for i in range(len(es)) if es[i] > mp.mpf('1e-16')*smax]
        S = vecs*cols(Us, good); S = mp.matrix([[S[a, b]/mp.sqrt(es[good[b]]) for b in range(len(good))] for a in range(M)])
        k = S.cols; Q = complement(S)
        WS = W*S; WSS = S.T*WS; R = Q.T*W*Q; WQS = Q.T*WS
        er = eigsy_sorted(R)[0]
        bil = mp.sqrt(max(eigsy_sorted(WS.T*WS)[0]))/nW
        ess = eigsy_sorted(WSS)[0]
        r.update(dimS=k, bilinear_resid_rel=mp.nstr(bil, 4), W_on_S=[mp.nstr(x, 4) for x in ess],
                 lam_min_restr=mp.nstr(er[0], 8), low_restr=[mp.nstr(x, 5) for x in er[:6]], n_neg_restr=sum(1 for x in er if x < -1e-7*nW),
                 interlace_upper=mp.nstr(evW[k], 8), lam_k1_W=mp.nstr(evW[k], 8))
        if all(x != 0 for x in ess):
            Sc = R - WQS*mp.inverse(WSS)*WQS.T; Sc = (Sc + Sc.T)/2; eq = eigsy_sorted(Sc)[0]
            r.update(lam_min_quot=mp.nstr(eq[0], 8), low_quot=[mp.nstr(x, 5) for x in eq[:6]], n_neg_quot=sum(1 for x in eq if x < -1e-7*nW),
                     n_neg_WSS=sum(1 for x in ess if x < 0))
        cosv = eigsy_sorted((S.T*cols(UW, list(range(k)))) * (S.T*cols(UW, list(range(k)))).T)[0]
        r['sin_max_angle_S_vs_bottom_k_eigvecs'] = mp.nstr(mp.sqrt(max(mp.mpf(0), 1 - cosv[0])), 4)
        # energy fraction of S inside the bottom-K eigenspace, K = 2k (does S at least live in the near-null cluster?)
        K2 = min(M, 2*k); c2 = eigsy_sorted((S.T*cols(UW, list(range(K2)))) * (S.T*cols(UW, list(range(K2)))).T)[0]
        r['sin_max_angle_S_vs_bottom_2k'] = mp.nstr(mp.sqrt(max(mp.mpf(0), 1 - c2[0])), 4)
        r['sec'] = round(time.time() - t0); rows.append(r)
        print(json.dumps({x: r[x] for x in ['model', 'L', 'F', 'M', 'tau', 'eps', 'dimS', 'bilinear_resid_rel', 'lam_min_restr', 'lam_min_quot', 'interlace_upper',
                                             'sin_max_angle_S_vs_bottom_k_eigvecs', 'n_neg_W', 'n_neg_restr', 'sec'] if x in r}), flush=True)
    return rows
if __name__ == '__main__':
    model, L, F, tau = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), mp.mpf(sys.argv[4])
    rows = run(model, L, F, tau)
    json.dump(rows, open(os.path.join(HERE, 'mp80', f'gap_{model}_{L}_{F}_{sys.argv[4]}' + ('_control' if os.environ.get('CONTROL') else '') + '.json'), 'w'), indent=1)
