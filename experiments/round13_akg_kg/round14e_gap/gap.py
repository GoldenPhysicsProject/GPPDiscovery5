"""Round 14E: Weil form on the complement of the window-truncated range of the Connes-Consani map E_F (Connes 1999 quotient).
Window: g supported in u = log x in [-L/2, L/2]; Galerkin: M Dirichlet sines (round14b weil_galerkin; W_F from the model's own
integer data: Dirichlet coefficients, gamma factor, poles).  E_F H(u) = e^{u/2} sum_n a_F(n) H(n e^u),  (E_F H)^(t) = F(1/2+it) H^(t).
Family: H = bump((v-c)/w) cos(k(v-c)+phi) in v = log y, centres c in the window, k <= KF*F; pole constraints int h e^{s v}dv = 0
at F's poles.  Window coefficients a_i = <phi_i, E_F H> by Parseval from F on the critical line (cline_*.npy; data, no zeros).
Window-truncated range S_L(eps): span of a(H) over H in the family with ||P_M P_window E_F H||^2 >= (1-eps)||E_F H||^2
(generalised eigenproblem, prolate-type concentration).  Report: radical residual of W on S, lambda_min of W on S^perp."""
import sys, os, math, json, numpy as np, scipy.linalg as sl
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'round14b_euler')); sys.path.insert(0, os.path.join(HERE, '..'))
import weil_galerkin as wg, models14b as md

def phihat(L, M, t):  # int phi_j(u) e^{itu} du, closed form (same as round14d pencil.py)
    k = (np.arange(1, M+1)*np.pi/L)[:, None]; j = np.arange(1, M+1)[:, None]; t = t[None, :]
    den = k**2 - t**2; den = np.where(np.abs(den) < 1e-9, 1e-9, den)
    return np.exp(-1j*t*L/2)*math.sqrt(2/L)*k*(1 - (-1.0)**j*np.exp(1j*t*L))/den
_line = {}
def line(name):
    if name not in _line: _line[name] = np.load(os.path.join(HERE, f'cline_{name}.npy'))[:, ::TSUB]
    return _line[name]
TSUB = int(os.environ.get('TSUB', '4'))   # t-step 0.01*TSUB; Re(integrand) is even in t, so trapezoid on [0,T] is spectrally accurate
def bump(s): return np.where(np.abs(s) < 1, np.exp(-1/np.maximum(1-s*s, 1e-300)), 0.0)
SG = np.linspace(-1, 1, 801)[1:-1]; DS = SG[1]-SG[0]
def family(L, F, KF=0.75, wmax=1.5):
    """returns list of (c, w, k, phi)."""
    w = min(wmax, L/2*0.98); cs = np.arange(-L/2 + w, L/2 - w + 1e-9, w/2) if L/2 - w > -L/2 + w else np.array([0.0])
    dk = math.pi/(2*w); ks = np.round(np.arange(0, KF*F + 1e-9, dk), 2)
    return [(c, w, k, ph) for c in cs for k in ks for ph in ((0.0,) if k == 0 else (0.0, math.pi/2))], w
def Hhat(fam, t):
    """H^(t) = int h(v) e^{(1/2+it)v} dv = e^{(1/2+it)c} [e^{i phi} K_w(t+k) + e^{-i phi} K_w(t-k)]/2,
    K_w(tau) = w int bump(s) e^{w s/2} e^{i tau w s} ds  (k on the t-grid, so shifts are exact index shifts)."""
    dt = t[1]-t[0]; w = fam[0][1]; kmax = max(f[2] for f in fam)
    nk = int(round(kmax/dt)) + 2; tau = dt*np.arange(-nk - len(t), nk + len(t) + 1); off = nk + len(t)
    K = np.zeros(len(tau), complex); hb = bump(SG)*np.exp(w*SG/2)*w*DS
    for a in range(0, len(tau), 4000): K[a:a+4000] = np.exp(1j*np.outer(tau[a:a+4000], w*SG)) @ hb
    idx = np.arange(len(t)); out = np.zeros((len(fam), len(t)), complex)
    for j, (c, w_, k, ph) in enumerate(fam):
        ik = int(round(k/dt))
        out[j] = np.exp((0.5+1j*t)*c)*0.5*(np.exp(1j*ph)*K[off+idx+ik] + np.exp(-1j*ph)*K[off+idx-ik])
    return out
def Hmellin(fam, s):
    return np.array([np.sum(bump(SG)*np.cos(k*(w*SG)+ph)*np.exp(s*(c+w*SG)))*w*DS for (c, w, k, ph) in fam])
def wmatrix(mo, L, M):
    a, p, q = wg.operators(mo, L, M); W = a+p-q; return 0.5*(W+W.T)
WIDTHS = (1.5, 1.0, 0.6); KF = 1.2
def radical_space(name, mo, L, M, F, eps, widths=None, kf=None):
    t, Fl = line(name); t = t.real
    fams = [family(L, F, KF=kf or KF, wmax=wm)[0] for wm in (widths or WIDTHS)]
    fam = [f for fa in fams for f in fa]; w = sorted(set(f[1] for f in fam))
    Hh = np.vstack([Hhat(fa, t) for fa in fams])
    # pole constraints: H^(s) = 0 at s = 1/2 + a, a in poles
    Cn = np.array([Hmellin(fam, 0.5 + a) for a, _ in mo.get('poles', [])])
    B = np.eye(len(fam))
    if len(Cn): B = sl.null_space(Cn)
    G = Hh*Fl[None, :]                       # (E_F H)^ on t >= 0
    wt = np.full(len(t), t[1]-t[0]); wt[0] *= 0.5; wt[-1] *= 0.5
    Gtot = (G.conj()*wt) @ G.T; Gtot = (Gtot.real)/math.pi    # ||E H||^2 Gram (real parts; H real)
    P = phihat(L, M, t)
    A = ((P*wt) @ G.conj().T).real/math.pi                    # a_i = <phi_i, E_F H_j>
    Gt = B.T @ Gtot @ B; At = A @ B
    ev, U = np.linalg.eigh(Gt); keep = ev > 1e-13*ev.max(); Z = U[:, keep]/np.sqrt(ev[keep])   # Gtot-whitening
    rho, V = np.linalg.eigh(Z.T @ (At.T @ At) @ Z)
    sel = rho >= 1 - eps
    vecs = At @ Z @ V[:, sel]                 # window coefficient vectors of concentrated E_F H
    if vecs.shape[1] == 0: return np.zeros((M, 0)), dict(nfam=len(fam), w=w, nconc=0, rho_top=float(rho.max()) if len(rho) else None)
    Uq, sv, _ = np.linalg.svd(vecs, full_matrices=False); S = Uq[:, sv > 1e-8*sv.max()]
    return S, dict(nfam=len(fam), w=w, nconc=int(sel.sum()), dimS=int(S.shape[1]), rho_top=float(rho.max()), rho_cut=float(rho[sel].min()))
def analyse(name, mo, L, M, F, eps, W=None, S=None, info=None, **kw):
    W = wmatrix(mo, L, M) if W is None else W
    if S is None: S, info = radical_space(name, mo, L, M, F, eps, **kw)
    nW = np.abs(np.linalg.eigvalsh(W)).max(); evW = np.linalg.eigvalsh(W)
    r = dict(L=L, M=M, F=F, eps=eps, **info, lam_min_W=float(evW[0]), n_tiny_W=int((np.abs(evW) < 1e-8*nW).sum()), normW=float(nW),
             n_neg_W=int((evW < -1e-7*nW).sum()))
    if S.shape[1]:
        ws = np.linalg.eigvalsh(S.T @ W @ S); r['radical_resid_rel'] = float(np.abs(ws).max()/nW)
        r['bilinear_resid_rel'] = float(np.linalg.norm(W @ S, 2)/nW)   # sup |W(v,phi)|/(|v||phi||W|), v in S, phi in V_M
        Q = sl.null_space(S.T)
    else:
        r['radical_resid_rel'] = 0.0; Q = np.eye(M)
    r['dim_perp'] = int(Q.shape[1])
    if Q.shape[1]:
        e = np.linalg.eigvalsh(Q.T @ W @ Q); r['lam_min_perp'] = float(e[0]); r['lam_min_perp_rel'] = float(e[0]/nW)
        r['n_neg_perp'] = int((e < -1e-7*nW).sum()); r['lam2_perp'] = float(e[1]) if len(e) > 1 else None
        r['low_perp'] = [float(x) for x in e[:8]]
    r['low_W'] = [float(x) for x in evW[:8]]
    return r
