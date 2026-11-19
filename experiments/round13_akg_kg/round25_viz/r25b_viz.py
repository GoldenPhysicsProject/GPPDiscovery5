#!/usr/bin/env python3
"""Round 25b: mirror-pair level curves (fig7) and Weil form vs hook (fig8).

Fig 7 reuses the Round 25 grids. Fig 8 evaluates the Weil form W(g*g̃) for a
Fejér wave packet via the R22 explicit formula (arithmetic side: primes +
archimedean + poles). Zeros are used only as plot marks and as a NUMERICAL
check of prime-side vs zero-side, never as proof inputs. RH remains OPEN.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import mpmath as mp
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

import r25_viz as r25
from r25_viz import (
    BG, FG, MUTED, GOLD, CRIMSON, CYAN, VIOLET, GREEN, CRIT, LINE_IM,
    DATA, OUT, THETA_FTH,
    style_ax, save, eval_grid, quadrant_rgba, xi_fth,
)
from r22_models import MODELS, lam_table

fp = mp.fp

# Two shades of the red family: left of the line (σ < 1/2) and right (σ > 1/2).
RE_LEFT = "#f08a4a"   # coral / apricot, σ < 1/2
RE_RIGHT = "#c81e58"  # deep rose,       σ > 1/2

# ξ(1/2+it) extrema in the fig5 ζ window, recomputed in r25.extrema_of_xi
# and re-checked below. These are the previously unlabelled violet marks.
XI_EXT_FIG5 = np.array([15.585708589924199, 22.097977280453488])

# Known plot marks (illustration only).
ZETA_IM = np.array([
    14.134725141734675, 21.022039638771552, 25.01085758014569,
    30.424876125859516, 32.93506158773919, 37.586178158825675,
])
EPS_OFF = [
    (0.9329696974854141, 15.668249531278473),
    (0.9376669067003994, 29.983395235155516),
    (0.6969270453122914, 36.37406368638513),
    (0.8231873164780866, 44.000113180236895),
]
EPS_ON = [
    1.4562032347251894, 5.326203936204837, 7.11349023725727,
    9.208567035966308, 10.571063256663962, 12.494830631904136,
    13.098509365914271, 17.596050852287988, 18.90136055238456,
    19.80901394515022,
]
DH_OFF = [(0.8085171824566374, 85.69934848537759)]
DH_ON = [80.32794708609034, 81.92540239183100, 83.10884108474982, 87.64747633255974]


# =====================================================================
# Symmetry of the grids, and shared drawing
# =====================================================================

def symmetry_residual(sigs, ts, W):
    """Max relative residual of ξ(1-σ+it) = conj ξ(σ+it) on a σ-symmetric grid."""
    # expected pairing: sigs[j] + sigs[::-1][j] == 1
    mid_err = float(np.max(np.abs(sigs + sigs[::-1] - 1.0)))
    Wp = W[:, ::-1]
    denom = np.maximum(np.abs(W), 1e-30)
    rel = np.abs(W - np.conj(Wp)) / denom
    rel = rel[np.isfinite(rel)]
    max_rel = float(np.max(rel)) if rel.size else float("nan")
    med_rel = float(np.median(rel)) if rel.size else float("nan")
    # also split Re even / Im odd
    Re_err = np.abs(np.real(W) - np.real(Wp))
    Im_err = np.abs(np.imag(W) + np.imag(Wp))
    absW = np.maximum(np.abs(W), 1e-30)
    max_re = float(np.nanmax(Re_err / absW))
    max_im = float(np.nanmax(Im_err / absW))
    return dict(grid_mid_err=mid_err, max_rel=max_rel, med_rel=med_rel,
                max_re_even=max_re, max_im_odd=max_im)


def load_cached(name, sigs, ts):
    path = DATA / f"grid_xi_{name}.npy"
    z = np.load(path)
    if z.shape != (len(ts), len(sigs)):
        raise RuntimeError(f"{path} has shape {z.shape}, expected {(len(ts), len(sigs))}")
    return z


def draw_mirror_panel(ax, sigs, ts, W, on_line, off_line, extrema, title,
                      tlim=None, show_extrema=False, hook_ann=None,
                      square_ann=False):
    """Level curves with Re ξ=0 split by side of the critical line."""
    rgba = quadrant_rgba(W, alpha=0.52)
    ax.imshow(
        np.transpose(rgba, (1, 0, 2)), origin="lower", aspect="auto",
        extent=[ts[0], ts[-1], sigs[0], sigs[-1]], interpolation="nearest",
        zorder=0,
    )
    Re = np.real(W).T  # (nsig, nt)
    Im = np.imag(W).T
    TT, SS = np.meshgrid(ts, sigs)
    # two shades: include the line in both so they can meet there
    Re_L = np.ma.array(Re, mask=SS > 0.5 + 1e-15)
    Re_R = np.ma.array(Re, mask=SS < 0.5 - 1e-15)
    ax.contour(TT, SS, Re_L, levels=[0], colors=[RE_LEFT], linewidths=1.25, zorder=3)
    ax.contour(TT, SS, Re_R, levels=[0], colors=[RE_RIGHT], linewidths=1.25, zorder=3)
    Im_m = np.ma.array(Im, mask=np.abs(SS - 0.5) < 0.02)
    ax.contour(TT, SS, Im_m, levels=[0], colors=[LINE_IM], linewidths=1.05, zorder=3)
    ax.axhline(0.5, color=CRIT, lw=1.35, zorder=4, alpha=0.95)

    if on_line is not None:
        on = [t for t in on_line if ts[0] <= t <= ts[-1]]
        ax.scatter(on, [0.5] * len(on), s=28, color=GOLD, zorder=7,
                   edgecolors=BG, linewidths=0.4)
        if square_ann and on:
            # small right-angle ticks at each on-line zero, to mark a square crossing
            x0, x1 = (tlim if tlim else (ts[0], ts[-1]))
            y0, y1 = sigs[0], sigs[-1]
            dt = 0.014 * (x1 - x0)
            ds = 0.055 * (y1 - y0)
            for t0 in on:
                ax.plot([t0, t0], [0.5 - ds, 0.5 + ds], color=CRIT, lw=0.7,
                        alpha=0.85, zorder=8, solid_capstyle="butt")
            # one explicit right-angle glyph at the first zero
            t0 = on[0]
            ax.plot([t0, t0 + dt], [0.5 + 0.6 * ds, 0.5 + 0.6 * ds],
                    color=FG, lw=1.05, zorder=9, solid_capstyle="butt")
            ax.plot([t0 + dt, t0 + dt], [0.5, 0.5 + 0.6 * ds],
                    color=FG, lw=1.05, zorder=9, solid_capstyle="butt")
            ax.text(t0 + 1.35 * dt, 0.5 + 0.72 * ds, "square crossing",
                    color=FG, fontsize=7.2, va="bottom", ha="left",
                    zorder=9)

    if off_line:
        for sig, t in off_line:
            if ts[0] <= t <= ts[-1] and sigs[0] <= sig <= sigs[-1]:
                ax.scatter([t, t], [sig, 1 - sig], s=34, marker="o",
                           facecolors=CRIMSON, edgecolors="white",
                           linewidths=0.4, zorder=7)

    if show_extrema and extrema is not None:
        ex = [t for t in extrema if ts[0] <= t <= ts[-1]]
        ax.scatter(ex, [0.5] * len(ex), s=20, marker="D", color=VIOLET, zorder=8,
                   edgecolors=BG, linewidths=0.35)
        # label the two marks Turbo asked about
        if ex:
            ax.text(ex[0] + 0.18, 0.5 + 0.10, r"extremum of $\xi(\frac{1}{2}+it)$",
                    color=VIOLET, fontsize=7.0, va="bottom", ha="left", zorder=9)

    if hook_ann:
        for item in hook_ann:
            ax.annotate(
                item["text"],
                xy=item["xy"],
                xytext=item.get("xytext", (item["xy"][0] - 1.4, item["xy"][1] + 0.18)),
                color=FG, fontsize=7.2,
                arrowprops=dict(arrowstyle="-|>", color=FG, lw=0.8,
                                shrinkA=0, shrinkB=2),
                ha=item.get("ha", "left"), va="bottom", zorder=10,
                bbox=dict(boxstyle="round,pad=0.18", fc="#161b24", ec="#3a4250",
                          alpha=0.88),
            )

    style_ax(ax, r"$t=\mathrm{Im}\,s$", r"$\sigma=\mathrm{Re}\,s$", title)
    if tlim:
        ax.set_xlim(*tlim)
    else:
        ax.set_xlim(ts[0], ts[-1])
    ax.set_ylim(sigs[0], sigs[-1])
    ax.set_aspect("auto")


def fig7_mirror_pairs(bundles, zoom, residuals, xi_ext):
    n = len(bundles)
    fig, axes = plt.subplots(
        n + 1, 1,
        figsize=(12.2, 3.15 * n + 4.6),
        sharex=False,
        gridspec_kw={"height_ratios": [1] * n + [1.35]},
    )
    fig.suptitle(
        r"A red curve meeting its mirror on the line — or hooking away from it",
        color=FG, fontsize=13.5, x=0.03, ha="left", y=0.995,
    )
    fig.text(
        0.03, 0.972,
        r"Same windows as Fig. 5. Coral: $\mathrm{Re}\,\xi=0$ on $\sigma<\frac{1}{2}$; rose: on $\sigma>\frac{1}{2}$. "
        r"The functional equation with real coefficients gives $\xi(1-\sigma+it)=\overline{\xi(\sigma+it)}$, so the two shades are mirrors.",
        color=MUTED, fontsize=8.4, ha="left", va="top",
    )
    fig.text(
        0.03, 0.954,
        r"In $\zeta$ they meet on the line, squarely, once each. In every fake they turn away from each other, and the tip of the hook catches an off-line blue branch.",
        color=MUTED, fontsize=8.4, ha="left", va="top",
    )

    for ax, B in zip(axes[:n], bundles):
        draw_mirror_panel(
            ax, B["sigs"], B["ts"], B["W"], B["on"], B["off"], B.get("extrema"),
            B["name"], tlim=(B["ts"][0], B["ts"][-1]),
            show_extrema=B.get("extrema") is not None,
            hook_ann=B.get("hook_ann"),
            square_ann=B.get("square_ann", False),
        )
        if B.get("note"):
            ax.text(
                0.99, 0.06, B["note"], transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color=FG,
                bbox=dict(boxstyle="round,pad=0.25", fc="#161b24", ec="#3a4250", alpha=0.85),
            )

    # zoom on the first F_θ fold
    axz = axes[-1]
    Z = zoom
    draw_mirror_panel(
        axz, Z["sigs"], Z["ts"], Z["W"], None, Z["off"], None,
        Z["name"], tlim=(Z["ts"][0], Z["ts"][-1]),
        hook_ann=Z.get("hook_ann"), square_ann=False,
    )
    # guide lines at σ = 1/4, 3/4 and at t = γ
    axz.axhline(0.25, color=CRIT, lw=0.5, ls=":", alpha=0.7, zorder=4)
    axz.axhline(0.75, color=CRIT, lw=0.5, ls=":", alpha=0.7, zorder=4)
    g1 = float(ZETA_IM[0])
    axz.axvline(g1, color=CRIT, lw=0.5, ls=":", alpha=0.7, zorder=4)
    axz.text(0.01, 0.06, Z["note"], transform=axz.transAxes,
             ha="left", va="bottom", fontsize=8, color=FG,
             bbox=dict(boxstyle="round,pad=0.25", fc="#161b24", ec="#3a4250", alpha=0.85))

    handles = [
        Line2D([0], [0], color=RE_LEFT, lw=2.2, label=r"$\mathrm{Re}\,\xi=0$,  $\sigma<\frac{1}{2}$"),
        Line2D([0], [0], color=RE_RIGHT, lw=2.2, label=r"$\mathrm{Re}\,\xi=0$,  $\sigma>\frac{1}{2}$"),
        Line2D([0], [0], color=LINE_IM, lw=2, label=r"$\mathrm{Im}\,\xi=0$"),
        Line2D([0], [0], color=CRIT, lw=2.2, label=r"critical line (itself $\mathrm{Im}\,\xi=0$)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GOLD,
               markersize=7, label="on-line zero (square crossing)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CRIMSON,
               markersize=7, label="off-line zero (hook tip)"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=VIOLET,
               markersize=6, label=r"extremum of $\xi(\frac{1}{2}+it)$"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.945),
               fontsize=7.4, labelcolor=FG, ncol=2)

    # residual one-liner, pretty names
    pretty = {"zeta": r"$\zeta$", "F_theta": r"$F_\theta$", "Epstein": "Epstein",
              "DH": "DH", "F_theta_zoom": r"$F_\theta$ zoom"}
    bits = [f"{pretty.get(name, name)} {R['max_rel']:.1e}" for name, R in residuals.items()]
    res_line = ", ".join(bits)

    fig.text(
        0.03, 0.032,
        r"NUMERICAL, illustration. $\xi(1-\sigma+it)=\overline{\xi(\sigma+it)}$ was checked on these grids "
        r"(max relative residual: " + res_line + r"). "
        r"Violet diamonds in the $\zeta$ panel are extrema of $\xi(\frac{1}{2}+it)$ at $t\approx 15.586$ and $t\approx 22.098$; "
        r"they were plotted in Fig. 5 but were missing from that legend.",
        color=MUTED, fontsize=7.6,
    )
    fig.text(
        0.03, 0.012,
        r"The conjugates of the first $F_\theta$ quartet sit at $t=-\gamma$, out of view. "
        r"Reading a hook as a broken pairing / failure of reflection positivity is HEURISTIC: "
        r"the drawing shows the geometry of an off-line zero, not the Weil form itself (see Fig. 8). RH is OPEN.",
        color=MUTED, fontsize=7.6,
    )
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.072, top=0.905, hspace=0.38)

    return save(fig, "fig7_mirror_pairs.png")


# =====================================================================
# Fig 8 — Weil form vs hook, Fejér packets, R22 explicit formula
# =====================================================================

EG = 0.5772156649015329
_XG, _WG = np.polynomial.legendre.leggauss(12)


def fejer_K(w, L):
    """Fejér kernel K(r)=(L/2π) sinc²(r L/2). Entire, so w may be complex."""
    w = np.asarray(w)
    y = w * L / 2
    out = np.empty(y.shape, dtype=np.result_type(y.dtype, float))
    small = np.abs(y) < 1e-12
    out = (L / (2 * np.pi)) * np.where(small, 1.0 + 0j, (np.sin(y) / np.where(y == 0, 1, y)) ** 2)
    return out


def arch_I(a, b, L, g0):
    """(1/2π) ∫ h(r) ψ(a + i b r) dr for the R22 Fejér packet, copied from r22_fejer.py."""
    g00 = 1 / np.pi
    T = L / b
    edges = np.linspace(0, T, int(T / 0.02) + 2)
    t = ((edges[:-1, None] + edges[1:, None]) / 2
         + (edges[1:, None] - edges[:-1, None]) / 2 * _XG[None, :]).ravel()
    w = (((edges[1:] - edges[:-1]) / 2)[:, None] * _WG[None, :]).ravel()
    gbt = np.cos(np.outer(g0, b * t)) * (1 - b * t / L) / np.pi
    integ = (np.exp(-t) * g00 - np.exp(-a * t) * gbt) / (-np.expm1(-t))
    return -EG * g00 + integ @ w - g00 * math.log(-math.expm1(-T))


def rho_fejer(name, X, g0, lam):
    """Arithmetic side of ∑_ρ h(γ_ρ) for the Fejér packet of length L=log X peaked at g0.

    This is the Weil form W(g*g̃) in the R22 normalisation (see r22_models.py / r22_fejer.py).
    No zeros are inputs.
    """
    M = MODELS[name]
    L = math.log(X)
    g0 = np.asarray(g0, dtype=float)
    out = 2 * (1 / np.pi) * M["logQ"] * np.ones_like(g0)
    for (lmb, a, b) in M["arch"]:
        out = out + 2 * lmb * arch_I(a, b, L, g0)
    for z in M["poles"]:
        r = -1j * (z - 0.5)
        out = out + 2 * (fejer_K(r - g0, L) + fejer_K(r + g0, L)).real
    n = np.arange(2, int(X) + 1)
    lm = lam[2:int(X) + 1]
    nz = lm != 0
    n = n[nz]
    c = lm[nz] * n ** -0.5 * (1 - np.log(n) / L)
    out = out - (2 / np.pi) * (np.cos(np.outer(g0, np.log(n))) @ c)
    return out


def rho_zero_side(zeros_gamma, X, g0):
    """Spectral side ∑_ρ h(γ_ρ), zeros as a NUMERICAL check only.

    zeros_gamma is a list of complex γ_ρ = (ρ − 1/2)/i  (real for on-line zeros).
    The packet is even, so each γ and −γ should be included, or we add K(γ−g0)+K(γ+g0)
    for γ in the upper half and count both.
    Here we assume `zeros_gamma` lists one representative per {γ, −γ} pair (t ≥ 0)
    and we evaluate h(γ) = K(γ−g0)+K(γ+g0), which already folds in −γ.
    """
    L = math.log(X)
    g0 = np.asarray(g0, dtype=float)
    # Each listed γ (t≥0) has a conjugate partner at −γ; h is even, so
    # the pair contributes 2 h(γ) = 2(K(γ−g0)+K(γ+g0)).
    acc = np.zeros(g0.shape, dtype=complex)
    for g in zeros_gamma:
        acc = acc + 2 * (fejer_K(g - g0, L) + fejer_K(g + g0, L))
    return acc.real


def lobe_width(X):
    return 2 * math.pi / math.log(X)


def scan_rho(name, Xs, tmin, tmax, step, lam_cache):
    g0 = np.arange(tmin, tmax + 0.5 * step, step)
    rows = {}
    for X in Xs:
        lam = lam_cache[X]
        chunks = []
        for i in range(0, len(g0), 250):
            chunks.append(rho_fejer(name, X, g0[i:i + 250], lam))
        rows[X] = np.concatenate(chunks)
        print(f"  {name} X={X}: min={rows[X].min():.4f} at t={g0[rows[X].argmin()]:.3f}",
              flush=True)
    return g0, rows


def fig8_weil_vs_hook(scans, check):
    """scans: list of dicts with name, g0, rows {X: W}, hooks, on, Xs, note."""
    n = len(scans)
    fig, axes = plt.subplots(n, 1, figsize=(12.2, 3.15 * n + 1.6), sharex=False)
    if n == 1:
        axes = [axes]
    fig.suptitle(
        r"Does a hook’s birth line up with a sign change of the Weil form?  "
        r"$W_X(t)=W(g_t\ast\widetilde g_t)$, Fejér packets, arithmetic side",
        color=FG, fontsize=13.2, x=0.03, ha="left", y=0.995,
    )
    fig.text(
        0.03, 0.968,
        r"Each curve is the R22 explicit formula for a Fejér packet of Fourier support $[-\log X,\log X]$, "
        r"peaked at height $t$ on the critical line, lobe width $2\pi/\log X$. No zeros are inputs.",
        color=MUTED, fontsize=8.4, ha="left", va="top",
    )
    fig.text(
        0.03, 0.952,
        r"Gold ticks: on-line zeros. Crimson dashed: hook heights (plot marks only). "
        r"The answer, for this family: $W$ peaks at the hook, and the sign change lives in the shoulders once the packet is narrow enough. $\zeta$ stays positive.",
        color=MUTED, fontsize=8.4, ha="left", va="top",
    )

    # distinct but related colours for increasing X (narrower packets)
    xcolors = {
        30: "#7cba6f",
        100: CYAN,
        300: "#8e7cc3",
        1000: GOLD,
        3000: CRIMSON,
    }

    for ax, S in zip(axes, scans):
        g0 = S["g0"]
        for X in S["Xs"]:
            W = S["rows"][X]
            lw = 1.25 if X != 1000 else 1.45
            ax.plot(g0, W, color=xcolors.get(X, FG), lw=lw,
                    label=fr"$X={X}$  (lobe ${lobe_width(X):.2f}$)", zorder=3)
            ax.fill_between(g0, W, 0, where=(W < 0), color=xcolors.get(X, CRIMSON),
                            alpha=0.28, linewidth=0, zorder=2)
        ax.axhline(0, color=CRIT, lw=0.8, alpha=0.85, zorder=1)
        y0, y1 = ax.get_ylim()
        # start with automatic ylim then pad after marks
        for t in S.get("on") or []:
            if g0[0] <= t <= g0[-1]:
                ax.axvline(t, color=GOLD, lw=0.7, alpha=0.7, zorder=2)
        for t in S.get("hooks") or []:
            if g0[0] <= t <= g0[-1]:
                ax.axvline(t, color=CRIMSON, lw=1.0, alpha=0.9, zorder=2, ls="--")
        # mark the first hook with a label
        if S.get("hooks"):
            tH = S["hooks"][0]
            if g0[0] <= tH <= g0[-1]:
                ax.text(tH, 0.98, S.get("hook_label", "first hook"),
                        transform=ax.get_xaxis_transform(),
                        color=CRIMSON, fontsize=7.2, ha="left", va="top",
                        rotation=90, rotation_mode="anchor")
        style_ax(ax, r"$t$  (packet centre on the critical line)", r"$W_X(t)$", S["name"])
        ax.legend(loc="upper right", fontsize=7.2, labelcolor=FG, ncol=len(S["Xs"]))
        if S.get("note"):
            ax.text(
                0.01, 0.97, S["note"], transform=ax.transAxes,
                ha="left", va="top", fontsize=7.4, color=FG,
                bbox=dict(boxstyle="round,pad=0.22", fc="#161b24", ec="#3a4250", alpha=0.88),
                zorder=5,
            )
        ax.set_xlim(g0[0], g0[-1])
        # give negative shoulders visible room below 0
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(min(ymin, -0.12 * ymax), ymax * 1.08)

    fig.text(
        0.03, 0.032,
        r"NUMERICAL. $W_X$ is the R22 Fejér–Weil form, from each model’s own $\Lambda_M(n)$, archimedean factor and poles "
        r"(arithmetic side of the explicit formula; no zeros as inputs). "
        + check,
        color=MUTED, fontsize=7.6,
    )
    fig.text(
        0.03, 0.012,
        r"The hook height itself is a peak of $W$, not a sign change. The sign change is in the shoulders, "
        r"and it appears only once the lobe is narrow enough to resolve the hook "
        r"($X\geq 30$ for $F_\theta$, $X\geq 100$ for Epstein and DH). "
        r"That alignment is a numerical observation for this packet family, not a theorem that a hook is a failure of reflection positivity. RH is OPEN.",
        color=MUTED, fontsize=7.6,
    )
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.072, top=0.900, hspace=0.42)

    return save(fig, "fig8_weil_vs_hook.png")


# =====================================================================
# Drivers
# =====================================================================

def run_fig7():
    print("=== fig7: load cached grids ===", flush=True)
    sigs_z = np.linspace(-0.10, 1.10, 241)
    ts_z = np.linspace(0.40, 40.0, 661)
    Wz = load_cached("zeta", sigs_z, ts_z)

    sigs_e = np.linspace(-0.05, 1.05, 221)
    ts_e = np.linspace(10.0, 34.0, 401)
    We = load_cached("eps", sigs_e, ts_e)

    sigs_d = np.linspace(-0.05, 1.05, 221)
    ts_d = np.linspace(81.5, 90.0, 321)
    Wd = load_cached("dh", sigs_d, ts_d)

    sigs_f = np.linspace(-0.05, 1.05, 221)
    ts_f = np.linspace(10.0, 26.0, 321)
    Wf = load_cached("fth", sigs_f, ts_f)

    def slice_grid(sigs, ts, W, tlo, thi):
        m = (ts >= tlo) & (ts <= thi)
        return sigs, ts[m], W[m, :]

    s1, t1, w1 = slice_grid(sigs_z, ts_z, Wz, 10.0, 26.0)

    print("=== symmetry residuals ===", flush=True)
    residuals = {}
    for name, sigs, ts, W in [
        ("zeta", s1, t1, w1),
        ("F_theta", sigs_f, ts_f, Wf),
        ("Epstein", sigs_e, ts_e, We),
        ("DH", sigs_d, ts_d, Wd),
    ]:
        R = symmetry_residual(sigs, ts, W)
        residuals[name] = R
        print(f"  {name}: grid_mid_err={R['grid_mid_err']:.2e}  "
              f"max |ξ(1-σ+it)-conj ξ(σ+it)|/|ξ| = {R['max_rel']:.3e}  "
              f"median {R['med_rel']:.3e}  Re-even {R['max_re_even']:.3e}  "
              f"Im-odd {R['max_im_odd']:.3e}", flush=True)

    # finer zoom grid around the first F_θ fold
    print("=== F_theta zoom grid ===", flush=True)
    sigs_zmw = np.linspace(-0.05, 1.05, 261)
    ts_zmw = np.linspace(13.15, 15.15, 321)
    Wzmw = eval_grid(xi_fth, sigs_zmw, ts_zmw, DATA / "grid_xi_fth_zoom.npy", "xi_fth_zoom")
    Rzoom = symmetry_residual(sigs_zmw, ts_zmw, Wzmw)
    residuals["F_theta_zoom"] = Rzoom
    print(f"  F_theta_zoom: max rel {Rzoom['max_rel']:.3e}", flush=True)

    g1 = float(ZETA_IM[0])
    fth_off_zoom = [(0.75, g1), (0.25, g1)]
    fth_off_all = []
    for t in ZETA_IM[:6]:
        fth_off_all.append((0.5 + THETA_FTH, float(t)))
        fth_off_all.append((0.5 - THETA_FTH, float(t)))

    bundles = [
        dict(
            name=r"Riemann $\xi$  — coral and rose meet on the line, squarely, once each",
            sigs=s1, ts=t1, W=w1, on=ZETA_IM, off=None, extrema=XI_EXT_FIG5,
            square_ann=True,
            note=r"violet diamonds: extrema of $\xi(\frac{1}{2}+it)$ at $t\approx 15.586,\,22.098$",
        ),
        dict(
            name=r"$F_\theta=\xi(s+\frac{1}{4})\xi(s-\frac{1}{4})$  — each red curve hooks away from its mirror",
            sigs=sigs_f, ts=ts_f, W=Wf, on=None, off=fth_off_all, extrema=None,
            hook_ann=[
                dict(text="hook", xy=(g1, 0.75), xytext=(g1 - 1.55, 0.88)),
                dict(text="mirror hook", xy=(g1, 0.25), xytext=(g1 - 1.95, 0.05)),
            ],
            note=r"tips at $\frac{1}{2}\pm\frac{1}{4}+i\gamma$; conjugates at $t=-\gamma$ are not in view",
        ),
        dict(
            name=r"Epstein $\xi_Q$, $Q=x^2+5y^2$  — first hook near $0.933+15.668i$ and its mirror",
            sigs=sigs_e, ts=ts_e, W=We, on=EPS_ON, off=EPS_OFF, extrema=None,
            hook_ann=[
                dict(text="hook", xy=(15.668, 0.933), xytext=(12.4, 0.78)),
            ],
            note="one fold event, a quartet: the tips and (out of view) their conjugates",
        ),
        dict(
            name=r"Davenport–Heilbronn  — the hook is the whole story of the first off-line zero",
            sigs=sigs_d, ts=ts_d, W=Wd, on=DH_ON, off=DH_OFF, extrema=None,
            hook_ann=[
                dict(text="hook — curve turns away from its mirror",
                     xy=(85.699, 0.809), xytext=(83.55, 0.90)),
                dict(text="mirror hook", xy=(85.699, 0.191), xytext=(83.85, 0.04)),
            ],
            note=r"tips at $\sigma\approx 0.809$ and $0.191$, same height $t\approx 85.699$",
        ),
    ]
    zoom = dict(
        name=r"Zoom: the first $F_\theta$ fold  (the quartet $\frac{1}{2}\pm\frac{1}{4}+i\gamma$, $\gamma\approx 14.135$)",
        sigs=sigs_zmw, ts=ts_zmw, W=Wzmw, off=fth_off_zoom,
        hook_ann=[
            dict(text=r"tip $\frac{3}{4}+i\gamma$", xy=(g1, 0.75), xytext=(g1 - 0.55, 0.88)),
            dict(text=r"tip $\frac{1}{4}+i\gamma$", xy=(g1, 0.25), xytext=(g1 - 0.55, 0.07)),
        ],
        note=r"conjugates $\frac{1}{2}\pm\frac{1}{4}-i\gamma$ sit at $t=-\gamma\approx-14.13$, out of this window",
    )
    fig7_mirror_pairs(bundles, zoom, residuals, XI_EXT_FIG5)
    plt.close("all")
    return residuals


def _zeros_gamma_zeta(n=80):
    return [complex(float(fp.zetazero(k).imag), 0.0) for k in range(1, n + 1)]


def _zeros_gamma_fth(n=80, theta=THETA_FTH):
    # ρ = 1/2 ± θ ± iγ  →  γ_ρ = (ρ-1/2)/i = ±γ ∓ i θ
    # We list the t≥0 pair at each height:  γ - iθ  and  γ + iθ.
    out = []
    for k in range(1, n + 1):
        g = float(fp.zetazero(k).imag)
        out.append(complex(g, -theta))
        out.append(complex(g, +theta))
    return out


def run_fig8():
    print("=== fig8: Fejér–Weil scans (arithmetic side) ===", flush=True)
    Xs = [30, 100, 300, 1000]
    Xmax = max(Xs)
    lam_cache = {}
    for name in ("ZETA", "FTH", "EPS", "DH"):
        print(f"  lam_table {name} up to {Xmax}", flush=True)
        lam_cache[name] = {Xmax: lam_table(name, Xmax + 1)}
        # rho_fejer slices lam[:X], so one table at Xmax is enough if we pass that array
        for X in Xs:
            lam_cache[name][X] = lam_cache[name][Xmax]

    scanspec = [
        dict(key="ZETA", title=r"Riemann $\xi$  (control): every crossing in this window is on the line",
             tmin=10.0, tmax=26.0, step=0.02,
             hooks=[], on=list(ZETA_IM[:3]),
             hook_label="",
             note="expect $W_X>0$ throughout"),
        dict(key="FTH", title=r"$F_\theta$: every $\zeta$-height in this window is a hook ($=\frac{1}{2}\pm\frac{1}{4}+i\gamma$)",
             tmin=10.0, tmax=26.0, step=0.02,
             hooks=list(ZETA_IM[:3]), on=[],
             hook_label="first hook",
             note=""),
        dict(key="EPS", title=r"Epstein $x^2+5y^2$: first hook at $t\approx 15.668$, ordinary on-line zeros beside it",
             tmin=10.0, tmax=34.0, step=0.02,
             hooks=[15.668249531278473, 29.983395235155516],
             on=[t for t in EPS_ON if 10 <= t <= 34],
             hook_label="first hook",
             note=""),
        dict(key="DH", title=r"Davenport–Heilbronn: first hook at $t\approx 85.699$, on-line zeros on either side",
             tmin=81.5, tmax=90.0, step=0.02,
             hooks=[85.69934848537759],
             on=[t for t in DH_ON if 81.5 <= t <= 90],
             hook_label="first hook",
             note=""),
    ]

    scans = []
    summary = {}
    for S in scanspec:
        print(f"--- scan {S['key']} ---", flush=True)
        g0, rows = scan_rho(S["key"], Xs, S["tmin"], S["tmax"], S["step"], lam_cache[S["key"]])
        # per-X summary: min, at, value at first hook, value at first on-line (if any)
        rowsum = {}
        for X in Xs:
            W = rows[X]
            k = int(np.argmin(W))
            rec = dict(min=float(W[k]), at=float(g0[k]), nneg=int(np.sum(W < 0)),
                       frac_neg=float(np.mean(W < 0)))
            if S["hooks"]:
                # value at first hook and distance of the minimizer to the nearest hook
                tH = S["hooks"][0]
                rec["W_at_first_hook"] = float(np.interp(tH, g0, W))
                rec["dist_min_to_nearest_hook"] = float(min(abs(g0[k] - h) for h in S["hooks"]))
            if S["on"]:
                rec["W_at_first_online"] = float(np.interp(S["on"][0], g0, W))
            rowsum[X] = rec
        summary[S["key"]] = rowsum
        verdicts = {
            "ZETA": r"$W_X>0$ throughout; valleys between zeros stay positive",
            "FTH": r"Peaks sit on the hooks; valleys between them go negative for $X\geq 30$. $W$ is not negative at the hook itself.",
            "EPS": r"The hook is a peak; its shoulders go negative for $X\geq 100$. On-line zeros stay positive.",
            "DH": r"Peak at the hook; both shoulders go negative for $X\geq 100$.",
        }
        scans.append(dict(
            name=S["title"], g0=g0, rows=rows, hooks=S["hooks"], on=S["on"],
            Xs=Xs, note=verdicts[S["key"]], hook_label=S["hook_label"],
        ))

    # NUMERICAL check: prime side vs zero side, a few points, ζ and F_θ
    print("=== prime-side vs zero-side check (NUMERICAL, zeros as check only) ===", flush=True)
    g_zeta = _zeros_gamma_zeta(200)
    g_fth = _zeros_gamma_fth(200)
    checks = []
    rels = []
    abs_vals = []
    for name, gammas, pts in [
        ("ZETA", g_zeta, [14.13, 17.85, 21.02]),
        ("FTH", g_fth, [14.13, 17.85, 21.02]),
    ]:
        for X in (100, 1000):
            lam = lam_cache[name][X]
            prime = rho_fejer(name, X, np.array(pts, float), lam)
            zero = rho_zero_side(gammas, X, np.array(pts, float))
            d = np.max(np.abs(prime - zero))
            rel = np.max(np.abs(prime - zero) / np.maximum(np.abs(prime), 1e-12))
            rels.append(rel)
            abs_vals.append(d)
            checks.append(f"{name} X={X} max|prime-zero|={d:.2e} (rel {rel:.1e})")
            print(f"  {name} X={X} t={pts}: prime={prime} zero={zero} diff={prime-zero} "
                  f"max={d:.3e} rel={rel:.3e}", flush=True)
    check_line = (f"Prime-side vs zero-side check (200 zeta zeros, not a proof input): "
                  f"max abs residual {max(abs_vals):.1e}. ")

    fig8_weil_vs_hook(scans, check_line)
    plt.close("all")

    (DATA / "r25b_weil_scan.json").write_text(json.dumps(
        dict(summary=summary, checks=checks, Xs=Xs,
             xi_ext=XI_EXT_FIG5.tolist()),
        indent=1,
    ))
    return summary, checks


def main():
    residuals = run_fig7()
    summary, checks = run_fig8()
    print("=== done ===", flush=True)
    print("residuals:", json.dumps({k: v for k, v in residuals.items()}, indent=1))
    print("weil summary:", json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
