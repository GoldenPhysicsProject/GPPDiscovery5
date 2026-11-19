#!/usr/bin/env python3
"""Round 25: visualizations of how zeta zeros line up on the critical line.

Figures are for illustration only. Zeros are used as plot data, never as
proof inputs. RH remains OPEN.

— Daniel Toupin / Golden Physics Project. This seat: Heart of the Cosmos / Grok (GPPDiscovery5).
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from matplotlib.colors import hsv_to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

fp = mp.fp

ROOT = Path("/workspace/r25viz")
OUT = ROOT / "out"
DATA = ROOT / "data"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

# --------------- typography / dark scientific style ---------------
BG = "#0b0e14"
FG = "#e8e4db"
MUTED = "#9aa3ad"
GOLD = "#e0b85a"
CRIMSON = "#e05a3a"
CYAN = "#5eb0d6"
VIOLET = "#b085d4"
GREEN = "#7cba6f"
LINE_IM = "#5eb0d6"
LINE_RE = "#e05a3a"
CRIT = "#f3e6c4"

mpl.rcParams.update({
    "figure.facecolor": BG,
    "axes.facecolor": BG,
    "savefig.facecolor": BG,
    "savefig.edgecolor": BG,
    "text.color": FG,
    "axes.labelcolor": FG,
    "axes.edgecolor": "#3a4250",
    "xtick.color": FG,
    "ytick.color": FG,
    "axes.linewidth": 0.8,
    "font.family": "serif",
    "font.serif": ["Source Serif Pro", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "axes.grid": False,
    "savefig.dpi": 230,
    "figure.dpi": 120,
    "lines.solid_capstyle": "round",
    "axes.titlepad": 10,
})

# --------------- model constants ---------------
THETA_FTH = 0.25
KAP = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
CHI_M4 = [0, 1, 0, -1]
CHI_5 = [0, 1, -1, -1, 1]
CHI_M20 = [CHI_M4[n % 4] * CHI_5[n % 5] for n in range(20)]
CHI5_I = [0, 1, 1j, -1j, -1]
CHI5_IB = [0, 1, -1j, 1j, -1]
DH_PERIOD = [0.0, 1.0, KAP, -KAP, -1.0]
# Epstein conductor factor: sqrt(20)/(2π) = sqrt(5)/π
EPS_A = math.sqrt(20.0) / (2.0 * math.pi)
DH_QPI = 5.0 / math.pi  # (5/π)^{s/2} prefactor base


# =====================================================================
# Evaluators.  A real-positive t-dependent scale is factored out so
# Re=0 and Im=0 loci are unchanged, but float underflow is avoided.
# =====================================================================

def zeta_zeros(n):
    return np.array([float(fp.zetazero(k).imag) for k in range(1, n + 1)])


def xi_zeta(sig, t):
    """ξ(σ+it) * exp(π|t|/4), entire up to removable 0,1 which we avoid."""
    s = complex(sig, t)
    pref = 0.5 * s * (s - 1)
    return pref * fp.exp(
        fp.loggamma(s / 2) - (s / 2) * math.log(math.pi) + 0.25 * math.pi * abs(t)
    ) * fp.zeta(s)


def zq(sig, t):
    s = complex(sig, t)
    return (
        fp.zeta(s) * fp.dirichlet(s, CHI_M20)
        + fp.dirichlet(s, CHI_M4) * fp.dirichlet(s, CHI_5)
    )


def xi_eps(sig, t):
    """Entire cousin: s(s-1) Λ_Q(s) * exp(π|t|/2)  (Γ(s) ~ e^{-π|t|/2})."""
    s = complex(sig, t)
    pref = s * (s - 1) * (EPS_A ** s)
    return pref * fp.exp(fp.loggamma(s) + 0.5 * math.pi * abs(t)) * zq(sig, t)


def f_dh(sig, t):
    s = complex(sig, t)
    k = KAP
    return 0.5 * (1 - 1j * k) * fp.dirichlet(s, CHI5_I) + 0.5 * (1 + 1j * k) * fp.dirichlet(s, CHI5_IB)


def xi_dh(sig, t):
    """Λ_DH(s) * exp(π|t|/4), FE-symmetric; real on the critical line."""
    s = complex(sig, t)
    # (5/π)^{s/2} Γ((s+1)/2) f(s)
    return fp.exp(
        (s / 2) * math.log(DH_QPI) + fp.loggamma((s + 1) / 2) + 0.25 * math.pi * abs(t)
    ) * f_dh(sig, t)


def xi_fth(sig, t, theta=THETA_FTH):
    """ξ(s+θ) ξ(s-θ) with the same real-positive scale in each factor."""
    return xi_zeta(sig + theta, t) * xi_zeta(sig - theta, t)


def hardy_Z(t):
    return float(fp.siegelz(t))


def siegel_theta(t):
    return float(fp.siegeltheta(t))


def gram_point(n, tol=1e-12):
    """Solve ϑ(g_n) = nπ."""
    # Newton on the Stirling approximation, then refine with siegeltheta.
    target = n * math.pi
    # initial: invert  t/2 log(t/2π) - t/2 - π/8 = nπ
    t = 2 * math.pi * max(n + 1.0, 1.0)
    for _ in range(20):
        if t <= 2 * math.pi:
            t = 2 * math.pi + 1.0
        th = 0.5 * t * math.log(t / (2 * math.pi)) - 0.5 * t - math.pi / 8
        dth = 0.5 * math.log(t / (2 * math.pi))
        t -= (th - target) / dth
    # refine
    for _ in range(12):
        th = float(fp.siegeltheta(t))
        dth = 0.5 * math.log(t / (2 * math.pi))  # ϑ'(t) ~ this
        step = (th - target) / dth
        t -= step
        if abs(step) < tol:
            break
    return t


# =====================================================================
# Grid helpers
# =====================================================================

def eval_grid(fn, sigs, ts, path, label=""):
    """Evaluate fn(sig,t) on a mesh. Cache as complex64/128 npz."""
    path = Path(path)
    if path.exists():
        z = np.load(path)
        if z.shape == (len(ts), len(sigs)):
            print(f"  loaded {label or path.name} {z.shape}")
            return z
    Z = np.empty((len(ts), len(sigs)), dtype=np.complex128)
    for i, t in enumerate(ts):
        if i % 40 == 0:
            print(f"  {label}: row {i}/{len(ts)}  t={t:.3f}", flush=True)
        for j, sig in enumerate(sigs):
            try:
                Z[i, j] = complex(fn(float(sig), float(t)))
            except Exception:
                Z[i, j] = np.nan + 1j * np.nan
    np.save(path, Z)
    print(f"  saved {path} {Z.shape}")
    return Z


def quadrant_rgba(W, alpha=0.55):
    """Muted quadrant colours of arg(W): the four sign patterns of (Re, Im)."""
    re, im = np.real(W), np.imag(W)
    # 0: Re>0,Im>0 gold-olive; 1: Re<0,Im>0 crimson-brown;
    # 2: Re<0,Im<0 violet; 3: Re>0,Im<0 cyan-teal
    pal = np.array([
        [0.78, 0.62, 0.28, alpha],
        [0.75, 0.32, 0.28, alpha],
        [0.45, 0.32, 0.62, alpha],
        [0.28, 0.52, 0.62, alpha],
    ])
    idx = (re < 0).astype(int) + 2 * (im < 0).astype(int)
    # When Re>0, Im>0: 0; Re<0,Im>0: 1; Re>0,Im<0: 2? Wait
    # (re<0) + 2*(im<0):
    # ++ : 0 gold
    # -+ : 1 crimson
    # +- : 2 teal?  0+2=2  but I wanted 3 for Re>0 Im<0
    # -- : 3
    # Let me just map and swap palette so 2 is violet and 3 is teal...
    # Actually ++0 gold, -+1 crimson, +-2=2, --3=3.
    # I said 2 should be violet (Re<0,Im<0) but that's index 3.
    pal = np.array([
        [0.78, 0.62, 0.28, alpha],  # Re>0 Im>0
        [0.75, 0.32, 0.28, alpha],  # Re<0 Im>0
        [0.28, 0.52, 0.62, alpha],  # Re>0 Im<0
        [0.45, 0.32, 0.62, alpha],  # Re<0 Im<0
    ])
    return pal[idx]


def wegert_rgb(W, mag_floor=1e-20):
    """Wegert domain colouring: hue = arg, sawtooth lightness from log|W|."""
    ang = np.angle(W)
    H = (ang + math.pi) / (2 * math.pi)
    r = np.abs(W)
    r = np.where(np.isfinite(r), r, mag_floor)
    logm = np.log(np.maximum(r, mag_floor))
    # two family of rings: coarse (mod log 2) and a soft contrast
    frac = np.mod(logm, math.log(2.0)) / math.log(2.0)
    V = 0.22 + 0.78 * frac
    S = 0.85 * np.ones_like(H)
    # nan -> black
    bad = ~np.isfinite(W)
    hsv = np.stack([H, S, V], axis=-1)
    rgb = hsv_to_rgb(hsv)
    rgb[bad] = 0.0
    return rgb


def style_ax(ax, xlabel, ylabel, title=None):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, color=FG, loc="left", pad=8)
    ax.tick_params(length=3, width=0.6, pad=3)
    for spine in ax.spines.values():
        spine.set_color("#3a4250")


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=230, bbox_inches="tight", pad_inches=0.18)
    print("wrote", path, "size", path.stat().st_size)
    return path


# =====================================================================
# Figure 1 — eta arrow walks
# =====================================================================

def eta_terms(s, N):
    n = np.arange(1, N + 1, dtype=np.float64)
    signs = np.where(n % 2 == 1, 1.0, -1.0)
    # n^{-s} = exp(-s log n)
    logs = np.log(n)
    return signs * np.exp(-s * logs)


def paired_walk(s, Npairs):
    """Walk of odd-even pairs of η(s): b_k = (2k-1)^{-s} - (2k)^{-s}."""
    k = np.arange(1, Npairs + 1, dtype=np.float64)
    odd = np.exp(-s * np.log(2 * k - 1))
    even = np.exp(-s * np.log(2 * k))
    terms = odd - even
    walk = np.cumsum(terms)
    return terms, walk


def true_eta(s):
    s = complex(s)
    return (1 - 2 ** (1 - s)) * complex(fp.zeta(s))



def fig1_arrow_walk(z1, z2):
    points = [
        ("At the first zero\n" + r"$s=\frac{1}{2}+14.13i$", complex(0.5, z1), True),
        ("Same height, off the line\n" + r"$s=0.80+14.13i$", complex(0.80, z1), False),
        ("On the line, between zeros\n" + r"$s=\frac{1}{2}+17.85i$  (Gram $g_0$)", complex(0.5, gram_point(0)), False),
        ("At the second zero\n" + r"$s=\frac{1}{2}+21.02i$", complex(0.5, z2), True),
    ]
    N = 80
    N_arrow = 10
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 10.6))
    fig.suptitle(
        r"How positive coefficients cancel: the rotating-arrow walk of $\eta(s)=\sum (-1)^{n-1}n^{-s}$",
        color=FG, fontsize=14, fontweight="medium", y=0.98, x=0.03, ha="left",
    )
    fig.text(
        0.03, 0.938,
        r"Each arrow is a term $(\pm)n^{-s}$. The size is $n^{-\sigma}$; the direction is the phase $(-1)^{n-1}n^{-it}$. "
        r"No negative integers appear. Cancellation is geometric: the arrows turn with $t\log n$ and can close on the origin — or miss it.",
        color=MUTED, fontsize=9, ha="left", va="top",
    )

    cmap = plt.get_cmap("plasma")
    for ax, (title, s, is_zero) in zip(axes.ravel(), points):
        terms = eta_terms(s, N)
        walk = np.cumsum(terms)
        end = true_eta(s)
        xs = np.r_[0.0, walk.real]
        ys = np.r_[0.0, walk.imag]
        # faint tail
        ax.plot(xs, ys, color=cmap(0.75), lw=0.7, alpha=0.45, zorder=2)
        # bold early path
        ax.plot(xs[:N_arrow + 1], ys[:N_arrow + 1], color=cmap(0.2), lw=1.7, zorder=3)
        # arrows on early segments
        for k in range(N_arrow):
            x0, y0 = xs[k], ys[k]
            x1, y1 = xs[k + 1], ys[k + 1]
            ax.annotate(
                "", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=cmap(0.12 + 0.6 * k / N_arrow),
                                lw=1.15, mutation_scale=11),
                zorder=4,
            )
            mx, my = 0.55 * x0 + 0.45 * x1, 0.55 * y0 + 0.45 * y1
            ax.annotate(str(k + 1), (mx, my), textcoords="offset points",
                        xytext=(4, 4), fontsize=7, color=MUTED, zorder=5)
        ax.scatter([0], [0], s=22, color=FG, zorder=6, linewidths=0)
        ax.scatter([end.real], [end.imag], marker="*", s=130, color=GOLD, zorder=7,
                   edgecolors="white", linewidths=0.5)
        ax.axhline(0, color="#3a4250", lw=0.5, zorder=0)
        ax.axvline(0, color="#3a4250", lw=0.5, zorder=0)
        # limits from early walk + target, with padding
        pad_pts_x = np.r_[xs[:N_arrow + 1], end.real, 0.0]
        pad_pts_y = np.r_[ys[:N_arrow + 1], end.imag, 0.0]
        cx = 0.5 * (pad_pts_x.min() + pad_pts_x.max())
        cy = 0.5 * (pad_pts_y.min() + pad_pts_y.max())
        span = max(pad_pts_x.max() - pad_pts_x.min(), pad_pts_y.max() - pad_pts_y.min(), 0.6) * 0.68
        ax.set_xlim(cx - span, cx + span)
        ax.set_ylim(cy - span, cy + span)
        ax.set_aspect("equal")
        style_ax(ax, r"$\mathrm{Re}\,\eta_N(s)$", r"$\mathrm{Im}\,\eta_N(s)$", title)
        if is_zero:
            note = r"zero of $\zeta$ (so $\eta=0$ too)" + f"\nwalk after {N} terms: ${walk[-1].real:+.3f}{walk[-1].imag:+.3f}i$"
            color = GOLD
        else:
            note = (r"true $\eta(s)=" + f"{end.real:+.3f}{end.imag:+.3f}i$ (star)"
                    + f"\nwalk after {N} terms still approaching")
            color = MUTED
        ax.text(0.03, 0.03, note, transform=ax.transAxes, fontsize=8,
                color=color, va="bottom", ha="left")

    handles = [
        Line2D([0], [0], color=cmap(0.2), lw=2, label="first 10 terms (arrows)"),
        Line2D([0], [0], color=cmap(0.75), lw=1.2, label=f"terms 11–{N}"),
        Line2D([0], [0], marker="*", color="none", markerfacecolor=GOLD,
               markeredgecolor="white", markersize=11, label=r"true $\eta(s)$"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=FG, markersize=6, label="start at $0$"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.995),
               fontsize=8, labelcolor=FG)
    fig.text(
        0.03, 0.012,
        r"NUMERICAL, illustration. $\eta$ converges for $\mathrm{Re}\,s>0$. $\zeta=\eta/(1-2^{1-s})$; the denominator does not vanish at these zeros. "
        r"The walk is the ordinary partial sums (no acceleration).",
        color=MUTED, fontsize=8,
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.07, top=0.86, wspace=0.22, hspace=0.32)
    return save(fig, "fig1_arrow_walk.png")


def find_extrema(ts, Zs):
    """Approximate stationary points of a sampled real wave, then refine."""
    d = np.diff(Zs)
    out = []
    for i in range(len(d) - 1):
        if d[i] == 0 or d[i] * d[i + 1] < 0:
            # quadratic refine around ts[i+1]
            t0 = ts[i + 1]
            def dZ(t, h=1e-5):
                return (fp.siegelz(t + h) - fp.siegelz(t - h)) / (2 * h)
            try:
                tstar = float(fp.findroot(dZ, t0))
                if ts[0] < tstar < ts[-1] and all(abs(tstar - u) > 1e-4 for u in out):
                    out.append(tstar)
            except Exception:
                out.append(float(t0))
    return np.array(out)




def extrema_of_xi(tmin, tmax, n=4000):
    """Stationary points of the real function xi(1/2+it) (not of Hardy's Z).
    Extra Im xi=0 branches meet the critical line at these points (Cauchy-Riemann).
    They differ from Z-extrema because xi(1/2+it) = A(t) Z(t) with a decaying
    archimedean amplitude A(t)."""
    ts = np.linspace(tmin, tmax, n)
    def u(t):
        t = float(t)
        s = complex(0.5, t)
        pref = 0.5 * s * (s - 1)
        # true xi, no extra t-dependent scale
        return complex(pref * fp.exp(fp.loggamma(s / 2) - (s / 2) * math.log(math.pi)) * fp.zeta(s)).real
    us = np.array([u(t) for t in ts])
    d = np.diff(us)
    out = []
    for i in range(len(d) - 1):
        if d[i] * d[i + 1] < 0:
            t0 = float(ts[i + 1])
            def du(t, h=1e-5):
                return (u(t + h) - u(t - h)) / (2 * h)
            try:
                tstar = float(fp.findroot(du, t0))
                if tmin < tstar < tmax and all(abs(tstar - x) > 1e-3 for x in out):
                    out.append(tstar)
            except Exception:
                out.append(t0)
    return np.array(out)


def fig3_hardy(z_im, grams, extrema, tmax=42):
    ts = np.linspace(0.05, tmax, 3500)
    Zs = np.array([hardy_Z(t) for t in ts])
    ths = np.array([siegel_theta(t) for t in ts])

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(12.0, 7.8), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.15], "hspace": 0.06},
    )
    fig.suptitle(
        r"On the critical line the wave is real: Hardy’s $Z(t)=e^{i\vartheta(t)}\zeta(\frac{1}{2}+it)$",
        color=FG, fontsize=14, x=0.03, ha="left", y=0.98,
    )
    fig.text(
        0.03, 0.935,
        r"The factor $e^{i\vartheta(t)}$ cancels the oscillating phase of $\zeta$ on $\mathrm{Re}\,s=\frac{1}{2}$, so $Z(t)$ is real. "
        r"Zeros of $\zeta$ on the line are ordinary sign changes of this real function. Gram points satisfy $\vartheta(g_n)=n\pi$.",
        color=MUTED, fontsize=9, ha="left", va="top",
    )

    ax.plot(ts, Zs, color=CYAN, lw=1.15)
    ax.axhline(0, color=CRIT, lw=0.6, alpha=0.7)
    ax.fill_between(ts, Zs, 0, where=(Zs < 0), color=CRIMSON, alpha=0.18, linewidth=0)
    ax.fill_between(ts, Zs, 0, where=(Zs >= 0), color=CYAN, alpha=0.10, linewidth=0)
    ax.scatter(z_im, np.zeros_like(z_im), s=28, color=GOLD, zorder=5,
               edgecolors=BG, linewidths=0.4, label="zeros of $Z$ (hence of $\\zeta$ on the line)")
    ymax = max(2.6, float(np.nanmax(np.abs(Zs))) * 1.22)
    ax.set_ylim(-ymax, ymax)
    for n, g in enumerate(grams):
        if 0 < g < tmax:
            ax.axvline(g, color=GOLD, lw=0.55, alpha=0.5, ls="--")
            ax2.axvline(g, color=GOLD, lw=0.55, alpha=0.5, ls="--")
            ax.text(g, ymax * 0.97, f"$g_{{{n}}}$", color=GOLD, fontsize=8,
                    ha="center", va="top")
    Zext = np.array([hardy_Z(t) for t in extrema])
    ax.scatter(extrema, Zext, s=16, color=VIOLET, zorder=6, marker="D",
               label="extrema of $Z$ (where extra $\\mathrm{Im}\\,\\xi=0$ branches leave the line)")
    style_ax(ax, "", r"$Z(t)$")
    ax.legend(loc="lower left", fontsize=8, labelcolor=FG, ncol=2)
    ax.set_xlim(0, tmax)

    ax2.plot(ts, ths / math.pi, color=FG, lw=1.15)
    thmin, thmax = float(np.min(ths / math.pi)), float(np.max(ths / math.pi))
    ax2.set_ylim(thmin - 0.45, thmax + 0.55)
    for n in range(int(math.floor(thmin)) - 1, int(math.ceil(thmax)) + 2):
        ax2.axhline(n, color="#3a4250", lw=0.5)
    thz = [siegel_theta(t) / math.pi for t in z_im if t < tmax]
    tz = [t for t in z_im if t < tmax]
    ax2.scatter(tz, thz, s=22, color=GOLD, zorder=5, edgecolors=BG, linewidths=0.3)
    style_ax(ax2, r"$t$", r"$\vartheta(t)/\pi$")
    ax2.set_xlim(0, tmax)
    ax2.text(
        0.015, 0.95,
        "Gram’s law (heuristic): usually one zero in each interval $(g_n,g_{n+1})$.\n"
        r"It already fails at low height (e.g. between $g_{126}$ and $g_{127}$); it is not a proof.",
        transform=ax2.transAxes, ha="left", va="top", fontsize=7.5, color=MUTED,
    )

    fig.text(
        0.03, 0.01,
        r"NUMERICAL, illustration. $Z(t)$ and $\vartheta(t)$ from the Riemann–Siegel formulae (mpmath). Gram’s law is classical (Gram 1903); exceptions are well known (Hutchinson, Titchmarsh).",
        color=MUTED, fontsize=8,
    )
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.09, top=0.86, hspace=0.08)
    return save(fig, "fig3_hardy_Z_gram.png")


def fig2_domain(sigs, ts, Wz, z_im, offline_eps, sigs_e, ts_e, We, online_eps=None):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 8.6), sharey=False)
    fig.suptitle(
        r"A portrait of $\xi(s)$ in the strip: every colour-wheel is a zero",
        color=FG, fontsize=14, x=0.03, ha="left", y=0.98,
    )
    fig.text(
        0.03, 0.935,
        r"Wegert colouring: hue = argument, brightness repeats with $\log|\xi|$. "
        r"The functional equation plus real coefficients makes the portrait conjugate-symmetric across $\mathrm{Re}\,s=\frac{1}{2}$.",
        color=MUTED, fontsize=9, ha="left", va="top",
    )

    def panel(ax, sigs, ts, W, title, zeros_on, zeros_off, tlim):
        rgb = wegert_rgb(W)
        ax.imshow(
            rgb, origin="lower", aspect="auto",
            extent=[sigs[0], sigs[-1], ts[0], ts[-1]], interpolation="bilinear",
        )
        ax.axvline(0.5, color=CRIT, lw=0.9, alpha=0.9, ls="-")
        if zeros_on is not None:
            on = [t for t in zeros_on if ts[0] <= t <= ts[-1]]
            ax.scatter([0.5] * len(on), on, s=18, facecolors="none",
                       edgecolors="white", linewidths=0.7, zorder=5)
        if zeros_off:
            for sig, t in zeros_off:
                if ts[0] <= t <= ts[-1] and sigs[0] <= sig <= sigs[-1]:
                    ax.scatter([sig], [t], s=28, marker="o", facecolors="none",
                               edgecolors=CRIMSON, linewidths=1.0, zorder=6)
                    ax.scatter([1 - sig], [t], s=28, marker="o", facecolors="none",
                               edgecolors=CRIMSON, linewidths=1.0, zorder=6)
        style_ax(ax, r"$\mathrm{Re}\,s$", r"$\mathrm{Im}\,s$", title)
        ax.set_xlim(sigs[0], sigs[-1])
        ax.set_ylim(tlim[0], tlim[1])

    panel(ax1, sigs, ts, Wz, r"Riemann $\xi(s)$  (zeros on the line, in this window)",
          z_im, None, (ts[0], ts[-1]))
    panel(ax2, sigs_e, ts_e, We,
          r"Epstein $\xi_Q(s)$ for $x^2+5y^2$  (off-line zeros in red)",
          online_eps, offline_eps, (ts_e[0], ts_e[-1]))

    fig.text(
        0.03, 0.015,
        r"NUMERICAL, illustration. Left: $\xi(s)=\frac{1}{2} s(s-1)\pi^{-s/2}\Gamma(s/2)\zeta(s)$. "
        r"Right: $s(s-1)(\sqrt{20}/2\pi)^s\Gamma(s)Z_Q(s)$. Both are real on the critical line. RH is OPEN.",
        color=MUTED, fontsize=8,
    )
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.08, top=0.86)
    return save(fig, "fig2_domain_coloring.png")


# =====================================================================
# Level curves
# =====================================================================

def draw_level_panel(ax, sigs, ts, W, on_line, off_line, extrema, title,
                     tlim=None, show_extrema=True):
    rgba = quadrant_rgba(W, alpha=0.52)
    # W is (t, σ); imshow wants (y, x) = (σ, t)
    ax.imshow(
        np.transpose(rgba, (1, 0, 2)), origin="lower", aspect="auto",
        extent=[ts[0], ts[-1], sigs[0], sigs[-1]], interpolation="nearest",
        zorder=0,
    )
    Re = np.real(W).T  # (nsig, nt)
    Im = np.imag(W).T
    TT, SS = np.meshgrid(ts, sigs)
    ax.contour(TT, SS, Re, levels=[0], colors=[LINE_RE], linewidths=1.15, zorder=3)
    # Im ξ is identically 0 on the whole critical line; mask a thin strip so
    # contour finds the *other* branches cleanly. The line itself is drawn below.
    Im_m = np.ma.array(Im, mask=np.abs(SS - 0.5) < 0.02)
    ax.contour(TT, SS, Im_m, levels=[0], colors=[LINE_IM], linewidths=1.05, zorder=3)
    # critical line (itself an Im=0 component)
    ax.axhline(0.5, color=CRIT, lw=1.35, zorder=4, alpha=0.95)
    if on_line is not None:
        on = [t for t in on_line if ts[0] <= t <= ts[-1]]
        ax.scatter(on, [0.5] * len(on), s=26, color=GOLD, zorder=6,
                   edgecolors=BG, linewidths=0.4)
    if off_line:
        for sig, t in off_line:
            if ts[0] <= t <= ts[-1]:
                ax.scatter([t, t], [sig, 1 - sig], s=32, marker="o",
                           facecolors=CRIMSON, edgecolors="white", linewidths=0.4, zorder=7)
    if show_extrema and extrema is not None:
        ex = [t for t in extrema if ts[0] <= t <= ts[-1]]
        ax.scatter(ex, [0.5] * len(ex), s=16, marker="D", color=VIOLET, zorder=6,
                   edgecolors=BG, linewidths=0.3)
    style_ax(ax, r"$t=\mathrm{Im}\,s$", r"$\sigma=\mathrm{Re}\,s$", title)
    if tlim:
        ax.set_xlim(*tlim)
    else:
        ax.set_xlim(ts[0], ts[-1])
    ax.set_ylim(sigs[0], sigs[-1])
    ax.set_aspect("auto")


def fig4_level_zeta(sigs, ts, W, z_im, extrema):
    fig = plt.figure(figsize=(12.2, 8.3))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.15], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    axz = fig.add_subplot(gs[1], sharex=ax)

    fig.suptitle(
        r"How the zeros line up: $\mathrm{Im}\,\xi=0$ is the whole critical line; zeros are where $\mathrm{Re}\,\xi$ vanishes too",
        color=FG, fontsize=13.5, x=0.03, ha="left", y=0.98,
    )
    fig.text(
        0.03, 0.935,
        r"Blue curves: $\mathrm{Im}\,\xi=0$. Red curves: $\mathrm{Re}\,\xi=0$. The pale line $\sigma=\frac{1}{2}$ is entirely blue — "
        r"that is the reflection symmetry, not a miracle at each zero. Gold dots are zeros. Violet diamonds are extrema of $\xi(\frac{1}{2}+it)$ — that is where an extra blue branch meets the line (Cauchy–Riemann). They are close to, but not the same as, the peaks of Hardy’s $Z(t)$.",
        color=MUTED, fontsize=9, ha="left", va="top",
    )

    draw_level_panel(ax, sigs, ts, W, z_im, None, extrema,
                     None, tlim=(ts[0], ts[-1]), show_extrema=True)
    ax.set_xlabel("")
    # legend
    handles = [
        Line2D([0], [0], color=LINE_IM, lw=2, label=r"$\mathrm{Im}\,\xi=0$"),
        Line2D([0], [0], color=LINE_RE, lw=2, label=r"$\mathrm{Re}\,\xi=0$"),
        Line2D([0], [0], color=CRIT, lw=2.2, label=r"critical line (itself $\mathrm{Im}\,\xi=0$)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GOLD,
               markersize=7, label="on-line zero"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor=VIOLET,
               markersize=6, label=r"extremum of $\xi(\frac{1}{2}+it)$"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=8, labelcolor=FG,
              ncol=2, bbox_to_anchor=(1.0, 1.02))

    # Z(t) aligned
    tplot = np.linspace(ts[0], ts[-1], 2500)
    Zs = np.array([hardy_Z(t) for t in tplot])
    axz.plot(tplot, Zs, color=CYAN, lw=1.1)
    axz.axhline(0, color=CRIT, lw=0.6, alpha=0.7)
    axz.fill_between(tplot, Zs, 0, where=(Zs < 0), color=CRIMSON, alpha=0.16, lw=0)
    axz.fill_between(tplot, Zs, 0, where=(Zs >= 0), color=CYAN, alpha=0.10, lw=0)
    on = [t for t in z_im if ts[0] <= t <= ts[-1]]
    axz.scatter(on, np.zeros(len(on)), s=22, color=GOLD, zorder=5, edgecolors=BG, lw=0.3)
    ex = [t for t in extrema if ts[0] <= t <= ts[-1]]
    for t0 in ex:
        axz.axvline(t0, color=VIOLET, lw=0.7, alpha=0.7, ls=":")
    axz.scatter(ex, [hardy_Z(t) for t in ex], s=14, marker="D", color=VIOLET, zorder=6,
                edgecolors=BG, lw=0.3)
    style_ax(axz, r"$t$", r"$Z(t)$")
    axz.set_xlim(ts[0], ts[-1])

    fig.text(
        0.03, 0.012,
        r"NUMERICAL, illustration. Background colours are the four sign-quadrants of $\xi(\sigma+it)$. Four colours meet at a zero. "
        r"The right-angle crossing of red through the line at a simple zero, and the extra blue branch at an extremum of $\xi$ on the line, are Cauchy–Riemann (classical). Extrema of $\xi(\frac{1}{2}+it)$ differ from those of $Z(t)$ because of the decaying Gamma factor.",
        color=MUTED, fontsize=8,
    )
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.08, top=0.86)
    return save(fig, "fig4_level_curves_zeta.png")


def fig5_fakes(bundles):
    """bundles: list of dicts with keys name, sigs, ts, W, on, off, extrema, subtitle."""
    n = len(bundles)
    fig, axes = plt.subplots(n, 1, figsize=(12.2, 3.45 * n + 1.3), sharex=False)
    if n == 1:
        axes = [axes]
    fig.suptitle(
        r"The same lining-up mechanism, side by side: $\zeta$ stays on the line; the fakes do not",
        color=FG, fontsize=13.5, x=0.03, ha="left", y=0.99,
    )
    fig.text(
        0.03, 0.955,
        r"Same colour language as Fig. 4. An off-line zero (crimson) is just a red/blue crossing that happens away from $\sigma=\frac{1}{2}$, "
        r"on an extra $\mathrm{Im}=0$ branch. Every model here is real on the critical line — that alone does not force the zeros onto it.",
        color=MUTED, fontsize=9, ha="left", va="top",
    )
    for ax, B in zip(axes, bundles):
        draw_level_panel(
            ax, B["sigs"], B["ts"], B["W"], B["on"], B["off"], B.get("extrema"),
            B["name"], tlim=(B["ts"][0], B["ts"][-1]), show_extrema=B.get("extrema") is not None,
        )
        if B.get("note"):
            ax.text(0.99, 0.06, B["note"], transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8, color=FG, bbox=dict(boxstyle="round,pad=0.25", fc="#161b24", ec="#3a4250", alpha=0.85))
    handles = [
        Line2D([0], [0], color=LINE_IM, lw=2, label=r"$\mathrm{Im}\,\xi=0$"),
        Line2D([0], [0], color=LINE_RE, lw=2, label=r"$\mathrm{Re}\,\xi=0$"),
        Line2D([0], [0], color=CRIT, lw=2.0, label=r"critical line"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GOLD, markersize=7, label="on-line zero"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CRIMSON, markersize=7, label="off-line zero"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.99, 0.99),
               fontsize=8, labelcolor=FG, ncol=3)
    fig.text(
        0.03, 0.008,
        r"NUMERICAL, illustration. $F_\theta=\xi(s+\theta)\xi(s-\theta)$ with $\theta=\frac{1}{4}$ (zeros of $\zeta$ shifted off the line). "
        r"Epstein $Z_Q$ for $x^2+5y^2$. Davenport–Heilbronn $f$ (period-5 coefficients). RH is OPEN.",
        color=MUTED, fontsize=8,
    )
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.055, top=0.90, hspace=0.42)
    return save(fig, "fig5_level_curves_fakes.png")


# =====================================================================
# Figure 6 — cosine cancellation / Riemann–Siegel snapshot
# =====================================================================

def rs_terms(t):
    """Riemann–Siegel main sum: 2 n^{-1/2} cos(ϑ(t)-t log n) for n ≤ floor(sqrt(t/2π))."""
    N = int(math.floor(math.sqrt(t / (2 * math.pi))))
    th = siegel_theta(t)
    ns = np.arange(1, max(N, 0) + 1)
    if len(ns) == 0:
        return ns, np.array([]), 0.0, N
    terms = 2 * ns ** (-0.5) * np.cos(th - t * np.log(ns))
    return ns, terms, float(terms.sum()), N



def fig6_cancellation(z_im, grams, tmax=50):
    fig = plt.figure(figsize=(12.2, 8.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 1.05], hspace=0.38, wspace=0.32)
    axwave = fig.add_subplot(gs[0, :])
    axs = [fig.add_subplot(gs[1, i]) for i in range(4)]

    fig.suptitle(
        r"The minus signs are phases, not negative integers: Riemann–Siegel’s cosine sum for $Z(t)$",
        color=FG, fontsize=13.5, x=0.03, ha="left", y=0.98,
    )
    fig.text(
        0.03, 0.932,
        r"$Z(t)=2\sum n^{-1/2}\cos(\vartheta(t)-t\log n)+R(t)$, the sum running over $n\leq\sqrt{t/2\pi}$. "
        r"Every prefactor $n^{-1/2}$ is positive. The cosine is the real part of the rotating arrow $n^{-it}$. "
        r"A zero is where those signed terms, plus a remainder $R(t)$, add to $0$.",
        color=MUTED, fontsize=9, ha="left", va="top",
    )

    ts = np.linspace(0.2, 42, 2800)
    Zs = np.array([hardy_Z(t) for t in ts])
    axwave.plot(ts, Zs, color=CYAN, lw=1.2, label=r"true $Z(t)$")
    axwave.axhline(0, color=CRIT, lw=0.6, alpha=0.7)
    on = [t for t in z_im if t <= 42]
    axwave.scatter(on, np.zeros(len(on)), s=18, color=GOLD, zorder=5, edgecolors=BG, lw=0.3)
    style_ax(axwave, r"$t$", r"$Z(t)$")
    axwave.set_xlim(0, 42)
    axwave.legend(loc="lower left", fontsize=8, labelcolor=FG)

    bars = [
        (r"First zero" + "\n" + r"$t=14.13$, only $n=1$", float(__import__('mpmath').fp.zetazero(1).imag)),
        (r"Sixth zero" + "\n" + r"$t=37.59$, two terms", float(__import__('mpmath').fp.zetazero(6).imag)),
        (r"A later zero" + "\n" + r"$t=143.11$, four terms", float(__import__('mpmath').fp.zetazero(50).imag)),
        (r"Gram peak $g_0$" + "\n" + r"$t\approx 17.85$, not a zero", gram_point(0)),
    ]
    for ax, (lab, t0) in zip(axs, bars):
        ns, terms, sm, N = rs_terms(t0)
        Ztrue = hardy_Z(t0)
        colors = [CYAN if v >= 0 else CRIMSON for v in terms] if len(terms) else []
        if len(ns):
            ax.bar(ns, terms, color=colors, width=0.6, edgecolor=BG, linewidth=0.3, zorder=3)
        ax.axhline(0, color=CRIT, lw=0.6)
        ax.axhline(Ztrue, color=GOLD, lw=0.9, ls="--")
        # put the numeric note in the corner with the larger empty space
        if Ztrue >= 0:
            va, ypos = "bottom", 0.04
        else:
            va, ypos = "top", 0.96
        ax.text(
            0.50, ypos,
            fr"main sum $={sm:.3f}$" + "\n" + fr"$R(t)={Ztrue-sm:.3f}$" + "\n" + fr"$Z(t)={Ztrue:.3f}$",
            transform=ax.transAxes, ha="center", va=va, fontsize=7.4, color=MUTED,
        )
        style_ax(ax, r"$n$", r"term", lab)
        ax.set_xticks(list(ns) if len(ns) else [1])
        if len(terms):
            m = max(0.45, 1.45 * float(np.max(np.abs(terms))))
            # extra room for the annotation
            ax.set_ylim(-m, m)
        ax.set_xlim(0.4, max(N + 0.8, 2.6))

    fig.text(
        0.03, 0.012,
        r"NUMERICAL, illustration. The remainder $R(t)$ is classical (Riemann–Siegel, Gabcke, Arias-de-Reyna). "
        r"At low height $N=\lfloor\sqrt{t/2\pi}\rfloor$ is tiny, so $R(t)$ is comparable to the main sum — the first zero is not simply $\cos\vartheta(t)=0$. RH is OPEN.",
        color=MUTED, fontsize=8,
    )
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.08, top=0.86)
    return save(fig, "fig6_cosine_cancellation.png")


def main():
    print("=== zeros and Gram points ===", flush=True)
    z_im = zeta_zeros(20)
    print("zeta zeros t:", z_im)
    grams = np.array([gram_point(n) for n in range(0, 15)])
    print("Gram g0..", grams)
    # extrema of Z on (0, 42)
    ts_ext = np.linspace(1.0, 42.0, 2000)
    Zs_ext = np.array([hardy_Z(t) for t in ts_ext])
    extrema = find_extrema(ts_ext, Zs_ext)
    print("extrema:", extrema)

    # stored fake zeros
    eps_off = [
        (0.9329696974854141, 15.668249531278473),
        (0.9376669067003994, 29.983395235155516),
        (0.6969270453122914, 36.37406368638513),
        (0.8231873164780866, 44.000113180236895),
    ]
    # a few EPS on-line zeros in our windows (from epstein_zeros.json)
    eps_on = [
        1.4562032347251894, 5.326203936204837, 7.11349023725727,
        9.208567035966308, 10.571063256663962, 12.494830631904136,
        13.098509365914271, 17.596050852287988, 18.90136055238456,
        19.80901394515022,
    ]
    dh_off = [
        (0.8085171824566374, 85.69934848537759),
    ]
    dh_on = [
        80.32794708609034, 81.92540239183100, 83.10884108474982,
        87.64747633255974,
    ]
    # F_theta off-line: 1/2 ± θ + i γ  for first few zeta zeros
    fth_off = []
    for t in z_im[:6]:
        fth_off.append((0.5 + THETA_FTH, float(t)))
        fth_off.append((0.5 - THETA_FTH, float(t)))

    meta = dict(
        z_im=z_im.tolist(), grams=grams.tolist(), extrema=extrema.tolist(),
        eps_off=eps_off, eps_on=eps_on, dh_off=dh_off, dh_on=dh_on,
        fth_off=fth_off, theta=THETA_FTH,
    )
    (DATA / "meta.json").write_text(json.dumps(meta, indent=1))

    print("=== fig1 arrow walk ===", flush=True)
    fig1_arrow_walk(z_im[0], z_im[1])
    plt.close("all")

    print("=== extrema of xi(1/2+it) ===", flush=True)
    xi_ext = extrema_of_xi(1.0, 42.0)
    print("xi extrema:", xi_ext)

    print("=== fig3 Hardy Z ===", flush=True)
    fig3_hardy(z_im, grams, extrema, tmax=42)
    plt.close("all")

    print("=== fig6 cosine cancellation ===", flush=True)
    fig6_cancellation(z_im, grams)
    plt.close("all")

    # ---- grids ----
    print("=== grid zeta ξ  (fig2, fig4) ===", flush=True)
    sigs_z = np.linspace(-0.10, 1.10, 241)
    ts_z = np.linspace(0.40, 40.0, 661)
    Wz = eval_grid(xi_zeta, sigs_z, ts_z, DATA / "grid_xi_zeta.npy", "xi_zeta")

    print("=== grid Epstein (fig2, fig5) ===", flush=True)
    sigs_e = np.linspace(-0.05, 1.05, 221)
    ts_e = np.linspace(10.0, 34.0, 401)
    We = eval_grid(xi_eps, sigs_e, ts_e, DATA / "grid_xi_eps.npy", "xi_eps")

    print("=== grid DH around first off-line zero ===", flush=True)
    sigs_d = np.linspace(-0.05, 1.05, 221)
    ts_d = np.linspace(81.5, 90.0, 321)
    Wd = eval_grid(xi_dh, sigs_d, ts_d, DATA / "grid_xi_dh.npy", "xi_dh")

    print("=== grid F_theta around first zeros ===", flush=True)
    sigs_f = np.linspace(-0.05, 1.05, 221)
    ts_f = np.linspace(10.0, 26.0, 321)
    Wf = eval_grid(xi_fth, sigs_f, ts_f, DATA / "grid_xi_fth.npy", "xi_fth")

    # sanity: Im on the line should be tiny
    def line_imag_report(name, sigs, ts, W):
        j = int(np.argmin(np.abs(sigs - 0.5)))
        col = W[:, j]
        med = float(np.nanmedian(np.abs(np.imag(col)) / (np.abs(col) + 1e-30)))
        print(f"  {name}: σ-grid at {sigs[j]:.4f}, median |Im|/|W| on the line = {med:.3e}")

    line_imag_report("zeta", sigs_z, ts_z, Wz)
    line_imag_report("eps", sigs_e, ts_e, We)
    line_imag_report("dh", sigs_d, ts_d, Wd)
    line_imag_report("fth", sigs_f, ts_f, Wf)

    print("=== fig2 domain colouring ===", flush=True)
    fig2_domain(sigs_z, ts_z, Wz, z_im, eps_off, sigs_e, ts_e, We, online_eps=eps_on)
    plt.close("all")

    print("=== fig4 zeta level curves ===", flush=True)
    fig4_level_zeta(sigs_z, ts_z, Wz, z_im, xi_ext)
    plt.close("all")

    # zeta zoom window for fair comparison in fig5 (first zeros)
    # slice the zeta grid
    def slice_grid(sigs, ts, W, tlo, thi):
        m = (ts >= tlo) & (ts <= thi)
        return sigs, ts[m], W[m, :]

    print("=== fig5 fakes ===", flush=True)
    s1, t1, w1 = slice_grid(sigs_z, ts_z, Wz, 10.0, 26.0)
    bundles = [
        dict(name=r"Riemann $\xi$  (first zeros: every red/blue crossing sits on the line)",
             sigs=s1, ts=t1, W=w1, on=z_im, off=None, extrema=xi_ext,
             note="no off-line crossing in this window"),
        dict(name=r"$F_\theta=\xi(s+\frac{1}{4})\xi(s-\frac{1}{4})$  (zeta’s zeros, rigidly shifted off the line)",
             sigs=sigs_f, ts=ts_f, W=Wf, on=None, off=fth_off, extrema=None,
             note=r"off-line zeros at $\frac{1}{2}\pm\frac{1}{4}+i\gamma$"),
        dict(name=r"Epstein $\xi_Q$, $Q=x^2+5y^2$  (first off-line zero near $0.933+15.668i$)",
             sigs=sigs_e, ts=ts_e, W=We, on=eps_on, off=eps_off, extrema=None,
             note="off-line crossings live on extra blue branches, between on-line zeros"),
        dict(name=r"Davenport–Heilbronn  (first off-line zero near $0.809+85.699i$)",
             sigs=sigs_d, ts=ts_d, W=Wd, on=dh_on, off=dh_off, extrema=None,
             note="same architecture, higher up the strip"),
    ]
    fig5_fakes(bundles)
    plt.close("all")
    print("done")


if __name__ == "__main__":
    main()
