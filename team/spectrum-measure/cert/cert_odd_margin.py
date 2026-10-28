"""CERT (Arb ball arithmetic, python-flint): odd-sector H-margin of the archimedean form R at 2L = log 2 + ell.
Spectrum & Measure seat, GPP RH team. Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN.

Certifies a rigorous lower bound for
   inf over odd xi in C_c^inf(-L, L), ||xi||_H = 1, of  R(xi) - c2 * [ (1+1/eta)(2 ell T/pi) nu(<T) + (1+eta) nu(>=T) ]
(ell > 0), or of R(xi) (ell = 0, = m0^odd), where R(xi) = int Omega dnu, Omega(t) = Re psi(1/4 + it/2) - log pi,
nu = H-energy distribution on t >= 0: dnu = (1/pi)(1/4 + t^2)|xi^(t)|^2 dt.
Method (see CERT_p2_constants.md): nu([0,T]) <= D(T) = int_0^T dens, dens(t) = (1/pi)(1/4+t^2) sup|xi^(t)|^2/||xi||_H^2
(closed form via the Dirichlet Green function of -d^2 + 1/4 on (-L, L)), then the bathtub bound for increasing integrands.
Semantics: the step minorant uses the TRUE values omega~(e_i) (enclosed by balls); the printed ball encloses the
computed bound for this grid (a rigorous lower bound for the infimum), not the infimum itself; runs at different
precisions use slightly different (rigorous) upper sums for D, so their bounds may differ in late digits.
All quantities are arb balls; every ball is checked finite; integrals are bounded by sup over subintervals (ball evaluation).
"""
import sys
from flint import arb, acb, ctx

def chk(x, what):
    if not x.is_finite():
        raise SystemExit("NON-FINITE BALL: %s = %s" % (what, x))
    return x

def dens_odd(t, h):
    # (1/pi)[ h - h*sinc(2th) - 2 sin(th) ( (1/2)coth(h/2) sin(th) - t cos(th) ) / (t^2 + 1/4) ]
    half = arb(1)/2
    s = (t*h).sin(); c = (t*h).cos()
    coth = (h/2).cosh()/(h/2).sinh()
    v = (h - h*(2*t*h).sinc() - 2*s*(half*coth*s - t*c)/(t*t + arb(1)/4))/arb.pi()
    return chk(v, "dens(%s)" % t)

def Omega(t):
    v = (acb(arb(1)/4, t/2).digamma()).real - arb.pi().log()
    return chk(v, "Omega(%s)" % t)

def certify(ell, T, eta, Tmax, N, lohi=None):
    """ell, T, eta, Tmax exact rationals as strings; N subintervals of [0, Tmax]. Returns arb lower bound (ball)."""
    a = arb(2).log(); c2 = a/arb(2).sqrt()
    ell = arb(ell); T = arb(T); eta = arb(eta); Tmax = arb(Tmax)
    h = (a + ell)/2   # L
    dt = Tmax/N
    edges = [dt*i for i in range(N+1)]
    # D at right edges: sum of sup(dens) over [e_i, e_{i+1}] times dt (upper endpoint used)
    D = []; acc = arb(0)
    for i in range(N):
        box = arb((edges[i]+edges[i+1])/2, dt/2)   # ball covering [e_i, e_{i+1}]
        d = dens_odd(box, h)
        up = d.upper()
        if up < 0: up = arb(0)
        acc = acc + up*dt
        D.append(acc.upper())
    if lohi is not None:
        lo = arb(lohi[0]); hi = arb(lohi[1])   # Celestial C.13.H band-split: (mu_lo, mu_hi)
    elif ell == 0:
        lo = arb(0); hi = arb(0)
    else:
        lo = c2*(1 + 1/eta)*2*ell*T/arb.pi(); hi = c2*(1 + eta)
    OmT = Omega(T)
    def wt(t):  # increasing minorant omega~
        if t < T: 
            x = Omega(t) - lo; y = OmT - hi
            # min of balls: exact if disjoint, else take lower envelope safely
            if x < y: return x
            if y < x: return y
            return x.union(y)  # overlap: a ball containing both contains the true min (keeps exact-step semantics)
        return Omega(t) - hi
    # make sure T is an edge (choose N so that T/dt integer)
    total = wt(edges[0])
    prev = total
    for i in range(N):
        cur = wt(edges[i+1])
        inc = cur - prev
        f = 1 - D[i]
        f = f.lower()
        if f < 0: f = arb(0)
        # increments are >= 0 in truth; using the ball keeps rigour (lower endpoint of final sum)
        total = total + inc*f
        prev = cur
    return chk(total, "total"), D

if __name__ == "__main__":
    ctx.prec = int(sys.argv[1]) if len(sys.argv) > 1 else 80
    print("prec bits:", ctx.prec)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 2**16
    for (ell, T, eta, Tmax) in [("0", "23", "1", "64"), ("0.025", "23", "1", "64"), ("0.028", "22", "1", "64")]:
        tot, D = certify(ell, T, eta, Tmax, N)
        print("ell=%s T=%s eta=%s Tmax=%s N=%d : lower bound ball %s ; certified lower endpoint %s ; D(Tmax) upper %s"
              % (ell, T, eta, Tmax, N, tot, tot.lower(), D[-1]))
    # Celestial Conj C.13.H at L0 = (1/2)log 2, odd: R >= mu_hi nu(>=T) + mu_lo nu(<T); certified slack (lower bound of
    # inf [R - mu_hi nu_> - mu_lo nu_<]/||xi||_H^2) must be >= 0.
    for (T, mulo, muhi) in [("12.875", "0.25", "0.6"), ("16", "0.3", "0.7")]:
        tot, D = certify("0", T, "1", "64", N, lohi=(mulo, muhi))
        print("C.13.H  T=%s mu_lo=%s mu_hi=%s : certified slack lower endpoint %s" % (T, mulo, muhi, tot.lower()))
    a = arb(2).log(); c2 = a/arb(2).sqrt()
    print("c2 =", c2, "; Omega(12.8) =", Omega(arb("12.8")), "; Omega(12.875) =", Omega(arb("12.875")), "; Omega(16) =", Omega(arb(16)))
