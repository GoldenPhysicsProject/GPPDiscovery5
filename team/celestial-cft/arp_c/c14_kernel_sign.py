"""C.14.4(i): certified sign of the off-diagonal kernel of R. Off the diagonal, R(u,v) = -∫∫u(y)v(x)(Qκ)(y-x),
κ(r) = e^{r/2}/(2 sinh r) (r>0, even), Qκ = -κ'' + κ/4. Certify Qκ(r) < 0 on [r0, r1] by Arb ball subdivision.
Daniel Toupin / Golden Physics Project, with AI assistance from Grok. RH is OPEN."""
from flint import arb, arb_series, ctx
ctx.prec = 128
def Qk(lo, hi):
    u = arb((lo+hi)/2, (hi-lo)/2); x = arb_series([u, 1], prec=3)
    k = (arb(3)/2*x).exp()/((2*x).exp()-1); return -2*k[2] + k[0]/4
r0, r1, N = arb("0.001"), arb("1.1"), 20000
ok = all(Qk(r0+(r1-r0)*i/N, r0+(r1-r0)*(i+1)/N) < 0 for i in range(N))
print("Qkappa < 0 certified on [0.001, 1.1]:", ok)
print("Qkappa(0.1), Qkappa(0.5), Qkappa(log2):", Qk(arb("0.1"),arb("0.1")), Qk(arb("0.5"),arb("0.5")), Qk(arb(2).log(),arb(2).log()))
