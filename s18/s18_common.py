"""S18 shared helpers (reuse S17: Tarb, consts, E_comb, build). Daniel Toupin / Golden Physics Project, AI assistance: Grok."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "s17"))
from s17_build import Tarb, consts, build, forms, vM, arb, acb, fmpq, ctx, arb_mat
from s17_bounds import E_comb
from flint import dirichlet_char
from math import gcd
ctx.prec = 256
HSTAR = arb(3).log()/2
LOG2S2 = arb(2).log()/arb(2).sqrt()
def T_eps(h, e):            # Codex T_eps(h) = 2[atanh(e^{-h/2}) +- atan(e^{-h/2})] ; S17 Tarb(L) uses r=e^{-L/2}
    return Tarb(h, e)
def primitive_chars(f):
    """all primitive characters of conductor f: (conrey label, parity, order, eta(2) as acb)"""
    out = []
    for l in range(1, f+1):
        if gcd(l, f) != 1: continue
        c = dirichlet_char(f, l)
        if c.is_primitive(): out.append((l, c.parity(), c.order(), c(2), c))
    return out
