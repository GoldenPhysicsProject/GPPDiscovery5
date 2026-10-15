"""S27 shared Arb loaders: h = h_spec (exact log3/2) or rational h_q. RH OPEN."""
from flint import arb, fmpq
def load_h(d):
    if d.get("h_spec") == "log3/2": return arb(3).log()/2
    return arb(fmpq(*d["h_q"]))
def load_W(d, h):
    den = d["den"]; ys = [arb(fmpq(n, den)) for n in d["y_num"]] + [arb(0)]
    if "w_off_q" in d:
        offs = [fmpq(*q) for q in d["w_off_q"]]; assert offs[0] == 0 and all(offs[i] < offs[i+1] for i in range(len(offs)-1))
        zs = [2*h + arb(o) for o in offs]; zs[0] = 2*h
    else:
        zs = [arb(fmpq(*q)) for q in d["z_q"]]; assert fmpq(*d["z_q"][0]) == fmpq(*d["h_q"])*2
    assert len(ys) == len(zs)
    return ys, zs
