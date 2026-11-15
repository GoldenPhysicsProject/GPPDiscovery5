import mpmath as mp
mp.mp.dps = 30
def arch_multi(ds, K=4000, m=40):
    # Gauss-Legendre with m nodes on each period [2πk, 2π(k+1)], k < K; tail beyond 2πK is < 1e-11.
    xs, ws = zip(*[(mp.mpf(a), mp.mpf(b)) for a, b in mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(6, mp.mp.prec)][:0]) if False else (None, None)
    gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
    nodes = gl.calc_nodes(6, mp.mp.prec)   # degree 6 -> 3*2^5 = 96 nodes on [-1,1]
    tot = [mp.mpf(0) for _ in ds]
    h = mp.pi  # half-width
    lp = mp.log(mp.pi)
    for k in range(K):
        c = (2*k+1)*mp.pi
        for (u, w) in nodes:
            t = c + h*u
            psi2 = (mp.sin(t/2)/(t/2))**4
            g = psi2*(mp.re(mp.digamma(mp.mpc(mp.mpf('0.25'), t/2))) - lp)*w*h
            for i, d in enumerate(ds):
                tot[i] += g*mp.cos(t*d)
    return [v/mp.pi for v in tot]
if __name__ == "__main__":
    import sys, time
    t0=time.time()
    K = int(sys.argv[1]) if len(sys.argv)>1 else 200
    v = arch_multi([mp.mpf(0)], K=K)
    print(K, mp.nstr(v[0], 20), mp.nstr(v[0] - mp.mpf('-1.54293305945377535652'), 5), round(time.time()-t0,1))
