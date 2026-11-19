"""F_{p^2} and hyperelliptic function-field arithmetic for Round 16 (Bombieri-Stepanov executed).
F_q, q = p^2 = F_p[r]/(r^2 - nr).  Elements: ints e = u0 + p*u1.  Polynomials over F_q: numpy int64 arrays shape (2, n)
(component 0 = u0 coeffs, 1 = u1 coeffs, index = power of x).  Functions on y^2 = f(x) (deg f = 2g+1, one point O at
infinity, ord_O x = -2, ord_O y = -(2g+1)) as pairs (a, b) meaning a(x) + b(x) y."""
import numpy as np

class Fq:
    def __init__(self, p):
        self.p = p; self.q = p*p
        self.nr = next(x for x in range(2, p) if pow(x, (p-1)//2, p) == p-1)
        q = self.q; u0 = np.arange(q) % p; u1 = np.arange(q) // p
        self.M = ((u0[:, None]*u0[None, :] + self.nr*u1[:, None]*u1[None, :]) % p
                  + p*((u0[:, None]*u1[None, :] + u1[:, None]*u0[None, :]) % p))
        self.A = ((u0[:, None]+u0[None, :]) % p) + p*((u1[:, None]+u1[None, :]) % p)
        self.NEG = ((-u0) % p) + p*((-u1) % p)
        self.INV = np.zeros(q, dtype=np.int64)
        for a in range(1, q): self.INV[a] = int(np.nonzero(self.M[a] == 1)[0][0])
        self.FROB = np.array([self.pw(a, p) for a in range(q)])          # a -> a^p (conjugation)
        sq = set(int(self.M[a, a]) for a in range(1, q)); self.nonsq = next(a for a in range(1, q) if a not in sq)
        self.CHI = np.array([0] + [1 if a in sq else -1 for a in range(1, q)])
    def pw(self, a, n):
        r = 1
        while n:
            if n & 1: r = int(self.M[r, a])
            a = int(self.M[a, a]); n >>= 1
        return r
    # polynomial layer
    def P(self, coeffs):                 # list of field ints -> (2, n) array
        c = np.array(coeffs, dtype=np.int64); return np.stack([c % self.p, c // self.p])
    def trim(self, a):
        nz = np.nonzero(a[0] | a[1])[0]; return a[:, :nz[-1]+1] if len(nz) else a[:, :0]
    def add(self, a, b):
        n = max(a.shape[1], b.shape[1]); out = np.zeros((2, n), dtype=np.int64)
        out[:, :a.shape[1]] += a; out[:, :b.shape[1]] += b; return self.trim(out % self.p)
    def smul(self, c, a):                # field scalar times poly
        c0, c1 = c % self.p, c // self.p
        return self.trim(np.stack([(c0*a[0] + self.nr*c1*a[1]) % self.p, (c0*a[1] + c1*a[0]) % self.p]))
    def mul(self, a, b):
        if a.shape[1] == 0 or b.shape[1] == 0: return a[:, :0]
        cv = lambda u, v: np.convolve(u, v) % self.p if min(len(u), len(v)) < 64 else self._fconv(u, v)
        r0 = (cv(a[0], b[0]) + self.nr*cv(a[1], b[1])) % self.p; r1 = (cv(a[0], b[1]) + cv(a[1], b[0])) % self.p
        return self.trim(np.stack([r0, r1]))
    def _fconv(self, u, v):              # exact integer convolution via split FFT-free blocks (values < p, safe in int64)
        return np.convolve(u.astype(np.int64), v.astype(np.int64)) % self.p
    def shiftx(self, a, k):
        return np.concatenate([np.zeros((2, k), dtype=np.int64), a], axis=1) if a.shape[1] else a
    def deg(self, a): return a.shape[1] - 1
    def ev(self, a, x):                  # Horner at field element x
        r = 0
        for k in range(a.shape[1]-1, -1, -1): r = int(self.A[self.M[r, x], a[0, k] + self.p*a[1, k]])
        return r
    def frobpoly(self, a, e=1):          # coefficients -> c^(p^e)  (x untouched)
        c = a[0] + self.p*a[1]
        for _ in range(e): c = self.FROB[c]
        return np.stack([c % self.p, c // self.p])

class Curve:
    """y^2 = f(x), f given as field ints low->high, deg f = 2g+1."""
    def __init__(self, F, f):
        self.F = F; self.f = F.trim(F.P(f)); self.d = F.deg(self.f); self.g = (self.d - 1)//2
        self._fp = {0: F.P([1])}
    def fpow(self, k):
        if k not in self._fp:
            h = k//2; A = self.fpow(h); B = self.F.mul(A, A)
            self._fp[k] = self.F.mul(B, self.f) if k % 2 else B
        return self._fp[k]
    def mul(self, u, v):
        F = self.F; a1, b1 = u; a2, b2 = v
        return (F.add(F.mul(a1, a2), F.mul(F.mul(b1, b2), self.f)), F.add(F.mul(a1, b2), F.mul(a2, b1)))
    def add(self, u, v): return (self.F.add(u[0], v[0]), self.F.add(u[1], v[1]))
    def smul(self, c, u): return (self.F.smul(c, u[0]), self.F.smul(c, u[1]))
    def pole(self, u):
        a, b = u; o = []
        if a.shape[1]: o.append(2*self.F.deg(a))
        if b.shape[1]: o.append(2*self.F.deg(b) + self.d)
        return max(o) if o else None
    def mono(self, i, j):                # x^i y^j, j in {0,1}
        one = self.F.shiftx(self.F.P([1]), i); z = one[:, :0]
        return (one, z) if j == 0 else (z, one)
    def basis(self, m):                  # L(mO): monomials with 2i + j d <= m (distinct pole orders)
        return [(i, j) for j in (0, 1) for i in range(m//2 + 1) if 2*i + j*self.d <= m]
    def monopow(self, i, j, Q):          # (x^i y^j)^Q for Q a power of p, via additivity: x^{iQ} y^{jQ}, y^Q = y f^{(Q-1)/2}
        F = self.F; X = F.shiftx(F.P([1]), i*Q); z = X[:, :0]
        return (X, z) if j == 0 else (z, F.mul(X, self.fpow((Q-1)//2)))
    def ev(self, u, x, y):
        F = self.F; return int(F.A[F.ev(u[0], x), F.M[F.ev(u[1], x), y]])
    def points(self):                    # affine F_q-points
        F = self.F; pts = []
        for x in range(F.q):
            r = F.ev(self.f, x)
            for y in range(F.q):
                if F.M[y, y] == r: pts.append((x, y))
        return pts
