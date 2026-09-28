"""
Self-similar blow-up profiles of the 1D Hou–Luo model (boundary model of 2D Boussinesq / 3D axisymmetric Euler):
    ω_t + u ω_x = θ_x,   θ_t + u θ_x = 0,   u_x = H ω      (Hf(x) = (1/π) p.v.∫ f(y)/(x−y) dy).
Ansatz as for 2D Boussinesq: ω = (1−t)^{-1}Ω(ξ), θ = (1−t)^{λ−1}Θ(ξ), u = (1−t)^λ U(ξ), ξ = x/(1−t)^{1+λ}:
    Ω + ((1+λ)ξ + U)Ω' = Θ',   (1−λ)Θ + ((1+λ)ξ + U)Θ' = 0,   U' = HΩ    (Θ even, Ω and U odd).
The local structure at the origin is identical to 2D: U ≈ −Aξ, Θ ≈ ξ^m, m = (λ−1)/(1+λ−A); smooth ⇔ m = 2.
Here A = −HΩ(0) = (2/π)∫_0^∞ Ω dξ/ξ.

Log variable η = ln ξ, D(η) = 1+λ+U/ξ (characteristic speed; D → ε = 1+λ−A at the origin):
    d ln Θ/dη = (λ−1)/D,       D dΩ/dη + Ω = e^{−η} dΘ/dη,
    U/ξ = e^{cη} (M ∗ e^{−cη}Ω),  M̂(k) = cot(π(c+ik)/2)/(c+ik+1)   (Mellin symbol of Ω ↦ U/ξ), −1/(1+λ) < c < m−1.
Unknown q = U/ξ on the grid; residual q − Q[q] (march Θ, Ω from the exact local solution, then Mellin multiplier).
Quadratures: 8th-order cumulative integration for ln Θ; 2-stage Gauss–Legendre (A-stable, order 4) for Ω with
8-point Lagrange interpolation of q and ln Θ to the Gauss points.
"""
import numpy as np
from numba import njit

SQ3 = np.sqrt(3.0)


def _lag_w(xs, x0, der=0):
    n = len(xs); w = np.zeros(n)
    for i in range(n):
        oth = [xs[k] for k in range(n) if k != i]
        den = np.prod([xs[i] - o for o in oth])
        if der == 0:
            w[i] = np.prod([x0 - o for o in oth]) / den
        else:
            w[i] = sum(np.prod([x0 - o for kk, o in enumerate(oth) if kk != k]) for k in range(len(oth))) / den
    return w


class _Stencils:
    """precomputed 8-point Lagrange weights (values at offsets c, interval integrals), shifted at the ends"""
    def __init__(self, N, cs):
        xs = np.arange(8.0)
        self.N = N
        self.off = {}
        for c in cs:
            wi = _lag_w(xs, 3 + c)
            ends = {}
            for j in list(range(0, 3)) + list(range(N - 4, N - 1)):
                k0 = min(max(j - 3, 0), N - 8)
                ends[j] = (k0, _lag_w(xs, j - k0 + c))
            self.off[c] = (wi, ends)
        gx, gw = np.polynomial.legendre.leggauss(6)
        def iw(pos):
            t = pos + 0.5 * (gx + 1)
            return sum(gw[q] * 0.5 * _lag_w(xs, t[q]) for q in range(6))
        self.iwi = iw(3.0)
        self.iends = {}
        for j in list(range(0, 3)) + list(range(N - 4, N - 1)):
            k0 = min(max(j - 3, 0), N - 8)
            self.iends[j] = (k0, iw(float(j - k0)))

    def offset(self, X, c):
        N = self.N; wi, ends = self.off[c]
        V = np.empty(N - 1)
        V[3:N - 4] = sum(wi[k] * X[k:N - 7 + k] for k in range(8))
        for j, (k0, w) in ends.items():
            V[j] = w @ X[k0:k0 + 8]
        return V

    def cumint(self, X, h, i0):
        N = self.N
        mid = np.empty(N - 1)
        mid[3:N - 4] = sum(self.iwi[k] * X[k:N - 7 + k] for k in range(8))
        for j, (k0, w) in self.iends.items():
            mid[j] = w @ X[k0:k0 + 8]
        Cm = np.concatenate([[0.0], np.cumsum(mid * h)])
        return Cm - Cm[i0]


@njit(cache=True)
def _gl2_omega(h, i0, Om0, a1, a2, f1, f2):
    """Ω' = a(η)(f(η) − Ω): 2-stage Gauss–Legendre, coefficients at the Gauss points of each interval"""
    N = len(a1) + 1
    Om = np.zeros(N)
    Om[i0] = Om0
    A11, A12, A21, A22 = 0.25, 0.25 - np.sqrt(3.0) / 6, 0.25 + np.sqrt(3.0) / 6, 0.25
    for j in range(i0, N - 1):
        y = Om[j]
        # K_i = a_i (f_i − y − h Σ A_ik K_k)
        m11 = 1 + h * a1[j] * A11; m12 = h * a1[j] * A12
        m21 = h * a2[j] * A21;     m22 = 1 + h * a2[j] * A22
        r1 = a1[j] * (f1[j] - y); r2 = a2[j] * (f2[j] - y)
        det = m11 * m22 - m12 * m21
        K1 = (r1 * m22 - m12 * r2) / det
        K2 = (m11 * r2 - m21 * r1) / det
        Om[j + 1] = y + 0.5 * h * (K1 + K2)
    return Om


class HL:
    def __init__(self, lam, L1=40.0, L2=160.0, N=8192, c=-0.25, eta_start=-30.0, sign=1.0):
        self.lam, self.c, self.sign = lam, c, sign
        self.N = N; self.h = (L1 + L2) / N
        self.eta = -L1 + self.h * np.arange(N)
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.h)
        z = c + 1j * k
        self.M = 1.0 / np.tan(np.pi * z / 2) / (z + 1)
        self.Ec = np.exp(c * self.eta); self.Emc = np.exp(-c * self.eta)
        self.i0 = int(np.argmin(np.abs(self.eta - eta_start)))
        self.nA = 40
        self.c1, self.c2 = 0.5 - SQ3 / 6, 0.5 + SQ3 / 6
        self.st = _Stencils(N, (self.c1, self.c2))

    def strain(self, q):
        return float(-np.mean(q[self.i0:self.i0 + self.nA]))

    def march(self, q):
        lam, h, eta, i0 = self.lam, self.h, self.eta, self.i0
        d = lam - 1
        A = self.strain(q)
        eps = 1 + lam - A
        m = d / eps
        D = (1 + lam) + q
        g = d / D - m
        g[:i0] = 0.0
        lnTh = m * eta + self.st.cumint(g, h, i0)
        lnTh[:i0] = m * eta[:i0]
        C = m / (A - 1)
        Om0 = self.sign * C * np.exp((m - 1) * eta[i0])
        # Gauss-point values
        e1 = eta[:-1] + self.c1 * h; e2 = eta[:-1] + self.c2 * h
        D1 = (1 + lam) + self.st.offset(q, self.c1); D2 = (1 + lam) + self.st.offset(q, self.c2)
        T1 = np.exp(self.st.offset(lnTh, self.c1)); T2 = np.exp(self.st.offset(lnTh, self.c2))
        f1 = self.sign * np.exp(-e1) * (d / D1) * T1; f2 = self.sign * np.exp(-e2) * (d / D2) * T2
        Om = _gl2_omega(h, i0, Om0, 1.0 / D1, 1.0 / D2, f1, f2)
        Om[:i0] = self.sign * C * np.exp((m - 1) * eta[:i0])
        Th = self.sign * np.exp(lnTh)
        return dict(A=A, m=m, eps=eps, Th=Th, Om=Om, D=D, Dmin=float(D.min()))

    def Q(self, q):
        r = self.march(q)
        qn = self.Ec * np.real(np.fft.ifft(self.M * np.fft.fft(self.Emc * r['Om'])))
        return qn, r

    def residual(self, q):
        qn, r = self.Q(q)
        R = q - qn
        R[:self.i0] = 0.0          # below the start the exact local structure (constant strain) is imposed
        return R, r

    def full(self, q):
        """extend q below the start as the constant strain"""
        q = q.copy(); q[:self.i0] = q[self.i0]
        return q

    def guess(self, A0):
        xi = np.exp(self.eta)
        return -A0 / (1 + xi ** 2) ** (0.5 / (1 + self.lam))


def newton(S, q, tol=1e-12, maxit=30, verbose=True, pert=1e-6):
    from scipy.sparse.linalg import LinearOperator, gmres
    q = S.full(q)
    R, info = S.residual(q)
    nrm = np.abs(R).max()
    for it in range(maxit):
        if verbose:
            print(f"   HL NK it {it}: |R|={nrm:.3e} A={info['A']:.12f} m={info['m']:.12f} Dmin={info['Dmin']:.5f}", flush=True)
        if nrm < tol:
            return q, info, True

        def mv(v):
            e = pert / max(np.abs(v).max(), 1e-300)
            vv = S.full(v) if False else v.copy()
            vv[:S.i0] = vv[S.i0]
            Rp, _ = S.residual(q + e * vv); Rm, _ = S.residual(q - e * vv)
            out = (Rp - Rm) / (2 * e)
            out[:S.i0] = v[:S.i0]           # identity on the frozen part
            return out
        J = LinearOperator((S.N, S.N), matvec=mv, dtype=float)
        dx, gi = gmres(J, -R, rtol=1e-11, atol=0.0, restart=100, maxiter=5)
        dx[:S.i0] = dx[S.i0]
        t = 1.0
        while t > 1e-3:
            Rn, infon = S.residual(q + t * dx)
            if infon['Dmin'] > 0 and np.all(np.isfinite(Rn)) and np.abs(Rn).max() < (1 - 0.1 * t) * nrm:
                break
            t *= 0.5
        if t <= 1e-3:
            return q, info, False
        q, R, info = q + t * dx, Rn, infon
        nrm = np.abs(R).max()
    return q, info, nrm < tol
