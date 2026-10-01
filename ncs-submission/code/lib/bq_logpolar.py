"""
Self-similar blow-up profiles of the 2D Boussinesq equations with boundary (the 2D analogue of 3D
axisymmetric Euler with boundary, Hou–Luo / Chen–Hou scenario), in LOG-POLAR coordinates.

Physical problem (upper half-plane y2 ≥ 0, no-penetration on y2 = 0):
    ω_t + u·∇ω = ∂_x θ,   θ_t + u·∇θ = 0,   u = (ψ_y, −ψ_x),  −Δψ = ω,  ψ = 0 on y2 = 0.
Self-similar ansatz (T = 1):  ω = (1−t)^{-1} Ω(y),  θ = (1−t)^{λ−1} Θ(y),  u = (1−t)^λ U(y),  y = x/(1−t)^{1+λ}:
    Ω + ((1+λ)y + U)·∇Ω = ∂_{y1}Θ,     (1−λ)Θ + ((1+λ)y + U)·∇Θ = 0,     U = (Ψ_2, −Ψ_1), −ΔΨ = Ω.
Symmetry: Θ even in y1, Ω and Ψ odd in y1  →  quarter plane, polar angle β ∈ [0, π/2] from the boundary.

Local structure at the stagnation point y = 0 (strain U ≈ (−A y1, A y2), A = −∂_1U_1(0)):
    Θ ≈ c|y1|^m,  m = (λ−1)/(1+λ−A)   ("least singular" family; smooth profiles ⇔ m = 2 ⇔ λ = 2A − 3).
This is the 2D analogue of the CCF local exponent p(λ) = λ/(1+λ+h1): the strain A is a nonlocal functional of Ω.

Log-polar variables s = ln r, β:  V·∇ = (V_r/r)∂_s + w ∂_β,  V_r/r = (1+λ) + U_r/r,  w = U_β/r,
    U_r/r = e^{−2s} Ψ_β,   w = −e^{−2s} Ψ_s,   Ψ_ss + Ψ_ββ = −e^{2s} Ω.
The axis singularity of the non-smooth family members is factored out exactly:
    Θ = (cos β)^m Θ̂,  Ω = (cos β)^{m−1} Ω̂   (Θ̂, Ω̂ smooth), giving the marching equations
    Θ̂_s = [(λ−1)Θ̂ + m (w tanβ) Θ̂ − w Θ̂_β] / (V_r/r)
    Ω̂_s = [S − Ω̂ + (m−1)(w tanβ) Ω̂ − w Ω̂_β] / (V_r/r),   S = e^{−s}[cos²β Θ̂_s − sinβ cosβ Θ̂_β + m sin²β Θ̂].
All characteristics leave the origin (V_r > 0: no sonic point), so Θ̂, Ω̂ are obtained by marching in s from the
exact local solution  Θ̂ = e^{ms},  Ω̂ = m/(A−1) e^{(m−1)s}  (normalisation c = 1).
Biot–Savart: Ψ = e^{as}Ψ̃, (∂_s + a)²Ψ̃ + Ψ̃_ββ = −e^{(2−a)s}Ω, FFT in s, Chebyshev collocation in β (Ψ = 0 at β = 0,
π/2), with (1+2λ)/(1+λ) < a < 2 (excludes the harmonic strain r² sin 2β at infinity and r^{-2} at the origin).
"""
import numpy as np
from numba import njit


def cheb(N):
    """Trefethen's Chebyshev differentiation matrix on x_j = cos(πj/N), j = 0..N."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1); c[0] = c[-1] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


@njit(cache=True)
def _march(s, h, Th0, Om0, Ur, w, wt, Urh, wh, wth, Db, cb, sb, lam, m, i0):
    """RK4 march in s of (Θ̂, Ω̂) from index i0 to the end; U fields at grid points and half points."""
    Ns, Nb1 = Ur.shape
    Th = np.zeros((Ns, Nb1)); Om = np.zeros((Ns, Nb1))
    Th[i0] = Th0; Om[i0] = Om0
    cb2 = cb * cb; sbcb = sb * cb; sb2 = sb * sb

    def rhs(sv, th, om, ur, ww, wtt):
        vr = (1.0 + lam) + ur
        thb = Db @ th
        ths = ((lam - 1.0) * th + m * wtt * th - ww * thb) / vr
        S = np.exp(-sv) * (cb2 * ths - sbcb * thb + m * sb2 * th)
        omb = Db @ om
        oms = (S - om + (m - 1.0) * wtt * om - ww * omb) / vr
        return ths, oms

    for j in range(i0, Ns - 1):
        th = Th[j]; om = Om[j]; sv = s[j]
        k1t, k1o = rhs(sv, th, om, Ur[j], w[j], wt[j])
        k2t, k2o = rhs(sv + 0.5 * h, th + 0.5 * h * k1t, om + 0.5 * h * k1o, Urh[j], wh[j], wth[j])
        k3t, k3o = rhs(sv + 0.5 * h, th + 0.5 * h * k2t, om + 0.5 * h * k2o, Urh[j], wh[j], wth[j])
        k4t, k4o = rhs(sv + h, th + h * k3t, om + h * k3o, Ur[j + 1], w[j + 1], wt[j + 1])
        Th[j + 1] = th + h / 6.0 * (k1t + 2 * k2t + 2 * k3t + k4t)
        Om[j + 1] = om + h / 6.0 * (k1o + 2 * k2o + 2 * k3o + k4o)
    return Th, Om


class BQLogPolar:
    def __init__(self, lam, s_min=-120.0, s_max=100.0, hs=0.025, Nb=32, a=None, s_start=-20.0, sign=-1.0):
        self.lam = lam
        self.sign = sign        # Θ ≈ sign·|y1|^m near the origin (sign = −1 gives compressive strain A > 0)
        if a is None:
            lo = (1 + 2 * lam) / (1 + lam)
            a = 0.5 * (lo + 2.0)
        self.a = a
        self.Ns = int(round((s_max - s_min) / hs)); self.Ns += self.Ns % 2
        self.hs = hs
        self.s = s_min + hs * np.arange(self.Ns)
        self.i0 = int(np.argmin(np.abs(self.s - s_start)))
        D, x = cheb(Nb)
        self.Nb = Nb
        self.beta = (np.pi / 4) * (1 - x)                      # 0 … π/2
        self.Db = -(4.0 / np.pi) * D
        self.Dbb = self.Db @ self.Db
        self.cb, self.sb = np.cos(self.beta), np.sin(self.beta)
        self.cb[-1] = 0.0
        k = 2 * np.pi * np.fft.fftfreq(self.Ns, d=hs)
        self.k = k
        # Poisson operators per Fourier mode (interior Chebyshev points, Dirichlet ends)
        L0 = self.Dbb[1:-1, 1:-1]
        n = Nb - 1
        self.Linv = np.linalg.inv(L0[None, :, :] + ((1j * k + a) ** 2)[:, None, None] * np.eye(n)[None, :, :])
        self.ea2 = np.exp((a - 2) * self.s)
        self.e2a = np.exp((2 - a) * self.s)

    # ---------------- Biot–Savart ----------------
    def poisson(self, Om):
        """Ψ̃ (Ns×(Nb+1)) from Ω on the grid: (∂_s+a)²Ψ̃ + Ψ̃_ββ = −e^{(2−a)s}Ω, Ψ̃ = 0 at β = 0, π/2."""
        F = -self.e2a[:, None] * Om
        Fh = np.fft.fft(F[:, 1:-1], axis=0)
        Ph = np.einsum('kij,kj->ki', self.Linv, Fh)
        P = np.zeros_like(Om)
        P[:, 1:-1] = np.real(np.fft.ifft(Ph, axis=0))
        return P

    def velocity(self, P):
        """U_r/r, w = U_β/r, w·tanβ on grid and at half points s_j + h/2 (spectral interpolation)."""
        Pb = P @ self.Db.T
        Ph = np.fft.fft(P, axis=0)
        Ps = np.real(np.fft.ifft(1j * self.k[:, None] * Ph, axis=0))
        shift = np.exp(1j * self.k * self.hs / 2)[:, None]
        P_h = np.real(np.fft.ifft(Ph * shift, axis=0))
        Pb_h = P_h @ self.Db.T
        Ps_h = np.real(np.fft.ifft(1j * self.k[:, None] * Ph * shift, axis=0))
        eah = np.exp((self.a - 2) * (self.s + self.hs / 2))
        Ur = self.ea2[:, None] * Pb
        w = -self.ea2[:, None] * (self.a * P + Ps)
        Urh = eah[:, None] * Pb_h
        wh = -eah[:, None] * (self.a * P_h + Ps_h)
        return Ur, w, self._wtan(w), Urh, wh, self._wtan(wh)

    def _wtan(self, w):
        t = np.empty_like(w)
        t[:, :-1] = w[:, :-1] * (self.sb[:-1] / self.cb[:-1])
        t[:, -1] = -(w @ self.Db.T)[:, -1]                    # limit β → π/2 (w = 0 there)
        return t

    # ---------------- transport ----------------
    def strain(self, Ur):
        """A = −∂_1U_1(0) from U_r/r ≈ −A cos 2β at small r (β = 0 value, averaged over a window)."""
        j = slice(self.i0, self.i0 + 40)
        return float(-np.mean(Ur[j, 0]))

    def march(self, P, return_all=False, A_fixed=None):
        lam = self.lam
        Ur, w, wt, Urh, wh, wth = self.velocity(P)
        A = self.strain(Ur) if A_fixed is None else A_fixed
        m = (lam - 1.0) / (1.0 + lam - A)
        C = m / (A - 1.0)
        s0 = self.s[self.i0]
        Th0 = np.full(self.Nb + 1, self.sign * np.exp(m * s0))
        Om0 = np.full(self.Nb + 1, self.sign * C * np.exp((m - 1) * s0))
        Th, Om = _march(self.s, self.hs, Th0, Om0, Ur, w, wt, Urh, wh, wth, self.Db, self.cb, self.sb, lam, m, self.i0)
        # below s_start: exact local solution
        sl = self.s[:self.i0]
        Th[:self.i0] = self.sign * np.exp(m * sl)[:, None]
        Om[:self.i0] = self.sign * C * np.exp((m - 1) * sl)[:, None]
        cpow = np.where(self.cb > 0, self.cb, 0.0)
        Omega = Om * cpow[None, :] ** (m - 1)
        out = dict(A=A, m=m, Th=Th, Om=Om, Omega=Omega, vrmin=float(((1 + lam) + Ur).min()))
        if return_all:
            out.update(Ur=Ur, w=w)
        return out

    def T(self, P, A_fixed=None):
        r = self.march(P, A_fixed=A_fixed)
        return self.poisson(r['Omega']), r
