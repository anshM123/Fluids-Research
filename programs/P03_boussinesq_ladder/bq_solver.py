"""
Improved log-polar solver for self-similar 2D Boussinesq blow-up profiles (see bq_logpolar.py for the formulation).

Changes with respect to bq_logpolar.BQLogPolar:
  1. Velocity from X = Ψ/r² with 8th-order finite differences in s (U_r/r = X_β, w = −(2X + X_s)).
     The FFT derivative of Ψ̃ = e^{(2−a)s}X used before assumes periodicity in s. Newton directions that do not
     decay at large s (e.g. a uniform strain) grow like e^{(2−a)s} ~ e^{17} at s = 100, and the Gibbs ringing
     corrupted the Jacobian. This caused the residual floor ~3e-6 in the strain mode.
  2. Hybrid march: 2-stage Gauss–Legendre (implicit, A-stable, order 4) for s < s_sw, RK4 beyond.
     Near the stagnation point the angular transport rate w/(V_r/r) ~ 2Aβ/ε, ε = 1+λ−A, becomes stiff as λ → 1
     (ε = (λ−1)/2 on smooth profiles), and explicit RK4 is unstable for λ ≲ 1.3.
The Biot–Savart solve (FFT in s of the decaying forcing e^{(2−a)s}Ω, Chebyshev in β) is unchanged.
"""
import numpy as np
from numba import njit
from bq_logpolar import BQLogPolar

SQ3 = np.sqrt(3.0)
C1, C2 = 0.5 - SQ3 / 6, 0.5 + SQ3 / 6
A11, A12, A21, A22 = 0.25, 0.25 - SQ3 / 6, 0.25 + SQ3 / 6, 0.25


def _lag_w(xs, x0, der):
    n = len(xs); w = np.zeros(n)
    for i in range(n):
        oth = [xs[k] for k in range(n) if k != i]
        den = np.prod([xs[i] - o for o in oth])
        if der == 0:
            w[i] = np.prod([x0 - o for o in oth]) / den
        else:
            w[i] = sum(np.prod([x0 - o for kk, o in enumerate(oth) if kk != k]) for k in range(len(oth))) / den
    return w


def _offset_values(X, h, c):
    """values and s-derivatives at s_j + c h (j = 0..Ns−2) from 8-point Lagrange stencils (shifted at the ends)"""
    Ns = X.shape[0]
    V = np.empty((Ns - 1,) + X.shape[1:]); D = np.empty_like(V)
    xs = np.arange(8.0)
    # interior: stencil j−3..j+4, evaluation point 3 + c
    w0 = _lag_w(xs, 3 + c, 0); w1 = _lag_w(xs, 3 + c, 1) / h
    jlo, jhi = 3, Ns - 5                                  # intervals with a full centred stencil: j = 3 .. Ns−5
    V[jlo:jhi + 1] = 0.0; D[jlo:jhi + 1] = 0.0
    for k in range(8):
        seg = X[jlo - 3 + k: jhi - 3 + k + 1]
        V[jlo:jhi + 1] += w0[k] * seg; D[jlo:jhi + 1] += w1[k] * seg
    for j in list(range(0, jlo)) + list(range(jhi + 1, Ns - 1)):
        k0 = min(max(j - 3, 0), Ns - 8)
        wv = _lag_w(xs, j - k0 + c, 0); wd = _lag_w(xs, j - k0 + c, 1) / h
        V[j] = np.tensordot(wv, X[k0:k0 + 8], axes=(0, 0)); D[j] = np.tensordot(wd, X[k0:k0 + 8], axes=(0, 0))
    return V, D


def _grid_deriv(X, h):
    """8th-order s-derivative at the grid points (9-point centred, shifted one-sided at the ends)"""
    Ns = X.shape[0]
    D = np.empty_like(X)
    xs = np.arange(9.0)
    w = _lag_w(xs, 4.0, 1) / h
    D[4:Ns - 4] = 0.0
    for k in range(9):
        D[4:Ns - 4] += w[k] * X[k:Ns - 8 + k]
    for j in list(range(0, 4)) + list(range(Ns - 4, Ns)):
        k0 = min(max(j - 4, 0), Ns - 9)
        wd = _lag_w(xs, float(j - k0), 1) / h
        D[j] = np.tensordot(wd, X[k0:k0 + 9], axes=(0, 0))
    return D


@njit(cache=True)
def _march_hybrid(s, h, Th0, Om0, Ur, w, wt, Urh, wh, wth, Ur1, w1, wt1, Ur2, w2, wt2,
                  Db, cb, sb, lam, m, i0, jsw):
    Ns, n = Ur.shape
    Th = np.zeros((Ns, n)); Om = np.zeros((Ns, n))
    Th[i0] = Th0; Om[i0] = Om0
    cb2 = cb * cb; sbcb = sb * cb; sb2 = sb * sb
    I = np.eye(n)
    a11, a12, a21, a22 = 0.25, 0.25 - np.sqrt(3.0) / 6, 0.25 + np.sqrt(3.0) / 6, 0.25
    c1, c2 = 0.5 - np.sqrt(3.0) / 6, 0.5 + np.sqrt(3.0) / 6
    M = np.zeros((2 * n, 2 * n)); rhs = np.zeros(2 * n)

    def f_rk(sv, th, om, ur, ww, wtt):
        vr = (1.0 + lam) + ur
        thb = Db @ th
        ths = ((lam - 1.0) * th + m * wtt * th - ww * thb) / vr
        S = np.exp(-sv) * (cb2 * ths - sbcb * thb + m * sb2 * th)
        oms = (S - om + (m - 1.0) * wtt * om - ww * (Db @ om)) / vr
        return ths, oms

    for j in range(i0, Ns - 1):
        sv = s[j]; th = Th[j]; om = Om[j]
        if j < jsw:
            # ---- Gauss–Legendre 2-stage (implicit) ----
            vr1 = (1.0 + lam) + Ur1[j]; vr2 = (1.0 + lam) + Ur2[j]
            P1 = np.empty((n, n)); P2 = np.empty((n, n)); Q1 = np.empty((n, n)); Q2 = np.empty((n, n))
            for i in range(n):
                for k in range(n):
                    P1[i, k] = -w1[j, i] * Db[i, k] / vr1[i]
                    P2[i, k] = -w2[j, i] * Db[i, k] / vr2[i]
                Q1[i, :] = P1[i, :]; Q2[i, :] = P2[i, :]
                P1[i, i] += ((lam - 1.0) + m * wt1[j, i]) / vr1[i]
                P2[i, i] += ((lam - 1.0) + m * wt2[j, i]) / vr2[i]
                Q1[i, i] += (-1.0 + (m - 1.0) * wt1[j, i]) / vr1[i]
                Q2[i, i] += (-1.0 + (m - 1.0) * wt2[j, i]) / vr2[i]
            # Θ̂ stages
            M[:n, :n] = I - h * a11 * P1; M[:n, n:] = -h * a12 * P1
            M[n:, :n] = -h * a21 * P2;    M[n:, n:] = I - h * a22 * P2
            rhs[:n] = P1 @ th; rhs[n:] = P2 @ th
            K = np.linalg.solve(M, rhs)
            K1 = K[:n]; K2 = K[n:]
            Y1 = th + h * (a11 * K1 + a12 * K2); Y2 = th + h * (a21 * K1 + a22 * K2)
            S1 = np.exp(-(sv + c1 * h)) * (cb2 * K1 - sbcb * (Db @ Y1) + m * sb2 * Y1) / vr1
            S2 = np.exp(-(sv + c2 * h)) * (cb2 * K2 - sbcb * (Db @ Y2) + m * sb2 * Y2) / vr2
            # Ω̂ stages
            M[:n, :n] = I - h * a11 * Q1; M[:n, n:] = -h * a12 * Q1
            M[n:, :n] = -h * a21 * Q2;    M[n:, n:] = I - h * a22 * Q2
            rhs[:n] = Q1 @ om + S1; rhs[n:] = Q2 @ om + S2
            L = np.linalg.solve(M, rhs)
            Th[j + 1] = th + 0.5 * h * (K1 + K2)
            Om[j + 1] = om + 0.5 * h * (L[:n] + L[n:])
        else:
            k1t, k1o = f_rk(sv, th, om, Ur[j], w[j], wt[j])
            k2t, k2o = f_rk(sv + 0.5 * h, th + 0.5 * h * k1t, om + 0.5 * h * k1o, Urh[j], wh[j], wth[j])
            k3t, k3o = f_rk(sv + 0.5 * h, th + 0.5 * h * k2t, om + 0.5 * h * k2o, Urh[j], wh[j], wth[j])
            k4t, k4o = f_rk(sv + h, th + h * k3t, om + h * k3o, Ur[j + 1], w[j + 1], wt[j + 1])
            Th[j + 1] = th + h / 6.0 * (k1t + 2 * k2t + 2 * k3t + k4t)
            Om[j + 1] = om + h / 6.0 * (k1o + 2 * k2o + 2 * k3o + k4o)
    return Th, Om


class BQ(BQLogPolar):
    def __init__(self, lam, s_min=-120.0, s_max=100.0, hs=0.025, Nb=32, a=None, s_start=-20.0, sign=-1.0, s_sw=4.0):
        super().__init__(lam, s_min=s_min, s_max=s_max, hs=hs, Nb=Nb, a=a, s_start=s_start, sign=sign)
        self.s_sw = s_sw
        self.jsw = int(np.searchsorted(self.s, s_sw))

    def _vel(self, X, Xs):
        Ur = X @ self.Db.T
        w = -(2.0 * X + Xs)
        return Ur, w, self._wtan(w)

    def velocity(self, P):
        """all velocity fields needed by the hybrid march, from X = e^{(a−2)s}Ψ̃ with finite differences in s"""
        X = self.ea2[:, None] * P
        h = self.hs
        Xs = _grid_deriv(X, h)
        out = list(self._vel(X, Xs))
        for c in (0.5, C1, C2):
            V, D = _offset_values(X, h, c)
            out += list(self._vel(V, D))
        return out            # Ur, w, wt, (half) Urh, wh, wth, (c1) Ur1, w1, wt1, (c2) Ur2, w2, wt2

    def march(self, P, return_all=False, A_fixed=None):
        lam = self.lam
        F = self.velocity(P)
        Ur = F[0]
        A = self.strain(Ur) if A_fixed is None else A_fixed
        m = (lam - 1.0) / (1.0 + lam - A)
        C = m / (A - 1.0)
        s0 = self.s[self.i0]
        Th0 = np.full(self.Nb + 1, self.sign * np.exp(m * s0))
        Om0 = np.full(self.Nb + 1, self.sign * C * np.exp((m - 1) * s0))
        pad = lambda Z: np.vstack([Z, Z[-1:]])            # offset arrays have Ns−1 rows
        Th, Om = _march_hybrid(self.s, self.hs, Th0, Om0, F[0], F[1], F[2],
                               pad(F[3]), pad(F[4]), pad(F[5]), pad(F[6]), pad(F[7]), pad(F[8]),
                               pad(F[9]), pad(F[10]), pad(F[11]),
                               self.Db, self.cb, self.sb, lam, m, self.i0, self.jsw)
        sl = self.s[:self.i0]
        Th[:self.i0] = self.sign * np.exp(m * sl)[:, None]
        Om[:self.i0] = self.sign * C * np.exp((m - 1) * sl)[:, None]
        cpow = np.where(self.cb > 0, self.cb, 0.0)
        Omega = Om * cpow[None, :] ** (m - 1)
        out = dict(A=A, m=m, Th=Th, Om=Om, Omega=Omega, vrmin=float(((1 + lam) + Ur).min()))
        if return_all:
            out.update(Ur=Ur, w=F[1])
        return out
