"""
Linear stability of self-similar Boussinesq profiles via the eigen-condition ν(μ) = 1 for the linearised
"march + Biot–Savart" map T_μ (the method validated on CCF in P02, R011).

Self-similar time τ = −ln(1−t); perturbations (θ', ω', ψ') e^{μτ} of a profile (Θ, Ω, Ψ) satisfy
    (μ + 1 − λ) θ' + V·∇θ' + u'·∇Θ = 0,
    (μ + 1) ω' + V·∇ω' + u'·∇Ω = ∂₁θ',        u' = ∇^⊥ψ',  −Δψ' = ω',     V = (1+λ)y + U.
All characteristics of V leave the stagnation point. For a given velocity perturbation (X' = ψ'/r²), (θ', ω')
are obtained by marching outward from the origin with regular (vanishing) data. The Biot–Savart law then returns
a new X'. This defines a linear map T_μ, and μ is an eigenvalue of the linearised operator iff T_μ has an
eigenvalue ν = 1. The time-translation mode is μ = 1 exactly (trivial); unstable modes are the crossings
ν(μ) = 1 with Re μ > 0, μ ≠ 1.

Log-polar form (D = V_r/r = (1+λ) + X_β, w = V_β/r = −(2X + X_s), base fields unhatted):
    D θ'_s = −(μ+1−λ) θ' − w θ'_β − (U'_r/r) Θ_s − (U'_β/r) Θ_β
    D ω'_s = −(μ+1) ω' − w ω'_β − (U'_r/r) Ω_s − (U'_β/r) Ω_β + e^{−s}(cosβ θ'_s − sinβ θ'_β)
with U'_r/r = X'_β and U'_β/r = −(2X' + X'_s). Numerics are those of bq_solver: 8th-order FD in s for
velocities, Gauss–Legendre (implicit) march for s < s_sw and RK4 beyond, FFT/Chebyshev Biot–Savart.
"""
import numpy as np
from numba import njit
from scipy.sparse.linalg import LinearOperator, eigs
from bq_solver import BQ, _grid_deriv, _offset_values, C1, C2
from bq_newton import full


@njit(cache=True)
def _lin_march(s, h, mu, lam, Db, cb, sb, i0, jsw,
               D0, w0, D1, w1, D2, w2, Dh, wh, Dg1, wg1,
               Ft0, Fo0, Ft1, Fo1, Ft2, Fo2, Fth, Foh, Ftg1, Fog1):
    """linear march of (θ', ω') with base speeds (D, w) and forcings (F_θ, F_ω) given at the grid points (0), the
    Gauss points (1, 2) and the half points (h); the '1'-suffixed grid arrays (Dg1, wg1, Ftg1, Fog1) hold the
    grid values shifted by one (j+1) for the RK4 end stage."""
    Ns, n = D0.shape
    Th = np.zeros((Ns, n)); Om = np.zeros((Ns, n))
    I = np.eye(n)
    a11, a12, a21, a22 = 0.25, 0.25 - np.sqrt(3.0) / 6, 0.25 + np.sqrt(3.0) / 6, 0.25
    c1, c2 = 0.5 - np.sqrt(3.0) / 6, 0.5 + np.sqrt(3.0) / 6
    M = np.zeros((2 * n, 2 * n)); rhs = np.zeros(2 * n)

    def f_rk(sv, th, om, D, w, Ft, Fo):
        thb = Db @ th
        ths = (-(mu + 1.0 - lam) * th - w * thb + Ft) / D
        oms = (-(mu + 1.0) * om - w * (Db @ om) + Fo + np.exp(-sv) * (cb * ths - sb * thb)) / D
        return ths, oms

    for j in range(i0, Ns - 1):
        sv = s[j]; th = Th[j]; om = Om[j]
        if j < jsw:
            P1 = np.empty((n, n)); P2 = np.empty((n, n)); Q1 = np.empty((n, n)); Q2 = np.empty((n, n))
            for i in range(n):
                for k in range(n):
                    P1[i, k] = -w1[j, i] * Db[i, k] / D1[j, i]
                    P2[i, k] = -w2[j, i] * Db[i, k] / D2[j, i]
                Q1[i, :] = P1[i, :]; Q2[i, :] = P2[i, :]
                P1[i, i] += -(mu + 1.0 - lam) / D1[j, i]
                P2[i, i] += -(mu + 1.0 - lam) / D2[j, i]
                Q1[i, i] += -(mu + 1.0) / D1[j, i]
                Q2[i, i] += -(mu + 1.0) / D2[j, i]
            g1 = Ft1[j] / D1[j]; g2 = Ft2[j] / D2[j]
            M[:n, :n] = I - h * a11 * P1; M[:n, n:] = -h * a12 * P1
            M[n:, :n] = -h * a21 * P2;    M[n:, n:] = I - h * a22 * P2
            rhs[:n] = P1 @ th + g1; rhs[n:] = P2 @ th + g2
            K = np.linalg.solve(M, rhs)
            K1 = K[:n]; K2 = K[n:]
            Y1 = th + h * (a11 * K1 + a12 * K2); Y2 = th + h * (a21 * K1 + a22 * K2)
            e1 = np.exp(-(sv + c1 * h)); e2 = np.exp(-(sv + c2 * h))
            q1 = (Fo1[j] + e1 * (cb * K1 - sb * (Db @ Y1))) / D1[j]
            q2 = (Fo2[j] + e2 * (cb * K2 - sb * (Db @ Y2))) / D2[j]
            M[:n, :n] = I - h * a11 * Q1; M[:n, n:] = -h * a12 * Q1
            M[n:, :n] = -h * a21 * Q2;    M[n:, n:] = I - h * a22 * Q2
            rhs[:n] = Q1 @ om + q1; rhs[n:] = Q2 @ om + q2
            L = np.linalg.solve(M, rhs)
            Th[j + 1] = th + 0.5 * h * (K1 + K2)
            Om[j + 1] = om + 0.5 * h * (L[:n] + L[n:])
        else:
            k1t, k1o = f_rk(sv, th, om, D0[j], w0[j], Ft0[j], Fo0[j])
            k2t, k2o = f_rk(sv + 0.5 * h, th + 0.5 * h * k1t, om + 0.5 * h * k1o, Dh[j], wh[j], Fth[j], Foh[j])
            k3t, k3o = f_rk(sv + 0.5 * h, th + 0.5 * h * k2t, om + 0.5 * h * k2o, Dh[j], wh[j], Fth[j], Foh[j])
            k4t, k4o = f_rk(sv + h, th + h * k3t, om + h * k3o, Dg1[j], wg1[j], Ftg1[j], Fog1[j])
            Th[j + 1] = th + h / 6.0 * (k1t + 2 * k2t + 2 * k3t + k4t)
            Om[j + 1] = om + h / 6.0 * (k1o + 2 * k2o + 2 * k3o + k4o)
    return Th, Om


class BQStab:
    def __init__(self, lam, Y, Nb=32, hs=0.025, s_sw=12.0):
        self.B = B = BQ(lam, Nb=Nb, hs=hs, s_sw=s_sw)
        self.lam = lam
        X = full(B, Y)
        r = B.march(X / B.ea2[:, None], return_all=True)
        self.m, self.A = r['m'], r['A']
        cp = np.where(B.cb > 0, B.cb, 0.0)
        Th = r['Th'] * cp[None, :] ** self.m
        Om = r['Omega']
        self.Th, self.Om = Th, Om
        h = B.hs
        Db = B.Db
        # base speeds at grid / half / Gauss points
        V = B.velocity(X / B.ea2[:, None])          # Ur, w, wt, Urh, wh, wth, Ur1, w1, wt1, Ur2, w2, wt2
        pad = lambda Z: np.vstack([Z, Z[-1:]])
        self.D0 = (1 + lam) + V[0]; self.w0 = V[1]
        self.Dh = pad((1 + lam) + V[3]); self.wh = pad(V[4])
        self.D1 = pad((1 + lam) + V[6]); self.w1 = pad(V[7])
        self.D2 = pad((1 + lam) + V[9]); self.w2 = pad(V[10])
        self.Dg1 = np.vstack([self.D0[1:], self.D0[-1:]]); self.wg1 = np.vstack([self.w0[1:], self.w0[-1:]])
        # base gradients (unhatted) at grid / half / Gauss points
        def grads(F):
            Fs = _grid_deriv(F, h); Fb = F @ Db.T
            out = {'0': (Fs, Fb)}
            for key, c in (('h', 0.5), ('1', C1), ('2', C2)):
                Vv, Dd = _offset_values(F, h, c)
                out[key] = (pad(Dd), pad(Vv @ Db.T))
            return out
        self.gT = grads(Th); self.gO = grads(Om)
        self.i0, self.jsw = B.i0, B.jsw
        self.n_unk = (B.Ns - B.i0) * (B.Nb + 1)

    def _pert_vel(self, Xp):
        """U'_r/r and U'_β/r at grid, half and Gauss points from X' (full grid)"""
        B = self.B; h = B.hs; Db = B.Db
        pad = lambda Z: np.vstack([Z, Z[-1:]])
        Xs = _grid_deriv(Xp, h)
        out = {'0': (Xp @ Db.T, -(2 * Xp + Xs))}
        for key, c in (('h', 0.5), ('1', C1), ('2', C2)):
            Vv, Dd = _offset_values(Xp, h, c)
            out[key] = (pad(Vv @ Db.T), pad(-(2 * Vv + Dd)))
        return out

    def T(self, Yp, mu):
        """T_μ: X'[i0:] ↦ X'_new[i0:] (constant extension of X' below the start, as for the base state)"""
        B = self.B
        Xp = full(B, Yp.reshape(B.Ns - B.i0, B.Nb + 1))
        u = self._pert_vel(Xp)
        F = {}
        for key in ('0', 'h', '1', '2'):
            Ur, w = u[key]
            Ts, Tb = self.gT[key]; Os, Ob = self.gO[key]
            F[key] = (-(Ur * Ts + w * Tb), -(Ur * Os + w * Ob))
        Ftg1 = np.vstack([F['0'][0][1:], F['0'][0][-1:]]); Fog1 = np.vstack([F['0'][1][1:], F['0'][1][-1:]])
        th, om = _lin_march(B.s, B.hs, mu, self.lam, B.Db, B.cb, B.sb, self.i0, self.jsw,
                            self.D0, self.w0, self.D1, self.w1, self.D2, self.w2, self.Dh, self.wh,
                            self.Dg1, self.wg1,
                            F['0'][0], F['0'][1], F['1'][0], F['1'][1], F['2'][0], F['2'][1],
                            F['h'][0], F['h'][1], Ftg1, Fog1)
        Pn = B.poisson(om)
        return (B.ea2[:, None] * Pn)[B.i0:].ravel()

    def spectrum(self, mu, k=8, tol=1e-10, v0=None):
        op = LinearOperator((self.n_unk, self.n_unk), matvec=lambda v: self.T(v, mu), dtype=float)
        vals, vecs = eigs(op, k=k, which='LM', tol=tol, v0=v0, maxiter=2000)
        o = np.argsort(-np.abs(vals))
        return vals[o], vecs[:, o]

    def time_translation_mode(self):
        """X' of the trivial μ = 1 mode: θ' = (λ−1)Θ − (1+λ)Θ_s, ω' = −Ω − (1+λ)Ω_s → ψ' via Biot–Savart"""
        B = self.B
        om = -self.Om - (1 + self.lam) * self.gO['0'][0]
        om[:B.i0] = 0.0
        return (B.ea2[:, None] * B.poisson(om))[B.i0:].ravel()
