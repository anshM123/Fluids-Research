"""Self-similar blow-up profiles of 2D IPM (incompressible porous media) with boundary, in the log-polar framework of
the Boussinesq solver (bq_solver.BQ) — the third model for the blind test of the one-phase theory.

Physical problem (upper half-plane, no penetration on y2 = 0), Darcy's law with gravity:
    ρ_t + u·∇ρ = 0,   u = −∇p − ρ e₂,  ∇·u = 0   ⇒   ω = curl u = −∂_{x1}ρ,   u = (ψ_y, −ψ_x),  −Δψ = ω.
So IPM is the Boussinesq system with the vorticity equation Ω + V·∇Ω = ∂₁Θ replaced by its "stalled" form
Ω = −∂₁R (vorticity slaved instantly to the density gradient).
Self-similar ansatz:  ρ = (1−t)^{λ} R(y),  u = (1−t)^{λ} U(y),  ω = (1−t)^{−1} Ω(y),  y = x/(1−t)^{1+λ}:
    V·∇R = λ R,   V = (1+λ)y + U,   Ω = −∂_{y1}R,   U = (Ψ_2, −Ψ_1),  −ΔΨ = Ω.
The transported scalar has exponent λ (Boussinesq: λ − 1); the ladder accumulates where it vanishes, λ → 0.
Symmetry: R even in y1, Ω, Ψ odd. Local structure at the stagnation point (strain A):
    R ≈ c|y1|^m,  m = λ/(1+λ−A);  smooth ⇔ m = 2 ⇔ λ = 2A − 2.
Hat variables as in Boussinesq: R = (cos β)^m R̂,  Ω = (cos β)^{m−1} Ω̂, with
    R̂_s = [λR̂ + m (w tanβ) R̂ − w R̂_β]/(V_r/r),   Ω̂ = −e^{−s}[cos²β R̂_s − sinβ cosβ R̂_β + m sin²β R̂].
The march (implicit Gauss–Legendre near the origin, RK4 beyond) carries R̂ only; Ω̂ follows algebraically."""
import numpy as np
from numba import njit
from bq_solver import BQ, C1, C2


@njit(cache=True)
def _march_R(s, h, R0, Ur, w, wt, Urh, wh, wth, Ur1, w1, wt1, Ur2, w2, wt2, Db, lam_s, lam, m, i0, jsw):
    Ns, n = Ur.shape
    R = np.zeros((Ns, n)); R[i0] = R0
    I = np.eye(n)
    a11, a12, a21, a22 = 0.25, 0.25 - np.sqrt(3.0) / 6, 0.25 + np.sqrt(3.0) / 6, 0.25
    M = np.zeros((2 * n, 2 * n)); rhs = np.zeros(2 * n)

    def f(r, ur, ww, wtt):
        vr = (1.0 + lam) + ur
        return (lam_s * r + m * wtt * r - ww * (Db @ r)) / vr

    for j in range(i0, Ns - 1):
        r = R[j]
        if j < jsw:
            vr1 = (1.0 + lam) + Ur1[j]; vr2 = (1.0 + lam) + Ur2[j]
            P1 = np.empty((n, n)); P2 = np.empty((n, n))
            for i in range(n):
                for k in range(n):
                    P1[i, k] = -w1[j, i] * Db[i, k] / vr1[i]
                    P2[i, k] = -w2[j, i] * Db[i, k] / vr2[i]
                P1[i, i] += (lam_s + m * wt1[j, i]) / vr1[i]
                P2[i, i] += (lam_s + m * wt2[j, i]) / vr2[i]
            M[:n, :n] = I - h * a11 * P1; M[:n, n:] = -h * a12 * P1
            M[n:, :n] = -h * a21 * P2;    M[n:, n:] = I - h * a22 * P2
            rhs[:n] = P1 @ r; rhs[n:] = P2 @ r
            K = np.linalg.solve(M, rhs)
            R[j + 1] = r + 0.5 * h * (K[:n] + K[n:])
        else:
            k1 = f(r, Ur[j], w[j], wt[j])
            k2 = f(r + 0.5 * h * k1, Urh[j], wh[j], wth[j])
            k3 = f(r + 0.5 * h * k2, Urh[j], wh[j], wth[j])
            k4 = f(r + h * k3, Ur[j + 1], w[j + 1], wt[j + 1])
            R[j + 1] = r + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    return R


class IPM(BQ):
    """same grid, Biot–Savart and Newton interface as bq_solver.BQ; sign = +1 gives compressive strain A > 0"""
    def __init__(self, lam, s_min=-120.0, s_max=100.0, hs=0.025, Nb=32, a=None, s_start=-20.0, sign=1.0, s_sw=4.0):
        super().__init__(lam, s_min=s_min, s_max=s_max, hs=hs, Nb=Nb, a=a, s_start=s_start, sign=sign, s_sw=s_sw)

    def march(self, P, return_all=False, A_fixed=None):
        lam = self.lam
        F = self.velocity(P)
        Ur, w, wt = F[0], F[1], F[2]
        A = self.strain(Ur) if A_fixed is None else A_fixed
        m = lam / (1.0 + lam - A)
        s0 = self.s[self.i0]
        R0 = np.full(self.Nb + 1, self.sign * np.exp(m * s0))
        pad = lambda Z: np.vstack([Z, Z[-1:]])
        R = _march_R(self.s, self.hs, R0, Ur, w, wt, pad(F[3]), pad(F[4]), pad(F[5]), pad(F[6]), pad(F[7]),
                     pad(F[8]), pad(F[9]), pad(F[10]), pad(F[11]), self.Db, lam, lam, m, self.i0, self.jsw)
        sl = self.s[:self.i0]
        R[:self.i0] = self.sign * np.exp(m * sl)[:, None]
        # vorticity from the density gradient (R̂_s from the transport equation itself)
        vr = (1.0 + lam) + Ur
        Rb = R @ self.Db.T
        Rs = (lam * R + m * wt * R - w * Rb) / vr
        cb2 = self.cb ** 2; sbcb = self.sb * self.cb; sb2 = self.sb ** 2
        Om = -np.exp(-self.s)[:, None] * (cb2[None, :] * Rs - sbcb[None, :] * Rb + m * sb2[None, :] * R)
        Om[:self.i0] = -self.sign * m * np.exp((m - 1) * sl)[:, None]
        cpow = np.where(self.cb > 0, self.cb, 0.0)
        Omega = Om * cpow[None, :] ** (m - 1)
        out = dict(A=A, m=m, Th=R, Om=Om, Omega=Omega, vrmin=float(vr.min()))
        if return_all:
            out.update(Ur=Ur, w=w)
        return out


if __name__ == "__main__":
    import sys, time
    from bq_newton import newton, initial_guess
    lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    B = IPM(lam, hs=0.025, Nb=32, s_sw=12.0)
    Y = initial_guess(B, 1 + lam / 2)
    t0 = time.time()
    Y, info, ok = newton(B, Y, tol=1e-10, maxit=30, verbose=True, fd='central', pert=1e-6)
    print(f"IPM λ={lam}: ok={ok} A={info['A']:.10f} m={info['m']:.10f} (smooth ⇔ A = 1 + λ/2 = {1 + lam/2:.6f}) "
          f"min V_r/r = {info['vrmin']:.4f} ({time.time()-t0:.0f}s)")
    np.save(f"ipm_Y_lam{lam:.6f}.npy", Y)
