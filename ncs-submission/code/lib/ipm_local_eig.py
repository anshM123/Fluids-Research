"""
Local (WKB) dispersion relation for rapidly oscillating perturbations in the stalled boundary layer of a 2D IPM
self-similar profile (the IPM analogue of bq_local_eig.py).

Linearised steady operator about the profile:  −λr + V·∇r + u'·∇R = 0,  u' = (ψ'_y, −ψ'_x),  −Δψ' = −∂₁r.
Near a wall point (x, 0) where the self-similar radial speed D = V₁/x = ε D̂ is small, with Y = y/(εx),
r = r̃(Y) e^{iϕ}, ψ' = εx ψ̃(Y) e^{iϕ}, ϕ = ε⁻¹∫κ ds (s = ln x), at leading order in ε:
    (a) [iκ(D̂ + ĉY) + μ Y ∂_Y] r̃ = −G ∂_Yψ̃ + iκ R_y ψ̃
    (c) −∂_Y²ψ̃ + κ²ψ̃ = −iκ r̃,        ψ̃(0) = 0, no growing mode e^{κY}
with ĉ = ∂_yV₁ = −Ω_b = ∂_xR (IPM: Ω = −∂₁R), μ = ∂_yV₂ = 2(1+λ) − D − xD', G = ∂_xR, R_y = ∂_yR at the wall.
Compared with Boussinesq, the vorticity transport equation (b) is replaced by the slaving ω̃ = −iκ r̃.
Dropped at leading order: −λr (O(ε)), the Y² terms of V₁, slow x-variation, the log-polar metric.
"""
import numpy as np
from scipy.linalg import eig
from scipy.integrate import solve_ivp
from bq_local_eig import cheb


class LocalEigIPM:
    def __init__(self, N=120, Ymax=80.0):
        D, t = cheb(N)
        t = t[::-1]; D = D[::-1, ::-1]
        self.Y = Ymax * (1 + t) ** 2 / 4
        dYdt = Ymax * (1 + t) / 2; dYdt[0] = 1.0
        self.DY = np.diag(1 / dYdt) @ D
        self.DY[0, :] = (D @ D)[0, :] / (Ymax / 2)
        self.DY2 = self.DY @ self.DY
        self.YDY = np.diag((1 + t) / 2) @ D
        self.N = N

    def matrices(self, Dh, ch, mu, G, Ry):
        n = self.N + 1; I = np.eye(n); Z = np.zeros((n, n)); Y = np.diag(self.Y)
        M0 = np.block([[mu * self.YDY, G * self.DY], [Z, -self.DY2]])
        M1 = np.block([[1j * (Dh * I + ch * Y), -1j * Ry * I], [1j * I, Z]])
        M2 = np.block([[Z, Z], [Z, I]])
        for row in (n, 2 * n - 1):
            M0[row, :] = 0; M1[row, :] = 0; M2[row, :] = 0; M0[row, row] = 1.0
        return M0, M1, M2

    def roots(self, Dh, ch, mu, G, Ry):
        M0, M1, M2 = self.matrices(Dh, ch, mu, G, Ry)
        n = M0.shape[0]; I = np.eye(n); Z = np.zeros((n, n))
        w, V = eig(np.block([[Z, I], [-M0, -M1]]), np.block([[I, Z], [Z, M2]]))
        ok = np.isfinite(w)
        return w[ok], V[:n, ok]


def shoot(kap, Dh, ch, mu, G, Ry, Y0=1e-8, Y1=None, rtol=1e-10):
    """regular solution with ψ(0)=0, ψ'(0)=1; returns the coefficient of the growing mode e^{κY} at Y1"""
    if Y1 is None:
        Y1 = 25.0 / max(kap.real, 0.05)
    r0 = -G / (1j * kap * Dh)                  # iκD̂ r(0) = −G ψ'(0)

    def f(Y, u):
        r, ps, dps = u
        dr = (-G * dps + 1j * kap * Ry * ps - 1j * kap * (Dh + ch * Y) * r) / (mu * Y)
        return [dr, dps, kap ** 2 * ps + 1j * kap * r]
    u0 = np.array([r0, Y0, 1.0], dtype=complex)
    sol = solve_ivp(f, (Y0, Y1), u0, method='DOP853', rtol=rtol, atol=1e-14)
    r, ps, dps = sol.y[:, -1]
    return (dps + kap * ps) * np.exp(-kap * Y1) / 2


def local_root(Dh, ch, mu, G, Ry, guess, tol=1e-10, maxit=60):
    k0 = complex(guess); k1 = k0 * (1 + 1e-3) + 1e-4
    f0 = shoot(k0, Dh, ch, mu, G, Ry); f1 = shoot(k1, Dh, ch, mu, G, Ry)
    for it in range(maxit):
        if f1 == f0:
            return k1, False
        k2 = k1 - f1 * (k1 - k0) / (f1 - f0)
        if abs(k2 - k1) < tol * max(1, abs(k1)):
            return k2, True
        k0, f0, k1 = k1, f1, k2
        f1 = shoot(k1, Dh, ch, mu, G, Ry)
    return k1, False
