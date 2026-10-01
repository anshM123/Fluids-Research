"""
Local (WKB) dispersion relation for rapidly oscillating perturbations in the quasi-stagnant boundary region of a
2D Boussinesq self-similar profile near λ → 1.

Near a boundary point (x, 0) where the self-similar radial speed D = V₁/x = ε D̂ is small, perturbations
∝ exp(i∫K dx) with K = κ/(εx) and vertical scale y = εx·Y satisfy, at leading order in ε,
    (a) [iκ(D̂ + ĉY) + μ Y ∂_Y] θ̃ = −G ∂_Yψ̃ + iκ Θ_y ψ̃
    (b) [1 + iκ(D̂ + ĉY) + μ Y ∂_Y] ω̃ = iκ θ̃
    (c) −∂_Y²ψ̃ + κ²ψ̃ = ω̃,     ψ̃(0) = 0, ψ̃ → 0 (Y → ∞),
where ĉ = ∂_y(V₁/x)·(εx)/ε·(1/x)… = −Ω_b (boundary vorticity, sign convention of bq_logpolar: Ω_b < 0),
μ = ∂_yV₂ = 2(1+λ) − D − xD', G = ∂_xΘ(x, 0) and Θ_y = ∂_yΘ(x, 0).
(a),(b) are transport equations in Y with outflow at Y = ∞ (no boundary condition); regularity at Y = 0 is automatic.
The problem is a quadratic eigenvalue problem in κ; the relevant root continues the Hou–Luo root
κ_HL = (−1 + sqrt(1 + 4iΩ))/(2iD̂) (Ω in the HL sign convention).

The WKB phase Φ = (1/ε)∫κ ds (s = ln x) then predicts the crossing spacing Δz = π/(m Re εΦ) along the branch.
"""
import numpy as np
from scipy.linalg import eig


def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1); c[0] = c[-1] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


class LocalEig:
    def __init__(self, N=120, Ymax=80.0, Lmap=6.0):
        # algebraic map Y = Lmap (1+t)/(1−t+2Lmap/Ymax) … simpler: Y = Ymax (1+t)/2 with clustering via t = sin map
        D, t = cheb(N)
        t = t[::-1]; D = D[::-1, ::-1]                       # t from −1 to 1
        # map t ∈ [−1,1] → Y ∈ [0, Ymax] with clustering near 0: Y = Ymax (1+t)² / 4
        self.Y = Ymax * (1 + t) ** 2 / 4
        dYdt = Ymax * (1 + t) / 2
        dYdt[0] = 1.0                                         # placeholder, row 0 replaced below
        self.DY = np.diag(1 / dYdt) @ D
        self.DY[0, :] = (D @ D)[0, :] / (Ymax / 2)            # L'Hôpital at Y = 0: f_Y = f_tt / Y_tt
        self.DY2 = self.DY @ self.DY
        self.N = N
        # at Y = 0 the map is singular: use Y ∂_Y = ((1+t)/2) ∂_t, regular
        self.YDY = np.diag((1 + t) / 2) @ D
        self.t = t

    def matrices(self, Dh, ch, mu, G, Thy):
        """M(κ) = M0 + κ M1 + κ² M2 on v = (θ, ω, ψ)"""
        n = self.N + 1; I = np.eye(n); Z = np.zeros((n, n)); Y = np.diag(self.Y)
        M0 = np.block([[mu * self.YDY, Z, G * self.DY],
                       [Z, I + mu * self.YDY, Z],
                       [Z, -I, -self.DY2]])
        M1 = np.block([[1j * (Dh * I + ch * Y), Z, -1j * Thy * I],
                       [-1j * I, 1j * (Dh * I + ch * Y), Z],
                       [Z, Z, Z]])
        M2 = np.block([[Z, Z, Z], [Z, Z, Z], [Z, Z, I]])
        # boundary conditions for ψ: ψ(0) = 0, ψ(Ymax) = 0 (rows of the Poisson block)
        for row, col in ((2 * n, 2 * n), (3 * n - 1, 3 * n - 1)):
            M0[row, :] = 0; M1[row, :] = 0; M2[row, :] = 0
            M0[row, col] = 1.0
        return M0, M1, M2

    def roots(self, Dh, ch, mu, G, Thy):
        M0, M1, M2 = self.matrices(Dh, ch, mu, G, Thy)
        n = M0.shape[0]; I = np.eye(n); Z = np.zeros((n, n))
        A = np.block([[Z, I], [-M0, -M1]])
        B = np.block([[I, Z], [Z, M2]])
        w, V = eig(A, B)
        ok = np.isfinite(w)
        return w[ok], V[:n, ok]


def hl_root(Dh, Om_hl):
    return (-1 + np.sqrt(1 + 4j * Om_hl)) / (2j * Dh)


if __name__ == "__main__":
    import sys
    L = LocalEig(N=int(sys.argv[1]) if len(sys.argv) > 1 else 100)
    # sample parameters typical of the stagnant region: D̂ = 1, ĉ = 2x, G = −2x, μ = 4 at x = 0.3
    for x in (0.1, 0.2, 0.3):
        w, V = L.roots(1.0, 2 * x, 4.0, -2 * x, 0.0)
        sel = (w.real > 0) & (np.abs(w) < 20)
        ws = w[sel][np.argsort(np.abs(w[sel]))]
        print(f"x={x}: HL-like root {hl_root(1.0, 2*x):.4f}; smallest QEP roots (Re>0):", np.round(ws[:6], 4))


# ---------------------------------------------------------------------------------------------------------------
# Shooting formulation (avoids the spurious spectrum of the truncated collocation problem)
from scipy.integrate import solve_ivp


def shoot(kap, Dh, ch, mu, G, Thy, Y0=1e-8, Y1=None, rtol=1e-10):
    """integrate the regular solution with ψ(0)=0, ψ'(0)=1 upward; return the growing-mode coefficient at Y1"""
    if Y1 is None:
        Y1 = 25.0 / max(kap.real, 0.05)
    th0 = -G / (1j * kap * Dh)                  # iκD̂ θ(0) = −G ψ'(0)
    om0 = 1j * kap * th0 / (1 + 1j * kap * Dh)

    def f(Y, u):
        th, om, ps, dps = u
        a = 1j * kap * (Dh + ch * Y)
        dth = (-G * dps + 1j * kap * Thy * ps - a * th) / (mu * Y)
        dom = (1j * kap * th - (1 + a) * om) / (mu * Y)
        return [dth, dom, dps, kap ** 2 * ps - om]
    u0 = np.array([th0, om0, Y0, 1.0], dtype=complex)
    sol = solve_ivp(f, (Y0, Y1), u0, method='DOP853', rtol=rtol, atol=1e-14)
    th, om, ps, dps = sol.y[:, -1]
    return (dps + kap * ps) * np.exp(-kap * Y1) / 2


def local_root(Dh, ch, mu, G, Thy, guess, tol=1e-10, maxit=50):
    """complex secant iteration on the growing-mode coefficient"""
    k0 = complex(guess); k1 = k0 * (1 + 1e-3) + 1e-4
    f0 = shoot(k0, Dh, ch, mu, G, Thy); f1 = shoot(k1, Dh, ch, mu, G, Thy)
    for it in range(maxit):
        k2 = k1 - f1 * (k1 - k0) / (f1 - f0)
        if abs(k2 - k1) < tol * max(1, abs(k1)):
            return k2, True
        k0, f0, k1 = k1, f1, k2
        f1 = shoot(k1, Dh, ch, mu, G, Thy)
    return k1, False
