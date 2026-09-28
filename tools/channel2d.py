"""
channel2d.py — 2D wall-bounded incompressible Navier–Stokes (plane Poiseuille / plane Couette)
Fourier (x, periodic) × Legendre–Galerkin (y ∈ [−1,1]) streamfunction formulation.

    Δψ_t = (1/Re) Δ²ψ + (u ω_x + v ω_y),   ω = −Δψ,  u = ψ_y, v = −ψ_x,
    ψ = Ψ_base(y) + ψ̃,   ψ̃(±1) = ψ̃'(±1) = 0   (Shen 1994 basis; constant flux for Poiseuille)

Base flows (steady exact solutions): 'poiseuille' U = 1 − y² (Re on centreline velocity & half-gap,
linear threshold Re_c = 5772.22), 'couette' U = y (Re on wall speed & half-gap).
Time stepping: Spalart–Moser–Rogers low-storage RK3 (explicit NL) + Crank–Nicolson (viscous).
x-dealiasing: 3/2 rule; y-products: over-integrated Gauss–Legendre quadrature.
"""
import numpy as np
import numpy.polynomial.legendre as npleg
import scipy.fft as sfft

RK_G = (8 / 15, 5 / 12, 3 / 4)
RK_Z = (0.0, -17 / 60, -5 / 12)
RK_A = (4 / 15, 1 / 15, 1 / 6)     # alpha_k = beta_k (Crank–Nicolson split)


def shen_basis(M, y):
    """phi_j(y) and derivatives 0..3 at points y; phi_j = L_j - 2(2j+5)/(2j+7) L_{j+2} + (2j+3)/(2j+7) L_{j+4}."""
    out = []
    for d in range(4):
        P = np.zeros((len(y), M))
        for j in range(M):
            c = np.zeros(j + 5)
            c[j] = 1.0
            c[j + 2] = -2 * (2 * j + 5) / (2 * j + 7)
            c[j + 4] = (2 * j + 3) / (2 * j + 7)
            if d:
                c = npleg.legder(c, d)
            P[:, j] = npleg.legval(y, c)
        out.append(P)
    return out


class Channel2D:
    def __init__(self, Nx=128, M=80, Lx=2 * np.pi, Re=5000.0, base="poiseuille", dealias=1.5, workers=1):
        self.Nx, self.M, self.Lx, self.Re = Nx, M, Lx, Re
        self.workers = workers
        self.nk = Nx // 2 + 1
        self.k = 2 * np.pi / Lx * np.arange(self.nk)
        kmax_keep = Nx // 2          # keep modes |k|< Nx/2 (drop Nyquist)
        self.keep = np.arange(self.nk) < kmax_keep
        self.Nxp = int(np.ceil(Nx * dealias / 2) * 2)        # padded physical grid in x
        # quadrature for Galerkin matrices and products
        nq = int(np.ceil(1.5 * (M + 4))) + 2
        self.yq, self.wq = npleg.leggauss(nq)
        self.P0, self.P1, self.P2, self.P3 = shen_basis(M, self.yq)
        W = self.wq[:, None]
        self.C = self.P0.T @ (W * self.P0)
        self.B = self.P1.T @ (W * self.P1)
        self.A = self.P2.T @ (W * self.P2)
        self.PT = (self.P0 * W).T          # projection: F_i = sum_q w_q phi_i(y_q) NL(y_q)
        # base flow evaluated at quadrature points
        y = self.yq
        if base == "poiseuille":
            self.U, self.Uy, self.Uyy = 1 - y**2, -2 * y, -2 * np.ones_like(y)
        elif base == "couette":
            self.U, self.Uy, self.Uyy = y.copy(), np.ones_like(y), np.zeros_like(y)
        else:
            raise ValueError(base)
        self.base = base
        self.dt = None

    # ---------------------------------------------------------------------------------------------
    def set_dt(self, dt):
        """precompute per-k, per-substage inverse matrices for (Mk − β dt Lk) and explicit (Mk + α dt Lk)"""
        self.dt = dt
        M = self.M
        Mk = -(self.B[None] + (self.k**2)[:, None, None] * self.C[None])
        Lk = (self.A[None] + 2 * (self.k**2)[:, None, None] * self.B[None]
              + (self.k**4)[:, None, None] * self.C[None]) / self.Re
        self.Minv = np.linalg.inv(Mk)          # for diagnostics / projections
        self.Mk, self.Lk = Mk, Lk
        self.Linv, self.Rexp = [], []
        for s in range(3):
            lhs = Mk - RK_A[s] * dt * Lk
            self.Linv.append(np.linalg.inv(lhs))
            self.Rexp.append(Mk + RK_A[s] * dt * Lk)

    def _to_phys(self, fk):
        """fk: (nq, nk) complex Fourier coefficients (unnormalised, numpy convention for Nx) → (nq, Nxp) real"""
        pad = np.zeros((fk.shape[0], self.Nxp // 2 + 1), dtype=complex)
        pad[:, : self.nk] = fk * self.keep[None, :]
        return sfft.irfft(pad, n=self.Nxp, axis=1, workers=self.workers) * (self.Nxp / self.Nx)

    def _to_spec(self, f):
        fk = sfft.rfft(f, axis=1, workers=self.workers)[:, : self.nk] * (self.Nx / self.Nxp)
        return fk * self.keep[None, :]

    def fields(self, a):
        """a: (M, nk) coefficients of ψ̃ → physical u, v, ω, ω_x, ω_y at (yq, x_pad) incl. base flow"""
        k = self.k[None, :]
        psi = self.P0 @ a
        psi1 = self.P1 @ a
        psi2 = self.P2 @ a
        psi3 = self.P3 @ a
        uk = psi1
        vk = -1j * k * psi
        wk = -(psi2 - k**2 * psi)
        wyk = -(psi3 - k**2 * psi1)
        u = self._to_phys(uk) + self.U[:, None]
        v = self._to_phys(vk)
        wx = self._to_phys(1j * k * wk)
        wy = self._to_phys(wyk) - self.Uyy[:, None]      # base vorticity −U', its y-derivative −U''
        return u, v, wx, wy

    def rhs_nl(self, a):
        u, v, wx, wy = self.fields(a)
        NL = self._to_spec(u * wx + v * wy)     # (nq, nk)
        return self.PT @ NL                     # (M, nk)

    def step(self, a):
        Nprev = None
        for s in range(3):
            N = self.rhs_nl(a)
            rhs = np.einsum("kij,jk->ik", self.Rexp[s], a) + self.dt * RK_G[s] * N
            if s > 0:
                rhs = rhs + self.dt * RK_Z[s] * Nprev
            a = np.einsum("kij,jk->ik", self.Linv[s], rhs)
            a[:, ~self.keep] = 0
            a[:, 0] = a[:, 0].real            # k=0 mode real
            Nprev = N
        return a

    # --- diagnostics ------------------------------------------------------------------------------
    def perturbation_energy(self, a):
        """(1/(2·Lx·2)) ∫∫ |u'|² dx dy per unit area — using Fourier/Parseval; includes mean-flow distortion (k=0)"""
        k = self.k[None, :]
        uk = self.P1 @ a
        vk = -1j * k * (self.P0 @ a)
        e = (np.abs(uk) ** 2 + np.abs(vk) ** 2) / self.Nx**2
        wts = np.full(self.nk, 2.0)
        wts[0] = 1.0
        return 0.5 * np.sum(self.wq[:, None] * e * wts[None, :]) / 2.0

    def energy_by_k(self, a):
        k = self.k[None, :]
        uk = self.P1 @ a
        vk = -1j * k * (self.P0 @ a)
        e = (np.abs(uk) ** 2 + np.abs(vk) ** 2) / self.Nx**2
        wts = np.full(self.nk, 2.0)
        wts[0] = 1.0
        return 0.5 * np.sum(self.wq[:, None] * e, axis=0) * wts / 2.0

    def wall_shear(self, a):
        """mean wall shear (du/dy at y=-1 and +1, x-averaged) incl. base"""
        ym = np.array([-1.0, 1.0])
        P2 = shen_basis(self.M, ym)[2]
        dU = {"poiseuille": np.array([2.0, -2.0]), "couette": np.array([1.0, 1.0])}[self.base]
        return (P2 @ a[:, 0]).real / self.Nx + dU

    def random_perturbation(self, amp, kmax=None, seed=0):
        rng = np.random.default_rng(seed)
        a = (rng.standard_normal((self.M, self.nk)) + 1j * rng.standard_normal((self.M, self.nk)))
        decay = np.exp(-np.arange(self.M) / 6.0)[:, None]
        a = a * decay
        if kmax is not None:
            a[:, self.k > kmax] = 0
        a[:, 0] = 0
        a[:, ~self.keep] = 0
        e = self.perturbation_energy(a)
        return a * np.sqrt(amp / max(e, 1e-300))


def orr_sommerfeld_check(Re=10000.0, alpha=1.0, M=80):
    """Linear growth rate from the Galerkin operator vs Orr–Sommerfeld reference (Orszag 1971:
    Re=10000, α=1 → c = 0.23752649 + 0.00373967 i)."""
    ch = Channel2D(Nx=8, M=M, Lx=2 * np.pi / alpha, Re=Re)
    k = alpha
    Mk = -(ch.B + k**2 * ch.C)
    Lk = (ch.A + 2 * k**2 * ch.B + k**4 * ch.C) / Re
    # linearised advection about U: NL_lin = U ω'_x + v' Ω_y with Ω = −U', Ω_y = −U''
    P0, P1, P2 = ch.P0, ch.P1, ch.P2
    wk = -(P2 - k**2 * P0)                 # ω' in terms of coefficients
    vk = -1j * k * P0
    NLmat = ch.U[:, None] * (1j * k * wk) + vk * (-ch.Uyy[:, None])
    Nlin = ch.PT @ NLmat
    ev = np.linalg.eigvals(np.linalg.solve(Mk, Lk + Nlin))
    # eigen σ; phase speed c = i σ / k
    c = 1j * ev / k
    i = np.argmax(ev.real)
    return c[i], ev[i]


if __name__ == "__main__":
    c, s = orr_sommerfeld_check()
    print("OS check  c =", c, " (ref 0.23752649+0.00373967j)")
