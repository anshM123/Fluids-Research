"""
ps2d.py — minimal, verified 2D pseudo-spectral Navier–Stokes solver (vorticity form) on a periodic box.

    ω_t + u·∇ω = ν Δω − α ω − ν_h (−Δ)^p ω + F(x,y),      u = (ψ_y, −ψ_x),  Δψ = −ω

Domain: [0, 2π Lx) × [0, 2π Ly). Integrating-factor RK4 in time, 2/3-rule dealiasing.
Shared tool for several programs (2D turbulence probes, Kolmogorov flow, long-time Euler).
"""
import numpy as np
import scipy.fft as sfft

WORKERS = 4


class PS2D:
    def __init__(self, N, Ny=None, Lx=1.0, Ly=1.0, nu=1e-3, alpha=0.0, nuh=0.0, p=4,
                 forcing=None, dealias=True, workers=WORKERS):
        self.Nx = N
        self.Ny = Ny if Ny is not None else N
        self.Lx, self.Ly = Lx, Ly
        self.nu, self.alpha, self.nuh, self.p = nu, alpha, nuh, p
        self.workers = workers
        kx = sfft.fftfreq(self.Nx, 1.0 / self.Nx) / Lx
        ky = sfft.rfftfreq(self.Ny, 1.0 / self.Ny) / Ly
        self.KX, self.KY = np.meshgrid(kx, ky, indexing="ij")
        self.K2 = self.KX**2 + self.KY**2
        self.K2inv = np.where(self.K2 > 0, 1.0 / np.where(self.K2 > 0, self.K2, 1.0), 0.0)
        if dealias:
            kxmax = (self.Nx // 2) / Lx
            kymax = (self.Ny // 2) / Ly
            self.mask = ((np.abs(self.KX) < 2.0 / 3.0 * kxmax) & (np.abs(self.KY) < 2.0 / 3.0 * kymax)).astype(float)
        else:
            self.mask = np.ones_like(self.K2)
        self.L = -nu * self.K2 - alpha - nuh * self.K2**p      # linear operator (diagonal)
        x = np.arange(self.Nx) * 2 * np.pi * Lx / self.Nx
        y = np.arange(self.Ny) * 2 * np.pi * Ly / self.Ny
        self.X, self.Y = np.meshgrid(x, y, indexing="ij")
        self.Fh = None
        if forcing is not None:
            self.Fh = self.fft(forcing(self.X, self.Y)) * self.mask

    # --- transforms -----------------------------------------------------------------------------
    def fft(self, f):
        return sfft.rfft2(f, workers=self.workers)

    def ifft(self, fh):
        return sfft.irfft2(fh, s=(self.Nx, self.Ny), workers=self.workers)

    # --- physics --------------------------------------------------------------------------------
    def velocity_hat(self, wh):
        psih = wh * self.K2inv
        return 1j * self.KY * psih, -1j * self.KX * psih

    def nonlinear(self, wh):
        """returns −(u·∇ω)^ (dealiased) + forcing"""
        uh, vh = self.velocity_hat(wh)
        u, v = self.ifft(uh), self.ifft(vh)
        wx, wy = self.ifft(1j * self.KX * wh), self.ifft(1j * self.KY * wh)
        N = -self.fft(u * wx + v * wy) * self.mask
        if self.Fh is not None:
            N = N + self.Fh
        return N

    def step(self, wh, dt):
        """integrating-factor RK4"""
        E = np.exp(self.L * dt / 2)
        E2 = E * E
        k1 = self.nonlinear(wh)
        a = E * (wh + dt / 2 * k1)
        k2 = self.nonlinear(a)
        b = E * wh + dt / 2 * k2
        k3 = self.nonlinear(b)
        c = E2 * wh + dt * E * k3
        k4 = self.nonlinear(c)
        return E2 * wh + dt / 6 * (E2 * k1 + 2 * E * (k2 + k3) + k4)

    # --- diagnostics (per unit area, Parseval with rfft weights) -----------------------------------
    def _sum(self, q):
        w = np.full(q.shape, 2.0)
        w[:, 0] = 1.0
        if self.Ny % 2 == 0:
            w[:, -1] = 1.0
        return np.sum(w * q) / (self.Nx * self.Ny) ** 2

    def energy(self, wh):
        return 0.5 * self._sum(np.abs(wh) ** 2 * self.K2inv)

    def enstrophy(self, wh):
        return 0.5 * self._sum(np.abs(wh) ** 2)

    def palinstrophy(self, wh):
        return 0.5 * self._sum(self.K2 * np.abs(wh) ** 2)

    def cfl_dt(self, wh, cfl=0.5):
        uh, vh = self.velocity_hat(wh)
        umax = np.max(np.abs(self.ifft(uh))) + np.max(np.abs(self.ifft(vh))) + 1e-12
        dx = 2 * np.pi * min(self.Lx / self.Nx, self.Ly / self.Ny)
        return cfl * dx / umax


def _selftest():
    # Taylor–Green-type exact decaying solution ω = 2 sin x sin y e^{-2νt}
    s = PS2D(64, nu=0.05)
    w0 = 2 * np.sin(s.X) * np.sin(s.Y)
    wh = s.fft(w0)
    dt, T = 0.01, 2.0
    for _ in range(int(T / dt)):
        wh = s.step(wh, dt)
    err = np.max(np.abs(s.ifft(wh) - w0 * np.exp(-2 * 0.05 * T)))
    print("TG decay max error:", err)
    # Kolmogorov laminar state: F = ν n^3 ... check steady state ω_L = -(n) cos(n y)... use F such that ω_L = n cos(ny)
    n, nu = 4, 0.2
    s = PS2D(64, nu=nu, forcing=lambda X, Y: nu * n**2 * n * np.cos(n * Y))
    wh = s.fft(0 * s.X)
    for _ in range(4000):
        wh = s.step(wh, 0.01)
    print("Kolmogorov laminar error:", np.max(np.abs(s.ifft(wh) - n * np.cos(n * s.Y))))


if __name__ == "__main__":
    _selftest()
