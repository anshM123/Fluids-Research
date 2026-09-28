"""CCF self-similar profiles on a GEOMETRICALLY graded grid around a thin internal sonic layer.

The tanh map of ccf_mapped.py refines only a region of half-width ~σ ε^{3/2} around η_d, which for very thin
layers (ε ≲ 1e-2) covers only the core of the Lorentzian-like layer 1/den, den ≈ δ(1+(η-η_p)²/w²): the flanks
are then resolved by ~2-3 points per layer width, the discrete layer "pins" to grid points and Newton stalls.

Here  dη/ds = cosh u / (cosh u + k),  u = (s-s_d)/B,  k = h_s/h_c - 1  (h_c = spacing at the layer centre):
spacing h_c at the centre, spacing ≈ h_s |η-η_d| / B for h_c B/h_s ≲ |η-η_d| ≲ B (geometric grading:
a fixed number h_s/B of points per e-fold of distance from the layer), spacing h_s far away.
Closed form:  η = η_d + B [u - (2k/√(k²-1)) artanh(√((k-1)/(k+1)) tanh(u/2))].
Everything else (Sidi-Israeli alternating-point Hilbert rule, 8th-order cumulative quadrature, integral-form
equations, far-field normalisation) is inherited from CCFMapped.
"""
import numpy as np
from scipy.optimize import brentq
from ccf_mapped import CCFMapped, Kc


class CCFSinh(CCFMapped):
    def __init__(self, L1=30.0, L2=120.0, hs=0.03, eta_d=-1.0, hc=1e-4, B=1.0, c=0.7):
        self.c, self.hs, self.eta_d, self.hc, self.B = c, hs, eta_d, hc, B
        k = max(hs / hc - 1.0, 1.0 + 1e-9)
        r = np.sqrt((k - 1) / (k + 1))
        q = 2 * k / np.sqrt(k * k - 1)
        g = lambda s: eta_d + B * ((s - eta_d) / B - q * np.arctanh(r * np.tanh((s - eta_d) / (2 * B))))
        smin = brentq(lambda s: g(s) + L1, -L1 - 20 * B * np.log(2 * k) - 10, eta_d)
        smax = brentq(lambda s: g(s) - L2, eta_d, L2 + 20 * B * np.log(2 * k) + 10)
        N = int(np.ceil((smax - smin) / hs))
        N += N % 2
        self.N = N
        # place a grid point exactly at s = eta_d (layer centre)
        j0 = int(round((eta_d - smin) / hs))
        self.s = eta_d + hs * (np.arange(N) - j0)
        u = (hs / B) * (np.arange(N) - j0)
        # Local coordinate x = η − η_d, built by summing exact (8-point Gauss–Legendre) integrals of
        # dx/du = B cosh u/(cosh u + k) over each grid interval, outward from the centre.  The closed form
        # u − q artanh(r tanh(u/2)) loses ~log10(k) digits (both near the centre and where 1 − r tanh ≈ 1/k);
        # with k up to 1e9 that corrupted kernel arguments at the 1e−8 level.  Summation keeps every
        # difference x_i − x_j (the only thing the kernel sees) to ~1e−15 relative accuracy.
        gpf = lambda v: np.cosh(v) / (np.cosh(v) + k)
        xg, wg = np.polynomial.legendre.leggauss(8)
        m = max(j0, N - 1 - j0)
        a = (hs / B) * np.arange(m)
        mid, half = a + 0.5 * hs / B, 0.5 * hs / B
        inc = B * half * (gpf(mid[:, None] + half * xg[None, :]) @ wg)
        X = np.concatenate([[0.0], np.cumsum(inc)])
        idx = np.arange(N) - j0
        self.x = np.sign(idx) * X[np.abs(idx)]
        self.eta = eta_d + self.x
        self.s = eta_d + B * u
        self.gp = gpf(u)
        self.eps = 1.0 / (1.0 + k)          # centre spacing / hs (diagnostic compatibility)
        self.E = np.exp((c - 1) * self.eta)
        self.i0 = int(np.argmin(np.abs(self.eta)))
        self.iR = int(np.argmin(np.abs(self.eta - 50.0)))
        self.etaR = self.eta[self.iR]
        I = np.arange(N)
        D = self.x[:, None] - self.x[None, :]
        mask = ((I[:, None] - I[None, :]) % 2) == 1
        H = np.where(mask, Kc(np.where(mask, D, 1.0), c), 0.0)
        self.Hm = 2 * hs * H * self.gp[None, :]
        del D, H, mask
        nodes = np.arange(-3, 5)
        rhs = np.array([1.0 / (m + 1) for m in range(8)])
        self.w8 = np.linalg.solve(np.vander(nodes, 8, increasing=True).T, rhs)
        self.wl = [np.linalg.solve(np.vander(np.arange(-s, 8 - s), 8, increasing=True).T, rhs) for s in range(8)]


def make_sinh_grid(etad, w, hs=0.03, pts=24, B=1.0, c=0.7, L1=30.0, L2=120.0):
    """grid with ~pts points per layer width w at the layer centre etad"""
    return CCFSinh(L1, L2, hs=hs, eta_d=etad, hc=min(w / pts, hs / 2.0), B=B, c=c)


def spacing_profile(M, etad, w):
    """local spacing (in units of w) at distances 0, w, 3w, 10w from the layer centre"""
    out = []
    for dist in (0.0, 1.0, 3.0, 10.0):
        j = int(np.argmin(np.abs(M.eta - (etad + dist * w))))
        out.append(M.hs * M.gp[j] / w)
    return out
