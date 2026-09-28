"""
CCF self-similar profiles on an ADAPTIVE (mapped) grid, for profiles with thin internal sonic layers.

Map: η = g(s) = s − (1−ε) σ tanh((s−η_d)/σ), s uniform (spacing h_s)  → local spacing ε h_s at η_d.
Hilbert convolution (K_c * Ψ)(η_i) = ∫ K_c(η_i−g(s)) Ψ(g(s)) g'(s) ds, K_c(x) = e^{−cx}/(π sinh x) (p.v.),
discretised with the Sidi–Israeli alternating-point trapezoidal rule (spectrally accurate for p.v. kernels):
        (K*Ψ)_i ≈ 2 h_s Σ_{j−i odd} K_c(η_i−η_j) g'_j Ψ_j .
Cumulative integrals ∫ b dη = ∫ b g' ds with 8th-order quadrature in s.
Same integral-form equations, normalisation (far-field amplitude) and solvers as ccf_nk.py.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
from scipy.optimize import brentq


def Kc(x, c):
    out = np.empty_like(x)
    pos = x > 0
    neg = x < 0
    xp, xn = x[pos], x[neg]
    out[pos] = 2 * np.exp(-(1 + c) * xp) / (-np.expm1(-2 * xp))
    out[neg] = -2 * np.exp((1 - c) * xn) / (-np.expm1(2 * xn))
    out[~(pos | neg)] = 0.0
    return out / np.pi


class CCFMapped:
    def __init__(self, L1=30.0, L2=120.0, hs=0.03, eta_d=-1.0, eps=1.0, sigma=3.0, c=0.7, dense=True):
        self.c, self.hs, self.eta_d, self.eps, self.sigma = c, hs, eta_d, eps, sigma
        g = lambda s: s - (1 - eps) * sigma * np.tanh((s - eta_d) / sigma)
        smin = brentq(lambda s: g(s) + L1, -L1 - 10 * sigma - 10, eta_d)
        smax = brentq(lambda s: g(s) - L2, eta_d, L2 + 10 * sigma + 10)
        N = int(np.ceil((smax - smin) / hs))
        N += N % 2
        self.N = N
        self.s = smin + hs * np.arange(N)
        self.eta = g(self.s)
        self.gp = 1 - (1 - eps) / np.cosh((self.s - eta_d) / sigma) ** 2
        self.E = np.exp((c - 1) * self.eta)
        self.i0 = int(np.argmin(np.abs(self.eta)))
        self.iR = int(np.argmin(np.abs(self.eta - 50.0)))
        self.etaR = self.eta[self.iR]
        # dense alternating-point Hilbert matrix
        I = np.arange(N)
        D = self.eta[:, None] - self.eta[None, :]
        mask = ((I[:, None] - I[None, :]) % 2) == 1
        H = np.where(mask, Kc(np.where(mask, D, 1.0), c), 0.0)
        self.Hm = 2 * hs * H * self.gp[None, :]
        del D, H, mask
        nodes = np.arange(-3, 5)
        rhs = np.array([1.0 / (m + 1) for m in range(8)])
        self.w8 = np.linalg.solve(np.vander(nodes, 8, increasing=True).T, rhs)
        self.wl = [np.linalg.solve(np.vander(np.arange(-s, 8 - s), 8, increasing=True).T, rhs) for s in range(8)]

    def cumint(self, X):
        """∫_{η_{i0}}^{η_i} X dη  = ∫ X g' ds (8th order in s)"""
        Y = X * self.gp
        n, h, w = self.N, self.hs, self.w8
        mid = np.zeros(n - 1)
        for j, off in enumerate(range(-3, 5)):
            mid[3:n - 4] += w[j] * Y[3 + off:n - 4 + off]
        for i in list(range(0, 3)) + list(range(n - 4, n - 1)):
            s = i if i < 3 else 7 - (n - 1 - i)
            mid[i] = self.wl[s] @ Y[i - s:i - s + 8]
        C = np.zeros(n)
        C[1:] = np.cumsum(mid * h)
        return C - C[self.i0]

    def conv(self, f):
        return self.Hm @ f

    def G(self, phi):
        return self.E * self.conv(np.exp(phi))

    def h1(self, phi):
        # -(2/π) ∫ Ψ e^{(c-1)η} dη  (trapezoid in s is spectral for decaying integrands)
        return -(2 / np.pi) * self.hs * np.sum(np.exp(phi) * self.E * self.gp)

    def dh1(self, phi):
        return -(2 / np.pi) * self.hs * np.exp(phi) * self.E * self.gp

    def p_of(self, phi, lam):
        return lam / (1 + lam + self.h1(phi))

    def normval(self, phi, lam):
        return phi[self.iR] + (self.c - lam / (1 + lam)) * self.etaR

    def F(self, phi, lam, phi0):
        den = 1 + lam + self.G(phi)
        F = phi - phi[self.i0] - self.cumint(lam / den - self.c)
        F[self.i0] = self.normval(phi, lam) - phi0
        return F, den

    def Jv(self, v, phi, lam, den):
        dG = self.E * self.conv(np.exp(phi) * v)
        out = v - v[self.i0] - self.cumint(-lam / den**2 * dG)
        out[self.i0] = v[self.iR]
        return out

    def Fl(self, lam, den):
        out = -self.cumint(1 / den - lam / den**2)
        out[self.i0] = -self.etaR / (1 + lam) ** 2
        return out

    def _gmres(self, matvec, rhs, n, tol=1e-12):
        A = LinearOperator((n, n), matvec=matvec, dtype=float)
        x, info = gmres(A, rhs, rtol=tol, atol=0.0, restart=120, maxiter=40)
        return x, info

    def solve_fixed(self, phi, lam, phi0=None, tol=1e-12, maxit=30, verbose=False):
        if phi0 is None:
            phi0 = self.normval(phi, lam)
        for it in range(maxit):
            F, den = self.F(phi, lam, phi0)
            nrm = np.max(np.abs(F))
            if verbose:
                print(f"   it {it}: |F|={nrm:.2e} minden={den.min():.5f}", flush=True)
            if nrm < tol:
                return phi, True
            if den.min() <= 0:
                return phi, False
            d, info = self._gmres(lambda v: self.Jv(v, phi, lam, den), -F, self.N)
            t = 1.0
            while t > 1e-3:
                Fn, dn = self.F(phi + t * d, lam, phi0)
                if dn.min() > 0 and np.max(np.abs(Fn)) < (1 - 0.2 * t) * nrm:
                    break
                t *= 0.5
            phi = phi + t * d
        return phi, False

    def solve_p(self, phi, lam, p, phi0=None, tol=1e-12, maxit=30, verbose=False):
        N = self.N
        if phi0 is None:
            phi0 = self.normval(phi, lam)
        for it in range(maxit):
            F, den = self.F(phi, lam, phi0)
            g = self.h1(phi) - (lam / p - 1 - lam)
            nrm = max(np.max(np.abs(F)), abs(g))
            if verbose:
                print(f"   p-it {it}: |F|={nrm:.2e} lam={lam:.15f} minden={den.min():.6f}", flush=True)
            if nrm < tol:
                return phi, lam, True
            Fl = self.Fl(lam, den)
            dh = self.dh1(phi)
            sc = 1.0 / max(np.linalg.norm(dh), abs(1 / p - 1))

            def mv(x):
                v, s = x[:N], x[N]
                return np.concatenate([self.Jv(v, phi, lam, den) + s * Fl, [sc * (dh @ v - (1 / p - 1) * s)]])
            d, info = self._gmres(mv, -np.concatenate([F, [sc * g]]), N + 1)
            t = 1.0
            while t > 1e-3:
                pn, ln = phi + t * d[:N], lam + t * d[N]
                Fn, dn = self.F(pn, ln, phi0)
                gn = self.h1(pn) - (ln / p - 1 - ln)
                if dn.min() > 0 and max(np.max(np.abs(Fn)), abs(gn)) < (1 - 0.2 * t) * nrm:
                    break
                t *= 0.5
            phi, lam = pn, ln
        return phi, lam, False

    def from_other(self, other_eta, other_phi):
        return np.interp(self.eta, other_eta, other_phi)
