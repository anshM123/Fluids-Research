"""
Newton / pseudo-arclength continuation for CCF self-similar profiles in integral form.

Unknown φ(η) = ln Ψ(η), Ψ = e^{-cη} Θ(e^η).  Profile ODE (exact rewriting of
−λΘ + ((1+λ)ξ + HΘ) Θ' = 0):
        φ'(η) = b(η) = λ / (1 + λ + G(η)) − c,     G = HΘ/ξ = e^{(c−1)η} (K_c * e^φ)
Integral form:  F_i = φ_i − φ_{i0} − (Q b)_i = 0  (Q = cumulative trapezoid from η_{i0}); row i0 → φ_{i0} = φ0.
This square system has no spurious periodic constraint (index-1 issue of the differential form).
Local exponent at ξ→0:  p = λ/(1+λ+h1),  h1 = G(−∞) = −(2/π)∫ Ψ e^{(c−1)η} dη.
Smooth profiles: p ∈ {2,4,6,…}.
"""
import numpy as np
import scipy.linalg as sla


class CCFNewton:
    def __init__(self, L1=30.0, L2=150.0, N=4096, c=0.75):
        self.N, self.c = N, c
        self.h = (L1 + L2) / N
        self.eta = -L1 + self.h * np.arange(N)
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.h)
        self.mK = -1j * np.tanh(np.pi * (k - 1j * c) / 2)
        self.E = np.exp((c - 1) * self.eta)
        self.i0 = int(np.argmin(np.abs(self.eta)))
        # cumulative trapezoid matrix relative to i0
        N, h, i0 = self.N, self.h, self.i0
        I = np.eye(N)
        self.Km = np.real(np.fft.ifft(self.mK[:, None] * np.fft.fft(I, axis=0), axis=0))

    def cumtrap(self, X):
        """cumulative integral from index i0 along axis 0, 8th-order (interval integrals of the
        degree-7 Lagrange interpolant on nodes i-3..i+4; stencils shifted at the ends)."""
        h, i0 = self.h, self.i0
        n = X.shape[0]
        if not hasattr(self, '_w8'):
            nodes = np.arange(-3, 5)
            V = np.vander(nodes, 8, increasing=True).T          # V[m, j] = nodes_j^m
            rhs = np.array([1.0 / (m + 1) for m in range(8)])   # ∫_0^1 t^m dt
            self._w8 = np.linalg.solve(V, rhs)
            self._wl = [np.linalg.solve(np.vander(np.arange(-s, 8 - s), 8, increasing=True).T, rhs) for s in range(8)]
        w = self._w8
        mid = np.zeros((n - 1,) + X.shape[1:])
        # interior intervals i = 3 .. n-5
        for j, off in enumerate(range(-3, 5)):
            mid[3:n - 4] += w[j] * X[3 + off:n - 4 + off]
        # boundary intervals with shifted stencils (nodes i-s .. i+7-s)
        for i in list(range(0, 3)) + list(range(n - 4, n - 1)):
            s = i if i < 3 else 8 - (n - 1 - i) - 1
            s = min(max(s, 0), 7)
            ww = self._wl[s]
            mid[i] = np.tensordot(ww, X[i - s:i - s + 8], axes=(0, 0))
        mid *= h
        C = np.zeros(X.shape)
        C[1:] = np.cumsum(mid, axis=0)
        return C - C[i0]

    def G(self, phi):
        return self.E * np.real(np.fft.ifft(self.mK * np.fft.fft(np.exp(phi))))

    def h1(self, phi):
        return -(2 / np.pi) * self.h * np.sum(np.exp(phi) * self.E)

    def F(self, phi, lam, phi0):
        G = self.G(phi)
        den = 1 + lam + G
        b = lam / den - self.c
        F = phi - phi[self.i0] - self.cumtrap(b)
        F[self.i0] = phi[self.i0] - phi0
        return F, den

    def jac(self, phi, lam, den):
        N = self.N
        Psi = np.exp(phi)
        dbdG = -lam / den**2
        dG = (self.E[:, None] * self.Km) * Psi[None, :]          # dG_i/dφ_j
        J = np.eye(N) - self.cumtrap(dbdG[:, None] * dG)
        J[:, self.i0] -= 1.0
        J[self.i0, :] = 0.0
        J[self.i0, self.i0] = 1.0
        # dF/dλ
        dbdl = 1 / den - lam / den**2
        Fl = -self.cumtrap(dbdl)
        Fl[self.i0] = 0.0
        return J, Fl

    def picard(self, phi, lam, relax=0.3, tol=1e-7, maxit=3000):
        i0 = self.i0
        for it in range(maxit):
            den = 1 + lam + self.G(phi)
            if den.min() <= 0:
                return phi, False
            new = self.cumtrap(lam / den - self.c)
            new = new - new[i0] + phi[i0]
            diff = np.max(np.abs(new - phi))
            phi = (1 - relax) * phi + relax * new
            if diff < tol:
                return phi, True
        return phi, False

    def solve_fixed(self, phi, lam, phi0=None, tol=1e-12, maxit=30, verbose=False):
        if phi0 is None:
            phi0 = phi[self.i0]
        phi, okp = self.picard(phi, lam)
        if verbose:
            print("   picard ok", okp)
        for it in range(maxit):
            F, den = self.F(phi, lam, phi0)
            nrm = np.max(np.abs(F))
            if verbose:
                print(f"   it {it} |F|={nrm:.2e} min(den)={den.min():.4f}")
            if nrm < tol:
                return phi, True
            if den.min() <= 0:
                return phi, False
            J, _ = self.jac(phi, lam, den)
            d = -sla.solve(J, F)
            t = 1.0
            while t > 1e-4:
                Fn, dn = self.F(phi + t * d, lam, phi0)
                if dn.min() > 0 and np.max(np.abs(Fn)) < (1 - 0.2 * t) * nrm:
                    break
                t *= 0.5
            phi = phi + t * d
        return phi, False

    def p_of(self, phi, lam):
        return lam / (1 + lam + self.h1(phi))

    def solve_p(self, phi, lam, p, phi0=None, tol=1e-12, maxit=40, verbose=False):
        """unknowns (φ, λ); extra equation p(φ,λ) = p (smoothness)."""
        N = self.N
        if phi0 is None:
            phi0 = phi[self.i0]
        for it in range(maxit):
            F, den = self.F(phi, lam, phi0)
            h1 = self.h1(phi)
            g = h1 - (lam / p - 1 - lam)
            nrm = max(np.max(np.abs(F)), abs(g))
            if verbose:
                print(f"   it {it} |F|={nrm:.2e} lam={lam:.15f}")
            if nrm < tol:
                return phi, lam, True
            J, Fl = self.jac(phi, lam, den)
            A = np.zeros((N + 1, N + 1))
            A[:N, :N] = J
            A[:N, N] = Fl
            A[N, :N] = -(2 / np.pi) * self.h * np.exp(phi) * self.E
            A[N, N] = -(1 / p - 1)
            d = sla.solve(A, -np.concatenate([F, [g]]))
            t = 1.0
            while t > 1e-4:
                pn, ln = phi + t * d[:N], lam + t * d[N]
                Fn, dn = self.F(pn, ln, phi0)
                gn = self.h1(pn) - (ln / p - 1 - ln)
                if dn.min() > 0 and max(np.max(np.abs(Fn)), abs(gn)) < (1 - 0.2 * t) * nrm:
                    break
                t *= 0.5
            phi, lam = pn, ln
        return phi, lam, False

    def guess(self, lam, p0):
        beta = lam / (1 + lam)
        Theta = np.exp(p0 * self.eta - 0.5 * (p0 - beta) * np.logaddexp(0, 2 * self.eta))
        Psi = np.exp(-self.c * self.eta) * Theta
        Psi *= (lam / p0 - 1 - lam) / (-(2 / np.pi) * self.h * np.sum(Psi * self.E))
        return np.log(Psi)
