"""Linear stability (instability order) of CCF self-similar profiles on the mapped grid.
μ δ = λ δ − d̄ δ_η − (Hδ/ξ) Θ̄_η ;  for Re μ>0 the smooth solution is δ = T_μ δ with
(T_μ δ)(η) = −∫_{−∞}^{η} exp((λ−μ)∫_{η'}^{η} ds/d̄) (Hδ/ξ)(η') Θ̄_η(η')/d̄(η') dη'.
Eigenvalues μ ⇔ 1 ∈ spec(T_μ). Scan real μ, track eigenvalues ν(μ) of T_μ, locate ν=1 crossings (secant)."""
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs
from ccf_mapped import CCFMapped
from numba import njit


@njit(cache=True)
def _recur(E, inc):
    n = inc.shape[0] + 1
    J = np.zeros(n, dtype=np.complex128)
    for i in range(1, n):
        J[i] = E[i - 1] * J[i - 1] + inc[i - 1]
    return J


class MappedStability:
    def __init__(self, M, phi, lam):
        self.M, self.lam = M, lam
        self.Theta = np.exp(phi + M.c * M.eta)
        self.d = 1 + lam + M.G(phi)
        self.w = lam * self.Theta / self.d**2
        self.A = M.cumint(1.0 / self.d)
        self.dA = np.diff(self.A)
        self.deta = np.diff(M.eta)
        self.pre = np.exp(-M.c * M.eta)
        self.post = np.exp((M.c - 1) * M.eta)

    def T(self, delta, mu):
        x = self.pre * delta
        Hx = self.M.Hm @ x.real + 1j * (self.M.Hm @ x.imag)
        f = self.post * Hx * self.w
        E = np.exp((self.lam - mu) * self.dA)
        inc = 0.5 * self.deta * (E * f[:-1] + f[1:])
        return -_recur(E.astype(np.complex128), inc.astype(np.complex128))

    def spectrum(self, mu, k=8):
        n = self.M.N
        op = LinearOperator((n, n), matvec=lambda v: self.T(v, mu), dtype=complex)
        vals = eigs(op, k=k, which='LM', tol=1e-8, ncv=40, maxiter=500, return_eigenvectors=False)
        return vals

    def crossings(self, mus, k=8):
        """return list of (mu, sorted real eigenvalues) for a real-μ scan"""
        out = []
        for mu in mus:
            v = self.spectrum(mu, k)
            vr = np.sort(v[np.abs(v.imag) < 1e-7 * np.maximum(1, np.abs(v))].real)[::-1]
            out.append((mu, vr))
        return out
