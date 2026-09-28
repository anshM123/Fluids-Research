"""
Linear stability of CCF self-similar profiles (instability order).

Self-similar time τ = −ln(T−t):  Θ_τ = λΘ − (1+λ)ξΘ_ξ − (HΘ)Θ_ξ.  Linearise at Θ̄ (η = ln ξ):
        μ δ = λ δ − d̄ δ_η − (Hδ/ξ) Θ̄_η ,   d̄ = 1+λ+HΘ̄/ξ > 0 ,  Θ̄_η = λΘ̄/d̄.
For Re μ > 0 the smooth (ξ→0 regular) solution is the particular solution
        δ(η) = −∫_{−∞}^{η} exp((λ−μ)∫_{η'}^{η} ds/d̄) · (Hδ/ξ)(η') Θ̄_η(η')/d̄(η') dη'  =: (T_μ δ)(η)
(the homogeneous solution ξ^{2(1−μ/λ)} is non-smooth and excluded).  Eigenvalues: μ with 1 ∈ spec(T_μ).
Trivial modes: μ = 1 (time translation, δ = λΘ̄Ḡ/d̄); μ = 0 (scaling).  Eigenfunctions vanish like ξ² at 0,
so Hδ/ξ = e^{(c−1)η} K_c*(e^{−cη} δ) (same multiplier as the profile).
T_μ is discretised with a 2nd-order exponential-trapezoid integrator; its dominant spectrum is found by
Arnoldi; the eigen-condition ν(μ)=1 is solved by secant iteration in μ.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigs


class CCFStability:
    def __init__(self, S, phi, lam):
        self.S, self.lam = S, lam
        self.Theta = np.exp(phi + S.c * S.eta)
        self.d = 1 + lam + S.G(phi)
        self.w = lam * self.Theta / self.d**2          # Θ̄_η / d̄
        # A(η) = ∫ ds/d̄  (8th-order cumulative)
        self.A = S.cumint(1.0 / self.d)
        self.dA = np.diff(self.A)

    def Hxi(self, delta):
        S = self.S
        return S.E * np.real(np.fft.ifft(S.mK * np.fft.fft(np.exp(-S.c * S.eta) * delta))) * np.exp(S.c * S.eta) * 0 \
            + np.exp((S.c - 1) * S.eta) * np.real(np.fft.ifft(S.mK * np.fft.fft(np.exp(-S.c * S.eta) * delta)))

    def T(self, delta, mu):
        """δ ↦ T_μ δ  via J' = (λ−μ)J/d̄ + f, J(−∞)=0, δ = −J  (exponential trapezoid)"""
        f = self.Hxi(delta) * self.w
        a = (self.lam - mu) * self.dA                  # (λ−μ) ΔA_i
        E = np.exp(a)
        h = self.S.h
        J = np.zeros(len(f), dtype=complex)
        # J_i = E_i J_{i-1} + h/2 (E_i f_{i-1} + f_i)
        inc = 0.5 * h * (E * f[:-1] + f[1:])
        for i in range(1, len(f)):
            J[i] = E[i - 1] * J[i - 1] + inc[i - 1]
        return -J

    def spectrum(self, mu, k=6):
        n = self.S.N
        op = LinearOperator((n, n), matvec=lambda v: self.T(v, mu), dtype=complex)
        vals, vecs = eigs(op, k=k, which='LM', tol=1e-10, maxiter=5000)
        return vals, vecs
