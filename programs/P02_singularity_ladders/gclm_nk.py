"""
Self-similar blow-up profiles of the generalised Constantin–Lax–Majda (gCLM / Okamoto–Sakajo–Wunsch) model
        ω_t + a u ω_x = u_x ω ,   u_x = Hω ,
ω = (T−t)^{-1} Ω(ξ), ξ = x/(T−t)^{c_l},  u = (T−t)^{c_l−1} U(ξ), U' = HΩ, Ω odd (Ω<0 for ξ>0):
        (c_l ξ + a U) Ω' = (HΩ − 1) Ω .
Log form (η = ln ξ):  d ln|Ω|/dη = B(η) = (HΩ − 1) / (c_l + a U/ξ).
Local exponent at ξ→0:  |Ω| ~ ξ^q,  q = (h0 − 1)/(c_l + a h0),  h0 = HΩ(0) = −(2/π)∫_0^∞ Ω/y dy > 0.
Smooth odd profiles: q ∈ {1,3,5,…}. Far field |Ω| ~ A ξ^{−1/c_l}.
Odd-function Hilbert transform in log variables: HΩ(e^η) = (1/π) p.v.∫ Ω(e^μ) e^{−(η−μ)}/sinh(η−μ) dμ;
with |Ω| = e^{cη}Ψ (−1/c_l < c < 0): H|Ω|(e^η) = e^{cη} (K_{1+c} * Ψ), multiplier −i tanh(π(k − i(1+c))/2).
U/ξ = e^{−η} ∫_{−∞}^{η} HΩ(e^μ) e^{μ} dμ.
Unknowns: φ = ln Ψ (N values) and c_l (when q is prescribed). Normalisation: far-field amplitude.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres


class GCLM:
    def __init__(self, a=0.0, L1=40.0, L2=80.0, N=16384, c=-0.4):
        self.a, self.N, self.c = a, N, c
        self.h = (L1 + L2) / N
        self.eta = -L1 + self.h * np.arange(N)
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.h)
        self.mK = -1j * np.tanh(np.pi * (k - 1j * (1 + c)) / 2)
        self.i0 = int(np.argmin(np.abs(self.eta)))
        self.iR = int(np.argmin(np.abs(self.eta - 30.0)))
        self.etaR = self.eta[self.iR]
        nodes = np.arange(-3, 5)
        rhs = np.array([1.0 / (m + 1) for m in range(8)])
        self.w8 = np.linalg.solve(np.vander(nodes, 8, increasing=True).T, rhs)
        self.wl = [np.linalg.solve(np.vander(np.arange(-s, 8 - s), 8, increasing=True).T, rhs) for s in range(8)]

    def cumint(self, X, anchor=None):
        n, h, w = self.N, self.h, self.w8
        mid = np.zeros(n - 1)
        for j, off in enumerate(range(-3, 5)):
            mid[3:n - 4] += w[j] * X[3 + off:n - 4 + off]
        for i in list(range(0, 3)) + list(range(n - 4, n - 1)):
            s = i if i < 3 else 7 - (n - 1 - i)
            mid[i] = self.wl[s] @ X[i - s:i - s + 8]
        C = np.zeros(n)
        C[1:] = np.cumsum(mid * h)
        return C if anchor is None else C - C[anchor]

    def fields(self, phi):
        """returns HΩ (signed, Ω=−e^{cη}Ψ) and U/ξ"""
        Psi = np.exp(phi)
        Habs = np.exp(self.c * self.eta) * np.real(np.fft.ifft(self.mK * np.fft.fft(Psi)))
        HO = -Habs
        Uxi = np.exp(-self.eta) * self.cumint(HO * np.exp(self.eta))      # ∫_{−∞}^{η} (from left end)
        return HO, Uxi

    def normval(self, phi, cl):
        return phi[self.iR] + (self.c + 1.0 / cl) * self.etaR

    def F(self, phi, cl, target):
        HO, Uxi = self.fields(phi)
        D = cl + self.a * Uxi
        B = (HO - 1.0) / D
        F = phi - phi[self.i0] - self.cumint(B - self.c, anchor=self.i0)
        F[self.i0] = self.normval(phi, cl) - target
        return F, D, HO, Uxi

    def h0(self, phi):
        # h0 = HΩ(0) = −(2/π)∫_0^∞ Ω/y dy = (2/π)∫ Ψ e^{cη} dη
        return (2 / np.pi) * self.h * np.sum(np.exp(phi + self.c * self.eta))

    def q_of(self, phi, cl):
        h0 = self.h0(phi)
        return (h0 - 1) / (cl + self.a * h0)

    # --- Newton–Krylov with finite-difference Jacobian-vector products (robust, simple) ---------------
    def _jv(self, v, phi, cl, target, F0, eps=1e-7):
        nv = np.linalg.norm(v)
        if nv == 0:
            return np.zeros_like(v)
        e = eps * max(1.0, np.linalg.norm(phi)) / nv
        F1, *_ = self.F(phi + e * v, cl, target)
        return (F1 - F0) / e

    def solve_q(self, phi, cl, q, target=None, tol=1e-10, maxit=40, verbose=False):
        N = self.N
        if target is None:
            target = self.normval(phi, cl)
        for it in range(maxit):
            F0, D, HO, Uxi = self.F(phi, cl, target)
            g = self.q_of(phi, cl) - q
            nrm = max(np.max(np.abs(F0)), abs(g))
            if verbose:
                print(f"   it {it}: |F|={nrm:.2e} c_l={cl:.12f} minD={D.min():.4f}", flush=True)
            if nrm < tol:
                return phi, cl, True
            dcl = 1e-7
            Fc, *_ = self.F(phi, cl + dcl, target)
            Fcl = (Fc - F0) / dcl
            gcl = (self.q_of(phi, cl + dcl) - self.q_of(phi, cl)) / dcl

            def mv(x):
                v, s = x[:N], x[N]
                jv = self._jv(v, phi, cl, target, F0) + s * Fcl
                ev = 1e-7 * max(1.0, np.linalg.norm(phi)) / max(np.linalg.norm(v), 1e-300)
                gq = (self.q_of(phi + ev * v, cl) - self.q_of(phi, cl)) / ev if np.linalg.norm(v) > 0 else 0.0
                return np.concatenate([jv, [gq + gcl * s]])
            A = LinearOperator((N + 1, N + 1), matvec=mv, dtype=float)
            d, info = gmres(A, -np.concatenate([F0, [g]]), rtol=1e-10, atol=0.0, restart=100, maxiter=20)
            t = 1.0
            while t > 1e-3:
                pn, cn = phi + t * d[:N], cl + t * d[N]
                Fn, Dn, *_ = self.F(pn, cn, target)
                gn = self.q_of(pn, cn) - q
                if Dn.min() > 0 and max(np.max(np.abs(Fn)), abs(gn)) < (1 - 0.2 * t) * nrm:
                    break
                t *= 0.5
            phi, cl = pn, cn
        return phi, cl, False

    def guess(self, cl):
        # |Ω| ≈ ξ/(1+ξ²)^{(1+1/cl)/2}
        th = self.eta - 0.5 * (1 + 1 / cl) * np.logaddexp(0, 2 * self.eta)
        return th - self.c * self.eta


def solve_fixed(G, phi, cl, target=None, tol=1e-10, maxit=40, verbose=False):
    """Newton–Krylov at fixed c_l (normalisation by far-field amplitude)."""
    N = G.N
    if target is None:
        target = G.normval(phi, cl)
    for it in range(maxit):
        F0, D, HO, Uxi = G.F(phi, cl, target)
        nrm = np.max(np.abs(F0))
        if verbose:
            print(f"   fx it {it}: |F|={nrm:.2e} minD={D.min():.4f}", flush=True)
        if nrm < tol:
            return phi, True
        A = LinearOperator((N, N), matvec=lambda v: G._jv(v, phi, cl, target, F0), dtype=float)
        d, info = gmres(A, -F0, rtol=1e-10, atol=0.0, restart=100, maxiter=20)
        t = 1.0
        while t > 1e-3:
            Fn, Dn, *_ = G.F(phi + t * d, cl, target)
            if Dn.min() > 0 and np.max(np.abs(Fn)) < (1 - 0.2 * t) * nrm:
                break
            t *= 0.5
        phi = phi + t * d
    return phi, False
