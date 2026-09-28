"""
Fractional CCF family  θ_t + (HΛ^s θ) θ_x = 0,  0 ≤ s < 1   (s=0: CCF;  s→1: HΛ = −∂_x, θ_t = θ_x², i.e.
Burgers for u = −2θ_x, whose smooth self-similar blow-up profiles form the explicit infinite ladder
λ_i = 1 + 1/i, i = 1, 2, …  in the θ-scaling below).

Self-similar ansatz θ = (T−t)^λ Θ(ξ), ξ = x/(T−t)^b, b = (1+λ)/(1+s):   −λΘ + (bξ + V)Θ' = 0,  V = HΛ^sΘ.
Log variables η = ln ξ, Θ = e^{cη}Ψ, φ = ln Ψ:   φ' = λ/den − c,  den = b + V/ξ,
V/ξ = e^{(c−1−s)η} (m_s ∗ Ψ)(η) with the Mellin multiplier (Θ = |x|^a ↦ V/ξ = m_s(a)|x|^{a−1−s}, a = c+ik)
        m_s(a) = 2^s Γ((1+a)/2) Γ((1+s−a)/2) / [Γ(−a/2) Γ(1+(a−s)/2)]          (m_0 = −tan(πa/2): CCF,  m_1 = −a).
Admissible weight: β < c < 1+s with far-field exponent β = λ/b.
Local exponent at the origin p = λ/(b + h_s), h_s = lim_{ξ→0} V/ξ = κ_s ∫_0^∞ Θ ξ^{−2−s} dξ,
        κ_s = (1+s)/(Γ(−s) sin(πs/2))   (κ_0 = −2/π).
Smooth (analytic) profiles: p = 2.
"""
import numpy as np
from scipy.special import loggamma, gamma
from scipy.sparse.linalg import LinearOperator, gmres
from ccf_nk import CCFNK


def mellin_multiplier(s, a):
    if s == 0.0:
        return -np.tan(np.pi * a / 2)
    return np.exp(s * np.log(2.0) + loggamma((1 + a) / 2) + loggamma((1 + s - a) / 2)
                  - loggamma(-a / 2) - loggamma(1 + (a - s) / 2))


class FCCF(CCFNK):
    def __init__(self, s=0.0, L1=30.0, L2=120.0, N=16384, c=0.7):
        super().__init__(L1, L2, N, c)
        self.s = s
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.h)
        self.mK = mellin_multiplier(s, c + 1j * k)
        self.E = np.exp((c - 1 - s) * self.eta)
        self.kap = -2 / np.pi if s == 0.0 else (1 + s) / (gamma(-s) * np.sin(np.pi * s / 2))

    def b(self, lam):
        return (1 + lam) / (1 + self.s)

    def h1(self, phi):
        return self.kap * self.h * np.sum(np.exp(phi) * self.E)

    def dh1(self, phi):
        return self.kap * self.h * np.exp(phi) * self.E

    def p_of(self, phi, lam):
        return lam / (self.b(lam) + self.h1(phi))

    def normval(self, phi, lam):
        return phi[self.iR] + (self.c - lam / self.b(lam)) * self.etaR

    def F(self, phi, lam, phi0):
        den = self.b(lam) + self.G(phi)
        F = phi - phi[self.i0] - self.cumint(lam / den - self.c)
        F[self.i0] = self.normval(phi, lam) - phi0
        return F, den

    def Fl(self, lam, den):
        out = -self.cumint(1 / den - lam / (1 + self.s) / den**2)
        out[self.i0] = -(1 + self.s) / (1 + lam) ** 2 * self.etaR
        return out

    def solve_p(self, phi, lam, p, phi0=None, tol=1e-12, maxit=30, verbose=False):
        """smoothness-constrained solve: unknowns (φ, λ), constraint p(φ,λ) = p  ⇔  h_s(φ) = λ/p − b(λ)"""
        N, db = self.N, 1 / (1 + self.s)
        if phi0 is None:
            phi0 = self.normval(phi, lam)
        for it in range(maxit):
            F, den = self.F(phi, lam, phi0)
            if den.min() <= 0:
                return phi, lam, False
            g = self.h1(phi) - (lam / p - self.b(lam))
            nrm = max(np.max(np.abs(F)), abs(g))
            if verbose:
                print(f"   p-it {it}: |F|={nrm:.2e} lam={lam:.15f} minden={den.min():.4f}", flush=True)
            if nrm < tol:
                return phi, lam, True
            Fl = self.Fl(lam, den)
            dh = self.dh1(phi)
            gl = -(1 / p - db)
            sc = 1.0 / max(np.linalg.norm(dh), abs(gl))

            def mv(x):
                v, sl = x[:N], x[N]
                return np.concatenate([self.Jv(v, phi, lam, den) + sl * Fl, [sc * (dh @ v + gl * sl)]])
            A = LinearOperator((N + 1, N + 1), matvec=mv, dtype=float)
            d, info = gmres(A, -np.concatenate([F, [sc * g]]), rtol=1e-11, atol=0.0, restart=80, maxiter=40)
            t = 1.0
            while t > 1e-3:
                pn, ln = phi + t * d[:N], lam + t * d[N]
                Fn, dn = self.F(pn, ln, phi0)
                gn = self.h1(pn) - (ln / p - self.b(ln))
                if dn.min() > 0 and max(np.max(np.abs(Fn)), abs(gn)) < (1 - 0.2 * t) * nrm:
                    break
                t *= 0.5
            phi, lam = pn, ln
        return phi, lam, False


if __name__ == "__main__":
    # validation: s=0 reproduces CCF λ0, λ1; multiplier limits
    a = 0.7 + 1j * np.linspace(-5, 5, 11)
    print("m_0 vs -tan:", np.abs(mellin_multiplier(1e-12, a) + np.tan(np.pi * a / 2)).max())
    print("m_1 vs -a  :", np.abs(mellin_multiplier(1.0 - 1e-12, a) + a).max())
    for f, lam_ref in (("ladder_F16k_up_cross2.npy", 1.180777662899), ("ladder_F16k_up_cross1.npy", 0.60573370)):
        d = np.load(f)
        for s in (0.0, 1e-9):
            S = FCCF(s=s, N=16384)
            phi, lam, ok = S.solve_p(d[1:], d[0], 2.0)
            print(f"{f} s={s}: lam={lam:.12f} (ref {lam_ref}) ok={ok}")
