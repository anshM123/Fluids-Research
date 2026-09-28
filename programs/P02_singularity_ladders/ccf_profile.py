"""
Self-similar blow-up profiles of the Córdoba–Córdoba–Fontelos (CCF) equation
        θ_t + (Hθ) θ_x = 0,      H f(x) = (1/π) p.v.∫ f(y)/(x−y) dy,
θ(x,t) = (T−t)^λ Θ(ξ),  ξ = x/(T−t)^{1+λ}:
        −λ Θ + (1+λ) ξ Θ' + (HΘ) Θ' = 0,   Θ even, Θ(0)=0,  Θ ~ C|ξ|^{β}, β = λ/(1+λ).

Local analysis at ξ→0 (HΘ ≈ h1 ξ):  Θ ~ a ξ^p with p = λ / (1+λ+h1),  h1 = −(2/π)∫_0^∞ Θ/ξ² dξ.
Smooth (analytic, even) profiles require p ∈ {2,4,6,...}.  Hypothesis: the n-th unstable profile has p = 2n+2.

Numerics: η = ln ξ, Θ = e^{cη} Ψ(η) with β < c < 1, so Ψ decays exponentially at both ends.
For even Θ:  HΘ(e^η) = e^{cη} (K_c * Ψ)(η),   K̂_c(k) = −i tanh(π(k − i c)/2)   (exact Fourier multiplier).
Residual (divided by e^{cη}):
        R = −λΨ + (1+λ)(cΨ+Ψ') + G (cΨ+Ψ') + s[(1−c)Ψ − Ψ'],   G = e^{(c−1)η}(K_c*Ψ) = HΘ/ξ
where s multiplies the generator of the scaling symmetry Θ → κΘ(ξ/κ) (s=0 at a solution).
Extra equations: normalisation Ψ(η0)=1; smoothness h1 = λ/p − 1 − λ.  Newton on (Ψ, λ, s).
"""
import numpy as np
import scipy.linalg as sla


class CCFProfile:
    def __init__(self, L1=30.0, L2=150.0, N=4096, c=0.75, Lb=8.0):
        self.L1, self.L2, self.N, self.c = L1, L2, N, c
        self.P = L1 + L2
        self.h = self.P / N
        self.eta = -L1 + self.h * np.arange(N)
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.h)
        self.k = k
        self.mK = -1j * np.tanh(np.pi * (k - 1j * c) / 2)
        self.mD = 1j * k
        if N % 2 == 0:
            self.mD[N // 2] = 0.0
        self.E = np.exp((c - 1) * self.eta)
        self.i0 = np.argmin(np.abs(self.eta - 0.0))
        # bordering direction (must lie outside range of dR/dPsi; the scaling generator does NOT)
        self.bvec = np.exp(-(self.eta - 1.0) ** 2)
        # buffer at right end: pin Psi=0 there (true Psi ~ e^{(beta-c)eta} negligible) to remove the
        # spurious periodic monodromy constraint of the first-order ODE
        self.buf = self.eta > (L2 - Lb)

    def conv(self, f, m):
        return np.real(np.fft.ifft(m * np.fft.fft(f)))

    def circulant(self, m):
        """dense matrix of the Fourier multiplier m (real part), columns = response to unit vectors"""
        I = np.eye(self.N)
        return np.real(np.fft.ifft(m[:, None] * np.fft.fft(I, axis=0), axis=0))

    def h1(self, Psi):
        return -(2 / np.pi) * self.h * np.sum(Psi * self.E)

    def residual(self, Psi, lam, s):
        c = self.c
        D = self.conv(Psi, self.mD)
        G = self.E * self.conv(Psi, self.mK)
        W = c * Psi + D
        R = -lam * Psi + (1 + lam) * W + G * W + s * self.bvec
        R = np.where(self.buf, Psi, R)
        return R, G, W, D

    def solve(self, Psi, lam, p, s=0.0, tol=1e-11, maxit=60, fix_lam=False, verbose=False, v0=None):
        """Newton for smooth profile with prescribed local exponent p (unknown λ), or fixed λ (unknown p)."""
        N, c = self.N, self.c
        if v0 is None:
            v0 = Psi[self.i0]
        Dm = self.circulant(self.mD)
        Km = self.circulant(self.mK)
        for it in range(maxit):
            R, G, W, D = self.residual(Psi, lam, s)
            h1 = self.h1(Psi)
            eq_norm = Psi[self.i0] - v0
            if fix_lam:
                eq_smooth = 0.0
            else:
                eq_smooth = h1 - (lam / p - 1 - lam)
            F = np.concatenate([R, [eq_norm, eq_smooth]])
            nrm = np.max(np.abs(F))
            if verbose:
                print(f"  it {it}: |F|={nrm:.3e}  lam={lam:.12f}  h1={h1:.10f}  p_loc={lam/(1+lam+h1):.8f}")
            if nrm < tol:
                return Psi, lam, s, True
            J = np.zeros((N + 2, N + 2))
            A = (-lam + (1 + lam) * c) * np.eye(N) + (1 + lam) * Dm
            A += G[:, None] * (c * np.eye(N) + Dm)
            A += (W * self.E)[:, None] * Km
            A[self.buf, :] = 0.0
            A[self.buf, self.buf] = 1.0
            J[:N, :N] = A
            J[:N, N] = np.where(self.buf, 0.0, -Psi + W)          # dR/dλ
            J[:N, N + 1] = np.where(self.buf, 0.0, self.bvec)     # dR/ds
            J[N, self.i0] = 1.0
            if fix_lam:
                J[N + 1, N] = 1.0                     # dλ = 0
            else:
                J[N + 1, :N] = -(2 / np.pi) * self.h * self.E
                J[N + 1, N] = -(1 / p - 1)
            dx = sla.solve(J, -F)
            # backtracking line search on max-norm of F
            t = 1.0
            for _ in range(12):
                Pn, ln, sn = Psi + t * dx[:N], lam + t * dx[N], s + t * dx[N + 1]
                Rn = self.residual(Pn, ln, sn)[0]
                en = 0.0 if fix_lam else self.h1(Pn) - (ln / p - 1 - ln)
                Fn = max(np.max(np.abs(Rn)), abs(Pn[self.i0] - v0), abs(en))
                if Fn < (1 - 0.3 * t) * nrm or t < 1e-3:
                    break
                t *= 0.5
            Psi, lam, s = Pn, ln, sn
        return Psi, lam, s, False

    def initial_guess(self, lam, p):
        beta = lam / (1 + lam)
        xi = np.exp(self.eta)
        Theta = np.exp(p * self.eta - 0.5 * (p - beta) * np.logaddexp(0, 2 * self.eta))
        Psi = np.exp(-self.c * self.eta) * Theta
        # rescale amplitude so that h1 matches the smoothness target for exponent p
        target = lam / p - 1 - lam
        return Psi * (target / self.h1(Psi))


if __name__ == "__main__":
    import sys, time
    prof = CCFProfile(L1=30, L2=150, N=4096, c=0.75)
    for p, lam0 in [(2, 1.0), (4, 0.6), (6, 0.47)]:
        t0 = time.time()
        Psi0 = prof.initial_guess(lam0, p)
        Psi, lam, s, ok = prof.solve(Psi0, lam0, p, verbose=True)
        R, G, W, D = prof.residual(Psi, lam, s)
        print(f"p={p}: converged={ok} lambda={lam:.12f} s={s:.2e} h1={prof.h1(Psi):.10f} "
              f"G(left)={G[5]:.10f}  time={time.time()-t0:.1f}s")
        np.save(f"ccf_p{p}.npy", np.concatenate([[lam], Psi]))
