"""
Matrix-free Newton–Krylov solver for CCF self-similar profiles (integral form, see ccf_newton.py).
Jacobian = identity + (cumulative integral) ∘ (multiplication) ∘ (Hilbert convolution): a second-kind
Fredholm-type operator → GMRES converges in O(10–100) iterations independent of N.
Supports: fixed-λ solves, smoothness-constrained solves (λ unknown, p prescribed), pseudo-arclength.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres


class CCFNK:
    def __init__(self, L1=30.0, L2=120.0, N=16384, c=0.7):
        self.N, self.c, self.L1, self.L2 = N, c, L1, L2
        self.h = (L1 + L2) / N
        self.eta = -L1 + self.h * np.arange(N)
        k = 2 * np.pi * np.fft.fftfreq(N, d=self.h)
        self.mK = -1j * np.tanh(np.pi * (k - 1j * c) / 2)
        self.E = np.exp((c - 1) * self.eta)
        self.i0 = int(np.argmin(np.abs(self.eta)))
        self.iR = int(np.argmin(np.abs(self.eta - 50.0)))      # far-field normalisation point
        self.etaR = self.eta[self.iR]
        nodes = np.arange(-3, 5)
        rhs = np.array([1.0 / (m + 1) for m in range(8)])
        self.w8 = np.linalg.solve(np.vander(nodes, 8, increasing=True).T, rhs)
        self.wl = [np.linalg.solve(np.vander(np.arange(-s, 8 - s), 8, increasing=True).T, rhs) for s in range(8)]

    # ---- building blocks ----------------------------------------------------------------------
    def cumint(self, X):
        n, h, w = self.N, self.h, self.w8
        mid = np.zeros(n - 1)
        for j, off in enumerate(range(-3, 5)):
            mid[3:n - 4] += w[j] * X[3 + off:n - 4 + off]
        for i in list(range(0, 3)) + list(range(n - 4, n - 1)):
            s = i if i < 3 else 7 - (n - 1 - i)
            mid[i] = self.wl[s] @ X[i - s:i - s + 8]
        C = np.zeros(n)
        C[1:] = np.cumsum(mid * h)
        return C - C[self.i0]

    def conv(self, f):
        return np.real(np.fft.ifft(self.mK * np.fft.fft(f)))

    def G(self, phi):
        return self.E * self.conv(np.exp(phi))

    def h1(self, phi):
        return -(2 / np.pi) * self.h * np.sum(np.exp(phi) * self.E)

    def p_of(self, phi, lam):
        return lam / (1 + lam + self.h1(phi))

    def normval(self, phi, lam):
        """log far-field amplitude ln C (Θ ~ C ξ^β): invariant-free, monotone along the scaling orbit"""
        return phi[self.iR] + (self.c - lam / (1 + lam)) * self.etaR

    def F(self, phi, lam, phi0):
        den = 1 + lam + self.G(phi)
        F = phi - phi[self.i0] - self.cumint(lam / den - self.c)
        F[self.i0] = self.normval(phi, lam) - phi0
        return F, den

    def Jv(self, v, phi, lam, den):
        Psi = np.exp(phi)
        dG = self.E * self.conv(Psi * v)
        out = v - v[self.i0] - self.cumint(-lam / den**2 * dG)
        out[self.i0] = v[self.iR]
        return out

    def Fl(self, lam, den):
        out = -self.cumint(1 / den - lam / den**2)
        out[self.i0] = -self.etaR / (1 + lam) ** 2
        return out

    def dh1(self, phi):
        return -(2 / np.pi) * self.h * np.exp(phi) * self.E

    # ---- solvers --------------------------------------------------------------------------------
    def _gmres(self, matvec, rhs, n, tol=1e-13):
        A = LinearOperator((n, n), matvec=matvec, dtype=float)
        x, info = gmres(A, rhs, rtol=tol, atol=0.0, restart=80, maxiter=40)
        return x, info

    def solve_fixed(self, phi, lam, phi0=None, tol=1e-12, maxit=30, verbose=False):
        if phi0 is None:
            phi0 = self.normval(phi, lam)
        for it in range(maxit):
            F, den = self.F(phi, lam, phi0)
            nrm = np.max(np.abs(F))
            if verbose:
                print(f"   NK it {it}: |F|={nrm:.2e} minden={den.min():.4f}", flush=True)
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
                print(f"   NK-p it {it}: |F|={nrm:.2e} lam={lam:.15f}", flush=True)
            if nrm < tol:
                return phi, lam, True
            Fl = self.Fl(lam, den)
            dh = self.dh1(phi)
            sc = 1.0 / max(np.linalg.norm(dh), abs(1 / p - 1))

            def mv(x):
                v, s = x[:N], x[N]
                return np.concatenate([self.Jv(v, phi, lam, den) + s * Fl, [sc * (dh @ v - (1 / p - 1) * s)]])
            d, info = self._gmres(mv, -np.concatenate([F, [sc * g]]), N + 1, tol=1e-11)
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

    def guess(self, lam, p0):
        beta = lam / (1 + lam)
        Theta = np.exp(p0 * self.eta - 0.5 * (p0 - beta) * np.logaddexp(0, 2 * self.eta))
        Psi = np.exp(-self.c * self.eta) * Theta
        Psi *= (lam / p0 - 1 - lam) / self.h1(np.log(Psi))
        return np.log(Psi)

    def picard(self, phi, lam, relax=0.3, tol=1e-7, maxit=5000):
        for it in range(maxit):
            den = 1 + lam + self.G(phi)
            if den.min() <= 0:
                return phi, False
            new = self.cumint(lam / den - self.c)
            new = new - new[self.i0] + phi[self.i0]
            diff = np.max(np.abs(new - phi))
            phi = (1 - relax) * phi + relax * new
            if diff < tol:
                return phi, True
        return phi, False

    def resample(self, other, phi_other):
        """interpolate a solution from another grid (same c) onto this grid (tails extrapolated linearly)"""
        return np.interp(self.eta, other.eta, phi_other,
                         left=None, right=None) if True else None


def arclength(S, phi, lam, ds, nsteps, phi0=None, tol=1e-11, dsmax=0.08, verbose=True, core=10.0):
    N = S.N
    if phi0 is None:
        phi0 = S.normval(phi, lam)
    w = (np.abs(S.eta) < core) * S.h / (2 * core)
    F, den = S.F(phi, lam, phi0)
    Fl = S.Fl(lam, den)
    tphi, _ = S._gmres(lambda v: S.Jv(v, phi, lam, den), -Fl, N)
    tl = 1.0
    nrm = np.sqrt((w * tphi) @ tphi + tl**2)
    sg = np.sign(ds)
    tphi, tl, ds = tphi / nrm * sg, tl / nrm * sg, abs(ds)
    out = [(lam, S.p_of(phi, lam), den.min())]
    sols = [(lam, phi.copy())]
    step = 0
    while step < nsteps:
        ph, lm = phi + ds * tphi, lam + ds * tl
        ok = False
        for it in range(20):
            F, den = S.F(ph, lm, phi0)
            if den.min() <= 0:
                break
            g = (w * tphi) @ (ph - phi) + tl * (lm - lam) - ds
            if max(np.max(np.abs(F)), abs(g)) < tol:
                ok = True
                break
            Fl = S.Fl(lm, den)
            sc = 1.0 / max(np.linalg.norm(w * tphi), abs(tl))

            def mv(x, ph=ph, lm=lm, den=den, Fl=Fl, sc=sc):
                v, s = x[:N], x[N]
                return np.concatenate([S.Jv(v, ph, lm, den) + s * Fl, [sc * ((w * tphi) @ v + tl * s)]])
            d, info = S._gmres(mv, -np.concatenate([F, [sc * g]]), N + 1, tol=1e-10)
            ph, lm = ph + d[:N], lm + d[N]
        if not ok:
            ds *= 0.5
            if verbose:
                print(f"   corrector failed, ds -> {ds:.2e}", flush=True)
            if ds < 1e-6:
                break
            continue
        Fl = S.Fl(lm, den)

        sc = 1.0 / max(np.linalg.norm(w * tphi), abs(tl))

        def mvt(x, ph=ph, lm=lm, den=den, Fl=Fl, sc=sc):
            v, s = x[:N], x[N]
            return np.concatenate([S.Jv(v, ph, lm, den) + s * Fl, [sc * ((w * tphi) @ v + tl * s)]])
        t, info = S._gmres(mvt, np.concatenate([np.zeros(N), [sc]]), N + 1, tol=1e-10)
        nrm = np.sqrt((w * t[:N]) @ t[:N] + t[N] ** 2)
        tphi, tl = t[:N] / nrm, t[N] / nrm
        phi, lam = ph, lm
        p = S.p_of(phi, lam)
        out.append((lam, p, den.min()))
        sols.append((lam, phi.copy()))
        if verbose:
            print(f"step {step}: lam={lam:.10f} p={p:.10f} min_den={den.min():.6f} "
                  f"eta_min={S.eta[np.argmin(den)]:.3f} dlam/ds={tl:+.4f} ds={ds:.2e}", flush=True)
        if it < 5:
            ds = min(ds * 1.4, dsmax)
        step += 1
    return np.array(out), sols
