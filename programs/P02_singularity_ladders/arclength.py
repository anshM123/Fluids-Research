"""Pseudo-arclength continuation of the CCF self-similar profile family in (φ, λ).
Tracks p(λ) = λ/(1+λ+h1) (local exponent at the origin) and the sonic margin min(1+λ+G)."""
import numpy as np
import scipy.linalg as sla
from ccf_newton import CCFNewton


def weights(S):
    # φ-part of the arclength norm restricted to the profile core |η|<10 (tails of ln Ψ are irrelevant)
    return (np.abs(S.eta) < 10).astype(float) * S.h / 20.0


def continue_branch(S, phi, lam, ds, nsteps, phi0=None, tol=1e-11, verbose=True, lam_stop=(0.05, 3.0)):
    N = S.N
    if phi0 is None:
        phi0 = phi[S.i0]
    w = weights(S)
    # initial tangent
    F, den = S.F(phi, lam, phi0)
    J, Fl = S.jac(phi, lam, den)
    tphi = -sla.solve(J, Fl)
    tl = 1.0
    nrm = np.sqrt((w * tphi) @ tphi + tl**2)
    tphi, tl = tphi / nrm * np.sign(ds), tl / nrm * np.sign(ds)
    ds = abs(ds)
    out = [(lam, S.p_of(phi, lam), den.min())]
    sols = [(lam, phi.copy())]
    for step in range(nsteps):
        # predictor
        ph, lm = phi + ds * tphi, lam + ds * tl
        ok = False
        for it in range(25):
            F, den = S.F(ph, lm, phi0)
            if den.min() <= 0:
                break
            g = (w * tphi) @ (ph - phi) + tl * (lm - lam) - ds
            res = max(np.max(np.abs(F)), abs(g))
            if res < tol:
                ok = True
                break
            J, Fl = S.jac(ph, lm, den)
            A = np.zeros((N + 1, N + 1))
            A[:N, :N] = J
            A[:N, N] = Fl
            A[N, :N] = w * tphi
            A[N, N] = tl
            d = sla.solve(A, -np.concatenate([F, [g]]))
            ph, lm = ph + d[:N], lm + d[N]
        if not ok:
            ds *= 0.5
            if verbose:
                print(f"   step {step}: corrector failed, ds -> {ds:.2e}")
            if ds < 1e-5:
                break
            continue
        # new tangent (solve bordered system with previous tangent for orientation)
        J, Fl = S.jac(ph, lm, den)
        A = np.zeros((N + 1, N + 1))
        A[:N, :N] = J
        A[:N, N] = Fl
        A[N, :N] = w * tphi
        A[N, N] = tl
        t = sla.solve(A, np.concatenate([np.zeros(N), [1.0]]))
        nrm = np.sqrt((w * t[:N]) @ t[:N] + t[N] ** 2)
        tphi, tl = t[:N] / nrm, t[N] / nrm
        phi, lam = ph, lm
        p = S.p_of(phi, lam)
        out.append((lam, p, den.min()))
        sols.append((lam, phi.copy()))
        if verbose:
            print(f"step {step}: lam={lam:.8f} p={p:.8f} min_den={den.min():.5f} dlam/ds={tl:+.3f} ds={ds:.3e}", flush=True)
        if it < 4:
            ds = min(ds * 1.5, 0.02)
        if not (lam_stop[0] < lam < lam_stop[1]):
            break
    return np.array(out), sols


if __name__ == "__main__":
    import sys
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2048
    direction = float(sys.argv[2]) if len(sys.argv) > 2 else -1.0
    nsteps = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    S = CCFNewton(30, 150, N, 0.75)
    phi = S.guess(0.6, 2.0)
    phi, ok = S.solve_fixed(phi, 0.6)
    out, sols = continue_branch(S, phi, 0.6, direction * 0.01, nsteps)
    tag = "down" if direction < 0 else "up"
    np.save(f"arc_{tag}_N{N}.npy", out)
    np.save(f"arc_{tag}_N{N}_sols.npy", np.array([np.concatenate([[l], p]) for l, p in sols]))
