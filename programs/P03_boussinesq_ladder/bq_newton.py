"""Newton–Krylov solver for the least-singular self-similar Boussinesq profile at fixed λ:
unknown Ψ̃ (weighted streamfunction on the log-polar grid), residual R(Ψ̃) = Ψ̃ − T(Ψ̃), T = Biot–Savart ∘ march.
Line search keeps the self-similar flow reversal-free (min V_r/r > 0) and requires residual decrease."""
import numpy as np, time
from scipy.sparse.linalg import LinearOperator, gmres
from bq_logpolar import BQLogPolar


def full(B, Y):
    """X = Ψ/r² on the whole grid from the unknowns Y = X[i0:] (s ≥ s_start); below s_start the exact local
    structure is the constant strain X(s,β) = X(s_start,β) (corrections O(e^{(m−1)s}) are negligible)."""
    X = np.empty((B.Ns, B.Nb + 1))
    X[B.i0:] = Y
    X[:B.i0] = Y[0]
    return X


def residual(B, Y):
    """unknowns Y = X[i0:], X = Ψ/r² = e^{(a−2)s}Ψ̃ (velocity-gradient scale); residual on s ≥ s_start only"""
    X = full(B, Y)
    Pn, info = B.T(X / B.ea2[:, None])
    if not np.all(np.isfinite(Pn)) or info['vrmin'] <= 0.02 or info['m'] <= 1.0:
        return None, info
    return (X - B.ea2[:, None] * Pn)[B.i0:], info


def newton(B, P, tol=1e-10, maxit=30, gm_tol=1e-9, gm_max=80, verbose=True, t0=None):
    t0 = time.time() if t0 is None else t0
    R, info = residual(B, P)
    if R is None:
        raise RuntimeError(f"invalid start: A={info['A']:.4f} m={info['m']:.4f} vrmin={info['vrmin']:.3f}")
    nrm = np.abs(R).max()
    shape = P.shape
    for it in range(maxit):
        if verbose:
            print(f"   NK it {it}: |R|={nrm:.3e} A={info['A']:.10f} m={info['m']:.10f} vrmin={info['vrmin']:.4f} t={time.time()-t0:.0f}s",
                  flush=True)
        if nrm < tol:
            return P, info, True
        eps = 1e-7

        def mv(v):
            v = v.reshape(shape)
            Rp, _ = residual(B, P + eps * v)
            if Rp is None:
                Rp, _ = residual(B, P - eps * v)
                return ((R - Rp) / eps).ravel()
            return ((Rp - R) / eps).ravel()
        J = LinearOperator((P.size, P.size), matvec=mv, dtype=float)
        dx, gi = gmres(J, -R.ravel(), rtol=gm_tol, atol=0.0, restart=gm_max, maxiter=3)
        dx = dx.reshape(shape)
        t = 1.0
        while t > 1e-3:
            Rn, infon = residual(B, P + t * dx)
            if Rn is not None and np.abs(Rn).max() < (1 - 0.1 * t) * nrm:
                break
            t *= 0.5
        if t <= 1e-3:
            if verbose:
                print("   line search failed", flush=True)
            return P, info, False
        P, R, info = P + t * dx, Rn, infon
        nrm = np.abs(R).max()
    return P, info, nrm < tol


def initial_guess(B, A0):
    """Y = X[i0:], X = Ψ/r² for a localised strain of rate A0 decaying like r^{-1/(1+λ)}"""
    s = B.s[B.i0:, None]; b = B.beta[None, :]; r = np.exp(s)
    return -(A0 / 2) * np.sin(2 * b) / (1 + r ** 2) ** (1 / (2 * (1 + B.lam)))


if __name__ == "__main__":
    import sys
    lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.92
    B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
    P = initial_guess(B, (3 + lam) / 2)
    P, info, ok = newton(B, P)
    print(f"λ={lam}: ok={ok} A={info['A']:.10f} m={info['m']:.10f}  (smooth ⇔ m = 2 ⇔ A = (3+λ)/2 = {(3+lam)/2:.6f})")
    np.save(f"bq_X_lam{lam:.4f}.npy", P)
