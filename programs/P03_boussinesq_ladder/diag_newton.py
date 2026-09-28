import numpy as np, time
from scipy.sparse.linalg import LinearOperator, gmres
from bq_logpolar import BQLogPolar
from bq_newton import residual, initial_guess, full
lam = 1.92
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
X = initial_guess(B, (3 + lam) / 2)
R, info = residual(B, X)
shape = X.shape
for eps in (1e-5, 1e-7, 1e-9):
    v = np.random.default_rng(0).standard_normal(shape) * np.exp(-0.05 * B.s[B.i0:, None] ** 2 / 100)
    Rp, _ = residual(B, X + eps * v); Rm, _ = residual(B, X - eps * v)
    print(f"eps={eps:.0e}: |Jv| (fwd) {np.linalg.norm((Rp - R) / eps):.6e}  (central) {np.linalg.norm((Rp - Rm) / (2 * eps)):.6e}", flush=True)
eps = 1e-7
res_hist = []
def mv(v):
    v = v.reshape(shape)
    Rp, _ = residual(B, X + eps * v)
    return ((Rp - R) / eps).ravel()
J = LinearOperator((X.size, X.size), matvec=mv, dtype=float)
t0 = time.time()
cb = lambda rk: res_hist.append(rk)
dx, gi = gmres(J, -R.ravel(), rtol=1e-8, atol=0.0, restart=200, maxiter=2, callback=cb, callback_type='pr_norm')
print(f"GMRES info={gi}, its={len(res_hist)}, rel residuals: {[f'{h:.1e}' for h in res_hist[::10]]} ... {res_hist[-1]:.2e}  ({time.time()-t0:.0f}s)", flush=True)
dx = dx.reshape(shape)
print("|dx|max =", np.abs(dx).max(), " at s =", B.s[B.i0 + np.unravel_index(np.argmax(np.abs(dx)), shape)[0]])
for t in (1e-3, 1e-2, 0.1, 0.3, 1.0):
    Rn, infon = residual(B, X + t * dx)
    if Rn is None:
        print(f"t={t}: invalid (A={infon['A']:.4f} m={infon['m']:.4f} vrmin={infon['vrmin']:.4f})"); continue
    print(f"t={t}: |R|max={np.abs(Rn).max():.4e} |R|2={np.linalg.norm(Rn):.4e}  (linear model |R|2 ≈ {(1-t)*np.linalg.norm(R):.4e})  A={infon['A']:.6f} m={infon['m']:.6f} vrmin={infon['vrmin']:.4f}", flush=True)
