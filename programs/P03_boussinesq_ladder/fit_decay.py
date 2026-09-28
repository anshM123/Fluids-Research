"""Fit m(λ) − 2 = C z^{−p} e^{−γ z} cos(ω z + φ) on z ≥ zmin (scan arrays: λ, A, m, ...)."""
import numpy as np, sys
from scipy.optimize import least_squares
f = sys.argv[1]; zmin = float(sys.argv[2]); zmax = float(sys.argv[3]) if len(sys.argv) > 3 else 1e9
pfix = float(sys.argv[4]) if len(sys.argv) > 4 else None
D = np.load(f); lam, m = D[:, 0], D[:, 2]
z = 1 / (lam - 1); y = m - 2
o = np.argsort(z); z, y = z[o], y[o]
sel = (z >= zmin) & (z <= zmax); z, y = z[sel], y[sel]
def model(p, z):
    lnC, pp, g, w, ph = p
    if pfix is not None: pp = pfix
    return np.exp(lnC) * z ** (-pp) * np.exp(-g * z) * np.cos(w * z + ph)
def res(p):
    return (model(p, z) - y) / (np.abs(model(p, z)) + 1e-3 * np.abs(y).max() * 0 + 1e-12 + 0.05 * np.exp(p[0]) * z ** (-(pfix if pfix is not None else p[1])) * np.exp(-p[2] * z))
best = None
for w0 in (2.2, 2.4, 2.5, 2.1):
    for g0 in (1.0, 1.3, 1.6):
        p0 = [np.log(np.abs(y).max() * np.exp(g0 * z[0])), 1.0, g0, w0, 0.0]
        try:
            r = least_squares(res, p0, max_nfev=20000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
p = best.x
pp = pfix if pfix is not None else p[1]
print(f"fit on z∈[{z.min():.2f},{z.max():.2f}] ({len(z)} pts): p={pp:.3f} γ={p[2]:.5f} ω={p[3]:.5f} (half-period π/ω={np.pi/p[3]:.5f}, "
      f"ratio e^(γπ/ω)={np.exp(p[2]*np.pi/p[3]):.4f}) rel.rms={np.sqrt(np.mean(best.fun**2)):.2e}")
