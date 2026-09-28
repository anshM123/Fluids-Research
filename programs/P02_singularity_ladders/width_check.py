"""Layer-width constant w/δ² at high resolution (the curvature estimate is resolution dependent: 3-point formula).
Re-solve saved pinned states on grids with 24/48/96 points per width and report w/δ², to compare with the
universal inner-layer prediction w/δ² = 4.7405/k² (k² = 2λΘ_s/ξ_s)."""
import numpy as np, sys
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info
import pinned_delta as P

for f in sys.argv[1:]:
    d = np.load(f)
    eta0, phi0, lam0 = d[0], d[1], float(d[2][0])
    sp = np.diff(eta0)
    ic0 = 1 + int(np.argmin(sp[:-1] + sp[1:]))
    etac = eta0[ic0]
    spl = CubicSpline(eta0, phi0)
    delta = None
    for pts in (24, 48, 96):
        P.HS = 0.03
        M, ic = P.new_grid(etac, 24 * sp.min(), pts)
        phi = spl(M.eta)
        F, den = P.Fpin(M, phi, lam0, ic)
        if delta is None:
            delta = den[ic]
        phi, lam, ok, hist, _ = P.solve(M, phi, lam0, ic, delta)
        F, den = P.Fpin(M, phi, lam, ic)
        e, md, w = layer_info(M, den)
        Th = np.exp(phi + M.c * M.eta)
        k2 = 2 * lam * Th[ic] / np.exp(M.eta[ic])
        print(f"{f}: δ={delta:.4e} pts={pts}: w/δ²={w/delta**2:.5f}   k²={k2:.5f}   inner-layer prediction 4.7405/k²={4.7405/k2:.5f}",
              flush=True)
