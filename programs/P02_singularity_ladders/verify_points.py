"""Grid-refinement check at FIXED sonic depth δ: re-solve a pinned-branch state on grids with different
centre resolution (points per layer width) and far-field spacing h_s, and compare λ and p.
Usage: verify_points.py STATE.npy [STATE2.npy ...]"""
import numpy as np, sys
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info
import pinned_delta as P

for f in sys.argv[1:]:
    d = np.load(f)
    eta0, phi0, lam0 = d[0], d[1], float(d[2][0])
    sp = np.diff(eta0)
    ic0 = 1 + int(np.argmin(sp[:-1] + sp[1:]))            # centre point: smallest pair of adjacent spacings
    etac = eta0[ic0]
    spl = CubicSpline(eta0, phi0)
    res = []
    for pts, hs in ((24, 0.03), (36, 0.02), (48, 0.03), (24, 0.02)):
        P.HS = hs
        M, ic = P.new_grid(etac, 24 * sp.min(), pts)
        phi = spl(M.eta)
        F, den = P.Fpin(M, phi, lam0, ic)
        # target δ: the source state's value (evaluated on the first grid)
        if not res:
            delta = den[ic]
        phi, lam, ok, hist, _ = P.solve(M, phi, lam0, ic, delta)
        p = M.p_of(phi, lam)
        res.append((pts, hs, M.N, lam, p, ok))
        print(f"{f}: δ={delta:.6e} pts/width={pts} h_s={hs} N={M.N}: λ={lam:.13f} p={p:.13f} ok={ok}", flush=True)
    lams = np.array([r[3] for r in res]); ps = np.array([r[4] for r in res])
    print(f"   spread: Δλ={lams.max()-lams.min():.2e}  Δp={ps.max()-ps.min():.2e}", flush=True)
