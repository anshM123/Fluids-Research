"""Fill the branch between λ = 0.458 and 0.4626 (natural continuation in λ, dense Newton on the geometric grid)
so that the plotted branch is continuous from the mapped-grid map (λ ≥ 0.4625) to the arclength/pinned data."""
import numpy as np
import scipy.linalg as sla
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info
from dense_arclength import jac
from ccf_sinhgrid import make_sinh_grid

d = np.load("pdiag_start.npy"); eta0, phi0, lam = d[0], d[1], float(d[2][0])
rows = []
M = make_sinh_grid(-0.9632, 0.00162, hs=0.03, pts=24)
phi = CubicSpline(eta0, phi0)(M.eta)
target = M.normval(phi, lam)
for L in np.arange(0.4580, 0.46301, 0.0005):
    for it in range(30):
        F, den = M.F(phi, L, target)
        if np.abs(F).max() < 1e-11:
            break
        J, Fl = jac(M, phi, L, den)
        phi = phi - sla.solve(J, F, check_finite=False)
    F, den = M.F(phi, L, target)
    e, md, w = layer_info(M, den)
    p = M.p_of(phi, L)
    rows.append((L, p, md, e, w))
    print(f"λ={L:.4f} p={p:.10f} δ={md:.5e} layer={e:.5f} w={w:.3e} |F|={np.abs(F).max():.1e}", flush=True)
    # re-centre the grid on the layer
    Mn = make_sinh_grid(e, w, hs=0.03, pts=24)
    phi = CubicSpline(M.eta, phi)(Mn.eta); M = Mn
np.save("gap_branch.npy", np.array(rows))
