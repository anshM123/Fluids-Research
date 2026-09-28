"""Map the CCF self-similar solution branch p(λ) with the adaptive mapped solver (natural continuation in λ),
from the stable profile λ0 through λ1 down to λ2 and beyond (until the fold), saving diagnostics."""
import numpy as np, sys, time
from ccf_mapped import CCFMapped
from mapped_continuation import layer_info, make_grid
from ccf_nk import CCFNK

S0 = CCFNK(30, 120, 16384, 0.7)
d = np.load('ladder_F16k_up_cross1.npy')                 # λ1 profile (well resolved)
M = make_grid(-1.45, 0.26, hs=0.03)
phi = M.from_other(S0.eta, d[1:]); phi, lam, ok = M.solve_p(phi, d[0], 2.0)
target = M.normval(phi, lam)
start_phi, start_lam, startM = phi.copy(), lam, M
rows = []


def record(M, phi, lam):
    F, den = M.F(phi, lam, target)
    e, md, w = layer_info(M, den)
    rows.append((lam, M.p_of(phi, lam), md, e, w, M.h1(phi)))


for direction, lams in [(+1, np.arange(0.61, 2.001, 0.01)), (-1, np.arange(0.600, 0.4575, -0.0025))]:
    M, phi, lam = startM, start_phi.copy(), start_lam
    for lam_new in lams:
        phi_new, ok = M.solve_fixed(phi.copy(), lam_new, phi0=target)
        if not ok:
            print(f"fail at {lam_new:.4f}", flush=True)
            break
        phi, lam = phi_new, lam_new
        record(M, phi, lam)
        F, den = M.F(phi, lam, target)
        e, md, w = layer_info(M, den)
        if w / (M.hs * M.eps) < 15 or (w / (M.hs * M.eps) > 80 and M.eps < 1):
            Mn = make_grid(e, w, hs=0.03)
            phi = Mn.from_other(M.eta, phi); M = Mn
            phi, ok = M.solve_fixed(phi, lam, phi0=target)
        print(f"lam={lam:.4f} p={rows[-1][1]:.8f} min_den={md:.4e} layer={e:.3f} w={w:.2e} N={M.N}", flush=True)
rows = np.array(sorted(rows))
np.save("branch_map.npy", rows)
print("saved", rows.shape)
