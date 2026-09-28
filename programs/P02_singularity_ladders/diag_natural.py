import numpy as np, time
from ccf_mapped import CCFMapped
from mapped_continuation import layer_info, make_grid
d = np.load('mapped_l2_hs0.02_eps0.025.npy'); lam = 0.4713242277712
M = make_grid(-0.9977, 0.009, hs=0.03)
phi = M.from_other(d[0], d[1]); phi, lam, ok = M.solve_p(phi, lam, 2.0)
target = M.normval(phi, lam)
saved = []
for lam_new in np.arange(0.470, 0.4539, -0.001):
    phi, ok = M.solve_fixed(phi, lam_new, phi0=target)
    F, den = M.F(phi, lam_new, target)
    e, md, w = layer_info(M, den)
    # count local minima of den below 0.3
    locmin = np.where((den[1:-1] < den[:-2]) & (den[1:-1] < den[2:]) & (den[1:-1] < 0.3))[0] + 1
    print(f"lam={lam_new:.4f} ok={ok} p={M.p_of(phi,lam_new):.8f} min_den={md:.3e} at {e:.4f} w={w:.2e} N={M.N} eps={M.eps:.1e} | den local minima: " +
          ", ".join(f"({M.eta[j]:.3f},{den[j]:.3e})" for j in locmin[:6]), flush=True)
    if not ok: break
    saved.append(np.vstack([M.eta, phi, np.full(M.N, lam_new)]))
    if w / (M.hs * M.eps) < 15 or w / (M.hs*M.eps) > 60:
        Mn = make_grid(e, w, hs=0.03); phi = Mn.from_other(M.eta, phi); M = Mn
        phi, ok = M.solve_fixed(phi, lam_new, phi0=target)
np.save("diag_natural.npy", np.array(saved, dtype=object), allow_pickle=True)
