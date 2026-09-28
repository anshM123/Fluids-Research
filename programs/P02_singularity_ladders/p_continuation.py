"""Continuation of the CCF profile family using the local exponent p as parameter (λ unknown)."""
import numpy as np, sys, time, json
from ccf_mapped import CCFMapped
from mapped_continuation import layer_info, make_grid
d = np.load('mapped_l2_hs0.02_eps0.025.npy'); lam = 0.4713242277712
M = make_grid(-0.9977, 0.009, hs=0.03)
phi = M.from_other(d[0], d[1]); phi, lam, ok = M.solve_p(phi, lam, 2.0)
target = M.normval(phi, lam)
# natural continuation in λ down to 0.459 (p max region), then switch to p-parameter
for lam_new in np.arange(0.470, 0.4575, -0.001):
    phi, ok = M.solve_fixed(phi, lam_new, phi0=target); lam = lam_new
    F, den = M.F(phi, lam, target); e, md, w = layer_info(M, den)
    if w / (M.hs * M.eps) < 15:
        Mn = make_grid(e, w, hs=0.03); phi = Mn.from_other(M.eta, phi); M = Mn
        phi, ok = M.solve_fixed(phi, lam, phi0=target)
p = M.p_of(phi, lam)
print(f"reached lam={lam:.6f} p={p:.8f}", flush=True)
pts = np.concatenate([np.arange(p - 0.0001, 2.0105, -0.0001), np.arange(2.0105, 1.99, -0.00025)])
out = []
for pt in pts:
    ph_try, l_try, ok = M.solve_p(phi.copy(), lam, pt, phi0=target, maxit=40)
    if not ok:
        # try a regrid then retry
        F, den = M.F(phi, lam, target); e, md, w = layer_info(M, den)
        Mn = make_grid(e, w * 0.7, hs=0.03); phi2 = Mn.from_other(M.eta, phi)
        ph_try, l_try, ok = Mn.solve_p(phi2, lam, pt, phi0=target, maxit=40)
        if ok: M = Mn
    if not ok:
        print(f"p={pt:.5f}: FAILED (possible fold in p)", flush=True); break
    phi, lam = ph_try, l_try
    F, den = M.F(phi, lam, target); e, md, w = layer_info(M, den)
    out.append((pt, lam, md, e, w, M.N, M.eps))
    print(f"p={pt:.5f}: lam={lam:.10f} min_den={md:.4e} layer={e:.4f} w={w:.2e} pts/w={w/(M.hs*M.eps):.1f} N={M.N}", flush=True)
    np.save(f"pcont_sol_p{pt:.5f}.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
    if w / (M.hs * M.eps) < 15 or w / (M.hs * M.eps) > 60:
        Mn = make_grid(e, w, hs=0.03); phi = Mn.from_other(M.eta, phi); M = Mn
        phi, ok = M.solve_fixed(phi, lam, phi0=target)
np.save("pcont_branch.npy", np.array(out))
