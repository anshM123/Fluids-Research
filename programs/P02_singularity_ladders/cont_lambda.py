import numpy as np, time, sys
from ccf_newton import CCFNewton
S = CCFNewton(30,150,2048,0.75)
phi = S.guess(0.6,2.0); phi, ok = S.solve_fixed(phi, 0.6)
start = phi.copy()
out = []
for lams in [np.arange(0.59, 0.30, -0.01), np.arange(0.61, 1.2, 0.01)]:
    phi = start.copy()
    for lam in lams:
        phi_new, ok = S.solve_fixed(phi.copy(), lam) if False else (None, False)
        # Newton only (no Picard) from previous solution
        p0 = phi[S.i0]; ph = phi.copy(); conv=False
        for it in range(30):
            F, den = S.F(ph, lam, p0)
            if den.min() <= 0: break
            if np.max(np.abs(F)) < 1e-12: conv=True; break
            J,_ = S.jac(ph, lam, den)
            import scipy.linalg as sla
            ph = ph - sla.solve(J, F)
        if not conv:
            print(f"lam={lam:.3f}: Newton failed (min den {den.min():.3e})", flush=True); break
        phi = ph
        F, den = S.F(phi, lam, p0)
        p = S.p_of(phi, lam)
        out.append((lam, p, den.min(), S.eta[np.argmin(den)]))
        print(f"lam={lam:.3f}  p={p:.8f}  min(1+lam+G)={den.min():.5f} at eta={S.eta[np.argmin(den)]:.2f}", flush=True)
np.save("cont_lambda_branch1.npy", np.array(sorted(out)))
