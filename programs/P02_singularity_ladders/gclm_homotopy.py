"""Homotopy in a for the smooth (q=1) gCLM self-similar profile, from the exact CLM profile at a=0."""
import numpy as np, time
from gclm_nk import GCLM
G = GCLM(a=0.0, L1=40, L2=80, N=8192, c=-0.4)
phi = G.guess(1.0); cl = 1.0
res = []
for a in np.round(np.arange(0.0, 0.96, 0.025), 4):
    G.a = a
    t0 = time.time()
    phi_n, cl_n, ok = G.solve_q(phi.copy(), cl, 1.0, maxit=30)
    if not ok:
        print(f"a={a}: FAILED (cl={cl_n:.6f})", flush=True); break
    phi, cl = phi_n, cl_n
    F, D, HO, Uxi = G.F(phi, cl, G.normval(phi, cl))
    res.append((a, cl, G.h0(phi), D.min()))
    print(f"a={a:.3f}: c_l={cl:.10f} h0={G.h0(phi):.8f} (check (1+c_l)/(1-a)={(1+cl)/(1-a):.8f}) minD={D.min():.4f} ({time.time()-t0:.0f}s)", flush=True)
    np.save(f"gclm_q1_a{a:.3f}.npy", np.concatenate([[cl], phi]))
np.save("gclm_homotopy.npy", np.array(res))
