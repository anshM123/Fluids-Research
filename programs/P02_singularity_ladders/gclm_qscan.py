"""For fixed a, continue the gCLM self-similar family in c_l (natural continuation, fixed-c_l Newton–Krylov)
starting from the smooth q=1 profile, and record the local exponent q(c_l) and sonic margin min D."""
import numpy as np, sys, time
from gclm_nk import GCLM, solve_fixed

a = float(sys.argv[1])
d = np.load(f"gclm_q1_a{a:.3f}.npy")
cl0, phi0 = d[0], d[1:]
G = GCLM(a=a, L1=40, L2=80, N=8192, c=-0.4)
rows = []
for direction in (-1, +1):
    phi, cl = phi0.copy(), cl0
    step = 0.01 * cl0
    for k in range(200):
        cl_new = cl + direction * step
        if cl_new < 0.03:
            break
        ph, ok = solve_fixed(G, phi.copy(), cl_new)
        if not ok:
            step *= 0.5
            if step < 1e-5 * cl0:
                print(f"a={a}: stop at c_l={cl:.6f} (dir {direction})", flush=True)
                break
            continue
        phi, cl = ph, cl_new
        F, D, HO, Uxi = G.F(phi, cl, G.normval(phi, cl))
        q = G.q_of(phi, cl)
        rows.append((cl, q, G.h0(phi), D.min(), G.eta[np.argmin(D)]))
        if k % 5 == 0:
            print(f"a={a}: c_l={cl:.6f} q={q:.8f} h0={G.h0(phi):.6f} minD={D.min():.4e} at eta={G.eta[np.argmin(D)]:.2f}", flush=True)
        step = min(step * 1.2, 0.03 * cl0)
rows = np.array(sorted(rows))
np.save(f"gclm_qscan_a{a:.3f}.npy", rows)
print("done", rows.shape)
