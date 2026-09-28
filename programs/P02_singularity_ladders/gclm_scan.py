import numpy as np, sys, time
from gclm_nk import GCLM, solve_fixed
a = float(sys.argv[1])
G = GCLM(a=a, L1=40, L2=80, N=8192, c=-0.4)
phi = G.guess(1.0); phi = phi + np.log(((1+1.0)/(1-a if a<1 else 0.5))/G.h0(phi))
# first get a q=1 profile at this a (c_l unknown)
phi, cl, ok = G.solve_q(phi, 1.0, 1.0)
print(f"a={a}: q=1 profile c_l={cl:.10f} ok={ok}", flush=True)
target = G.normval(phi, cl)
start_phi, start_cl = phi.copy(), cl
for direction in (+1, -1):
    phi, c0 = start_phi.copy(), start_cl
    for k in range(1, 60):
        cl_new = c0 + direction * 0.02 * k
        if cl_new <= 0.05: break
        ph, ok = solve_fixed(G, phi.copy(), cl_new, target=G.normval(phi, cl_new) if False else None)
        if not ok:
            print(f"   c_l={cl_new:.3f}: fail", flush=True); break
        phi = ph
        F, D, HO, Uxi = G.F(phi, cl_new, G.normval(phi, cl_new))
        print(f"   c_l={cl_new:.3f} q={G.q_of(phi, cl_new):.8f} h0={G.h0(phi):.6f} minD={D.min():.4f} at eta={G.eta[np.argmin(D)]:.2f}", flush=True)
