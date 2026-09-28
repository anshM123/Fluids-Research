import numpy as np, time
from gclm_nk import GCLM
for a in [0.0]:
    G = GCLM(a=a, L1=40, L2=80, N=8192, c=-0.4)
    for cl0 in [1.0]:
        phi = G.guess(cl0)
        # scale amplitude so that h0 is consistent with q=1: h0 = (1+cl)/(1-a)
        h0t = (1 + cl0) / (1 - a); phi = phi + np.log(h0t / G.h0(phi))
        t0 = time.time()
        phi, cl, ok = G.solve_q(phi, cl0, 1.0, verbose=True)
        F, D, HO, Uxi = G.F(phi, cl, G.normval(phi, cl))
        print(f"a={a}: ok={ok} c_l={cl:.12f} q={G.q_of(phi,cl):.10f} minD={D.min():.4f} h0={G.h0(phi):.8f} ({time.time()-t0:.1f}s)")
        Om = -np.exp(phi + G.c*G.eta); xi = np.exp(G.eta)
        for x in [0.01, 0.1, 1, 10, 100]:
            i = np.argmin(abs(xi-x)); print(f"   xi={xi[i]:.3g} Omega={Om[i]:.6e}  Omega*(1+xi^2)/xi={Om[i]*(1+xi[i]**2)/xi[i]:.6f}")
