import numpy as np
from ccf_newton import CCFNewton
S = CCFNewton(30,150,2048,0.75)
phi = S.guess(0.6, 2.0)
F, den = S.F(phi, 0.6, phi[S.i0]); print("init |F|", np.abs(F).max(), "min den", den.min())
phi, ok = S.solve_fixed(phi, 0.6, verbose=True)
print("fixed ok", ok, "p", S.p_of(phi, 0.6))
# FD check of dF/dlam and jac
F0, den = S.F(phi, 0.6, phi[S.i0]); J, Fl = S.jac(phi, 0.6, den)
F1,_ = S.F(phi, 0.6+1e-7, phi[S.i0]); print("dF/dlam FD err", np.abs((F1-F0)/1e-7 - Fl).max(), np.abs(Fl).max())
e = np.zeros(S.N); e[700]=1e-7; F2,_ = S.F(phi+e, 0.6, phi[S.i0]); print("J col FD err", np.abs((F2-F0)/1e-7 - J[:,700]).max())
print("cond J", np.linalg.cond(J))
