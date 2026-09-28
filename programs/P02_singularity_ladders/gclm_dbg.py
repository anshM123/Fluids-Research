import numpy as np
from gclm_nk import GCLM
a=0.5
G = GCLM(a=a, L1=40, L2=80, N=8192, c=-0.4)
phi = G.guess(1.0); phi = phi + np.log(((1+1.0)/(1-a))/G.h0(phi))
phi, cl, ok = G.solve_q(phi, 1.0, 1.0, verbose=True, maxit=25)
F, D, HO, Uxi = G.F(phi, cl, G.normval(phi, cl))
print("ok", ok, "cl", cl, "|F|", abs(F).max(), "q", G.q_of(phi,cl), "h0", G.h0(phi), "minD", D.min())
np.save("gclm_a05_try.npy", np.concatenate([[cl], phi]))
