import numpy as np, time
from ccf_newton import CCFNewton
S = CCFNewton(30,150,2048,0.75)
t0=time.time(); phi = S.guess(0.6,2.0); phi, ok = S.solve_fixed(phi, 0.6, verbose=True); print('fixed', ok, S.p_of(phi,0.6), time.time()-t0, flush=True)
phi, lam, ok2 = S.solve_p(phi, 0.6, 2.0, verbose=True); print('p=2', ok2, repr(lam), time.time()-t0)
