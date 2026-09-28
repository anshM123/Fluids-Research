import numpy as np, sys, time
from ccf_nk import CCFNK, arclength
N = int(sys.argv[1]); d = float(sys.argv[2]); nst = int(sys.argv[3]); c = float(sys.argv[4]) if len(sys.argv)>4 else 0.7
S = CCFNK(30, 120, N, c)
phi = S.guess(0.6, 2.0); phi, ok = S.picard(phi, 0.6, tol=1e-6); phi, ok = S.solve_fixed(phi, 0.6)
print("start", ok, S.p_of(phi, 0.6), flush=True)
t0 = time.time()
out, sols = arclength(S, phi, 0.6, d*0.02, nst, dsmax=0.1)
tag = f"N{N}_c{c}_{'down' if d<0 else 'up'}"
np.save(f"nkarc_{tag}.npy", out)
np.save(f"nkarcsols_{tag}.npy", np.array([np.concatenate([[l], p]) for l, p in sols[::5]]))
print("done", time.time()-t0)
