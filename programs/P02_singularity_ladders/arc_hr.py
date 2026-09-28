import numpy as np, sys
from ccf_newton import CCFNewton
from arclength import continue_branch
N = int(sys.argv[1]); c = float(sys.argv[2]); L1 = float(sys.argv[3]); L2 = float(sys.argv[4]); d = float(sys.argv[5]); nst = int(sys.argv[6])
S = CCFNewton(L1, L2, N, c)
phi = S.guess(0.6, 2.0); phi, ok = S.solve_fixed(phi, 0.6)
print("start ok", ok, "p", S.p_of(phi, 0.6), flush=True)
out, sols = continue_branch(S, phi, 0.6, d*0.01, nst)
tag = f"N{N}_c{c}_L{L1:.0f}_{L2:.0f}_{'down' if d<0 else 'up'}"
np.save(f"arc_{tag}.npy", out)
np.save(f"arcsols_{tag}.npy", np.array([np.concatenate([[l], p]) for l, p in sols]))
