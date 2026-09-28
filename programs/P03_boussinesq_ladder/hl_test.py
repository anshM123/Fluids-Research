import numpy as np, sys, time
from hl_solver import HL, newton
lam = float(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 8192
S = HL(lam, N=N)
q = S.guess((3 + lam) / 2)
t = time.time(); R, r = S.residual(q); print("residual eval", time.time() - t, "s; |R|", np.abs(R).max(), "A", r['A'], "m", r['m'])
q, info, ok = newton(S, q, tol=1e-11, maxit=30)
print("ok", ok, "A", info['A'], "m", info['m'], "time", time.time() - t)
np.save(f"hl_q_lam{lam:.4f}.npy", q)
