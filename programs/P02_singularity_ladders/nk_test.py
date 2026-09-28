import numpy as np, time
from ccf_nk import CCFNK
for N in [4096, 16384]:
    t0 = time.time()
    S = CCFNK(30, 120, N, 0.7)
    phi = S.guess(0.6, 2.0)
    phi, ok = S.picard(phi, 0.6, tol=1e-6)
    phi, ok = S.solve_fixed(phi, 0.6)
    phi, lam, ok2 = S.solve_p(phi, 0.6, 2.0, verbose=False)
    print(f"N={N}: fixed ok={ok} p2 ok={ok2} lambda1={lam:.14f}  time {time.time()-t0:.1f}s", flush=True)
