import numpy as np, time, sys
from ccf_newton import CCFNewton
for (L1, L2, N, c) in [(30,150,2048,0.75),(30,150,4096,0.75),(40,200,6144,0.75),(30,150,4096,0.70)]:
    t0=time.time()
    S = CCFNewton(L1, L2, N, c)
    phi = S.guess(0.6, 2.0)
    phi, ok = S.solve_fixed(phi, 0.6)
    phi, lam, ok2 = S.solve_p(phi, 0.6, 2.0)
    print(f"L1={L1} L2={L2} N={N} c={c}: fixed ok={ok}, p=2 ok={ok2}, lambda={lam:.13f}  ({time.time()-t0:.0f}s)", flush=True)
