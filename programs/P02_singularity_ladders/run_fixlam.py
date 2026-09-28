import numpy as np, time, sys
from ccf_profile import CCFProfile
N = int(sys.argv[1]); lams = [float(x) for x in sys.argv[2].split(',')]
prof = CCFProfile(L1=30, L2=150, N=N, c=0.75)
Psi = prof.initial_guess(lams[0], 2)
for lam in lams:
    Psi, l, s, ok = prof.solve(Psi, lam, 2, fix_lam=True, maxit=40)
    h1 = prof.h1(Psi)
    print(f"N={N} lam={lam:.4f}: ok={ok} s={s:.3e} h1={h1:.10f} p_loc={lam/(1+lam+h1):.10f}", flush=True)
