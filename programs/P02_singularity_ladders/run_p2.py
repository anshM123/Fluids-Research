import numpy as np, time, sys
from ccf_profile import CCFProfile
N = int(sys.argv[1]) if len(sys.argv)>1 else 1024
prof = CCFProfile(L1=30, L2=150, N=N, c=0.75)
for p, lam0 in [(2, 1.0)]:
    Psi0 = prof.initial_guess(lam0, p)
    t0=time.time()
    Psi, lam, s, ok = prof.solve(Psi0, lam0, p, verbose=True)
    print(f"p={p} N={N}: ok={ok} lam={lam:.14f} s={s:.2e} time={time.time()-t0:.1f}s")
    np.save(f"ccf_p{p}_N{N}.npy", np.concatenate([[lam], Psi]))
