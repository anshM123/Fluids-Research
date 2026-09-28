"""Uniqueness probe: at fixed λ, solve the CCF profile equation from many different initial guesses
(local exponent p0 at the origin, amplitude factor A) and compare the converged local exponent p with the
principal branch value p_branch(λ).  A different p would reveal a solution on another (disconnected) branch."""
import numpy as np, sys, time, json
from ccf_nk import CCFNK

bm = np.load("branch_map.npy")
bm = bm[np.argsort(bm[:, 0])]
S = CCFNK(30.0, 120.0, 16384, 0.7)
lams = [float(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else [0.5, 0.55, 0.6, 0.7, 0.85, 1.0, 1.3, 1.8]
out = []
t0 = time.time()
for lam in lams:
    pb = float(np.interp(lam, bm[:, 0], bm[:, 1]))
    found = []
    for p0 in (1.2, 1.6, 2.0, 2.4, 3.0, 4.0):
        for A in (0.3, 1.0, 3.0):
            phi = S.guess(lam, p0) + np.log(A)
            phi, okp = S.picard(phi, lam, relax=0.3, tol=1e-6, maxit=400)
            if not np.all(np.isfinite(phi)):
                found.append((p0, A, None)); continue
            ph, ok = S.solve_fixed(phi, lam, maxit=40)
            if ok:
                p = S.p_of(ph, lam)
                F, den = S.F(ph, lam, S.normval(ph, lam))
                found.append((p0, A, float(p), float(den.min())))
            else:
                found.append((p0, A, None))
    conv = [f for f in found if f[2] is not None]
    ps = sorted(set(round(f[2], 7) for f in conv))
    print(f"λ={lam:.3f}: branch p={pb:.7f}; {len(conv)}/{len(found)} starts converged; distinct p: {ps}  "
          f"t={time.time()-t0:.0f}s", flush=True)
    out.append(dict(lam=lam, p_branch=pb, converged=len(conv), starts=len(found), distinct_p=ps))
    json.dump(out, open("multistart.json", "w"), indent=1)
