import numpy as np, time
from scipy.sparse.linalg import LinearOperator, gmres
from ccf_nk import CCFNK
S = CCFNK(30, 120, 8192, 0.7)
phi = S.guess(0.6, 2.0); phi, ok = S.picard(phi, 0.6, tol=1e-6)
F, den = S.F(phi, 0.6, phi[S.i0])
for tol in [1e-6, 1e-10]:
    cnt = [0]
    A = LinearOperator((S.N, S.N), matvec=lambda v: S.Jv(v, phi, 0.6, den), dtype=float)
    t0 = time.time()
    x, info = gmres(A, -F, rtol=tol, atol=0, restart=400, maxiter=10, callback=lambda r: cnt.__setitem__(0, cnt[0]+1), callback_type='pr_norm')
    print(f"tol={tol}: info={info} iterations={cnt[0]} time={time.time()-t0:.2f}s  resid={np.linalg.norm(A.matvec(x)+F)/np.linalg.norm(F):.2e}")
