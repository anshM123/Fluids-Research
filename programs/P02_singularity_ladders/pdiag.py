import numpy as np, time
from ccf_mapped import CCFMapped
from mapped_continuation import layer_info, make_grid
d = np.load('mapped_l2_hs0.02_eps0.025.npy'); lam = 0.4713242277712
M = make_grid(-0.9977, 0.009, hs=0.03)
phi = M.from_other(d[0], d[1]); phi, lam, ok = M.solve_p(phi, lam, 2.0)
target = M.normval(phi, lam)
t0=time.time()
for lam_new in np.arange(0.470, 0.4575, -0.001):
    phi, ok = M.solve_fixed(phi, lam_new, phi0=target); lam = lam_new
    F, den = M.F(phi, lam, target); e, md, w = layer_info(M, den)
    if w / (M.hs * M.eps) < 15:
        Mn = make_grid(e, w, hs=0.03); phi = Mn.from_other(M.eta, phi); M = Mn
        phi, ok = M.solve_fixed(phi, lam, phi0=target)
print("natural done", lam, M.p_of(phi, lam), time.time()-t0, flush=True)
np.save("pdiag_start.npy", np.vstack([M.eta, phi, np.full(M.N, lam)]))
# one p-step with verbose + GMRES info
p0 = M.p_of(phi, lam)
from scipy.sparse.linalg import LinearOperator, gmres
N = M.N
F, den = M.F(phi, lam, target)
Fl = M.Fl(lam, den); dh = M.dh1(phi); pt = p0 - 1e-4
g = M.h1(phi) - (lam/pt - 1 - lam)
sc = 1.0/max(np.linalg.norm(dh), abs(1/pt-1))
cnt=[0]
def mv(x):
    cnt[0]+=1; v, s = x[:N], x[N]
    return np.concatenate([M.Jv(v, phi, lam, den) + s*Fl, [sc*(dh@v - (1/pt-1)*s)]])
A = LinearOperator((N+1,N+1), matvec=mv, dtype=float)
t1=time.time(); x, info = gmres(A, -np.concatenate([F,[sc*g]]), rtol=1e-11, atol=0, restart=120, maxiter=40)
print("gmres info", info, "iters", cnt[0], "time", time.time()-t1, "dlam", x[N], flush=True)
ph2, l2, ok2 = M.solve_p(phi.copy(), lam, pt, phi0=target, verbose=True)
print("solve_p ok", ok2, l2, flush=True)
