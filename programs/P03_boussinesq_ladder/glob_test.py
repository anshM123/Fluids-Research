import numpy as np, sys, time
from bq_global import BQGlobal
f, lam0 = sys.argv[1], float(sys.argv[2])
smin, smax, hs, Nb = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
src_hs = float(sys.argv[7]) if len(sys.argv) > 7 else 0.025
t = time.time()
G = BQGlobal(s_min=smin, s_max=smax, hs=hs, Nb=Nb)
U = G.initial_from_march(f, lam0, hs_src=src_hs)
print(f"setup {time.time()-t:.1f}s, unknowns {len(U)}", flush=True)
t = time.time(); J = G.jacobian(U); print(f"jacobian {time.time()-t:.1f}s nnz={J.nnz}", flush=True)
# finite-difference check of the Jacobian on a random direction
rng = np.random.default_rng(1); v = rng.standard_normal(len(U)) * 1e-3
v[-1] = 1e-4
e = 1e-6
Jv = J @ v; FDv = (G.residual(U + e * v) - G.residual(U - e * v)) / (2 * e)
print("Jacobian check: rel err", np.abs(Jv - FDv).max() / np.abs(FDv).max(), flush=True)
U, ok = G.newton(U, tol=1e-11, maxit=12)
print(f"RESULT global: λ = {U[-1]:.10f} ok={ok}  (march input λ = {lam0})  time {time.time()-t:.0f}s", flush=True)
np.save(f"glob_lam{U[-1]:.6f}_h{hs}_Nb{Nb}.npy", U)
