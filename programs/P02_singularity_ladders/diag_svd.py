"""Smallest singular values/vectors of the fixed-λ Jacobian J and the δ-bordered Jacobian at λ=0.458."""
import numpy as np, time
import scipy.linalg as sla
from mapped_continuation import layer_info, make_grid
from dense_arclength import jac
from ccf_sinhgrid import make_sinh_grid

d = np.load('pdiag_start.npy'); eta0, phi0, lam = d[0], d[1], float(d[2][0])
Mref = make_grid(-0.963, 0.0015, hs=0.03, pts=24)
target = Mref.normval(Mref.from_other(eta0, phi0), lam)
M = make_sinh_grid(-0.9632, 0.00162, hs=0.03, pts=24, B=1.0)
phi = M.from_other(eta0, phi0)
for it in range(25):
    F, den = M.F(phi, lam, target)
    if np.abs(F).max() < 1e-12: break
    J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
F, den = M.F(phi, lam, target)
N = M.N; ip = int(np.argmin(den))
J, Fl = jac(M, phi, lam, den)
t0 = time.time()
U, S, Vt = sla.svd(J, lapack_driver='gesdd')
print("J smallest sv:", S[-6:], "largest", S[0], f"({time.time()-t0:.0f}s)", flush=True)
A = np.zeros((N + 1, N + 1)); A[:N, :N] = J; A[:N, N] = Fl
A[N, :N] = M.E[ip] * M.Hm[ip, :] * np.exp(phi); A[N, N] = 1.0
S2 = sla.svdvals(A)
print("A(delta-bordered) smallest sv:", S2[-6:], flush=True)
# where do the two smallest right singular vectors live?
for k in (1, 2, 3):
    v = Vt[-k]; j = np.argmax(np.abs(v))
    wid = np.sum(np.abs(v) > 0.1 * np.abs(v).max())
    print(f"v_{k}: sigma={S[-k]:.3e} peak at eta={M.eta[j]:.5f} (layer {M.eta[ip]:.5f}), #pts>10% = {wid}, "
          f"Fl-component u.Fl={U[:, -k] @ Fl:.3e}, bordering-row.v={(A[N,:N] @ v):.3e}", flush=True)
np.save("diag_svd_vecs.npy", np.vstack([M.eta, Vt[-1], Vt[-2], Vt[-3], U[:, -1], U[:, -2]]))
