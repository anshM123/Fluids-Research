import numpy as np
from ccf_mapped import CCFMapped
from mapped_continuation import make_grid
from dense_arclength import jac, cumint_mat
d = np.load('pdiag_start.npy'); eta0, phi0, lam = d[0], d[1], float(d[2][0])
M = make_grid(-0.963, 0.002, hs=0.03)
phi = M.from_other(eta0, phi0); target = M.normval(phi, lam)
F0, den = M.F(phi, lam, target)
J, Fl = jac(M, phi, lam, den)
rng = np.random.default_rng(0)
for trial in range(3):
    v = rng.standard_normal(M.N) * np.exp(-(M.eta/5)**2)
    e = 1e-7
    F1, _ = M.F(phi + e*v, lam, target)
    fd = (F1 - F0)/e
    an = J @ v
    print("rel err J v:", np.abs(fd-an).max()/np.abs(fd).max())
# cumint consistency
X = rng.standard_normal((M.N, 2))
print("cumint_mat vs cumint:", np.abs(cumint_mat(M, X)[:,0] - M.cumint(X[:,0])).max())
F2, _ = M.F(phi, lam+1e-7, target); print("Fl rel err:", np.abs((F2-F0)/1e-7 - Fl).max()/np.abs(Fl).max())
