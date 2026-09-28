import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mapped_continuation import make_grid, layer_info
from ccf_mapped import CCFMapped
d = np.load('pdiag_start.npy'); eta0, phi0, lam = d[0], d[1], float(d[2][0])
M = make_grid(-0.963, 0.0015, hs=0.03, pts=24)
phi = M.from_other(eta0, phi0); target = M.normval(phi, lam)
import scipy.linalg as sla
from dense_arclength import jac
for it in range(20):
    F, den = M.F(phi, lam, target)
    if np.abs(F).max() < 1e-11: break
    J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
F, den = M.F(phi, lam, target)
b = lam/den - M.c
# local minima of den in |eta|<15
m = (np.abs(M.eta) < 15)
idx = np.where((den[1:-1] < den[:-2]) & (den[1:-1] < den[2:]))[0] + 1
idx = [i for i in idx if abs(M.eta[i]) < 15]
print("local minima of den:", [(round(M.eta[i],4), float('%.4e'%den[i])) for i in idx])
idx2 = np.where((den[1:-1] > den[:-2]) & (den[1:-1] > den[2:]))[0] + 1
print("local maxima of den:", [(round(M.eta[i],4), float('%.4e'%den[i])) for i in idx2 if abs(M.eta[i])<15])
Th = np.exp(phi + M.c*M.eta)
fig, ax = plt.subplots(1,2, figsize=(12,4))
mm = (M.eta > -3) & (M.eta < 1)
ax[0].plot(M.eta[mm], den[mm], '.-', ms=2); ax[0].set_ylabel('sonic factor'); ax[0].set_yscale('log')
ax[1].plot(M.eta[mm], np.log(Th[mm]), '.-', ms=2); ax[1].set_ylabel('ln Theta')
for a in ax: a.set_xlabel('eta = ln xi')
plt.tight_layout(); plt.savefig('fig_inspect_l0458.png', dpi=110)
# spacing near minimum
j = np.argmin(den); print("grid spacing at layer:", M.eta[j+1]-M.eta[j], " den curvature width:", layer_info(M, den)[2])
