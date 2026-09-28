"""Radial-velocity dip along the boundary: V_r/r(r, β=0) = 1+λ+U_r/r; location/value of its minimum vs ε = 1+λ−A."""
import numpy as np, sys
from bq_solver import BQ
from bq_newton import full
print(f"{'lam':>8} {'m':>10} {'eps':>8} {'vrmin_b':>8} {'ratio':>7} {'r_dip':>8} {'vrmin_all':>9} {'at r':>7} {'beta':>6}")
for f in sys.argv[1:]:
    lam = float(f.split('lam')[1][:6])
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
    Ur = X @ B.Db.T
    vr = (1 + lam) + Ur
    A = -np.mean(Ur[B.i0:B.i0 + 40, 0]); eps = 1 + lam - A; m = (lam - 1) / eps
    k = np.argmin(vr[:, 0]); i, j = np.unravel_index(np.argmin(vr), vr.shape)
    print(f"{lam:8.4f} {m:10.6f} {eps:8.5f} {vr[k,0]:8.5f} {vr[k,0]/eps:7.4f} {np.exp(B.s[k]):8.4f} {vr[i,j]:9.5f} {np.exp(B.s[i]):7.4f} {B.beta[j]:6.3f}")
