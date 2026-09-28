"""Checks of bq_solver.BQ against bq_logpolar.BQLogPolar on the converged λ = 1.92 state."""
import numpy as np, time
from bq_logpolar import BQLogPolar
from bq_solver import BQ
from bq_newton import full, residual
lam = 1.92
Bo = BQLogPolar(lam); Bn = BQ(lam)
Y = np.load('bq_X_lam1.9200.npy'); X = full(Bo, Y); P = X / Bo.ea2[:, None]
Vo = Bo.velocity(P); Vn = Bn.velocity(P)
sel = slice(Bo.i0, Bo.Ns - 400)
for name, a, b in zip(['Ur', 'w', 'wt', 'Urh', 'wh', 'wth'], Vo[:6], Vn[:6]):
    k = min(a.shape[0], b.shape[0])
    print(f"{name}: max|old-new| on s in [-20, 90] = {np.abs(a[:k][sel]-b[:k][sel]).max():.2e}  (max|.|={np.abs(a[sel]).max():.2f})")
t = time.time(); ro = Bo.march(P); to = time.time() - t
Bn.march(P); t = time.time(); rn = Bn.march(P); tn = time.time() - t
print(f"march: A {ro['A']:.12f} vs {rn['A']:.12f}; times old {to*1e3:.0f} ms new {tn*1e3:.0f} ms")
for key in ('Th', 'Om'):
    d = np.abs(ro[key] - rn[key]) / (np.abs(ro[key]) + 1e-300)
    print(key, "max rel diff (s in [-20, 60]):", d[Bo.i0:int(np.searchsorted(Bo.s, 60))].max())
Ro, _ = residual(Bo, Y); Rn, _ = residual(Bn, Y)
print("residual old", np.abs(Ro).max(), "new", np.abs(Rn).max())
