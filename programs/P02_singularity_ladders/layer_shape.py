"""Shape of the sonic factor den(η) near its minimum: parabolic (δ + a x²) vs hyperbolic (√(δ²+k²x²)) vs asymmetric."""
import numpy as np, sys
from ccf_sinhgrid import make_sinh_grid
from mapped_continuation import layer_info
import scipy.linalg as sla
from dense_arclength import jac
f = sys.argv[1]
d = np.load(f); eta0, phi0, lam = d[0], d[1], float(d[2][0])
# use the saved grid directly: rebuild sinh grid at the saved layer and interpolate (cubic)
from scipy.interpolate import CubicSpline
from mapped_continuation import make_grid
M0 = make_grid(-0.96, 0.0015, hs=0.03, pts=24)
target = M0.normval(M0.from_other(eta0, phi0), lam)
F, den = M0.F(M0.from_other(eta0, phi0), lam, target)
e, md, w = layer_info(M0, den)
M = make_sinh_grid(e, w, hs=0.03, pts=24)
phi = CubicSpline(eta0, phi0)(M.eta)
for it in range(20):
    F, den = M.F(phi, lam, target)
    if np.abs(F).max() < 1e-11: break
    J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
F, den = M.F(phi, lam, target)
j = int(np.argmin(den)); x = M.eta - M.eta[j]; dm = den[j]
print(f"lam={lam:.10f} p={M.p_of(phi, lam):.10f} delta={dm:.5e} at eta={M.eta[j]:.6f} |F|={np.abs(F).max():.1e}")
for X in [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256]:
    xx = X * w
    dl = np.interp(-xx, x, den); dr = np.interp(xx, x, den)
    print(f"  x=±{X:6.2f}w ({xx:.2e}): den_left/δ={dl/dm:10.4f} den_right/δ={dr/dm:10.4f}   parab 1+X²={1+X*X:10.2f}  hyperb √(1+X²·2)={np.sqrt(1+2*X*X):8.3f}")
# jump of ln Θ across the layer and the profile Θ around it
Th = np.exp(phi + M.c * M.eta)
for X in [-64, -16, -4, 0, 4, 16, 64]:
    k = int(np.argmin(np.abs(x - X * w)))
    print(f"  x={X:4d}w: lnΘ={np.log(Th[k]):.5f}  φ'={(np.gradient(phi, M.eta)[k]):.3f}")
