"""m-noise versus origin truncation at fixed λ: converge at λ0, then at λ0 + k·dλ (k = 1..K) from the same state,
for one s_start; the scatter about a straight line measures the noise.  usage: ipm_noise2.py STATE LAM0 DLAM K SS"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton, residual
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, lam0, dl, K, ss = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
t0 = time.time(); rows = []
Y = transfer(BQLogPolar(lam0, hs=0.0125, Nb=32, s_start=-20.0), np.load(f), IPM(lam0, hs=0.0125, Nb=32, s_sw=12.0, s_start=ss))
for k in range(K + 1):
    lam = lam0 + k * dl
    B = IPM(lam, hs=0.0125, Nb=32, s_sw=12.0, s_start=ss)
    Yk, info, ok = newton(B, Y, tol=1e-13, maxit=10, verbose=False, fd='central', pert=1e-6)
    R, _ = residual(B, Yk)
    if k == 0:
        Y = Yk
    rows.append((lam, info['m'] - 2))
    print(f"ss={ss} λ={lam:.9f} |R|={np.abs(R).max():.1e} m-2={info['m']-2:+.6e} t={time.time()-t0:.0f}s", flush=True)
r = np.array(rows); p = np.polyfit(r[:, 0] - lam0, r[:, 1], 1); res = r[:, 1] - np.polyval(p, r[:, 0] - lam0)
print(f"ss={ss}: slope {p[0]:.3e}, intercept {p[1]:+.4e}, scatter rms {np.sqrt(np.mean(res**2)):.2e}")
