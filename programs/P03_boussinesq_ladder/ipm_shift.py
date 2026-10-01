"""Grid-shift averaging against front–grid locking: solve at fixed λ on the log-polar grid shifted by a fraction
δ of a cell (s_min → s_min + δ h_s); the locking error is periodic in the front's sub-cell position, so averaging m
over δ = 0, ¼, ½, ¾ cancels its first three harmonics.  usage: ipm_shift.py STATE Z DELTA [hs] (STATE: unshifted
grid at the same hs, s_start −20, Nb 32)"""
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton, residual
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, z, dl = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
hs = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0125
lam = 1 / z; t0 = time.time()
B = IPM(lam, hs=hs, Nb=32, s_sw=12.0, s_start=-20.0, s_min=-120.0 + dl * hs)
import os
hs_in = float(os.environ.get('HS_IN', hs))                       # source grid (unshifted) may differ from the target
Y = transfer(BQLogPolar(lam, hs=hs_in, Nb=32, s_start=-20.0), np.load(f), B)
Y, info, ok = newton(B, Y, tol=7e-13, maxit=12, verbose=False, fd='central', pert=1e-6)
R, _ = residual(B, Y)
print(f"z={z:.4f} shift={dl:.3f} hs={hs} src={os.path.basename(f)}: |R|={np.abs(R).max():.1e} m-2={info['m']-2:+.6e} t={time.time()-t0:.0f}s", flush=True)
np.save(f"ipm_shift_z{z:.4f}_d{dl:.3f}_hs{hs}.npy", Y)
