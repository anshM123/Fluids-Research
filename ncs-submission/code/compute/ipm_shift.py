"""Grid-shift averaging against front–grid locking: solve at fixed λ on the log-polar grid shifted by a fraction
δ of a cell (s_min → s_min + δ h_s); the locking error is periodic in the front's sub-cell position, so averaging m
over δ = 0, ¼, ½, ¾ cancels its first three harmonics.  usage: ipm_shift.py STATE Z DELTA [hs] (STATE: unshifted
grid at the same hs, s_start −20, Nb 32)"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, time, os
from ipm_solver import IPM
from bq_newton import newton, residual
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, z, dl = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
hs = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0125
lam = 1 / z; t0 = time.time()
NB = int(os.environ.get('NB', 32)); SS = float(os.environ.get('SS', -20.0)); SMAX = float(os.environ.get('SMAX', 100.0))
B = IPM(lam, hs=hs, Nb=NB, s_sw=12.0, s_start=SS, s_max=SMAX, s_min=-120.0 + dl * hs)
hs_in = float(os.environ.get('HS_IN', hs))                       # source grid (unshifted) may differ from the target
Y = transfer(BQLogPolar(lam, hs=hs_in, Nb=32, s_start=-20.0), np.load(f), B)
tol = float(os.environ.get('NK_TOL', 1.5e-12 if hs < 0.01 else 7e-13))     # the h_s = 0.00625 floor is 5e-13–1.4e-12
Y, info, ok = newton(B, Y, tol=tol, maxit=12, verbose=False, fd='central', pert=1e-6)
R, _ = residual(B, Y)
print(f"z={z:.4f} shift={dl:.3f} hs={hs} Nb={NB} s_start={SS} s_max={SMAX} src={os.path.basename(f)}: |R|={np.abs(R).max():.1e} m-2={info['m']-2:+.6e} t={time.time()-t0:.0f}s", flush=True)
np.save(f"ipm_shift_z{z:.4f}_d{dl:.3f}_hs{hs}_nb{NB}_ss{SS:g}_sm{SMAX:g}.npy", Y)
