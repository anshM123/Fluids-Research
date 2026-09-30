"""Angular-resolution check: m at fixed λ for Nb = 48 (vs 32), same hs and s_start.
usage: ipm_nbcheck.py STATE(Nb32, hs, s_start −20) LAM HS SS NB"""
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton, residual
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, lam, hs, ss, nb = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
t0 = time.time()
B = IPM(lam, hs=hs, Nb=nb, s_sw=12.0, s_start=ss)
Y = transfer(BQLogPolar(lam, hs=hs, Nb=32, s_start=-20.0), np.load(f), B)
Y, info, ok = newton(B, Y, tol=1e-13, maxit=12, verbose=True, fd='central', pert=1e-6)
R, _ = residual(B, Y)
print(f"Nb={nb} hs={hs} s_start={ss} λ={lam}: ok={ok} |R|={np.abs(R).max():.1e} A={info['A']:.14f} m-2={info['m']-2:+.6e} "
      f"vrmin={info['vrmin']:.5f} t={time.time()-t0:.0f}s", flush=True)
np.save(f"ipm_nbcheck_Nb{nb}_lam{lam:.7f}.npy", Y)
