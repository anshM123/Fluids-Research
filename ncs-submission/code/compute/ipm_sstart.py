"""Origin-truncation study: m at fixed λ for several s_start (the unknowns start at s_start; below it the
constant-strain local structure is imposed, error O(e^{s_start}) for m = 2).  Also reports where the residual
floor sits.  usage: ipm_sstart.py STATE LAM HS SS1,SS2,... [NB_IN HS_IN SS_IN]"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton, residual
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, lam, hs = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
sss = [float(v) for v in sys.argv[4].split(',')]
nbi = int(sys.argv[5]) if len(sys.argv) > 5 else 32
hsi = float(sys.argv[6]) if len(sys.argv) > 6 else hs
ssi = float(sys.argv[7]) if len(sys.argv) > 7 else -20.0
Y0 = np.load(f); t0 = time.time()
for ss in sss:
    B = IPM(lam, hs=hs, Nb=32, s_sw=12.0, s_start=ss)
    Y = transfer(BQLogPolar(lam, hs=hsi, Nb=nbi, s_start=ssi), Y0, B)
    Y, info, ok = newton(B, Y, tol=1e-12, maxit=10, verbose=False, fd='central', pert=1e-6)
    R, _ = residual(B, Y); iR = np.unravel_index(np.argmax(np.abs(R)), R.shape)
    print(f"s_start={ss}: ok={ok} |R|={np.abs(R).max():.2e} at s={B.s[B.i0 + iR[0]]:.2f}, j={iR[1]}  "
          f"A={info['A']:.14f} m-2={info['m']-2:+.6e}  t={time.time()-t0:.0f}s", flush=True)
