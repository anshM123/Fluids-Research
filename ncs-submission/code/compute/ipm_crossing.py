"""IPM: locate a smooth profile (m(λ) = λ/(1+λ−A) = 2) by secant/Illinois iteration in λ with the Newton–Krylov
solver at each λ (the IPM version of bq_crossing.py).
Usage: ipm_crossing.py Y_a.npy lam_a Y_b.npy lam_b [Nb hs s_start s_max] [NB_IN HS_IN]"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
from bq_solver import BQ

fa, la, fb, lb = sys.argv[1], float(sys.argv[2]), sys.argv[3], float(sys.argv[4])
Nb = int(sys.argv[5]) if len(sys.argv) > 5 else 32
hs = float(sys.argv[6]) if len(sys.argv) > 6 else 0.025
sst = float(sys.argv[7]) if len(sys.argv) > 7 else -20.0
smax = float(sys.argv[8]) if len(sys.argv) > 8 else 100.0
nbi = int(sys.argv[9]) if len(sys.argv) > 9 else 32
hsi = float(sys.argv[10]) if len(sys.argv) > 10 else 0.025
tag = f"Nb{Nb}_hs{hs}_ss{sst}_sm{smax}"
t0 = time.time()

def solve(lam, Yg, from_grid=None):
    B = IPM(lam, hs=hs, Nb=Nb, s_start=sst, s_max=smax, s_sw=12.0)
    if from_grid is not None and ((nbi, hsi) != (Nb, hs) or Yg.shape != (B.Ns - B.i0, B.Nb + 1)):
        Yg = transfer(from_grid, Yg, B)
    Y, info, ok = newton(B, Yg, tol=float(os.environ.get('NK_TOL', 1e-11)), maxit=10, verbose=False, fd='central', pert=1e-6)
    print(f"  λ={lam:.10f}: ok={ok} A={info['A']:.12f} m={info['m']:.12f} vrmin={info['vrmin']:.5f} t={time.time()-t0:.0f}s", flush=True)
    return Y, info['m']

import os
_ssi, _smi = float(os.environ.get('SS_IN', -20.0)), float(os.environ.get('SM_IN', 100.0))
gin = lambda lam: BQLogPolar(lam, hs=hsi, Nb=nbi, s_start=_ssi, s_max=_smi)
Ya, ma = solve(la, np.load(fa), gin(la))
Yb, mb = solve(lb, np.load(fb), gin(lb))
bracket = (ma - 2) * (mb - 2) < 0
states = {la: Ya, lb: Yb}
for it in range(12):
    if bracket:                                   # Illinois (safeguarded regula falsi)
        lc = lb - (mb - 2) * (lb - la) / (mb - ma)
        Yg = states[la] + (states[lb] - states[la]) * (lc - la) / (lb - la)
        Yc, mc = solve(lc, Yg)
        if (mc - 2) * (mb - 2) < 0:
            la, ma = lb, mb
        else:
            ma = 2 + (ma - 2) / 2
        states[lc] = Yc
        lb, mb = lc, mc
        states = {la: states[la], lb: states[lb]}
    else:
        lc = lb + (2 - mb) * (lb - la) / (mb - ma)
        Yc = Yb + (Yb - Ya) * (lc - lb) / (lb - la)
        Yc, mc = solve(lc, Yc)
        la, ma, Ya, lb, mb, Yb = lb, mb, Yb, lc, mc, Yc
    if abs(mc - 2) < 1e-12:
        break
Yb = states[lb] if bracket else Yb
print(f"CROSSING {tag}: λ = {lb:.10f} (m−2 = {mb-2:+.2e})", flush=True)
np.save(f"ipm_rung_{tag}_lam{lb:.8f}.npy", Yb)
