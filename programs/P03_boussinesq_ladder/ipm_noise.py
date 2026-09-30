"""Noise floor of m(λ) on a given grid: is the scatter of m at deep rungs set by the Newton tolerance or by the
discretization?  (1) solve at λ0 to tol 1e-11 and then continue Newton as far as the residual floor allows;
(2) solve at λ0 + k·dλ (k = 1..K) to the tightest tolerance and test m(λ) for linearity on that tiny interval.
usage: ipm_noise.py STATE LAM0 DLAM K [Nb hs NB_IN HS_IN]"""
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton, residual
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, lam0, dl, K = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
Nb = int(sys.argv[5]) if len(sys.argv) > 5 else 32
hs = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0125
nbi = int(sys.argv[7]) if len(sys.argv) > 7 else 32
hsi = float(sys.argv[8]) if len(sys.argv) > 8 else 0.025
t0 = time.time()
B = IPM(lam0, hs=hs, Nb=Nb, s_sw=12.0)
Y = np.load(f)
if (nbi, hsi) != (Nb, hs):
    Y = transfer(BQLogPolar(lam0, hs=hsi, Nb=nbi), Y, B)
Y, info, ok = newton(B, Y, tol=1e-11, maxit=12, verbose=True, t0=t0, fd='central', pert=1e-6)
m11 = info['m']
print(f"tol 1e-11: m-2 = {m11-2:+.4e} A = {info['A']:.14f}", flush=True)
for tol in (1e-12, 1e-13, 3e-14):
    Y2, info2, ok2 = newton(B, Y, tol=tol, maxit=6, verbose=True, t0=t0, fd='central', pert=1e-6)
    print(f"tol {tol:.0e}: ok={ok2} m-2 = {info2['m']-2:+.4e}  change vs tol 1e-11: {info2['m']-m11:+.2e}", flush=True)
    if ok2:
        Y, info = Y2, info2
    else:
        break
np.save(f"ipm_noise_lam{lam0:.7f}_hs{hs}.npy", Y)
rows = [(lam0, info['m'])]
for k in range(1, K + 1):
    lam = lam0 + k * dl
    Bk = IPM(lam, hs=hs, Nb=Nb, s_sw=12.0)
    Yk, ik, okk = newton(Bk, Y, tol=3e-13, maxit=10, verbose=False, t0=t0, fd='central', pert=1e-6)
    R, _ = residual(Bk, Yk)
    rows.append((lam, ik['m']))
    print(f"λ={lam:.10f} ok={okk} |R|={np.abs(R).max():.1e} m-2={ik['m']-2:+.6e} t={time.time()-t0:.0f}s", flush=True)
r = np.array(rows); p = np.polyfit(r[:, 0] - lam0, r[:, 1], 2)
res = r[:, 1] - np.polyval(p, r[:, 0] - lam0)
print(f"quadratic fit of m(λ) over {K*dl:.1e}: slope {p[1]:.4e}, residual rms {np.sqrt(np.mean(res**2)):.2e}, max {np.abs(res).max():.2e}")
