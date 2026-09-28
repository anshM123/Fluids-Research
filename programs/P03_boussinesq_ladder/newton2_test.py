"""Newton–Krylov with the improved solver (bq_solver.BQ) at λ = 1.92 from the old converged state."""
import numpy as np, sys, time
from bq_solver import BQ
from bq_newton import newton
lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.92
Nb = int(sys.argv[2]) if len(sys.argv) > 2 else 32
hs = float(sys.argv[3]) if len(sys.argv) > 3 else 0.025
B = BQ(lam, hs=hs, Nb=Nb)
Y = np.load('bq_X_lam1.9200.npy')
Y, info, ok = newton(B, Y, tol=1e-11, maxit=12, verbose=True)
print(f"RESULT lam={lam} Nb={Nb} hs={hs}: ok={ok} A={info['A']:.12f} m={info['m']:.12f}", flush=True)
np.save(f"Y2_lam{lam:.4f}_Nb{Nb}_hs{hs}.npy", Y)
