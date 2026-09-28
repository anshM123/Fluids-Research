import numpy as np, sys, time
from bq_solver import BQ
from bq_newton import newton
lam = 1.92
fd = sys.argv[1]; pert = float(sys.argv[2])
B = BQ(lam)
Y = np.load('Y2_lam1.9200_Nb32_hs0.025.npy')
Y, info, ok = newton(B, Y, tol=1e-12, maxit=6, verbose=True, fd=fd, pert=pert)
print(f"RESULT fd={fd} pert={pert}: ok={ok} A={info['A']:.13f} m={info['m']:.13f}", flush=True)
np.save(f"Y3_{fd}_{pert}.npy", Y)
