import numpy as np, sys, time
from bq_solver import BQ
from bq_newton import newton
lam = float(sys.argv[1]); f = sys.argv[2]
B = BQ(lam, s_sw=12.0)
Y = np.load(f)
Y, info, ok = newton(B, Y, tol=1e-10, maxit=10, verbose=True, fd='central', pert=1e-6)
print("ok", ok, info['A'], info['m'])
