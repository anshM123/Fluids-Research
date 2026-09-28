"""Resolution check of m(λ) at fixed λ with the production solver: transfer a (Nb=32, hs=0.025) state to (Nb, hs)."""
import numpy as np, sys, time
from bq_solver import BQ
from bq_newton import newton
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, lam = sys.argv[1], float(sys.argv[2])
for spec in sys.argv[3:]:
    Nb, hs = spec.split(':'); Nb = int(Nb); hs = float(hs)
    B = BQ(lam, hs=hs, Nb=Nb)
    Y = transfer(BQLogPolar(lam, hs=0.025, Nb=32), np.load(f), B)
    t = time.time()
    Y, info, ok = newton(B, Y, tol=1e-11, maxit=10, verbose=False, fd='central', pert=1e-6)
    print(f"λ={lam} Nb={Nb} hs={hs}: ok={ok} m−2={info['m']-2:+.10e} A={info['A']:.12f} vrmin={info['vrmin']:.5f} ({time.time()-t:.0f}s)", flush=True)
