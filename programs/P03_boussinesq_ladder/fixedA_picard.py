"""Inner problem: Picard iteration for the velocity field at a PRESCRIBED origin strain A (sets m and local data);
reports the strain A_out produced by the converged field.  Outer: A = A_out(A)."""
import numpy as np, sys, time
from bq_logpolar import BQLogPolar
from bq_newton import initial_guess
lam = float(sys.argv[1]); A = float(sys.argv[2]); alpha = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
X = initial_guess(B, A)
t0 = time.time()
for it in range(300):
    Pn, info = B.T(X / B.ea2[:, None], A_fixed=A)
    Xn = B.ea2[:, None] * Pn
    d = np.abs(Xn - X).max()
    Aout = B.strain(B.velocity(Pn)[0])
    if it % 10 == 0 or d < 1e-11:
        print(f"it {it}: |ΔX|={d:.2e} A_in={A:.6f} A_out={Aout:.10f} vrmin={info['vrmin']:.4f} t={time.time()-t0:.0f}s", flush=True)
    if d < 1e-11 or not np.isfinite(d):
        break
    X = (1 - alpha) * X + alpha * Xn
