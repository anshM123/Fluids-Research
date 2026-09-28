"""Large-λ end of the CCF branch by pseudo-arclength continuation (uniform FFT log grid, long domain).
Weight c must satisfy β < c < min(1, p); here c = 0.9, L2 = 500 (Ψ decays like e^{(β−c)η})."""
import numpy as np, sys
from ccf_nk import CCFNK, arclength

c = float(sys.argv[1]) if len(sys.argv) > 1 else 0.9
L2 = float(sys.argv[2]) if len(sys.argv) > 2 else 500.0
N = int(sys.argv[3]) if len(sys.argv) > 3 else 65536
d = np.load("ladder_F16k_up_cross2.npy")
lam = d[0]
eta_old = -30.0 + (150.0 / 16384) * np.arange(16384)
lt_old = d[1:] + 0.7 * eta_old
S = CCFNK(30.0, L2, N, c)
lt = np.interp(S.eta, eta_old, lt_old)
right = S.eta > eta_old[-1]
lt[right] = lt_old[-1] + (lam / (1 + lam)) * (S.eta[right] - eta_old[-1])
phi = lt - c * S.eta
phi, lam, ok = S.solve_p(phi, lam, 2.0)
print(f"λ0 on long grid (c={c}, L2={L2}, N={N}): {lam:.12f} ok={ok}", flush=True)
out, sols = arclength(S, phi, lam, 0.02, 400, dsmax=0.15, verbose=True)
np.save(f"large_lambda_arc_c{c}.npy", out)
