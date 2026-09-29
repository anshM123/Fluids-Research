"""Re-converge saved branch states at an exactly specified λ.
The branch scans save states under file names with λ rounded (HL: 6 decimals; 2D: 4 decimals), and a state used with
the rounded λ is inconsistent: m shifts by Δλ/ε and the residual is ~1e-5 (HL) or ~1e-4 (2D). Newton at the stated λ
removes this.
usage: python3 reconverge.py hl  STATE LAM N L1 L2 ETA_START OUTPREFIX
       python3 reconverge.py bq  STATE LAM HS OUTPREFIX"""
import numpy as np, sys
kind, f, lam = sys.argv[1], sys.argv[2], float(sys.argv[3])
if kind == 'hl':
    from hl_solver import HL, newton
    N, L1, L2, e0, pre = int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), sys.argv[8]
    S = HL(lam, N=N, L1=L1, L2=L2, eta_start=e0)
    q0 = S.full(np.load(f)); R0, i0 = S.residual(q0)
    q, info, ok = newton(S, q0, tol=1e-12, maxit=20, verbose=False)
    np.save(f"{pre}_lam{lam:.8f}.npy", q)
    print(f"HL λ={lam:.8f}: residual {np.abs(R0).max():.1e} → converged={ok}, m−2 {i0['m']-2:+.3e} → {info['m']-2:+.3e}", flush=True)
else:
    from bq_solver import BQ
    from bq_newton import newton
    hs, pre = float(sys.argv[4]), sys.argv[5]
    B = BQ(lam, hs=hs, Nb=32, s_sw=12.0)
    Y, info, ok = newton(B, np.load(f), tol=1e-10, maxit=16, verbose=False, fd='central', pert=1e-6)
    np.save(f"{pre}_lam{lam:.6f}.npy", Y)
    print(f"2D λ={lam:.6f}: converged={ok}, m−2 = {info['m']-2:+.4e}", flush=True)
