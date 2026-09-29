"""Transfer a Hou–Luo crossing state to a finer grid (N → 2N or 4N, same domain) and re-solve at the same λ."""
import numpy as np, sys, re
from hl_solver import HL, newton
f, N0, N1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
S0 = HL(lam, N=N0); S1 = HL(lam, N=N1)
q = np.interp(S1.eta, S0.eta, np.load(f))
q, info, ok = newton(S1, q, tol=1e-12, maxit=20, verbose=False)
np.save(f"hl_crossN{N1}_lam{lam:.8f}.npy", q)
print(f"λ={lam:.8f} N={N1}: ok={ok} m−2={info['m']-2:+.2e}", flush=True)
