"""Grid check for a deep Hou–Luo branch state: transfer an F state (N = 65536 on η ∈ [−25, 75]) to N1 points on the
same domain and re-solve at the same λ.
usage: python3 hl_refine_F.py STATE N1"""
import numpy as np, sys, re
from hl_solver import HL, newton
f, N1 = sys.argv[1], int(sys.argv[2])
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
S0 = HL(lam, N=65536, L1=25, L2=75, eta_start=-20); S1 = HL(lam, N=N1, L1=25, L2=75, eta_start=-20)
q = np.interp(S1.eta, S0.eta, S0.full(np.load(f)))
q, info, ok = newton(S1, q, tol=1e-12, maxit=20, verbose=False)
np.save(f"hl_q_G{N1}_lam{lam:.6f}.npy", q)
print(f"λ={lam:.6f} N={N1}: ok={ok} m={info['m']:.10f} (N=65536 state: m−2 as in its file)", flush=True)
