import numpy as np, sys, re, time
from hl_solver import HL, newton
f = sys.argv[1]; lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
q8 = np.load(f); S8 = HL(lam, N=8192)
for N in (8192, 16384, 32768, 65536):
    S = HL(lam, N=N)
    q0 = np.interp(S.eta, S8.eta, q8)
    t = time.time()
    q, info, ok = newton(S, q0, tol=1e-12, maxit=20, verbose=False)
    print(f"λ={lam:.6f} z={1/(lam-1):.3f} N={N}: ok={ok} m−2={info['m']-2:+.10e} A={info['A']:.12f}  ({time.time()-t:.0f}s)", flush=True)
