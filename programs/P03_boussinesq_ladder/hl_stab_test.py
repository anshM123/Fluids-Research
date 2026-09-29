import numpy as np, sys, time
from hl_stability import HLStab
f, lam, N = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
t = time.time(); St = HLStab(lam, np.load(f), N); print(f"setup {time.time()-t:.1f}s m={St.m:.10f}", flush=True)
v = St.time_translation().astype(complex)
t = time.time(); Tv = St.T(v, 1.0); print(f"T apply {time.time()-t:.2f}s", flush=True)
t = time.time(); Tv = St.T(v, 1.0); print(f"T apply {time.time()-t:.3f}s")
w = slice(0, int(0.8 * len(v)))
print("time-translation mode: |T_1 v − v|/|v| =", np.linalg.norm((Tv - v)[w]) / np.linalg.norm(v[w]), flush=True)
for mu in (1.0, 0.5, 0.2, 0.1):
    t = time.time(); vals, _ = St.spectrum(mu, k=8)
    print(f"μ={mu}: ν = {np.round(vals, 4)} ({time.time()-t:.1f}s)", flush=True)
