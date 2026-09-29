import numpy as np, time, sys
from bq_stability import BQStab
f = sys.argv[1]; lam = float(sys.argv[2]); hs = float(sys.argv[3]) if len(sys.argv) > 3 else 0.025
t = time.time()
S = BQStab(lam, np.load(f), hs=hs)
print(f"setup {time.time()-t:.1f}s  m={S.m:.10f} A={S.A:.10f}", flush=True)
v = S.time_translation_mode()
t = time.time(); Tv = S.T(v, 1.0); print(f"T apply {time.time()-t:.2f}s", flush=True)
print("time-translation mode: |T_1 v − v|/|v| =", np.linalg.norm(Tv - v) / np.linalg.norm(v), flush=True)
for mu in (1.0, 0.5):
    t = time.time(); vals, _ = S.spectrum(mu, k=6); print(f"mu={mu}: ν = {np.round(vals, 5)}  ({time.time()-t:.0f}s)", flush=True)
