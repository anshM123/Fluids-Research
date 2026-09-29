import numpy as np, time, sys
from bq_stability import BQStab, Tc, spectrum_c
f = sys.argv[1]; lam = float(sys.argv[2]); hs = float(sys.argv[3])
S = BQStab(lam, np.load(f), hs=hs)
rng = np.random.default_rng(0); v = rng.standard_normal(S.n_unk)
t = time.time(); a = S.T(v, 0.5); t1 = time.time() - t
t = time.time(); b = Tc(S, v.astype(complex), 0.5 + 0j); t2 = time.time() - t
t = time.time(); b = Tc(S, v.astype(complex), 0.5 + 0j); t3 = time.time() - t
print("real vs complex T at real μ: rel diff", np.abs(a - b).max() / np.abs(a).max(), f"times {t1:.2f} {t2:.2f} {t3:.2f}", flush=True)
w = rng.standard_normal(S.n_unk)
lin = Tc(S, v + 1j * w, 0.4 + 0.3j) - (Tc(S, v, 0.4 + 0.3j) + 1j * Tc(S, w, 0.4 + 0.3j))
print("complex linearity check", np.abs(lin).max(), flush=True)
for mu in (0.5, 0.5 + 0.2j):
    t = time.time(); vals, _ = spectrum_c(S, mu, k=8); print(mu, np.round(vals, 4), f"{time.time()-t:.0f}s", flush=True)
