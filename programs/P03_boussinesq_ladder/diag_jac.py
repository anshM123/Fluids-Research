import numpy as np
from bq_logpolar import BQLogPolar
from bq_newton import residual, initial_guess
lam = 1.92
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
X = initial_guess(B, (3 + lam) / 2)
R, info = residual(B, X)
eps = 1e-8
rng = np.random.default_rng(1)
for lo, hi in ((-120, -60), (-60, -21), (-21, -18), (-18, -10), (-10, -3), (-3, 3), (3, 10), (10, 40), (40, 100)):
    v = np.zeros(X.shape)
    sel = (B.s >= lo) & (B.s < hi)
    v[sel] = rng.standard_normal((sel.sum(), X.shape[1]))
    v /= np.linalg.norm(v)
    Rp, ip = residual(B, X + eps * v)
    Jv = (Rp - R) / eps
    k = np.unravel_index(np.argmax(np.abs(Jv)), Jv.shape)
    print(f"band s∈[{lo},{hi}): |Jv|2={np.linalg.norm(Jv):.3e}  max at s={B.s[k[0]]:.2f}, β-index {k[1]}; ΔA={(ip['A']-info['A'])/eps:.3e}", flush=True)
