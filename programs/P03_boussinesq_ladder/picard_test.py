import numpy as np, time, sys
from bq_logpolar import BQLogPolar
lam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.92
alpha = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
s = B.s[:, None]; b = B.beta[None, :]; r = np.exp(s)
A0 = (3 + lam) / 2
P = np.exp(-B.a*s) * (-(A0/2) * r**2 * np.sin(2*b) / (1 + r**2) ** (1/(2*(1+lam))))
t0 = time.time()
for it in range(200):
    Pn, info = B.T(P)
    d = np.abs(Pn - P).max() / max(np.abs(Pn).max(), 1e-300)
    if it % 5 == 0 or d < 1e-10:
        print(f"it {it}: A={info['A']:.8f} m={info['m']:.8f} min(V_r/r)={info['vrmin']:.4f} |ΔP|/|P|={d:.2e} max|P|={np.abs(Pn).max():.3e} t={time.time()-t0:.1f}s", flush=True)
    if not np.isfinite(d) or d < 1e-10:
        break
    P = (1 - alpha) * P + alpha * Pn
