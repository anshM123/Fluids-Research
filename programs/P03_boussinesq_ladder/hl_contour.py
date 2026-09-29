"""Hou–Luo: argument-principle count of zeros of det(I − T_μ) in [x_lo, x_hi] × [−y_hi, y_hi] (upper half traversed,
conjugate symmetry), top-k eigenvalues, adaptive refinement to |Δarg| < π/6.
usage: python3 hl_contour.py STATE N X_LO X_HI Y_HI [K=20] [ETA_START=-30]"""
import numpy as np, sys, re, time
from hl_stability import HLStab
f, N = sys.argv[1], int(sys.argv[2])
x_lo, x_hi, y_hi = float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
k = int(sys.argv[6]) if len(sys.argv) > 6 else 20
e0 = float(sys.argv[7]) if len(sys.argv) > 7 else -30.0
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
St = HLStab(lam, np.load(f), N, eta_start=e0)
cache = {}
def fval(mu):
    key = (round(mu.real, 12), round(mu.imag, 12))
    if key not in cache:
        vals, _ = St.spectrum(mu, k=k)
        cache[key] = (np.prod(1 - vals), vals)
    return cache[key]
unres = []
def seg(a, b, fa, fb, depth=0):
    d = np.angle(fb / fa)
    if abs(d) <= np.pi / 6 or depth >= 10:
        if abs(d) > np.pi / 6:
            unres.append((a, b, d, np.min(np.abs(fval(b)[1] - 1))))
        return d
    c = 0.5 * (a + b); fc = fval(c)[0]
    return seg(a, c, fa, fc, depth + 1) + seg(c, b, fc, fb, depth + 1)
corners = [complex(x_hi, 0), complex(x_hi, y_hi), complex(x_lo, y_hi), complex(x_lo, 0)]
total = 0.0
t = time.time()
for a, b in zip(corners[:-1], corners[1:]):
    n = max(2, int(np.ceil(abs(b - a) / 0.05)))
    pts = [a + (b - a) * s for s in np.linspace(0, 1, n + 1)]
    for p, q in zip(pts[:-1], pts[1:]):
        total += seg(p, q, fval(p)[0], fval(q)[0])
Nz = 2 * total / (2 * np.pi)
kmax = max(abs(v[1][-1]) for v in cache.values())
print(f"HL λ={lam:.8f} z={1/(lam-1):.3f} N={N} η0={e0} R=[{x_lo},{x_hi}]x[±{y_hi}] k={k}: zeros = {Nz:.4f} "
      f"(max|ν_k|={kmax:.3f}, unresolved={len(unres)}{', min|1−ν| at them: ' + str([round(u[3],3) for u in unres]) if unres else ''}; "
      f"{len(cache)} evals, {time.time()-t:.0f}s)", flush=True)
