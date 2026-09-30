"""Stage 3 (post-outcome holdout, registered before λ7 and λ8 are computed): predictions for λ7, λ8 from the fine-grid
ladder λ0–λ6 by fixed rules, and the successive-subset estimates of the accumulation point λ_c.
Rules (fixed here, no refitting after λ7 is known):
  H1  shifted linear law 1/(λ_n − λ_c) = a n + b, least squares on λ2–λ6 (the WKB regime, D0 ≤ 0.16);
  H1' the same on λ0–λ6;
  H2  geometric contraction of the spacings in z = 1/λ, ratio from the last three spacings;
  H3  3-point λ_c from (λ4, λ5, λ6) with the corresponding a, b (the deepest local estimate)."""
import numpy as np
from scipy.optimize import least_squares, brentq
lam = np.array([1.0285722975, 0.4721297348, 0.3149618108, 0.2415663353, 0.1987224523, 0.1706180880, 0.15092])
err = np.array([1e-9, 1e-9, 1e-9, 1e-9, 1e-8, 1e-8, 2e-5])
n = np.arange(len(lam))

def fit_shift(idx):
    l, k = lam[idx], n[idx]
    r = lambda p: (1/(l - p[0]) - (p[1]*k + p[2]))
    s = least_squares(r, [0.035, 1.3, 1.0], bounds=([-1, 0, -50], [l.min()-1e-4, 10, 50]))
    return s.x, np.abs(s.fun).max()

def three(l):
    g = lambda c: (1/(l[0]-c) - 1/(l[1]-c)) - (1/(l[1]-c) - 1/(l[2]-c))
    cs = np.linspace(-0.5, l.min()-1e-6, 40001); v = np.array([g(c) for c in cs])
    return [brentq(g, cs[i], cs[i+1]) for i in range(len(cs)-1) if np.sign(v[i]) != np.sign(v[i+1])]

print("successive 3-point λ_c (fine ladder):")
for k in range(len(lam) - 2):
    print(f"  ({k}:{k+2}) λ_c = {three(lam[k:k+3])[0]:.4f}")
print("successive 4-point λ_c:")
for k in range(len(lam) - 3):
    p, e = fit_shift(np.arange(k, k+4)); print(f"  ({k}:{k+3}) λ_c = {p[0]:.4f} a = {p[1]:.4f} max res {e:.1e}")
lo, hi = lam[6] - err[6], lam[6] + err[6]
for d in (-1, 1):
    l = lam.copy(); l[6] += d * err[6]
    print(f"  λ6 {d:+d}σ: 3-pt(4:6) λ_c = {three(l[4:7])[0]:.4f}")

preds = {}
for name, idx in (("H1  shifted law on λ2–λ6", np.arange(2, 7)), ("H1' shifted law on λ0–λ6", np.arange(0, 7))):
    p, e = fit_shift(idx)
    l7 = p[0] + 1/(7*p[1] + p[2]); l8 = p[0] + 1/(8*p[1] + p[2])
    preds[name] = (l7, l8, p)
z = 1/lam; dzs = np.diff(z); rho = (dzs[-1]/dzs[-2] * dzs[-2]/dzs[-3]) ** 0.5
z7 = z[6] + dzs[-1]*rho; z8 = z7 + dzs[-1]*rho**2
preds["H2  geometric (ratio %.4f)" % rho] = (1/z7, 1/z8, None)
c = three(lam[4:7])[0]; a = 1/(lam[5]-c) - 1/(lam[4]-c); b = 1/(lam[4]-c) - 4*a
preds["H3  3-point (λ4, λ5, λ6)"] = (c + 1/(7*a+b), c + 1/(8*a+b), (c, a, b))
print("\npredictions:")
for k, (l7, l8, p) in preds.items():
    extra = "" if p is None else f"   (λ_c = {p[0]:.4f}, a = {p[1]:.4f}, b = {p[2]:.4f})"
    print(f"  {k:32s} λ7 = {l7:.5f} (z7 = {1/l7:.4f})   λ8 = {l8:.5f} (z8 = {1/l8:.4f}){extra}")
