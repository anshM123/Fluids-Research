"""Quantitative check of the square-root cusp law on the most singular computed profile:
   den ≈ δ + k|η−η_s|^{1/2} on both flanks (w ≪ |η−η_s| ≪ 1) with the predicted k = √(2λΘ_s/ξ_s)
(from B² = 2λΘ_s), and ln Θ − ln Θ_s ≈ ±(2λ/k)|η−η_s|^{1/2}.  Saves cusp_layer_profile.npy for the figure."""
import numpy as np, sys
from scipy.interpolate import CubicSpline
from mapped_continuation import layer_info
from pinned_delta import Fpin, solve, new_grid

f = sys.argv[1] if len(sys.argv) > 1 else "pin_C_last.npy"
d = np.load(f); eta0, phi0, lam = d[0], d[1], float(d[2][0])
sp = np.diff(eta0); j = int(np.argmin(sp))
etac = eta0[j] if sp[j - 1] > sp[j + 1] else eta0[j + 1]
etac = eta0[int(np.argmin(np.abs(eta0 - etac)))]
M, ic = new_grid(etac, 24 * sp.min(), 24)
phi = CubicSpline(eta0, phi0)(M.eta)
F, den = Fpin(M, phi, lam, ic)
phi, lam, ok, hist, _ = solve(M, phi, lam, ic, den[ic])
F, den = Fpin(M, phi, lam, ic)
e, delta, w = layer_info(M, den)
Th = np.exp(phi + M.c * M.eta)
Ths, xis = Th[ic], np.exp(M.eta[ic])
k_pred = np.sqrt(2 * lam * Ths / xis)
x = M.eta - M.eta[ic]
print(f"re-solve ok={ok}: λ={lam:.12f} δ={delta:.4e} w={w:.3e} (w/δ²={w/delta**2:.3f}) Θ_s={Ths:.6f} ξ_s={xis:.6f}")
print(f"predicted k = √(2λΘ_s/ξ_s) = {k_pred:.6f}")
rows = []
for X in np.logspace(0, 4, 41):
    xx = X * w
    dl, dr = np.interp(-xx, x, den), np.interp(xx, x, den)
    rows.append((X, dl / delta, dr / delta))
np.save("cusp_layer_profile.npy", np.array(rows))
print("  |x|        k_left(x)=(den−δ)/√|x|   k_right(x)     [prediction k_pred]")
for xx in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
    if xx < 5 * w:
        continue
    dl, dr = np.interp(-xx, x, den), np.interp(xx, x, den)
    print(f"  {xx:.0e}    {(dl-delta)/np.sqrt(xx):.6f}              {(dr-delta)/np.sqrt(xx):.6f}")
# symmetric part isolates the cusp (the regular part of den has a linear term a1 x that is odd)
print("  symmetric part [(den(x)+den(−x))/2 − δ]/√|x|:")
for xx in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3):
    if xx < 5 * w:
        continue
    dl, dr = np.interp(-xx, x, den), np.interp(xx, x, den)
    print(f"  {xx:.0e}    {((dl+dr)/2-delta)/np.sqrt(xx):.6f}   (ratio to k_pred {((dl+dr)/2-delta)/np.sqrt(xx)/k_pred:.4f})")
lT = np.log(Th)
print("  ln Θ jump: antisymmetric part [lnΘ(x) − lnΘ(−x)]/2 / √|x| vs 2λ/k_pred =", f"{2*lam/k_pred:.6f}")
for xx in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3):
    if xx < 5 * w:
        continue
    a = (np.interp(xx, x, lT) - np.interp(-xx, x, lT)) / 2 / np.sqrt(xx)
    print(f"  {xx:.0e}    {a:.6f}")
