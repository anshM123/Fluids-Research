"""Extrema and crossings of m(λ) − 2 along the branch (all scan2_*.npy), as functions of z = 1/(λ−1)."""
import numpy as np, glob
rows = [np.load(f) for f in sorted(glob.glob("scan2_*.npy"))]
D = np.vstack([r for r in rows if r.ndim == 2 and len(r)])
D = D[np.argsort(D[:, 0])]
_, iu = np.unique(np.round(D[:, 0], 9), return_index=True); D = D[iu]
lam, A, m, vr = D.T
z = 1 / (lam - 1); f = m - 2
o = np.argsort(z); z, f, lam = z[o], f[o], lam[o]
print("crossings (cubic interpolation in z):")
cz = []
for i in range(len(z) - 1):
    if f[i] * f[i + 1] < 0:
        j0 = max(0, min(i - 1, len(z) - 4)); zz, ff = z[j0:j0 + 4], f[j0:j0 + 4]
        p = np.polyfit(zz - z[i], ff, 3); r = np.roots(p); r = r[np.isreal(r)].real + z[i]
        r = r[(r >= z[i]) & (r <= z[i + 1])]
        zc = r[0] if len(r) else z[i] - f[i] * (z[i + 1] - z[i]) / (f[i + 1] - f[i])
        cz.append(zc)
        print(f"  z = {zc:.5f}   λ = {1 + 1/zc:.7f}")
print("  spacings:", np.round(np.diff(cz), 4))
print("extrema (quadratic fit in z):")
ez = []
for i in range(1, len(z) - 1):
    if (f[i] - f[i - 1]) * (f[i + 1] - f[i]) < 0:
        p = np.polyfit(z[i - 1:i + 2] - z[i], f[i - 1:i + 2], 2)
        zx = -p[1] / (2 * p[0]); fx = np.polyval(p, zx)
        ez.append((zx + z[i], fx))
        print(f"  z = {zx + z[i]:.4f}  λ = {1 + 1/(zx + z[i]):.5f}  m−2 = {fx:+.4e}")
ez = np.array(ez)
if len(ez) > 1:
    print("  |extremum| ratios:", np.round(np.abs(ez[:-1, 1] / ez[1:, 1]), 3), " z-spacings:", np.round(np.diff(ez[:, 0]), 4))
