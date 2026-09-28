import numpy as np, sys
D = np.load(sys.argv[1] if len(sys.argv) > 1 else "hl_scan_D.npy")
lam, A, m, Dmin = D.T
z = 1 / (lam - 1); f = m - 2
o = np.argsort(z); z, f = z[o], f[o]
cz = []
for i in range(len(z) - 1):
    if f[i] * f[i + 1] < 0:
        j0 = max(0, min(i - 1, len(z) - 4)); p = np.polyfit(z[j0:j0 + 4] - z[i], f[j0:j0 + 4], 3)
        r = np.roots(p); r = r[np.isreal(r)].real + z[i]; r = r[(r >= z[i]) & (r <= z[i + 1])]
        cz.append(r[0] if len(r) else z[i] - f[i] * (z[i + 1] - z[i]) / (f[i + 1] - f[i]))
print("crossings z:", np.round(cz, 5)); print("spacings:", np.round(np.diff(cz), 5))
print("λ_n:", np.round(1 + 1 / np.array(cz), 7))
ez = []
for i in range(1, len(z) - 1):
    if (f[i] - f[i - 1]) * (f[i + 1] - f[i]) < 0:
        p = np.polyfit(z[i - 1:i + 2] - z[i], f[i - 1:i + 2], 2); zx = -p[1] / (2 * p[0])
        ez.append((zx + z[i], np.polyval(p, zx)))
ez = np.array(ez)
for zz, ff in ez: print(f"  extremum z={zz:.4f} m−2={ff:+.4e}")
if len(ez) > 1:
    print("ratios:", np.round(np.abs(ez[:-1, 1] / ez[1:, 1]), 4), " spacings:", np.round(np.diff(ez[:, 0]), 4))
