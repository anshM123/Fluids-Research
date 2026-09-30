"""Rung location from the fixed-z scans: local fits of m−2 around each sign change.
For each scan: a cubic in z through the 4–6 points nearest the sign change (least squares), root and a noise-based
error (scatter of the fit residuals / slope).  The s16 scan is also shown after subtracting the s16−s20 offset
measured at common points near the crossing."""
import numpy as np, sys
def load(tag):
    return np.load(f"ipm_scan_{tag}.npy")          # columns: z, λ, A, m−2, |R|, vrmin
def crossings(z, d, k=6, deg=3):
    out = []
    for i in range(len(z) - 1):
        if d[i] * d[i + 1] < 0:
            lo, hi = max(0, i - k // 2 + 1), min(len(z), i + k // 2 + 1)
            zz, dd = z[lo:hi], d[lo:hi]
            dg = min(deg, len(zz) - 2)
            c = np.polyfit(zz - z[i], dd, dg); r = np.roots(c).real + z[i]
            r = r[(r > z[i] - 1e-9) & (r < z[i + 1] + 1e-9)]
            res = dd - np.polyval(c, zz - z[i]); slope = np.polyval(np.polyder(c), r[0] - z[i]) if len(r) else np.nan
            noise = np.sqrt(np.sum(res ** 2) / max(len(zz) - dg - 1, 1)) if len(zz) > dg + 1 else np.nan
            out.append((r[0] if len(r) else np.nan, abs(noise / slope) if np.isfinite(noise) else np.nan, slope, len(zz)))
    return out
tags = sys.argv[1:] or ["s20", "s16"]
S = {t: load(t) for t in tags}
for t, a in S.items():
    for zc, ez, sl, npts in crossings(a[:, 0], a[:, 3]):
        print(f"{t}: crossing z = {zc:.4f} ± {ez:.4f} (noise-based)  λ = {1/zc:.6f}  slope dm/dz = {sl:.2e}  ({npts} pts)")
if "s20" in S and "s16" in S:
    a, b = S["s20"], S["s16"]; n = min(len(a), len(b))
    off = a[:n, 3] - b[:n, 3]
    print("offset s20 − s16 at common z:", " ".join(f"{z:.3f}:{o:+.2e}" for z, o in zip(a[:n, 0], off)))
