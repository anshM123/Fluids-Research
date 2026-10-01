"""λ₇ decision analysis (rule fixed before the crossing is known).
1. For each matched scan (p7: h_s 0.0125, dz 0.03; h7: h_s 0.00625, dz 0.06) find the sign change of m − 2 and fit
   a local polynomial (degree ≤ 3) through the points within ±0.1 in z of it; z₇ is its root.
2. Front–grid error: on the dense p7 scan, the rms residual of a smooth (cubic) fit over the whole scan, divided by
   the local slope, gives σ_z(grid) for h_s = 0.0125.  For h7 the same σ is scaled by the ratio of the two grids'
   scatter about their own smooth fits (if h7 has ≥ 5 points), else taken as σ(p7)/4 (fourth-order scaling of the
   locking error, stated as an assumption).
3. If |z₇(p7) − z₇(h7)| > the combined σ, Richardson-extrapolate with order 2 (z∞ = z_h7 + (z_h7 − z_p7)/3) and
   report it with the difference as its uncertainty; otherwise report z₇(h7).
4. Compare with the registered predictions (H1 c922751; H2, H3 c922751; H4 273ae28) in units of σ_z.
   Decision rule: a hypothesis is disfavoured when |z₇ − z_H| > 2σ_z, and the comparison is called decisive only
   if σ_z ≤ 0.003."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np

PRED = {"H1 (shifted law, λ2–λ6)": 7.3415, "H1' (shifted law, λ0–λ6)": 7.3404, "H2 (geometric)": 7.3343,
        "H3 (3-point)": 7.3342, "H4 (phase, coarse branch)": 7.352}


def load(tag):
    a = np.load(f"ipm_scan_{tag}.npy")
    o = np.argsort(a[:, 0]); return a[o, 0], a[o, 3]


def crossing(z, d, half=0.1):
    i = [k for k in range(len(z) - 1) if d[k] * d[k + 1] < 0]
    if not i:
        return None
    k = i[0]; zc = z[k] - d[k] * (z[k + 1] - z[k]) / (d[k + 1] - d[k])
    sel = np.abs(z - zc) <= half + 1e-9
    deg = min(3, sel.sum() - 1)
    c = np.polyfit(z[sel] - zc, d[sel], deg); r = np.roots(c).real + zc
    r = r[(r >= z[k] - 1e-9) & (r <= z[k + 1] + 1e-9)]
    root = r[0] if len(r) else zc
    slope = np.polyval(np.polyder(c), root - zc)
    return root, slope


def scatter(z, d):
    if len(z) < 6:
        return np.nan
    c = np.polyfit(z - z.mean(), d, 3); return np.sqrt(np.mean((d - np.polyval(c, z - z.mean())) ** 2))


if __name__ == "__main__":
    res = {}
    for tag in ("p7", "h7"):
        try:
            z, d = load(tag)
        except FileNotFoundError:
            continue
        cr = crossing(z, d); sc = scatter(z, d)
        res[tag] = (z, d, cr, sc)
        print(f"{tag}: {len(z)} points, z ∈ [{z.min():.3f}, {z.max():.3f}]; crossing: "
              f"{'none in range' if cr is None else f'z = {cr[0]:.4f}, slope {cr[1]:.3e}'}; scatter about cubic {sc:.2e}")
    if "p7" in res and res["p7"][2]:
        zp, sp = res["p7"][2]; sig_p = res["p7"][3] / abs(sp)
        print(f"σ_z(h_s 0.0125) from front–grid scatter: {sig_p:.4f}")
        if "h7" in res and res["h7"][2]:
            zh, sh = res["h7"][2]
            sig_h = (res["h7"][3] / abs(sh)) if np.isfinite(res["h7"][3]) else sig_p / 4
            comb = np.hypot(sig_p, sig_h)
            if abs(zp - zh) > comb:
                z7 = zh + (zh - zp) / 3; sz = max(abs(zh - zp), sig_h)
                print(f"grids differ by {zh - zp:+.4f} (> {comb:.4f}): Richardson z₇ = {z7:.4f} ± {sz:.4f}")
            else:
                z7, sz = zh, sig_h
                print(f"grids agree ({zh - zp:+.4f} ≤ {comb:.4f}): z₇ = {z7:.4f} ± {sz:.4f}")
        else:
            z7, sz = zp, sig_p
            print(f"only h_s 0.0125: z₇ = {z7:.4f} ± {sz:.4f}")
        print(f"λ₇ = {1/z7:.6f}; decisive: {sz <= 0.003}")
        for k, v in PRED.items():
            dev = (v - z7) / sz
            print(f"  {k:28s} z = {v:.4f}: deviation {v - z7:+.4f} = {dev:+.1f}σ {'(disfavoured)' if abs(dev) > 2 else ''}")
