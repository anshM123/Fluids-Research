"""Two checks of the IPM local wavenumber K = κD̂ used by the phase (NCS Fig. 4c,d).
1. K against the wall density gradient G along the wall of the λ₇ state (z = 7.346, h_s 0.00625) and a deep state
   (z = 8.426, h_s 0.0125): K ≈ G for small G and saturates for G ≳ 2, which makes Φ₀ a weighted travel time.
2. The exact scaling κ = K(ĉ, μ, G, R_y)/D̂: at three wall points of the λ₇ state (inner layer, dip, front side) the
   root is recomputed with D̂ replaced by 0.012 … 2 and (ĉ, μ, G, R_y) fixed; D̂κ must not change.
Output: ipm_K_checks.out and ipm_K_checks.npz."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np
from ipm_wkb import wall_data, wkb_phase3
from ipm_local_eig import local_root

out = {}
for tag, f, lam, hs in (("l7", "ipm_h7x_z7.3460.npy", 1 / 7.346, 0.00625), ("deep", "ipm_deep_e5_lam0.1186803.npy", 0.1186803, 0.0125)):
    W = wall_data(f, lam, 32, hs)
    o = wkb_phase3(f, lam=lam, Nb=32, hs=hs, Dcut=2.0, return_kappa=True)
    g = o['grid']; Dh = np.interp(g, W['s'], W['D']) / W['D0']
    out[f"{tag}_G"] = np.interp(g, W['s'], W['G']); out[f"{tag}_K"] = o['kap'] * Dh; out[f"{tag}_x"] = np.exp(g)
    out[f"{tag}_Dh"] = Dh
    print(f"{tag}: z = {1 / W['lam']:.3f}, {len(g)} wall points; max G = {out[f'{tag}_G'].max():.2f}, "
          f"Re K at max G = {out[f'{tag}_K'][np.argmax(out[f'{tag}_G'])].real:.3f}", flush=True)
    if tag == "l7":
        s = W['s']; idip = np.argmin(np.abs(g - o['s_dip']))
        pts = {"inner layer (x = 0.5)": np.argmin(abs(g - np.log(0.5))), "dip": idip,
               "front side (D̂ ≈ 1)": idip + np.argmin(abs(Dh[idip:] - 1.0))}
        print("D̂-invariance of K = κD̂ (fixed ĉ, μ, G, R_y):")
        rows = []
        for lab, i in pts.items():
            sv = g[i]
            p = [np.interp(sv, s, W['D']) / W['D0'], -np.interp(sv, s, W['Omb']), np.interp(sv, s, W['mu']),
                 np.interp(sv, s, W['G']), np.interp(sv, s, W['Ry'])]
            K0 = o['kap'][i] * p[0]
            line = f"  {lab:24s} (D̂ = {p[0]:.3f}, G = {p[3]:.2f}): K = {K0.real:.5f}{K0.imag:+.5f}i;"
            for Dt in (0.012, 0.03, 0.1, 0.3, 1.0, 2.0):
                q = list(p); q[0] = Dt
                r, ok = local_root(*q, K0 / Dt)
                rows.append((i, Dt, (r * Dt).real, (r * Dt).imag, ok))
                line += f"  D̂={Dt}: {abs(r * Dt - K0) / abs(K0):.1e}{'' if ok else '(!)'}"
            print(line, flush=True)
        out["inv"] = np.array(rows)
np.savez("ipm_K_checks.npz", **out)
