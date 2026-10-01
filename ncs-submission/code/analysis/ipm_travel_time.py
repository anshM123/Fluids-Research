"""Geometry of the IPM wall layer and its travel time T = ∫ ds / D̂ (D̂ = D/D₀), the purely geometric factor of the
phase: by the exact scaling κ = K(G, μ, R_y)/D̂, I = D₀ Re Φ₀ = ∫ Re K / D̂ ds with K bounded (K ≈ G for small G,
saturating at ~0.85–1.0 for G ≳ 2), so the phase can grow without bound only if T does.
For each state: T split as inner (x < 0.5), approach (to the dip core), dip core (contiguous D̂ below half depth, D̂ < (1 + D̂_min)/2), front
(core to the cut-off D̂ = 2); dip depth D̂_min (grid minimum and a parabolic three-point refinement), core width
w (s), number of grid cells across the core, dip and cut-off positions.  States as in ipm_phase3_table.py
(env WKB_NB / WKB_HS for files without grid tags).  Usage: ipm_travel_time.py STATES..."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, re, os
from scipy.integrate import cumulative_trapezoid
from ipm_wkb import wall_data

S_LO, DCUT = -10.0, 2.0
print("file                                z       T_tot   inner  approach  dipcore   front   D̂min_g  D̂min_p   w_core  cells   x_dip   x_cut")
for f in sys.argv[1:]:
    lam = None
    mm = re.search(r'lam([0-9.]+?)(?:_hs[0-9.]+)?\.npy', f); lam = float(mm.group(1)) if mm else None
    mz = re.search(r'_z([0-9.]+?)(?:_d[0-9.]+_hs.*)?\.npy', f)
    if mz: lam = 1 / float(mz.group(1))
    Nb = int(os.environ.get('WKB_NB', 0)) or None; hs = float(os.environ.get('WKB_HS', 0)) or None
    W = wall_data(f, lam, Nb, hs)
    s, Dh = W['s'], W['D'] / W['D0']
    hgrid = s[1] - s[0]
    sel = (s > -3) & (s < 2); idip = np.where(sel)[0][np.argmin(Dh[sel])]
    icut = idip + np.argmax(Dh[idip:] > DCUT)
    s_cut = np.interp(DCUT, Dh[icut - 1:icut + 1], s[icut - 1:icut + 1])
    a, b, c = Dh[idip - 1], Dh[idip], Dh[idip + 1]                       # parabolic refinement of the minimum
    den = a - 2 * b + c; off = 0.5 * (a - c) / den if den > 0 else 0.0
    Dmin_p = b - 0.25 * (a - c) * off
    lo = idip
    while lo > 0 and Dh[lo - 1] < 0.5 * (1 + Dh[idip]): lo -= 1
    hi = idip
    while hi < len(s) - 1 and Dh[hi + 1] < 0.5 * (1 + Dh[idip]): hi += 1
    i0 = np.searchsorted(s, S_LO)
    g = np.concatenate([s[i0:icut], [s_cut]]); y = np.concatenate([1 / Dh[i0:icut], [1 / DCUT]])
    cum = np.concatenate([[0.0], cumulative_trapezoid(y, g)])
    c05 = np.interp(np.log(0.5), g, cum)
    clo, chi = np.interp(s[lo], g, cum), np.interp(s[hi], g, cum)
    parts = (c05, clo - c05, chi - clo, cum[-1] - chi)
    print(f"{os.path.basename(f)[-34:]:34s} {1/W['lam']:7.4f}  {cum[-1]:7.4f}  " + "  ".join(f"{p:7.4f}" for p in parts)
          + f"  {Dh[idip]:.4f}  {Dmin_p:.4f}  {s[hi]-s[lo]:7.4f}  {hi-lo+1:4d}   {np.exp(s[idip] + off * hgrid):.4f}  {np.exp(s_cut):.4f}", flush=True)
