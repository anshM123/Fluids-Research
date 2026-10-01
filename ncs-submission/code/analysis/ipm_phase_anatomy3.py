"""Where the IPM phase accumulates, with the front-safe tracker (wkb_phase3).  I = ∫ Re κ ds (κ in units of 1/D₀,
so I = D₀ Re Φ₀ = ∫ Re K / D̂ ds) is split by an exact cumulative partition into
  inner     x < 0.5                                   (s < ln 0.5)
  approach  x ≥ 0.5 up to the dip core
  dip core  the contiguous region around the dip below half depth, D̂ < (1 + D̂_min)/2
  front     from the dip core to the cut-off D̂ = 2
and the core width w (in s) and the mean of Re K over the core are reported.  Usage: ipm_phase_anatomy3.py STATES..."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, re, os
from scipy.integrate import cumulative_trapezoid
from ipm_wkb import wkb_phase3, wall_data

print("file                                z        D0      I_tot    inner  approach  dipcore    front   w_core  ReK_core  D̂min   x_dip   x_cut")
for f in sys.argv[1:]:
    lam = None
    mm = re.search(r'lam([0-9.]+?)(?:_hs[0-9.]+)?\.npy', f); lam = float(mm.group(1)) if mm else None
    mz = re.search(r'_z([0-9.]+?)(?:_d[0-9.]+_hs.*)?\.npy', f)
    if mz: lam = 1 / float(mz.group(1))
    Nb = int(os.environ.get('WKB_NB', 0)) or None; hs = float(os.environ.get('WKB_HS', 0)) or None
    o = wkb_phase3(f, lam=lam, Nb=Nb, hs=hs, Dcut=2.0, return_kappa=True)
    W = wall_data(f, o['lam'], Nb, hs)
    g, k = o['grid'], o['kap']
    Dh = np.interp(g, W['s'], W['D']) / o['D0']
    cum = np.concatenate([[0.0], cumulative_trapezoid(k.real, g)])
    idip = int(np.argmin(Dh))
    lo = idip
    while lo > 0 and Dh[lo - 1] < 0.5 * (1 + Dh[idip]): lo -= 1
    hi = idip
    while hi < len(g) - 1 and Dh[hi + 1] < 0.5 * (1 + Dh[idip]): hi += 1
    s05 = np.log(0.5)
    c05 = np.interp(s05, g, cum)
    parts = (c05, cum[lo] - c05, cum[hi] - cum[lo], cum[-1] - cum[hi])
    K = k * Dh
    print(f"{os.path.basename(f)[-34:]:34s} {1/o['lam']:7.4f} {o['D0']:.5f}  {cum[-1]:.5f}  " + "  ".join(f"{p:7.4f}" for p in parts)
          + f"  {g[hi]-g[lo]:7.4f}  {np.mean(K[lo:hi+1].real):7.3f}  {Dh[idip]:.4f}  {np.exp(g[idip]):.4f}  {np.exp(g[-1]):.4f}", flush=True)
