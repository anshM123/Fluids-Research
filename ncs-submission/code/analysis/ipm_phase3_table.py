"""Phase table with the front-safe tracker (wkb_phase3) on any list of states; env WKB_NB/WKB_HS/WKB_SS for files
without grid tags.  Prints z, D0, D̂_dip, s_cut, Re/Im Φ, I = ∫κ ds and the front failures."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, re, os
from ipm_wkb import wkb_phase3
for f in sys.argv[1:]:
    lam = None
    if not f.startswith("ipm_rung_"):
        mm = re.search(r'lam([0-9.]+?)(?:_hs[0-9.]+)?\.npy', f); lam = float(mm.group(1)) if mm else None
        mz = re.search(r'_z([0-9.]+?)\.npy', f)
        if mz: lam = 1 / float(mz.group(1))
    Nb = int(os.environ.get('WKB_NB', 0)) or None; hs = float(os.environ.get('WKB_HS', 0)) or None
    o = wkb_phase3(f, lam=lam, Nb=Nb, hs=hs, Dcut=2.0)
    print(f"{f[-34:]:34s} z={1/o['lam']:7.4f} D0={o['D0']:.5f} D̂dip={o['Dh_dip']:.4f} s_cut={o['s_cut']:6.3f} "
          f"Φ={o['Phi'].real:8.4f}{o['Phi'].imag:+8.4f}i I={o['I'].real:.5f} front_fail={o['nfail_front']}", flush=True)
