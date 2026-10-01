"""Diagnostic of the phase-tracker failure on deep IPM states (NCS Fig. 6): K = κD̂ along the wall from the original
tracker (`wkb_phase`, which holds κ across a failed point) and from the repaired one (`wkb_phase3`, K-continuation
beyond the dip), together with I = ∫Re κ ds from both.  For each state the front side (dip → cut-off) is printed
point by point; the arrays are saved to ipm_tracker_diag.npz for the figure.
Usage: WKB_NB=32 WKB_HS=0.0125 ipm_tracker_diag.py STATES..."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, re, os
from ipm_wkb import wkb_phase, wkb_phase3, wall_data

save = {}
for f in sys.argv[1:]:
    mm = re.search(r'lam([0-9.]+?)(?:_hs[0-9.]+)?\.npy', f); lam = float(mm.group(1)) if mm else None
    mz = re.search(r'_z([0-9.]+?)\.npy', f)
    if mz: lam = 1 / float(mz.group(1))
    Nb = int(os.environ.get('WKB_NB', 0)) or None; hs = float(os.environ.get('WKB_HS', 0)) or None
    W = wall_data(f, lam, Nb, hs)
    o1 = wkb_phase(f, lam=lam, Nb=Nb, hs=hs, Dcut=2.0, return_kappa=True)
    o3 = wkb_phase3(f, lam=lam, Nb=Nb, hs=hs, Dcut=2.0, return_kappa=True)
    Dh = lambda g: np.interp(g, W['s'], W['D']) / W['D0']
    z = 1 / W['lam']
    print(f"z = {z:.4f}: I(original) = {np.trapezoid(o1['kap'].real, o1['grid']):.5f}  (nfail {o1['nfail']}),  "
          f"I(repaired) = {o3['I'].real:.5f}  (front_fail {o3['nfail_front']});  dip s = {o3['s_dip']:.3f}, cut s = {o3['s_cut']:.3f}")
    print("     s       D̂      K original          K repaired")
    for g in np.linspace(o3['s_dip'] - 0.1, o3['s_cut'], 9):
        k1 = np.interp(g, o1['grid'], o1['kap'].real) + 1j * np.interp(g, o1['grid'], o1['kap'].imag)
        k3 = np.interp(g, o3['grid'], o3['kap'].real) + 1j * np.interp(g, o3['grid'], o3['kap'].imag)
        print(f"  {g:+.3f}  {Dh(g):6.3f}  {k1.real * Dh(g):7.3f}{k1.imag * Dh(g):+7.3f}i   {k3.real * Dh(g):7.3f}{k3.imag * Dh(g):+7.3f}i")
    key = f"z{z:.3f}"
    save[key + "_g1"], save[key + "_K1"] = o1['grid'], o1['kap'] * Dh(o1['grid'])
    save[key + "_g3"], save[key + "_K3"] = o3['grid'], o3['kap'] * Dh(o3['grid'])
    save[key + "_I"] = np.array([np.trapezoid(o1['kap'].real, o1['grid']), o3['I'].real])
np.savez("ipm_tracker_diag.npz", **save)
