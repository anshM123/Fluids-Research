"""Where along the wall the IPM phase accumulates: I = ∫Re κ ds split into bins of x = e^s (fixed physical
positions), plus the dip region; for a list of states (resolved fine-grid profiles)."""
import numpy as np, sys, re, os
from ipm_wkb import wkb_phase
edges = [0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 5.0]
print("z       D0      I_total  " + "  ".join(f"x<{e2:<4}" if e1 == 0 else f"{e1}-{e2}" for e1, e2 in zip(edges[:-1], edges[1:])))
for f in sys.argv[1:]:
    lam = None
    if not f.startswith("ipm_rung_"):
        mm = re.search(r'lam([0-9.]+?)(?:_hs[0-9.]+)?\.npy', f); lam = float(mm.group(1)) if mm else None
        mz = re.search(r'_z([0-9.]+?)\.npy', f)
        if mz: lam = 1 / float(mz.group(1))
    Nb = int(os.environ.get('WKB_NB', 0)) or None; hs = float(os.environ.get('WKB_HS', 0)) or None
    o = wkb_phase(f, lam=lam, Nb=Nb, hs=hs, Dcut=2.0, return_kappa=True)
    x = np.exp(o['grid']); k = o['kap'].real
    parts = []
    for e1, e2 in zip(edges[:-1], edges[1:]):
        sel = (x >= e1) & (x < e2)
        parts.append(np.trapezoid(k[sel], o['grid'][sel]) if sel.sum() > 1 else 0.0)
    print(f"{1/o['lam']:6.3f}  {o['D0']:.5f}  {np.trapezoid(k, o['grid']):.4f}   " + "  ".join(f"{p:7.4f}" for p in parts), flush=True)
