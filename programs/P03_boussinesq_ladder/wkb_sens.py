"""Sensitivity of the WKB phase εΦ = ∫κ ds (stagnation point → front) to the front cut-off D/ε > Dcut, the quadrature
step and the profile resolution (hs); also the sign of Re κ and Im κ along the region (C = 2 Re a₀ > 0 needs Re κ > 0)."""
import numpy as np, sys
from bq_wkb2d import wkb_phase
for f in sys.argv[1:]:
    for Dcut in (2.0, 3.0, 5.0, 8.0):
        for ds in ((0.025, 0.0125) if Dcut == 3.0 else (0.025,)):
            lam, m, eps, s_cut, Phi, nf, grid, kap = wkb_phase(f, ds=ds, Dcut=Dcut, return_kappa=True)
            z = 1 / (lam - 1)
            print(f"{f[:40]}: z={z:.5f} Dcut={Dcut} ds={ds}: s_cut={s_cut:.4f} εΦ={Phi.real:.5f}{Phi.imag:+.5f}i "
                  f"ReΦ=2z·Re εΦ={2*z*Phi.real:.4f} | min Re κ={kap.real.min():.3e} max Im κ={kap.imag.max():.3e} (filled {nf})", flush=True)
