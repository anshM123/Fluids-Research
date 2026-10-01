"""Where does the drift of δ_n = Re Φ₀(λ_n) − nπ at n = 6 come from?  Diagnostic, not a recalibration: the stage-4
predictions keep the registered cut-off D/D₀ = 2.  For the six fine-grid rungs, Re Φ₀ with the repaired tracker
(wkb_phase3) at front cut-offs D/D₀ = 1.5, 2, 3 and 5, and the offsets δ_n for each; a drift that changes with the
cut-off comes from the front side of the layer, one that does not comes from the inner layer or the dip."""
import numpy as np, os
from ipm_wkb import wkb_phase3

R = [(1, "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.47212973.npy", None),
     (2, "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.31496181.npy", None),
     (3, "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.24156634.npy", None),
     (4, "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.19872245.npy", None),
     (5, "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.17061809.npy", None),
     (6, "ipm_noise_lam0.1509200_hs0.0125.npy", 0.15092)]
CUTS = (1.5, 2.0, 3.0, 5.0)
print("n   " + "  ".join(f"δ(cut {c:g})  fail" for c in CUTS), flush=True)
D = {c: [] for c in CUTS}
for n, f, lam in R:
    kw = dict(Nb=32, hs=0.0125) if lam else {}
    line = f"{n}  "
    for c in CUTS:
        o = wkb_phase3(f, lam=lam, Dcut=c, **kw)
        d = o['Phi'].real - n * np.pi; D[c].append(d)
        line += f"  {d:9.4f}  {o['nfail_front']:3d}"
    print(line, flush=True)
for c in CUTS:
    d = np.array(D[c])
    print(f"cut {c:g}: mean δ(1–5) = {d[:5].mean():.4f} ± {d[:5].std(ddof=1):.4f};  δ₆ − mean = {d[5] - d[:5].mean():+.4f};  "
          f"spread over 1–6 = {np.ptp(d):.4f}")
