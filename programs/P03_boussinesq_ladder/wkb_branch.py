"""WKB phase Re Φ(z) along the whole branch (saved continuation states, Nb 32, hs 0.025, and the crossing states),
at a fixed front cut-off, for a fit Re Φ = C z + Φ₀ (+ Φ₁/z): C = 2 Re a₀ sets the asymptotic spacing π/C.
Output: wkb_branch_D{Dcut}.npy with rows (λ, z, Re εΦ, Im εΦ, m)."""
import numpy as np, sys, glob, warnings
warnings.filterwarnings('ignore')
from bq_wkb2d import wkb_phase
Dcut = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
files = sorted(glob.glob("Y2_[BCD]_lam*.npy") + glob.glob("Y2_C_lam*.npy") + glob.glob("Y2_A_lam1.46*.npy")
               + glob.glob("Ycross_Nb32_hs0.0125_ss-20.0_sm100.0_lam*.npy")
               + ["Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.184253.npy", "Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.252349.npy",
                  "Ycross_Nb32_hs0.025_ss-20.0_sm100.0_lam1.399096.npy"])
files = list(dict.fromkeys(files))
rows = []
for f in files:
    try:
        lam, m, eps, s_cut, Phi, nf = wkb_phase(f, Dcut=Dcut)
    except Exception as e:
        print("skip", f, e, flush=True); continue
    z = 1 / (lam - 1)
    rows.append((lam, z, Phi.real, Phi.imag, m))
    print(f"{f}: z={z:.4f} m={m:.8f} εΦ={Phi.real:.5f}{Phi.imag:+.5f}i ReΦ={m*z*Phi.real:.4f}", flush=True)
    np.save(f"wkb_branch_D{Dcut:g}.npy", np.array(sorted(rows)))
