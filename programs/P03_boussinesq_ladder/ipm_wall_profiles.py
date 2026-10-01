"""Wall speed D̂ = D/D₀ and the local wavenumber density Re K/D̂ along the wall for a few IPM states (NCS Fig. 1):
rungs 1, 3, 5, the λ₇ state (h_s 0.00625) and two deep states.  Saved to ipm_wall_profiles.npz (x = e^s, D̂, and
for the λ₇ state the repaired phase density on its tracking grid)."""
import numpy as np, os
from ipm_wkb import wall_data, wkb_phase3

S = [("n = 1", "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.47212973.npy", None, None, None),
     ("n = 3", "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.24156634.npy", None, None, None),
     ("n = 5", "ipm_rung_Nb32_hs0.0125_ss-20.0_sm100.0_lam0.17061809.npy", None, None, None),
     ("z = 7.346", "ipm_h7x_z7.3460.npy", 1 / 7.346, 32, 0.00625),
     ("z = 8.426", "ipm_deep_e5_lam0.1186803.npy", 0.1186803, 32, 0.0125),
     ("z = 9.506", "ipm_deep_e5_lam0.1051967.npy", 0.1051967, 32, 0.0125)]
out = {}
for k, (lab, f, lam, nb, hs) in enumerate(S):
    W = wall_data(f, lam, nb, hs)
    sel = (W['s'] > np.log(0.05)) & (W['s'] < np.log(3.0))
    out[f"x{k}"] = np.exp(W['s'][sel]); out[f"D{k}"] = W['D'][sel] / W['D0']; out[f"lab{k}"] = lab
    out[f"z{k}"] = 1 / W['lam']
    print(lab, f"z = {1 / W['lam']:.4f}", flush=True)
o = wkb_phase3("ipm_h7x_z7.3460.npy", lam=1 / 7.346, Nb=32, hs=0.00625, Dcut=2.0, return_kappa=True)
out["phase_x"] = np.exp(o['grid']); out["phase_k"] = o['kap']; out["phase_D0"] = o['D0']
np.savez("ipm_wall_profiles.npz", **out)
print("saved ipm_wall_profiles.npz")
