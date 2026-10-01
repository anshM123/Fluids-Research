"""Phase stability under the same numerical changes that move the defect (NCS Fig. 4b): the repaired phase
(wkb_phase3) on every variant state at z = 7.346 (h_s = 0.00625): four grid shifts, Nb 48 and 64, origin truncation
s_start −22 and −24, far field s_max 130 and 160, and the combined variant.  For each: m − 2 of the state (march),
Re Φ₀ and I; and the implied shifts of z₇, Δm/|dm/dz| (defect, slope −2.5e-8 per unit z at the crossing) and
ΔReΦ/(dReΦ/dz) (phase, 4.54 per unit z), relative to the unshifted production state."""
import numpy as np, os
from ipm_wkb import wkb_phase3

Z = 7.346; HS = 0.00625; DMDZ = 2.5e-8; DPDZ = 4.54
V = [("production (Nb 32, s_start −20, s_max 100)", "ipm_shift_z7.3460_d0.000_hs0.00625.npy", 32, -20, 100, 0.0),
     ("grid shift 1/4", "ipm_shift_z7.3460_d0.250_hs0.00625.npy", 32, -20, 100, 0.25),
     ("grid shift 1/2", "ipm_shift_z7.3460_d0.500_hs0.00625.npy", 32, -20, 100, 0.5),
     ("grid shift 3/4", "ipm_shift_z7.3460_d0.750_hs0.00625.npy", 32, -20, 100, 0.75),
     ("Nb 48", "ipm_shift_z7.3460_d0.000_hs0.00625_nb48_ss-20_sm100.npy", 48, -20, 100, 0.0),
     ("Nb 64", "ipm_shift_z7.3460_d0.000_hs0.00625_nb64_ss-20_sm100.npy", 64, -20, 100, 0.0),
     ("s_start −22", "ipm_shift_z7.3460_d0.000_hs0.00625_nb32_ss-22_sm100.npy", 32, -22, 100, 0.0),
     ("s_start −24", "ipm_shift_z7.3460_d0.000_hs0.00625_nb32_ss-24_sm100.npy", 32, -24, 100, 0.0),
     ("s_max 130", "ipm_shift_z7.3460_d0.000_hs0.00625_nb32_ss-20_sm130.npy", 32, -20, 130, 0.0),
     ("s_max 160", "ipm_shift_z7.3460_d0.000_hs0.00625_nb32_ss-20_sm160.npy", 32, -20, 160, 0.0),
     ("Nb 48 + s_start −22 + s_max 130", "ipm_shift_z7.3460_d0.000_hs0.00625_nb48_ss-22_sm130.npy", 48, -22, 130, 0.0)]
print(f"{'variant':36s}  {'m − 2':>11s}  {'Re Φ₀':>9s}  {'I':>8s}  {'Δz (defect)':>12s}  {'Δz (phase)':>11s}  front_fail")
base = None
for lab, f, nb, ss, sm, dl in V:
    if not os.path.exists(f):
        print(f"{lab:36s}  missing {f}"); continue
    os.environ['WKB_SS'] = str(ss); os.environ['WKB_SMAX'] = str(sm); os.environ['WKB_SMIN'] = str(-120.0 + dl * HS)
    o = wkb_phase3(f, lam=1 / Z, Nb=nb, hs=HS, Dcut=2.0)
    m2, P, I = o['m'] - 2, o['Phi'].real, o['I'].real
    if base is None: base = (m2, P)
    print(f"{lab:36s}  {m2:+.3e}  {P:9.4f}  {I:.5f}  {(m2 - base[0]) / DMDZ:+12.4f}  {(P - base[1]) / DPDZ:+11.5f}  {o['nfail_front']}", flush=True)
