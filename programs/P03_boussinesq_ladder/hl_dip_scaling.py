"""Hou–Luo: scaling of the dip/front structure of the stalled layer as ε → 0 (states of scans E and F, z = 4–40):
dip depth D̂_min, dip half-width, dip–front distance (front: first D̂ = 3 after the dip), and the local root
κ ~ D̂^{-3/2} at the dip. Power-law fits over the deepest states."""
import numpy as np, glob, re
from hl_solver import HL
rows = []
for pat, N, L1, L2, e0 in (("hl_q_E_lam*.npy", 32768, 40, 160, -30), ("hl_q_F_lam*.npy", 65536, 25, 75, -20)):
    for f in sorted(glob.glob(pat)):
        lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
        S = HL(lam, N=N, L1=L1, L2=L2, eta_start=e0); r = S.march(S.full(np.load(f)))
        eps, eta = r['eps'], S.eta; Dh = r['D'] / eps
        sel = (eta > -3) & (eta < 2); kd = np.argmax(sel) + np.argmin(Dh[sel]); dmin = Dh[kd]
        half = 0.5 * (1 + dmin)
        left = kd - np.argmax(Dh[kd::-1] > half); right = kd + np.argmax(Dh[kd:] > half)
        kc = kd + np.argmax(Dh[kd:] > 3)
        rows.append((eps, 1 / (lam - 1), eta[kd], dmin, eta[right] - eta[left], eta[kc] - eta[kd], eta[kc]))
rows = np.array(sorted(rows, key=lambda r: -r[0]))
print(" ε        z      η_dip    D̂_min   halfwidth  dip→front  η_front")
for r in rows:
    print(f"{r[0]:.5f} {r[1]:7.3f} {r[2]:8.4f} {r[3]:.4f} {r[4]:9.4f} {r[5]:9.4f} {r[6]:8.4f}")
deep = rows[rows[:, 0] < 0.03]
for j, name in ((3, "D̂_min"), (4, "half-width"), (5, "dip→front")):
    p = np.polyfit(np.log(deep[:, 0]), np.log(deep[:, j]), 1)
    print(f"{name} ∝ ε^{p[0]:.3f}  (fit over ε < 0.03, {len(deep)} states)")
np.save("hl_dip_scaling.npy", rows)
