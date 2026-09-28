"""m(λ) − 2 along the Boussinesq branch vs z = 1/(λ−1), from all scan2_*.npy files (+ old scans)."""
import numpy as np, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows = []
for f in sorted(glob.glob("scan2_*.npy")):
    d = np.load(f)
    if d.ndim == 2 and len(d):
        rows.append(d)
D = np.vstack(rows); D = D[np.argsort(D[:, 0])]
lam, A, m, vrmin = D.T
z = 1 / (lam - 1)
eps = 1 + lam - A
fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
ax[0].plot(lam, m, 'k.-', ms=3); ax[0].axhline(2, color='r', lw=0.7); ax[0].set_xlabel('λ'); ax[0].set_ylabel('m')
ax[1].plot(z, m - 2, 'k.-', ms=3); ax[1].axhline(0, color='r', lw=0.7); ax[1].set_yscale('symlog', linthresh=1e-6)
ax[1].set_xlabel('z = 1/(λ−1)'); ax[1].set_ylabel('m − 2 (symlog)')
for zz in 1.0863 + 1.4187 * np.arange(8):
    ax[1].axvline(zz, color='b', lw=0.5, ls=':')
ax[2].plot(z, vrmin / eps, 'k.-', ms=3); ax[2].set_xlabel('z'); ax[2].set_ylabel('min V_r/r on boundary / ε')
plt.tight_layout(); plt.savefig('branch_m.png', dpi=110)
sgn = np.sign(m - 2); idx = np.where(sgn[1:] * sgn[:-1] < 0)[0]
for i in idx:
    l0, l1, m0, m1 = lam[i], lam[i+1], m[i], m[i+1]
    lc = l0 + (2 - m0) * (l1 - l0) / (m1 - m0)
    print(f"crossing λ ≈ {lc:.7f}  z = {1/(lc-1):.5f}")
