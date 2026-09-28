"""Plot Θ, Ω of converged profiles (log-polar solution files) in physical coordinates."""
import numpy as np, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bq_logpolar import BQLogPolar
from bq_newton import full
files = sys.argv[1:]
fig, ax = plt.subplots(len(files), 3, figsize=(15, 4.2*len(files)))
ax = np.atleast_2d(ax)
for row, f in enumerate(files):
    lam = float(f.split('lam')[1][:6])
    B = BQLogPolar(lam)
    Y = np.load(f); X = full(B, Y); r = B.march(X / B.ea2[:, None], return_all=True)
    cpow = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cpow[None, :] ** r['m']
    Om = r['Omega']
    s = B.s; sel = (s > -6) & (s < 4)
    R = np.exp(s[sel])[:, None]; b = B.beta[None, :]
    y1 = R*np.cos(b); y2 = R*np.sin(b)
    for k, (F, name) in enumerate([(Th[sel], 'Θ'), (Om[sel], 'Ω')]):
        a = ax[row, k]
        lev = np.linspace(F.min(), F.max(), 41)
        cs = a.contourf(y1, y2, F, levels=lev, cmap='RdBu_r')
        plt.colorbar(cs, ax=a)
        a.set_xlim(0, 3); a.set_ylim(0, 3); a.set_aspect('equal')
        a.set_title(f"{name}, λ={lam:.4f}, m={r['m']:.4f}")
    a = ax[row, 2]
    for bi in (0, 4, 8, 16, 24, 30, 32):
        a.plot(s[s > -8], Th[s > -8, bi] / np.exp(r['m']*0) , label=f"β={B.beta[bi]:.2f}")
    a.set_xlim(-8, 20); a.set_title('Θ vs s=ln r (rays)'); a.legend(fontsize=6)
plt.tight_layout(); plt.savefig("viz_profiles.png", dpi=110)
print("saved")
