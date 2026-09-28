import numpy as np, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bq_solver import BQ
from bq_newton import full
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
for f in sys.argv[1:]:
    lam = float(f.split('lam')[1][:6])
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    vr = (1 + lam) + r['Ur']
    eps = 1 + lam - r['A']
    s = B.s; sel = (s > -5) & (s < 5)
    ax[0].plot(np.exp(s[sel]), vr[sel, 0], label=f"λ={lam:.3f}, ε={eps:.3f}")
    ax[1].plot(np.exp(s[sel]), vr[sel, 0] / eps, label=f"λ={lam:.3f}")
    cpow = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cpow[None, :] ** r['m']
    ax[2].loglog(np.exp(s[sel]), -Th[sel, 0], label=f"λ={lam:.3f}")
ax[0].set_xscale('log'); ax[1].set_xscale('log')
ax[0].set_xlabel('r (boundary)'); ax[0].set_ylabel('V_r/r on boundary'); ax[0].legend(fontsize=7)
ax[1].set_xlabel('r'); ax[1].set_ylabel('(V_r/r)/ε'); ax[2].set_xlabel('r'); ax[2].set_ylabel('−Θ(r, β=0)')
plt.tight_layout(); plt.savefig('viz_dip.png', dpi=100); print('ok')
