import numpy as np, sys, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bq_solver import BQ
from bq_newton import full
fig, ax = plt.subplots(2, 3, figsize=(17, 8))
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    A, m = r['A'], r['m']; eps = 1 + lam - A
    s = B.s; sel = (s > -8) & (s < 6)
    cpow = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cpow[None, :] ** m; Om = r['Omega']
    vr = (1 + lam) + r['Ur']
    lab = f"λ={lam:.3f} ε={eps:.3f}"
    ax[0, 0].plot(s[sel], (vr[sel, 0]) / eps, label=lab)
    ax[0, 1].plot(s[sel], Om[sel, 0] / np.exp(s[sel]), label=lab)        # Ω/r on boundary (≈ −C near origin)
    ax[0, 2].plot(s[sel], Th[sel, 0] / (-np.exp(2 * s[sel])), label=lab)  # Θ/(−r²) on boundary
    for bi, ls in ((2, '-'), (5, '--'), (8, ':')):
        ax[1, 0].plot(s[sel], vr[sel, bi] / eps, ls, label=f"{lab} β={B.beta[bi]:.3f}" if bi == 2 else None)
    # vorticity along the boundary relative to the local-balance prediction Ω ≈ ∂1Θ = e^{-s} ∂sΘ
    dTh = np.gradient(Th[:, 0], s)
    ax[1, 1].plot(s[sel], Om[sel, 0] / (np.exp(-s[sel]) * dTh[sel]), label=lab)
    ax[1, 2].plot(s[sel], r['w'][sel, 2], label=lab)
ax[0, 0].set_ylabel('(V_r/r)/ε on boundary'); ax[0, 1].set_ylabel('Ω/r on boundary'); ax[0, 2].set_ylabel('Θ/(−r²) on boundary')
ax[1, 0].set_ylabel('(V_r/r)/ε at β=beta_2,5,8'); ax[1, 1].set_ylabel('Ω/∂₁Θ on boundary'); ax[1, 2].set_ylabel('w at β_2')
for a in ax.flat: a.set_xlabel('s = ln r'); a.legend(fontsize=6)
ax[0, 0].set_ylim(0, 4); ax[1, 0].set_ylim(0, 6)
plt.tight_layout(); plt.savefig('viz_bnd.png', dpi=100); print('ok')
