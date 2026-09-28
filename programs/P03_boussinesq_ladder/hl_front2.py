import numpy as np, sys, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from hl_solver import HL
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    q = np.load(f); N = len(q)
    S = HL(lam, N=N); r = S.march(S.full(q))
    D, eps, eta = r['D'], r['eps'], S.eta
    x = np.exp(eta)
    # front position: where D = 2*Dmin-ish ... use max of dD/dx
    sel = (eta > -3) & (eta < 1)
    k = np.argmax(np.gradient(D, x)[sel]); xc = x[sel][k]
    ax[0].plot(x[sel], D[sel] / eps, label=f"z={1/(lam-1):.2f} ε={eps:.4f}")
    right = (x > xc) & (x < xc + 1)
    ax[1].loglog(x[right] - xc, D[right], label=f"z={1/(lam-1):.2f}")
    left = (x < xc) & (x > xc - 0.4)
    ax[2].loglog(xc - x[left], D[left] / eps, label=f"z={1/(lam-1):.2f}")
    kd = np.argmin(D[sel]); print(f"z={1/(lam-1):6.2f} ε={eps:.5f}: x_front(max dD/dx)={xc:.4f}  dip x={x[sel][kd]:.4f} D̂min={D[sel][kd]/eps:.4f}")
xs = np.logspace(-3, 0, 20)
ax[1].loglog(xs, 1.0 * xs ** 0.5, 'k--', label='∝ (x−x_c)^{1/2}')
ax[2].loglog(xs, 1.0 * xs ** 0.5, 'k--', label='∝ (x_c−x)^{1/2}')
ax[0].set_xlabel('x'); ax[0].set_ylabel('D/ε'); ax[0].set_ylim(0, 5); ax[0].legend(fontsize=7)
ax[1].set_xlabel('x − x_c'); ax[1].set_ylabel('D (outer side)'); ax[1].legend(fontsize=7)
ax[2].set_xlabel('x_c − x'); ax[2].set_ylabel('D/ε (inner side)'); ax[2].legend(fontsize=7)
plt.tight_layout(); plt.savefig('hl_front.png', dpi=100)
