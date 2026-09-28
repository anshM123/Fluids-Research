import numpy as np, sys, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bq_solver import BQ
from bq_newton import full
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
    Ur = X @ B.Db.T; D = (1 + lam) + Ur[:, 0]
    A = -np.mean(Ur[B.i0:B.i0 + 40, 0]); eps = 1 + lam - A
    x = np.exp(B.s); sel = (B.s > -3) & (B.s < 1.5)
    k = np.argmax(np.gradient(D, x)[sel]); xc = x[sel][k]
    kd = np.argmin(D[sel])
    right = (x > xc) & (x < xc + 2)
    ax[0].loglog(x[right] - xc, D[right], label=f"λ={lam:.3f} ε={eps:.3f}")
    left = (x < xc) & (x > 0.05)
    ax[1].loglog(xc - x[left], D[left] / eps, label=f"λ={lam:.3f}")
    print(f"λ={lam:.4f} z={1/(lam-1):.2f} ε={eps:.5f}: x_front={xc:.4f} dip at x={x[sel][kd]:.4f} D̂min={D[sel][kd]/eps:.4f}")
xs = np.logspace(-3, 0.3, 20)
ax[0].loglog(xs, 1.2 * xs ** 0.5, 'k--', label='∝ (x−x_c)^{1/2}')
ax[0].set_xlabel('x − x_c'); ax[0].set_ylabel('V_r/r on the boundary (outer side)'); ax[0].legend(fontsize=7)
ax[1].set_xlabel('x_c − x'); ax[1].set_ylabel('(V_r/r)/ε (inner side)'); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig('bq_front.png', dpi=100)
