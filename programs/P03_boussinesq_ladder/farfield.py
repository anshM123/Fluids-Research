import numpy as np, sys, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bq_solver import BQ
from bq_newton import full
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
for f, hs in ((sys.argv[1], 0.025), (sys.argv[2], 0.0125)):
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    B = BQ(lam, hs=hs); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    cp = np.where(B.cb > 0, B.cb, 0.0)
    Th = r['Th'] * cp[None, :] ** r['m']; Om = r['Omega']
    al = (lam - 1) / (1 + lam)
    for sv, ls in ((1, ':'), (3, '--'), (6, '-.'), (12, '-')):
        i = np.argmin(np.abs(B.s - sv))
        ax[0].plot(B.beta, Th[i] / np.exp(al * B.s[i]), ls, label=f"λ={lam:.3f} s={sv}")
        ax[1].plot(B.beta, Om[i] * np.exp(B.s[i] / (1 + lam)), ls, label=f"λ={lam:.3f} s={sv}")
        ax[2].plot(B.beta, X[i] * np.exp(B.s[i] / (1 + lam)), ls, label=f"λ={lam:.3f} s={sv}")
ax[0].set_title('Θ / r^α'); ax[1].set_title('Ω r^{1/(1+λ)}'); ax[2].set_title('X r^{1/(1+λ)}')
for a in ax: a.set_xlabel('β'); a.legend(fontsize=6)
plt.tight_layout(); plt.savefig('farfield.png', dpi=100)
