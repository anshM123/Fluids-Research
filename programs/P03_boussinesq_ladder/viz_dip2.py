import numpy as np, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bq_solver import BQ
from bq_newton import full
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for f in sys.argv[1:]:
    lam = float(f.split('lam')[1][:6])
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
    Ur = X @ B.Db.T
    vr = (1 + lam) + Ur
    A = -np.mean(Ur[B.i0:B.i0 + 40, 0]); eps = 1 + lam - A
    s = B.s; sel = (s > -6) & (s < 1)
    d = vr[sel, 0] - eps
    ax[0].loglog(np.exp(s[sel]), np.abs(d), label=f"λ={lam:.3f}")
    ax[1].plot(np.exp(s[sel]), d / np.exp(2*s[sel]), label=f"λ={lam:.3f}")
    # fit d = b3 r^2 + b5 r^4 + b7 r^6 on small r
    ss = (s > -5) & (s < -1.5)
    rr = np.exp(2*s[ss]); co = np.polyfit(rr, (vr[ss, 0] - eps) / rr, 3)
    print(f"λ={lam:.4f} ε={eps:.5f}: b3={co[-1]:+.4f} b5={co[-2]:+.4f} b7={co[-3]:+.4f} b9={co[-4]:+.4f}   b3^2/(4 ε b5)={co[-1]**2/(4*eps*co[-2]):.3f}")
ax[0].set_xlabel('r'); ax[0].set_ylabel('|V_r/r − ε| on boundary'); ax[0].legend(fontsize=7)
ax[1].set_xlabel('r'); ax[1].set_ylabel('(V_r/r − ε)/r²'); ax[1].set_xscale('log')
plt.tight_layout(); plt.savefig('viz_dip2.png', dpi=100)
