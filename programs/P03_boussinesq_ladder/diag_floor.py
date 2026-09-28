import numpy as np
from bq_logpolar import BQLogPolar
from bq_newton import residual
lam = 1.92
B = BQLogPolar(lam, s_min=-120, s_max=100, hs=0.025, Nb=32)
Y = np.load("bq_X_lam1.9200.npy")
R, info = residual(B, Y)
print("A m vrmin", info['A'], info['m'], info['vrmin'], " |R|max", np.abs(R).max())
ss = B.s[B.i0:]
for lo, hi in ((-20, -19.5), (-19.5, -15), (-15, -5), (-5, 5), (5, 20), (20, 60), (60, 100)):
    sel = (ss >= lo) & (ss < hi)
    Rs = np.abs(R[sel])
    k = np.unravel_index(np.argmax(Rs), Rs.shape)
    print(f"s∈[{lo},{hi}): max|R|={Rs.max():.2e} at β-index {k[1]} (β={B.beta[k[1]]:.3f}), s={ss[sel][k[0]]:.3f}; median {np.median(Rs):.1e}")
print("β profile of |R| max over s∈[-5,5]:", " ".join(f"{x:.0e}" for x in np.abs(R[(ss>-5)&(ss<5)]).max(axis=0)))
