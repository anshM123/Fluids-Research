import numpy as np
from bq_solver import BQ
from bq_newton import residual, full
B = BQ(1.92)
Y = np.load('Y2_lam1.9200_Nb32_hs0.025.npy')
R, info = residual(B, Y)
print("|R| =", np.abs(R).max(), "A =", info['A'])
s = B.s[B.i0:]
i, j = np.unravel_index(np.argmax(np.abs(R)), R.shape)
print("argmax at s =", s[i], "beta index", j)
for lo, hi in [(-20, -15), (-15, -5), (-5, 0), (0, 5), (5, 20), (20, 50), (50, 100)]:
    sel = (s >= lo) & (s < hi)
    print(f"  s in [{lo},{hi}): max|R| = {np.abs(R[sel]).max():.3e}")
# projection on the uniform strain mode sin 2b in s<-5
sel = s < -5
mode = -0.5*np.sin(2*B.beta)
coef = (R[sel] @ mode) / (mode @ mode)
print("strain-mode coefficient (s<-5): mean", coef.mean(), "std", coef.std())
# smoothness of R along random direction
rng = np.random.default_rng(0)
v = rng.standard_normal(Y.shape); v[:, 0] = 0; v[:, -1] = 0
v *= np.exp(-0.1*np.abs(s))[:, None]
v /= np.abs(v).max()
for eps in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
    Rp, _ = residual(B, Y + eps*v); Rm, _ = residual(B, Y - eps*v)
    Jc = (Rp - Rm)/(2*eps); Jf = (Rp - R)/eps
    print(f"eps={eps:.0e}: |Jc v|={np.abs(Jc).max():.6e} |Jf v - Jc v|={np.abs(Jf-Jc).max():.3e}  2nd-diff |Rp+Rm-2R|/eps^2={np.abs(Rp+Rm-2*R).max()/eps**2:.3e}")
