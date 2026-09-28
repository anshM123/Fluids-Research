import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres
from bq_solver import BQ
from bq_newton import residual, full
B = BQ(1.92)
Y = np.load('Y2_lam1.9200_Nb32_hs0.025.npy')
R, info = residual(B, Y)
s = B.s[B.i0:]
mode = -0.5*np.sin(2*B.beta)
def Jv(v, eps=1e-7):
    Rp, _ = residual(B, Y + eps*v); Rm, _ = residual(B, Y - eps*v)
    return (Rp - Rm)/(2*eps)
vA = np.tile(mode, (len(s), 1))
JA = Jv(vA)
def proj(Z, sel): return (Z[sel] @ mode)/(mode @ mode)
for lo, hi in [(-20, -15), (-10, -5), (-3, -1), (0, 2), (5, 10), (30, 40)]:
    sel = (s >= lo) & (s < hi)
    c = proj(JA, sel)
    print(f"J vA: s in [{lo},{hi}): strain coef mean {c.mean():+.6e}  max|JA| {np.abs(JA[sel]).max():.3e}")
# localized strain: taper to zero for s > -5
tap = 0.5*(1 - np.tanh((s + 8)/1.0))
vL = tap[:, None]*vA
JL = Jv(vL)
for lo, hi in [(-20, -15), (-10, -5), (-3, -1), (0, 2), (5, 10)]:
    sel = (s >= lo) & (s < hi)
    print(f"J vL: s in [{lo},{hi}): strain coef mean {proj(JL, sel).mean():+.6e}  max|JL| {np.abs(JL[sel]).max():.3e}")
# GMRES convergence history on J dx = -R
eps = 1e-7
def mv(v):
    v = v.reshape(Y.shape)
    Rp, _ = residual(B, Y + eps*v)
    return ((Rp - R)/eps).ravel()
J = LinearOperator((Y.size, Y.size), matvec=mv, dtype=float)
hist = []
dx, gi = gmres(J, -R.ravel(), rtol=1e-9, atol=0.0, restart=80, maxiter=3, callback=lambda r: hist.append(r), callback_type='pr_norm')
print("gmres info", gi, "n it", len(hist), "rel res history (every 10):", [f"{h:.2e}" for h in hist[::10]], f"last {hist[-1]:.2e}")
dx = dx.reshape(Y.shape)
Rn, infon = residual(B, Y + dx)
print("after full step |R| =", np.abs(Rn).max(), " dA =", infon['A'] - info['A'])
sel = s < -5
print("strain coef of R before", proj(R, sel).mean(), "after", proj(Rn, sel).mean())
print("strain coef of dx (s<-5)", proj(dx, sel).mean(), " dx at far field", np.abs(dx[s > 50]).max())
print("---- step analysis")
print("max|dx| =", np.abs(dx).max(), "at s =", s[np.unravel_index(np.argmax(np.abs(dx)), dx.shape)[0]])
from numpy.polynomial import chebyshev as Ch
xb = np.cos(np.pi*np.arange(B.Nb+1)/B.Nb)
cf = Ch.chebfit(xb, dx.T, B.Nb)
print("Chebyshev coef magnitude of dx by degree:", [f"{v:.1e}" for v in np.abs(cf).max(axis=1)[::4]])
d2 = np.abs(np.diff(dx, 2, axis=0)).max()
print("max second difference in s of dx:", d2, " vs max|dx|", np.abs(dx).max())
Jdx = (J @ dx.ravel()).reshape(Y.shape)
for t in (1e-3, 1e-2, 1e-1, 1.0):
    Rt, _ = residual(B, Y + t*dx)
    nl = Rt - R - t*Jdx
    print(f"t={t:.0e}: |R(Y+t dx)|={np.abs(Rt).max():.3e}  nonlinear part {np.abs(nl).max():.3e}  ratio/t^2 {np.abs(nl).max()/t**2:.3e}")
