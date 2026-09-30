"""P6 check (PREDICTIONS_IPM.md): where are the unstable eigenmodes of an IPM smooth profile localized?
For each unstable μ* of the profile, the eigenvector of T_μ* with ν ≈ 1 gives X'. The density perturbation r is
then marched as in IPMStab.Tc. We report where |r| peaks along the boundary and over the domain, against the dip of
the self-similar radial speed D = V_r/r along the boundary (the front of the stalled layer).
usage: python3 ipm_modes.py STATE LAM MU1 [MU2 ...]"""
import numpy as np, sys
from bq_newton import full
from ipm_solver import IPM
from ipm_stability import IPMStab, spectrum_c, _lin_march_r

f, lam = sys.argv[1], float(sys.argv[2]); mus = [float(v) for v in sys.argv[3:]]
ss = -30.0
Y = full(IPM(lam, hs=0.025, Nb=32), np.load(f))[IPM(lam, hs=0.025, Nb=32, s_start=ss).i0:]
S = IPMStab(lam, Y, hs=0.025, s_start=ss); B = S.B
jb = int(np.argmax(B.cb))                         # β = 0 (the boundary): cos β = 1
Dw = S.D0[:, jb]
sel = (B.s > -6) & (B.s < 4)
k_dip = np.argmin(np.where(sel, Dw, np.inf))
print(f"IPM λ={lam}: boundary radial speed D = V_r/r: D(origin) = {Dw[B.i0 + 5]:.4f}, dip D = {Dw[k_dip]:.4f} at "
      f"s = {B.s[k_dip]:.3f} (x = {np.exp(B.s[k_dip]):.3f}); D rises past 2×dip at s = "
      f"{B.s[k_dip + np.argmax(Dw[k_dip:] > 2 * Dw[k_dip])]:.3f}", flush=True)
for mu in mus:
    v, V = spectrum_c(S, complex(mu), k=6)
    i = int(np.argmin(np.abs(v - 1))); x = V[:, i]
    Yp = np.asarray(x, dtype=complex).reshape(B.Ns - B.i0, B.Nb + 1)
    Fr = S._forcing(full(B, Yp.real)); Fi = S._forcing(full(B, Yp.imag))
    F = {kk: Fr[kk] + 1j * Fi[kk] for kk in Fr}
    sh = lambda Z: np.vstack([Z[1:], Z[-1:]])
    r = _lin_march_r(B.s, B.hs, complex(mu), lam, B.Db, S.i0, S.jsw, S.D0, S.w0, S.D1, S.w1, S.D2, S.w2,
                     S.Dh, S.wh, S.Dg1, S.wg1, F['0'], F['1'], F['2'], F['h'], sh(F['0']))
    # r relative to the base density R (the modes are regular: r/R bounded at the origin)
    rw = np.abs(r[:, jb]); Rw = np.abs(S.Th[:, jb]) + 1e-300
    kk = sel & (B.s > B.s[B.i0] + 2)
    k1 = np.argmax(np.where(kk, rw / Rw, -1)); k2 = np.argmax(np.where(sel, np.abs(r).max(axis=1), -1))
    ss_ = [-20, -10, -5, -2, -1, -0.5, -0.25, 0, 0.5, 1, 2, 4]
    idx = [int(np.argmin(np.abs(B.s - q))) for q in ss_]
    nr = rw / rw[sel].max(); nR = Rw / Rw[sel].max()
    print("     s:     " + " ".join(f"{q:7.2f}" for q in ss_))
    print("     |r|:   " + " ".join(f"{nr[k]:7.3f}" for k in idx))
    print("     |R|:   " + " ".join(f"{nR[k]:7.3f}" for k in idx))
    print("     |r/R|: " + " ".join(f"{(rw[k] / Rw[k]) / (rw[k1] / Rw[k1]):7.3f}" for k in idx))
    print(f"  μ = {mu}: ν = {v[i]:.6f}; |r/R| on the boundary peaks at s = {B.s[k1]:.3f}; max_β |r| peaks at "
          f"s = {B.s[k2]:.3f}; |r/R|(origin side, s = {B.s[B.i0 + 200]:.1f}) / peak = {rw[B.i0 + 200] / Rw[B.i0 + 200] / (rw[k1] / Rw[k1]):.2e}",
          flush=True)
