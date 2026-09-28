"""Diagnose one sonic-depth continuation step (verbose Newton)."""
import numpy as np, sys, time
import scipy.linalg as sla
from mapped_continuation import layer_info, make_grid
from dense_arclength import jac

d = np.load('pdiag_start.npy'); eta0, phi0, lam = d[0], d[1], float(d[2][0])
M = make_grid(-0.963, 0.0015, hs=0.03, pts=24)
phi = M.from_other(eta0, phi0); target = M.normval(phi, lam)
for it in range(20):
    F, den = M.F(phi, lam, target)
    print("fixed it", it, np.abs(F).max(), flush=True)
    if np.max(np.abs(F)) < 1e-11:
        break
    J, Fl = jac(M, phi, lam, den); phi = phi - sla.solve(J, F)
F, den = M.F(phi, lam, target)
ip = int(np.argmin(den)); delta = den[ip]
N = M.N
print("N", N, "delta", delta, "ip", ip, "eta", M.eta[ip])
# tangent d(phi,lam)/d delta
J, Fl = jac(M, phi, lam, den)
A = np.zeros((N + 1, N + 1)); A[:N, :N] = J; A[:N, N] = Fl
A[N, :N] = M.E[ip] * M.Hm[ip, :] * np.exp(phi); A[N, N] = 1.0
print("cond est (1-norm) of bordered:", np.linalg.cond(A, 1), flush=True)
t = sla.solve(A, np.concatenate([np.zeros(N), [1.0]]))
print("dlam/ddelta =", t[N], " |dphi/ddelta|max =", np.abs(t[:N]).max(), flush=True)
# also d lam / d(normalization) sanity: solve J x = -Fl
x = sla.solve(J, -Fl)
dden_dlam_along = 1.0 + M.E[ip] * (M.Hm[ip, :] * np.exp(phi)) @ x
print("d den_p / d lam along branch =", dden_dlam_along, flush=True)
for frac in [0.99, 0.97, 0.93]:
    dt = delta * frac
    ph, lm = phi + (dt - delta) * t[:N], lam + (dt - delta) * t[N]
    print(f"--- target δ={dt:.5e} (frac {frac}), predictor lam={lm:.10f}")
    for it in range(12):
        F, den = M.F(ph, lm, target)
        gcon = den[ip] - dt
        print(f"  it {it}: |F|={np.abs(F).max():.3e} gcon={gcon:.3e} lam={lm:.12f} minden={den.min():.4e} argmin={np.argmin(den)-ip}", flush=True)
        if den.min() <= 0: break
        if np.abs(F).max() < 1e-10 and abs(gcon) < 1e-12: break
        J, Fl = jac(M, ph, lm, den)
        A = np.zeros((N + 1, N + 1)); A[:N, :N] = J; A[:N, N] = Fl
        A[N, :N] = M.E[ip] * M.Hm[ip, :] * np.exp(ph); A[N, N] = 1.0
        dd = sla.solve(A, -np.concatenate([F, [gcon]]))
        ph, lm = ph + dd[:N], lm + dd[N]
