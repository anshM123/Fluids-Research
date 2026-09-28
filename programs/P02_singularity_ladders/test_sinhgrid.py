"""Resolution test of the geometric (sinh) grid at λ=0.458 (thin sonic layer) + one sonic-depth Newton step."""
import numpy as np, sys, time
import scipy.linalg as sla
from mapped_continuation import layer_info, make_grid
from dense_arclength import jac
from ccf_sinhgrid import make_sinh_grid, spacing_profile

d = np.load('pdiag_start.npy'); eta0, phi0, lam = d[0], d[1], float(d[2][0])
Mref = make_grid(-0.963, 0.0015, hs=0.03, pts=24)
target = Mref.normval(Mref.from_other(eta0, phi0), lam)


def fixed_solve(M, phi, lam, target, verbose=False):
    for it in range(25):
        F, den = M.F(phi, lam, target)
        if verbose:
            print(f"    it {it} |F|={np.abs(F).max():.2e}", flush=True)
        if np.abs(F).max() < 1e-11:
            return phi, True
        J, Fl = jac(M, phi, lam, den)
        phi = phi - sla.solve(J, F)
    return phi, False


def delta_step(M, phi, lam, target, frac, verbose=True):
    F, den = M.F(phi, lam, target)
    ip = int(np.argmin(den)); delta = den[ip]; N = M.N
    J, Fl = jac(M, phi, lam, den)
    A = np.zeros((N + 1, N + 1)); A[:N, :N] = J; A[:N, N] = Fl
    A[N, :N] = M.E[ip] * M.Hm[ip, :] * np.exp(phi); A[N, N] = 1.0
    t = sla.solve(A, np.concatenate([np.zeros(N), [1.0]]))
    dt = delta * frac
    ph, lm = phi + (dt - delta) * t[:N], lam + (dt - delta) * t[N]
    for it in range(12):
        F, den = M.F(ph, lm, target)
        gcon = den[ip] - dt
        if verbose:
            print(f"    δ-step it {it}: |F|={np.abs(F).max():.3e} gcon={gcon:.2e} lam={lm:.12f} minden={den.min():.4e}", flush=True)
        if den.min() <= 0:
            return None
        if np.abs(F).max() < 1e-11 and abs(gcon) < 1e-13:
            return ph, lm, M.p_of(ph, lm), dt
        J, Fl = jac(M, ph, lm, den)
        A = np.zeros((N + 1, N + 1)); A[:N, :N] = J; A[:N, N] = Fl
        A[N, :N] = M.E[ip] * M.Hm[ip, :] * np.exp(ph); A[N, N] = 1.0
        dd = sla.solve(A, -np.concatenate([F, [gcon]]))
        ph, lm = ph + dd[:N], lm + dd[N]
    return None


cfgs = [(0.03, 24, 1.0), (0.03, 48, 1.0), (0.02, 36, 1.0), (0.03, 24, 0.5)]
if len(sys.argv) > 1:
    cfgs = cfgs[:int(sys.argv[1])]
for hs, pts, B in cfgs:
    t0 = time.time()
    M = make_sinh_grid(-0.9632, 0.00162, hs=hs, pts=pts, B=B)
    phi = M.from_other(eta0, phi0)
    phi, ok = fixed_solve(M, phi, lam, target)
    F, den = M.F(phi, lam, target)
    e, md, w = layer_info(M, den)
    sp = spacing_profile(M, e, w)
    print(f"hs={hs} pts={pts} B={B} N={M.N}: ok={ok} p={M.p_of(phi, lam):.12f} min_den={md:.8e} layer={e:.6f} w={w:.4e} "
          f"spacing/w at 0,w,3w,10w = {', '.join(f'{x:.3f}' for x in sp)} ({time.time()-t0:.0f}s)", flush=True)
    for frac in (0.99, 0.93):
        r = delta_step(M, phi, lam, target, frac, verbose=True)
        if r is None:
            print(f"   δ-step frac={frac} FAILED", flush=True)
        else:
            print(f"   δ-step frac={frac} converged: δ={r[3]:.6e} lam={r[1]:.12f} p={r[2]:.12f}", flush=True)
