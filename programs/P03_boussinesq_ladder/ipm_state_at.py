"""IPM branch states at prescribed λ (for spectra at branch points between rungs): continuation in z = 1/λ from a
start state in steps ≤ DZ, saving the converged state at every target.
usage: python3 ipm_state_at.py START.npy LAM_START DZ LAM_TARGET [LAM_TARGET ...]   → ipm_bp_lam{λ:.7f}.npy"""
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton

Y, lam_s, dz = np.load(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
targets = sorted((float(v) for v in sys.argv[4:]), key=lambda l: abs(1 / l - 1 / lam_s))
z = 1 / lam_s; t0 = time.time(); prev = None
for lt in targets:
    zt = 1 / lt
    while abs(zt - z) > 1e-12:
        zn = z + np.sign(zt - z) * min(dz, abs(zt - z))
        Yg = Y if prev is None else Y + (Y - prev[1]) * (zn - z) / (z - prev[0])
        B = IPM(1 / zn, hs=0.025, Nb=32, s_sw=12.0)
        Yn, info, ok = newton(B, Yg, tol=1e-10, maxit=16, verbose=False, fd='central', pert=1e-6)
        if not ok:
            raise SystemExit(f"no convergence at λ = {1 / zn}")
        prev = (z, Y); Y, z = Yn, zn
    np.save(f"ipm_bp_lam{lt:.7f}.npy", Y)
    print(f"λ = {lt:.10f} (z = {zt:.4f}): m = {info['m']:.10f}, A = {info['A']:.10f}, vrmin = {info['vrmin']:.4f} "
          f"({time.time() - t0:.0f}s)", flush=True)
