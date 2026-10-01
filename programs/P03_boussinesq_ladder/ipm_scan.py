"""Fixed-z scan of the IPM branch for rung location below the noise floor of a single solve.
Runs through z = Z0, Z0+DZ, …, Z1 (identical λ values for every s_start, so that runs at two origin cut-offs can
be combined point by point: m_∞ = m(s₁) − [m(s₂) − m(s₁)]/(e^{s₂−s₁} − 1), removing the O(e^{s_start}) bias).
Tolerance 1e-13 (the floor at s_start = −16/−14 is ≈ 1e-13).  Saves states every SAVE_EVERY points (not committed).
usage: ipm_scan.py START.npy LAM_START Z0 Z1 DZ SS TAG [hs]   (START on the hs grid with s_start = −20)"""
import numpy as np, sys, time, os
from ipm_solver import IPM
from bq_newton import newton, residual
from bq_regrid import transfer
from bq_logpolar import BQLogPolar
f, lam_s, z0, z1, dz, ss, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
hs = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0125
every = int(os.environ.get('SAVE_EVERY', 0))
t0 = time.time(); log = open(f"ipm_scan_{tag}.out", "a")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# ipm_scan {tag}: z {z0} -> {z1} step {dz}, s_start={ss}, hs={hs}, Nb={os.environ.get('NB', 32)}, "
    f"tol={os.environ.get('NK_TOL', 1e-13)}, start={f} (λ={lam_s})")
NB = int(os.environ.get('NB', 32))                                   # angular resolution of this run
NB_IN, HS_IN, SS_IN = int(os.environ.get('NB_IN', 32)), float(os.environ.get('HS_IN', hs)), float(os.environ.get('SS_IN', -20.0))
Y = transfer(BQLogPolar(lam_s, hs=HS_IN, Nb=NB_IN, s_start=SS_IN), np.load(f), IPM(lam_s, hs=hs, Nb=NB, s_sw=12.0, s_start=ss))
zs = np.arange(z0, z1 + 1e-9, dz); acc = []; rows = []
# resume after a container restart: skip the z already done and restart from the last saved state(s)
if os.path.exists(f"ipm_scan_{tag}.npy") and os.environ.get('RESUME', '1') == '1':
    done = np.load(f"ipm_scan_{tag}.npy"); rows = [tuple(r) for r in done]
    for zd in done[-2:, 0]:
        fz = f"ipm_scan_{tag}_z{zd:.4f}.npy"
        if os.path.exists(fz):
            acc.append((zd, np.load(fz)))
    if acc:
        out(f"# resumed after z = {done[-1, 0]:.4f} ({len(rows)} points done, {len(acc)} states reloaded)")
        zs = zs[zs > done[-1, 0] + 1e-9]
    else:
        rows = []
for i, z in enumerate(zs):
    lam = 1 / z
    B = IPM(lam, hs=hs, Nb=NB, s_sw=12.0, s_start=ss)
    if len(acc) >= 2:
        (za, Ya), (zb, Yb) = acc[-2], acc[-1]; Yg = Yb + (Yb - Ya) * (z - zb) / (zb - za)
    else:
        Yg = acc[-1][1] if acc else Y
    Y, info, ok = newton(B, Yg, tol=float(os.environ.get("NK_TOL", 1e-13)), maxit=12, verbose=False, fd="central", pert=1e-6)
    R, _ = residual(B, Y); nR = np.abs(R).max()
    if nR > 1e-10:
        out(f"z={z:.4f}: NOT CONVERGED |R|={nR:.1e}; stop"); break
    acc.append((z, Y)); acc = acc[-2:]
    rows.append((z, lam, info['A'], info['m'] - 2, nR, info['vrmin']))
    out(f"z={z:.4f} λ={lam:.9f} |R|={nR:.1e} A={info['A']:.14f} m-2={info['m']-2:+.6e} vrmin={info['vrmin']:.5f} t={time.time()-t0:.0f}s")
    np.save(f"ipm_scan_{tag}.npy", np.array(rows))
    if every and i % every == 0:
        np.save(f"ipm_scan_{tag}_z{z:.4f}.npy", Y)
