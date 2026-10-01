"""λ-continuation of the least-singular self-similar IPM profile (ipm_solver.IPM) downward from λ = 1, in uniform
steps of z = 1/λ (the rung law is linear in 1/λ), recording m(λ) = λ/(1+λ−A) and min V_r/r; smooth profiles ⇔ m = 2.
Same Newton–Krylov (central-difference matvecs) as the Boussinesq scans (bq_scan2.py).
usage: ipm_branch.py LAM_START LAM_END DZ TAG START_Y.npy [Nb] [hs] [s_sw]"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import newton

lam0, lam1, dz, tag, start = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], sys.argv[5]
Nb = int(sys.argv[6]) if len(sys.argv) > 6 else 32
hs = float(sys.argv[7]) if len(sys.argv) > 7 else 0.025
s_sw = float(sys.argv[8]) if len(sys.argv) > 8 else 12.0
acc, rows, t0 = [], [], time.time()
z, zmax, dzc = 1 / lam0, 1 / lam1, dz
log = open(f"ipm_branch_{tag}.log", "a")


def out(msg):
    print(msg, flush=True); log.write(msg + "\n"); log.flush()


out(f"# ipm_branch {tag}: λ {lam0} -> {lam1}, Δz = {dz}, Nb={Nb} hs={hs} s_sw={s_sw} start={start}")
sgn = 1.0 if dz > 0 else -1.0
while sgn * (z - zmax) <= 1e-12:
    lam = 1 / z
    B = IPM(lam, hs=hs, Nb=Nb, s_sw=s_sw)
    if len(acc) >= 2:
        (za, Ya), (zb, Yb) = acc[-2], acc[-1]
        Yg = Yb + (Yb - Ya) * (z - zb) / (zb - za)
    else:
        Yg = acc[-1][1] if acc else np.load(start)
    ok = False
    try:
        Y, info, ok = newton(B, Yg, tol=1e-10, maxit=16, verbose=False, t0=t0, fd='central', pert=1e-6)
    except RuntimeError as e:
        out(f"λ={lam:.7f}: invalid predictor ({e})")
    if not ok:
        if not acc:
            out("failed at the start"); break
        dzc *= 0.5
        out(f"λ={lam:.7f}: no convergence; Δz -> {dzc:.2e}")
        if abs(dzc) < 1e-3:
            break
        z = acc[-1][0] + dzc
        continue
    acc.append((z, Y)); acc = acc[-2:]
    rows.append((lam, info['A'], info['m'], info['vrmin']))
    out(f"λ={lam:.7f} z={z:.4f}: A={info['A']:.10f} m={info['m']:.10f} (m−2={info['m'] - 2:+.4e}) "
        f"vrmin={info['vrmin']:.4f} t={time.time() - t0:.0f}s")
    np.save(f"ipm_branch_{tag}.npy", np.array(rows))
    np.save(f"ipm_br_{tag}_last.npy", Y)
    if len(rows) > 1 and (rows[-2][2] - 2) * (rows[-1][2] - 2) < 0:
        l0, m0 = rows[-2][0], rows[-2][2]; l1, m1 = rows[-1][0], rows[-1][2]
        lc = l0 + (2 - m0) * (l1 - l0) / (m1 - m0)
        out(f"*** m=2 CROSSING near λ ≈ {lc:.7f} (1/λ ≈ {1 / lc:.4f})")
        np.save(f"ipm_br_{tag}_cross_lam{lam:.5f}.npy", Y)
    dzc = sgn * min(abs(dz), 1.5 * abs(dzc))
    z = z + dzc
