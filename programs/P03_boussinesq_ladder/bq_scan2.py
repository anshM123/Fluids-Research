"""λ-continuation of the least-singular self-similar Boussinesq profile with the improved solver (bq_solver.BQ:
finite-difference velocity, implicit Gauss–Legendre march near the stagnation point) and quadratically convergent
Newton–Krylov (central-difference matvecs). Records m(λ) = (λ−1)/(1+λ−A); smooth profiles ⇔ m = 2.
Usage: bq_scan2.py LAM_START LAM_END DLAM TAG START_Y.npy [Nb] [hs] [s_sw] [START_NB START_HS]"""
import numpy as np, sys, time
from bq_solver import BQ
from bq_newton import newton
from bq_regrid import transfer
from bq_logpolar import BQLogPolar

lam0, lam1, tag, start = float(sys.argv[1]), float(sys.argv[2]), sys.argv[4], sys.argv[5]
zstep = sys.argv[3].startswith('z')          # 'z0.1': uniform steps Δz in z = 1/(λ−1) (downward in λ)
dlam = -float(sys.argv[3][1:]) * (lam0 - 1) ** 2 if zstep else float(sys.argv[3])
Nb = int(sys.argv[6]) if len(sys.argv) > 6 else 32
hs = float(sys.argv[7]) if len(sys.argv) > 7 else 0.025
s_sw = float(sys.argv[8]) if len(sys.argv) > 8 else 12.0
snb = int(sys.argv[9]) if len(sys.argv) > 9 else Nb
shs = float(sys.argv[10]) if len(sys.argv) > 10 else hs
dmax = abs(dlam)
dz = float(sys.argv[3][1:]) if zstep else None
acc, rows, t0 = [], [], time.time()
lam = lam0
log = open(f"scan2_{tag}.log", "a")
def out(msg):
    print(msg, flush=True); log.write(msg + "\n"); log.flush()
out(f"# scan2 {tag}: λ {lam0} -> {lam1} step {dlam} Nb={Nb} hs={hs} s_sw={s_sw} start={start}")
while (dlam < 0 and lam >= lam1 - 1e-12) or (dlam > 0 and lam <= lam1 + 1e-12):
    B = BQ(lam, hs=hs, Nb=Nb, s_sw=s_sw)
    if len(acc) >= 2:
        (la, Ya), (lb, Yb) = acc[-2], acc[-1]
        Yg = Yb + (Yb - Ya) * (lam - lb) / (lb - la)
    elif len(acc) == 1:
        Yg = acc[-1][1]
    else:
        Yg = np.load(start)
        if (snb, shs) != (Nb, hs):
            Yg = transfer(BQLogPolar(lam, hs=shs, Nb=snb), Yg, B)
    ok = False
    try:
        Y, info, ok = newton(B, Yg, tol=1e-10, maxit=16, verbose=False, t0=t0, fd='central', pert=1e-6)
    except RuntimeError as e:
        out(f"λ={lam:.6f}: invalid predictor ({e})")
    if not ok:
        if not acc:
            out("failed at the start"); break
        dlam *= 0.5
        out(f"λ={lam:.6f}: no convergence; dλ -> {dlam:.2e}")
        if abs(dlam) < 1e-4:
            break
        lam = acc[-1][0] + dlam
        continue
    acc.append((lam, Y)); acc = acc[-2:]
    rows.append((lam, info['A'], info['m'], info['vrmin']))
    out(f"λ={lam:.6f}: A={info['A']:.10f} m={info['m']:.10f} (m−2={info['m']-2:+.4e}) vrmin={info['vrmin']:.4f} t={time.time()-t0:.0f}s")
    np.save(f"scan2_{tag}.npy", np.array(rows))
    np.save(f"Y2_{tag}_last.npy", Y)
    if len(rows) % 4 == 0:
        np.save(f"Y2_{tag}_lam{lam:.4f}.npy", Y)
    if len(rows) > 1 and (rows[-2][2] - 2) * (rows[-1][2] - 2) < 0:
        l0, m0 = rows[-2][0], rows[-2][2]; l1, m1 = rows[-1][0], rows[-1][2]
        out(f"*** m=2 CROSSING near λ ≈ {l0 + (2 - m0) * (l1 - l0) / (m1 - m0):.7f}")
        np.save(f"Y2_{tag}_cross_lam{lam:.4f}.npy", Y)
    if zstep:                                   # rescale the maximum step to the local Δz
        dmax = dz * (lam - 1) ** 2
    if abs(dlam) < dmax:
        dlam = np.sign(dlam) * min(dmax, 1.5 * abs(dlam))
    else:
        dlam = np.sign(dlam) * dmax
    lam = lam + dlam
