"""Continuation in λ of the least-singular self-similar Boussinesq family; records the local exponent
m(λ) = (λ−1)/(1+λ−A(λ)) at the stagnation point. Smooth self-similar profiles ⇔ m(λ) = 2.
Usage: bq_scan.py LAM_START LAM_END DLAM TAG [Nb] [hs] [start_Y.npy]"""
import numpy as np, sys, time
from bq_logpolar import BQLogPolar
from bq_newton import newton, initial_guess

lam0 = float(sys.argv[1]); lam1 = float(sys.argv[2]); dlam = float(sys.argv[3]); tag = sys.argv[4]
Nb = int(sys.argv[5]) if len(sys.argv) > 5 else 32
hs = float(sys.argv[6]) if len(sys.argv) > 6 else 0.025
start = sys.argv[7] if len(sys.argv) > 7 else None

acc = []                      # accepted (λ, Y)
rows, t0 = [], time.time()
lam = lam0
while (dlam < 0 and lam >= lam1 - 1e-12) or (dlam > 0 and lam <= lam1 + 1e-12):
    B = BQLogPolar(lam, s_min=-120, s_max=100, hs=hs, Nb=Nb)
    if len(acc) >= 2:
        (la, Ya), (lb, Yb) = acc[-2], acc[-1]
        Yg = Yb + (Yb - Ya) * (lam - lb) / (lb - la)
    elif len(acc) == 1:
        Yg = acc[-1][1]
    else:
        Yg = np.load(start) if start else initial_guess(B, (3 + lam) / 2)
    ok = False
    try:
        Y, info, ok = newton(B, Yg, tol=1e-5, maxit=12, verbose=False, t0=t0)
    except RuntimeError as e:
        print(f"λ={lam:.5f}: invalid predictor ({e})", flush=True)
    if not ok:
        if not acc:
            print("failed at the start", flush=True); break
        dlam *= 0.5
        print(f"λ={lam:.5f}: no convergence; dλ -> {dlam:.2e}", flush=True)
        if abs(dlam) < 2e-4:
            break
        lam = acc[-1][0] + dlam
        continue
    acc.append((lam, Y))
    acc = acc[-2:]
    rows.append((lam, info['A'], info['m'], info['vrmin']))
    print(f"λ={lam:.5f}: A={info['A']:.8f} m={info['m']:.8f} (m−2={info['m']-2:+.3e}) vrmin={info['vrmin']:.4f} "
          f"t={time.time()-t0:.0f}s", flush=True)
    np.save(f"scan_{tag}.npy", np.array(rows))
    if len(rows) % 5 == 0:
        np.save(f"Y_{tag}_lam{lam:.4f}.npy", Y)
    if len(rows) > 1 and (rows[-2][2] - 2) * (rows[-1][2] - 2) < 0:
        l0, m0 = rows[-2][0], rows[-2][2]; l1, m1 = rows[-1][0], rows[-1][2]
        print(f"*** m=2 CROSSING near λ ≈ {l0 + (2 - m0) * (l1 - l0) / (m1 - m0):.6f}", flush=True)
        np.save(f"Y_{tag}_cross_lam{lam:.4f}.npy", Y)
    lam = lam + dlam
