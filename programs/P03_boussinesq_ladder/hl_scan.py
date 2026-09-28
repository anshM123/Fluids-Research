"""λ-continuation of the least-singular Hou–Luo self-similar family; m(λ) − 2 vs z = 1/(λ−1).
Usage: hl_scan.py LAM_START LAM_END DZ TAG [N] [L1 L2]   (DZ > 0: steps uniform in z, downward in λ;
       a negative DZ is interpreted as a λ step)"""
import numpy as np, sys, time
from hl_solver import HL, newton
lam0, lam1, dz, tag = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
N = int(sys.argv[5]) if len(sys.argv) > 5 else 8192
L1 = float(sys.argv[6]) if len(sys.argv) > 6 else 40.0
L2 = float(sys.argv[7]) if len(sys.argv) > 7 else 160.0
start = sys.argv[8] if len(sys.argv) > 8 and sys.argv[8] != '-' else None      # optional start state (interpolated)
eta_start = float(sys.argv[9]) if len(sys.argv) > 9 else -30.0
sL1 = float(sys.argv[10]) if len(sys.argv) > 10 else 40.0        # domain of the start state
sL2 = float(sys.argv[11]) if len(sys.argv) > 11 else 160.0
log = open(f"hl_scan_{tag}.log", "a")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# HL scan {tag}: {lam0} -> {lam1}, dz={dz}, N={N}, L1={L1}, L2={L2}")
acc, rows, t0 = [], [], time.time()
lam = lam0; fac = 1.0
while (lam >= lam1 - 1e-12) if (dz > 0 or lam1 < lam0) else (lam <= lam1 + 1e-12):
    S = HL(lam, L1=L1, L2=L2, N=N, eta_start=eta_start)
    if len(acc) >= 2:
        (la, qa), (lb, qb) = acc[-2], acc[-1]
        qg = qb + (qb - qa) * (lam - lb) / (lb - la)
    elif acc:
        qg = acc[-1][1]
    elif start is not None:
        q8 = np.load(start); S8 = HL(lam, N=len(q8), L1=sL1, L2=sL2)
        qg = np.interp(S.eta, S8.eta, q8)
    else:
        qg = S.guess((3 + lam) / 2)
    q, info, ok = newton(S, qg, tol=1e-12, maxit=15, verbose=False)
    if not ok:
        fac *= 0.5
        out(f"λ={lam:.8f}: no convergence (|R| stalled); step factor -> {fac}")
        if fac < 1e-3 or not acc:
            break
        lb = acc[-1][0]
        lam = (1 + 1 / (1 / (lb - 1) + dz * fac)) if dz > 0 else lb + dz * fac
        continue
    acc.append((lam, q)); acc = acc[-2:]
    rows.append((lam, info['A'], info['m'], info['Dmin']))
    np.save(f"hl_scan_{tag}.npy", np.array(rows))
    out(f"λ={lam:.8f} z={1/(lam-1):.5f}: A={info['A']:.12f} m={info['m']:.12f} (m−2={info['m']-2:+.5e}) Dmin/ε={info['Dmin']/info['eps']:.5f} t={time.time()-t0:.0f}s")
    np.save(f"hl_q_{tag}_last.npy", q)
    if len(rows) > 1 and (rows[-2][2] - 2) * (rows[-1][2] - 2) < 0:
        out(f"*** m=2 CROSSING between λ={rows[-2][0]:.8f} and {rows[-1][0]:.8f}")
    if len(rows) % 10 == 0:
        np.save(f"hl_q_{tag}_lam{lam:.6f}.npy", q)
    fac = min(1.0, fac * 1.5)
    lam = (1 + 1 / (1 / (lam - 1) + dz * fac)) if dz > 0 else lam + abs(dz) * fac * (1 if lam1 > lam0 else -1)
