"""Smooth Hou–Luo profiles (m = 2) at each rung: regula falsi on m(λ) − 2 with Newton at each λ, starting from the
saved scan state closest in z. Brackets from the scan logs.
usage: python3 hl_cross.py N L1 L2 PREFIX OUTTAG  (states PREFIX_lam*.npy; brackets = all CROSSING lines of the log)"""
import numpy as np, sys, re, glob, time
from hl_solver import HL, newton
N, L1, L2, prefix, tag = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], sys.argv[5]
logname = sys.argv[6]
states = {float(re.search(r'lam([0-9.]+?)\.npy', f).group(1)): f for f in glob.glob(f"{prefix}_lam*.npy")}
br = [tuple(map(float, re.findall(r"([0-9]+\.[0-9]+)", l.split("between")[1]))) for l in open(logname) if "CROSSING" in l]
out = open(f"hl_cross_{tag}.log", "a")
def say(s):
    print(s, flush=True); out.write(s + "\n"); out.flush()
def solve(lam, qg):
    S = HL(lam, L1=L1, L2=L2, N=N)
    q, info, ok = newton(S, qg, tol=1e-12, maxit=20, verbose=False)
    return q, info['m'] - 2, ok
for (la, lb) in br:
    zc = 0.5 * (1 / (la - 1) + 1 / (lb - 1))
    ls = min(states, key=lambda l: abs(1 / (l - 1) - zc))
    q = np.load(states[ls])
    if len(q) != N:
        S0 = HL(ls, L1=40, L2=160, N=len(q)); S1 = HL(ls, L1=L1, L2=L2, N=N); q = np.interp(S1.eta, S0.eta, q)
    # walk from the state to the bracket in small z steps
    z0, za = 1 / (ls - 1), 1 / (la - 1)
    nstep = int(np.ceil(abs(za - z0) / 0.05))
    ok = True
    for zz in np.linspace(z0, za, nstep + 1)[1:]:
        q, fa, ok = solve(1 + 1 / zz, q)
        if not ok:
            break
    if not ok:
        say(f"bracket {la}-{lb}: walk failed"); continue
    qa = q; qb, fb, ok = solve(lb, q)
    if fa * fb > 0:
        say(f"bracket {la}-{lb}: no sign change ({fa:+.2e}, {fb:+.2e})"); continue
    for it in range(40):
        lc = lb - fb * (lb - la) / (fb - fa)
        qc, fc, ok = solve(lc, qa if abs(lc - la) < abs(lc - lb) else qb)
        if fc * fa < 0:
            lb, fb, qb = lc, fc, qc; fa *= 0.5
        else:
            la, fa, qa = lc, fc, qc; fb *= 0.5
        if abs(fc) < 1e-12 or abs(lb - la) < 1e-12:
            break
    np.save(f"hl_cross_{tag}_lam{lc:.8f}.npy", qc)
    say(f"CROSSING HL N={N}: λ = {lc:.10f}  z = {1/(lc-1):.6f}  m−2 = {fc:+.1e}")
