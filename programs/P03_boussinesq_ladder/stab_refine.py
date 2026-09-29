"""Refine the real unstable eigenvalues μ_k: in each bracket of the real-axis scan (stab_scan.py) where the number of
real ν > 1 increases as μ decreases, solve ν(μ) = 1 for the eigenvalue of T_μ closest to 1 (Brent). The bracket
containing μ = 1 (time translation) is refined too, as a check of the method (must return μ = 1).
usage: python3 stab_refine.py STATE LAM HS TAG [MU_MIN=0.06] [NB=32]"""
import numpy as np, sys, time
from scipy.optimize import brentq
from bq_stability import BQStab
f, lam, hs, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
mu_min = float(sys.argv[5]) if len(sys.argv) > 5 else 0.06
Nb = int(sys.argv[6]) if len(sys.argv) > 6 else 32
rows = np.load(f"stab_{tag}.npy")
S = BQStab(lam, np.load(f), hs=hs, Nb=Nb)
log = open(f"stabr_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# refine {tag}: λ={lam} hs={hs} Nb={Nb} m={S.m:.10f}")
def g(mu):
    vals, _ = S.spectrum(mu, k=10)
    re = vals[np.abs(vals.imag) < 1e-8].real
    return re[np.argmin(np.abs(re - 1))] - 1
brackets = [(rows[i + 1, 0], rows[i, 0]) for i in range(len(rows) - 1)
            if rows[i + 1, 1] > rows[i, 1] and rows[i + 1, 0] >= mu_min]
out(f"brackets: {[(round(a, 3), round(b, 3)) for a, b in brackets]}")
roots = []
for a, b in brackets:
    t = time.time()
    ga, gb = g(a), g(b)
    if ga * gb > 0:
        out(f"bracket ({a:.3f},{b:.3f}): no sign change of ν_closest − 1 ({ga:+.3e}, {gb:+.3e}); skipped")
        continue
    r = brentq(g, a, b, xtol=1e-7, rtol=1e-9, maxiter=40)
    roots.append(r)
    out(f"μ* = {r:.6f}   (bracket {a:.3f}–{b:.3f}, {time.time()-t:.0f}s)")
nontriv = [r for r in roots if abs(r - 1) > 1e-3]
out(f"RESULT {tag}: λ={lam}: {len(nontriv)} real unstable eigenvalues μ = {np.round(sorted(nontriv, reverse=True), 5).tolist()}"
    f"; trivial mode at μ = {[round(r, 6) for r in roots if abs(r - 1) <= 1e-3]}")
