"""Independent reproduction of a smooth profile with the global Newton solver (bq_global.py) for a sequence of
far-field truncations s_max, followed by geometric extrapolation s_max → ∞.
usage: python3 glob_run.py SRC_FILE LAM0 TAG HS NB SRC_HS SMIN SMAX1,SMAX2,..."""
import numpy as np, sys, time
from bq_global import BQGlobal
f, lam0, tag = sys.argv[1], float(sys.argv[2]), sys.argv[3]
hs, Nb, src_hs, smin = float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
smaxs = [float(v) for v in sys.argv[8].split(',')]
log = open(f"glob_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# global solver {tag}: march λ={lam0} src={f} hs={hs} Nb={Nb} s_min={smin}")
lams = []
for smax in smaxs:
    t = time.time()
    G = BQGlobal(s_min=smin, s_max=smax, hs=hs, Nb=Nb)
    U = G.initial_from_march(f, lam0, hs_src=src_hs)
    U, ok = G.newton(U, tol=1e-10, maxit=15, verbose=False)
    R = np.abs(G.residual(U)).max()
    lams.append(U[-1])
    out(f"s_max={smax:5.1f}: λ = {U[-1]:.10f}  |R|={R:.1e} ok={ok}  unknowns={len(U)}  ({time.time()-t:.0f}s)")
    np.save(f"glob_{tag}_smax{smax:g}.npy", U)
if len(lams) >= 3:
    a, b, c = lams[-3:]
    r = (c - b) / (b - a)
    ext = c + (c - b) * r / (1 - r)
    out(f"geometric extrapolation (last three, ratio {r:.3f}): λ_∞ = {ext:.10f}   (march λ = {lam0:.10f}, diff {ext-lam0:+.2e})")
