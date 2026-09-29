"""Real-μ scan near μ = 0 (0.06 → 0.005) to test whether eigenvalue crossings there depend on the truncation s_start
of the march (the smooth perturbation's r² coefficient has a 1/μ pole at μ = 0, regularised by the finite s_start).
usage: python3 stab_low.py STATE LAM HS TAG S_START [NB]"""
import numpy as np, sys, time
from bq_stability import BQStab
f, lam, hs, tag, ss = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], float(sys.argv[5])
Nb = int(sys.argv[6]) if len(sys.argv) > 6 else 32
S = BQStab(lam, np.load(f), hs=hs, Nb=Nb, s_start=ss)
log = open(f"stablow_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# low-μ scan {tag}: λ={lam} hs={hs} Nb={Nb} s_start={ss} m={S.m:.10f}")
for mu in np.arange(0.06, 0.004, -0.005):
    t = time.time()
    vals, _ = S.spectrum(mu, k=10)
    re = vals[np.abs(vals.imag) < 1e-8].real
    out(f"μ={mu:.3f}: N(ν>1)={int(np.sum(re > 1))} complex|ν|>1: {int(np.sum(np.abs(vals[np.abs(vals.imag) >= 1e-8]) > 1))}  "
        f"ν near 1: {np.array2string(vals[np.argsort(np.abs(vals - 1))[:3]], precision=4)}  ({time.time()-t:.0f}s)")
