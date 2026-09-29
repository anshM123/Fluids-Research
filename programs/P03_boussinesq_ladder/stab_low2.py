"""Low-μ real scan with the march started deeper (s_start = −30 instead of −20): genuine eigenvalues stay put, the
truncation-induced crossing (≈ 0.66/|s_start|) moves. The base state is the s_start = −20 solution extended by its
exact constant-strain structure below s = −20 (differences O(e^{−20})).
usage: python3 stab_low2.py STATE LAM HS TAG SS_NEW MU_HI MU_LO DMU"""
import numpy as np, sys, time
from bq_solver import BQ
from bq_newton import full
from bq_stability import BQStab
f, lam, hs, tag, ss = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], float(sys.argv[5])
mu_hi, mu_lo, dmu = float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])
B20 = BQ(lam, hs=hs, Nb=32); Bn = BQ(lam, hs=hs, Nb=32, s_start=ss)
Y = full(B20, np.load(f))[Bn.i0:]
S = BQStab(lam, Y, hs=hs, Nb=32, s_start=ss)
log = open(f"stablow_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# low-μ scan {tag}: λ={lam} hs={hs} s_start={ss} m={S.m:.10f}")
for mu in np.arange(mu_hi, mu_lo - 1e-9, -dmu):
    t = time.time()
    vals, _ = S.spectrum(mu, k=10)
    re = vals[np.abs(vals.imag) < 1e-8].real
    out(f"μ={mu:.3f}: N(ν>1)={int(np.sum(re > 1))} complex|ν|>1: {int(np.sum(np.abs(vals[np.abs(vals.imag) >= 1e-8]) > 1))}  "
        f"ν near 1: {np.array2string(vals[np.argsort(np.abs(vals - 1))[:2]], precision=4)}  ({time.time()-t:.0f}s)")
