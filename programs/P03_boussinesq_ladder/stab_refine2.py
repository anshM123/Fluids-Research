"""Refine one real eigenvalue μ* (ν(μ*) = 1 for the eigenvalue of T_μ closest to 1) in a given bracket, with an
optional deeper origin truncation s_start (base state extended by its exact constant-strain structure).
usage: python3 stab_refine2.py STATE LAM HS TAG MU_A MU_B [S_START=-20]"""
import numpy as np, sys, time
from scipy.optimize import brentq
from bq_solver import BQ
from bq_newton import full
from bq_stability import BQStab
f, lam, hs, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
a, b = float(sys.argv[5]), float(sys.argv[6])
ss = float(sys.argv[7]) if len(sys.argv) > 7 else -20.0
Y = np.load(f)
if ss != -20.0:
    B20 = BQ(lam, hs=hs, Nb=32); Bn = BQ(lam, hs=hs, Nb=32, s_start=ss)
    Y = full(B20, Y)[Bn.i0:]
S = BQStab(lam, Y, hs=hs, s_start=ss)
def g(mu):
    vals, _ = S.spectrum(mu, k=10)
    re = vals[np.abs(vals.imag) < 1e-8].real
    v = re[np.argmin(np.abs(re - 1))] - 1
    print(f"  μ={mu:.7f}: ν−1={v:+.3e}", flush=True)
    return v
t = time.time()
r = brentq(g, a, b, xtol=1e-7, rtol=1e-9)
print(f"RESULT {tag}: μ* = {r:.6f} (λ={lam}, hs={hs}, s_start={ss}, bracket {a}–{b}, {time.time()-t:.0f}s)", flush=True)
