"""Scan the real μ axis: eigenvalues ν of T_μ; N(μ) = #{real ν > 1}. Each downward jump of μ past a real
eigenvalue μ* of the linearised operator adds one to N. μ = 1 is the trivial time-translation mode."""
import numpy as np, sys, time
from bq_stability import BQStab
f, lam, hs, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
Nb = int(sys.argv[5]) if len(sys.argv) > 5 else 32
S = BQStab(lam, np.load(f), hs=hs, Nb=Nb)
log = open(f"stab_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# stability scan {tag}: λ={lam} m={S.m:.10f} hs={hs} Nb={Nb}")
mus = np.concatenate([np.arange(1.30, 0.20, -0.05), np.arange(0.20, 0.009, -0.01)])
rows = []
for mu in mus:
    t = time.time()
    vals, _ = S.spectrum(mu, k=10)
    re = vals[np.abs(vals.imag) < 1e-8].real
    cx = vals[np.abs(vals.imag) >= 1e-8]
    Nreal = int(np.sum(re > 1.0))
    Ncx = int(np.sum(np.abs(cx) > 1.0))
    rows.append((mu, Nreal, Ncx))
    out(f"μ={mu:.3f}: N(ν>1)={Nreal}  complex |ν|>1: {Ncx}  top ν: {np.array2string(vals[:6], precision=4, max_line_width=200)}  ({time.time()-t:.0f}s)")
np.save(f"stab_{tag}.npy", np.array(rows))
