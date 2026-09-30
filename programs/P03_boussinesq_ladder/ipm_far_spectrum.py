"""Spectral radius of T_μ far from the real axis (for the large-|Im μ| part of the index count): the eigenvalues of
T_μ of largest modulus along Re μ = X for Im μ = Y₁ … Y_k. Where ρ(T_μ) < 1, det(I − T_μ) ≠ 0 and μ is not an eigenvalue.
usage: python3 ipm_far_spectrum.py STATE LAM X Y1 [Y2 ...]   (model: ipm or bq, via env MODEL; default ipm)"""
import numpy as np, sys, os, time
from bq_newton import full
if os.environ.get('MODEL', 'ipm') == 'ipm':
    from ipm_solver import IPM as M
    from ipm_stability import IPMStab as St, spectrum_c
else:
    from bq_solver import BQ as M
    from bq_stability import BQStab as St, spectrum_c
f, lam, X = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
Ys = [float(v) for v in sys.argv[4:]]
ss = float(os.environ.get('S_START', -30)); hs = float(os.environ.get('HS', 0.025))
Y = np.load(f)
if ss != -20.0:
    Y = full(M(lam, hs=hs, Nb=32), Y)[M(lam, hs=hs, Nb=32, s_start=ss).i0:]
S = St(lam, Y, hs=hs, s_start=ss)
v0 = None
for y in Ys:
    t = time.time()
    vals, vecs = spectrum_c(S, complex(X, y), k=8, v0=v0); v0 = vecs[:, 0]
    print(f"μ = {X} + {y}i: |ν| top = {np.round(np.abs(vals[:4]), 3).tolist()}, min|1−ν| = {np.abs(1 - vals).min():.3f} "
          f"({time.time() - t:.0f}s)", flush=True)
