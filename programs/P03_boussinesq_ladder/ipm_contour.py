"""IPM version of stab_contour2.py (same algorithm; IPM profile and IPMStab). Right-half-plane eigenvalue count: argument principle for det(I − T_μ) on the rectangle
[x_lo, x_hi] × [−y_hi, y_hi] (upper half traversed, conjugate symmetry), with
  * k = 24 eigenvalues of T_μ per point (eigenvalue swaps at the truncation then involve |ν| ≲ 0.3 and cause only
    small phase errors),
  * adaptive refinement to |Δarg| < π/6 (depth ≤ 8), warm-started Arnoldi,
  * an optional deeper origin truncation s_start (base state extended by its exact constant-strain structure),
  * a diagnostic at any unresolved jump: the eigenvalue of T_μ closest to 1 there (a zero of det(I − T_μ) close to
    the contour shows up as |1 − ν| → 0).
usage: python3 stab_contour2.py STATE LAM HS TAG X_LO X_HI Y_HI [K=24] [S_START=-20]"""
import numpy as np, sys, time
from ipm_solver import IPM as BQ
from bq_newton import full
from ipm_stability import IPMStab as BQStab, spectrum_c
f, lam, hs, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
x_lo, x_hi, y_hi = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
k = int(sys.argv[8]) if len(sys.argv) > 8 else 16
ss = float(sys.argv[9]) if len(sys.argv) > 9 else -20.0
Y = np.load(f)
if ss != -20.0:
    B20 = BQ(lam, hs=hs, Nb=32); Bn = BQ(lam, hs=hs, Nb=32, s_start=ss)
    Y = full(B20, Y)[Bn.i0:]
S = BQStab(lam, Y, hs=hs, s_start=ss)
log = open(f"ipmc_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# contour v2 {tag}: λ={lam} m={S.m:.10f} hs={hs} s_start={ss}; R=[{x_lo},{x_hi}]x[-{y_hi},{y_hi}], k={k}")
cache = {}; last = {'v': None}
def fval(mu):
    key = (round(mu.real, 12), round(mu.imag, 12))
    if key not in cache:
        t = time.time()
        try:
            vals, vecs = spectrum_c(S, mu, k=k, v0=last['v'])
        except Exception:
            vals, vecs = spectrum_c(S, mu, k=k)
        last['v'] = vecs[:, 0]
        cache[key] = (np.prod(1 - vals), vals)
        near = vals[np.argmin(np.abs(vals - 1))]
        out(f"μ={mu.real:.5f}{mu.imag:+.5f}i f={cache[key][0]:.4e} |ν_k|={abs(vals[-1]):.3f} ν≈1: {near:.4f} "
            f"|ν|>0.5: {np.array2string(vals[np.abs(vals) > 0.5], precision=3, max_line_width=250)} ({time.time()-t:.0f}s)")
    return cache[key]
unresolved = []
def seg(a, b, fa, fb, depth=0):
    d = np.angle(fb / fa)
    if abs(d) <= np.pi / 6 or depth >= 8:
        if abs(d) > np.pi / 6:
            vals = fval(b)[1]
            unresolved.append((a, b, d, np.min(np.abs(vals - 1))))
            out(f"  UNRESOLVED Δarg={d:.3f} between {a:.6f} and {b:.6f}; min|1−ν| there = {np.min(np.abs(vals-1)):.2e}")
        return d
    c = 0.5 * (a + b); fc = fval(c)[0]
    return seg(a, c, fa, fc, depth + 1) + seg(c, b, fc, fb, depth + 1)
corners = [complex(x_hi, 0), complex(x_hi, y_hi), complex(x_lo, y_hi), complex(x_lo, 0)]
total = 0.0
for ie, (a, b) in enumerate(zip(corners[:-1], corners[1:])):
    n = max(2, int(np.ceil(abs(b - a) / (0.25 if ie < 2 else 0.1))))   # right/top edges are quiet (|ν| ≪ 1)
    pts = [a + (b - a) * t for t in np.linspace(0, 1, n + 1)]
    for p, q in zip(pts[:-1], pts[1:]):
        total += seg(p, q, fval(p)[0], fval(q)[0])
    out(f"edge {a} → {b}: accumulated Δarg/π = {total/np.pi:.4f}")
kmax = max(abs(v[1][-1]) for v in cache.values())
Nz = 2 * total / (2 * np.pi)
out(f"RESULT {tag}: zeros in R = {Nz:.4f} (nearest integer {int(round(Nz))}, deviation {abs(Nz-round(Nz)):.3f}); "
    f"max |ν_k| on ∂R = {kmax:.3f}; unresolved jumps: {len(unresolved)}; {len(cache)} evaluations")
