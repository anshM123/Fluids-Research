"""Argument-principle count of ALL eigenvalues μ (real and complex) of the linearised self-similar operator in the
rectangle R = [x_lo, x_hi] × [−y_hi, y_hi]:  #zeros of det(I − T_μ) in R = winding number of
f(μ) = Π_{i≤k} (1 − ν_i(μ)) along ∂R (ν_i = the k eigenvalues of T_μ of largest modulus; the remaining ones have
|ν| < |ν_k| ≪ 1 on ∂R and cannot wind around 1). T_μ is real, so f(μ̄) = conj f(μ) and only the upper half of ∂R
is traversed: winding = 2 Δarg_upper / 2π. Adaptive refinement keeps |Δarg| < π/4 between nodes.
usage: python3 stab_contour.py STATE LAM HS TAG X_LO X_HI Y_HI [K]"""
import numpy as np, sys, time
from bq_stability import BQStab, spectrum_c
f, lam, hs, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
x_lo, x_hi, y_hi = float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])
k = int(sys.argv[8]) if len(sys.argv) > 8 else 12
S = BQStab(lam, np.load(f), hs=hs)
log = open(f"stabc_{tag}.log", "w")
def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()
out(f"# contour count {tag}: λ={lam} m={S.m:.10f} hs={hs}; R=[{x_lo},{x_hi}]x[-{y_hi},{y_hi}], k={k}")
cache = {}
def fval(mu):
    key = (round(mu.real, 10), round(mu.imag, 10))
    if key not in cache:
        t = time.time()
        vals, _ = spectrum_c(S, mu, k=k)
        cache[key] = (np.prod(1 - vals), vals)
        big = vals[np.abs(vals) > 0.5]
        out(f"μ={mu.real:.4f}{mu.imag:+.4f}i  f={cache[key][0]:.4e}  |ν_k|={abs(vals[-1]):.3f}  ν(|ν|>0.5): "
            f"{np.array2string(big, precision=4, max_line_width=250)}  ({time.time()-t:.0f}s)")
    return cache[key]
def dang(a, b):
    return np.angle(b / a)
def seg(a, b, fa, fb, depth=0):
    """unwrapped Δarg f from a to b (fa, fb = f values), refining while |Δarg| > π/4"""
    d = dang(fa, fb)
    if abs(d) <= np.pi / 4 or depth >= 6:
        if abs(d) > np.pi / 4:
            out(f"  WARNING: unresolved Δarg={d:.3f} between {a} and {b}")
        return d
    c = 0.5 * (a + b); fc = fval(c)[0]
    return seg(a, c, fa, fc, depth + 1) + seg(c, b, fc, fb, depth + 1)
corners = [complex(x_hi, 0), complex(x_hi, y_hi), complex(x_lo, y_hi), complex(x_lo, 0)]
total = 0.0; kmax = 0.0
for a, b in zip(corners[:-1], corners[1:]):
    n = max(2, int(np.ceil(abs(b - a) / 0.1)))
    pts = [a + (b - a) * t for t in np.linspace(0, 1, n + 1)]
    for p, q in zip(pts[:-1], pts[1:]):
        total += seg(p, q, fval(p)[0], fval(q)[0])
    out(f"edge {a} → {b}: accumulated Δarg/π = {total/np.pi:.4f}")
kmax = max(abs(v[1][-1]) for v in cache.values())
Nz = 2 * total / (2 * np.pi)
out(f"RESULT {tag}: zeros of det(I−T_μ) in R = {Nz:.3f}   (max |ν_k| on ∂R = {kmax:.3f}; {len(cache)} evaluations)")
