"""Argument-principle count of eigenvalues of the linearised self-similar operator (zeros of det(I − T_μ)) inside
an arbitrary rectangle [x_lo, x_hi] × [y_lo, y_hi] of the μ-plane, traversed counter-clockwise. The same adaptive
refinement is used as in stab_contour2.py. The purpose is to close the strips above the symmetric contour boxes
(1.5 ≤ Im μ ≤ Y*), beyond which the spectral radius of T_μ is < 1.
usage: MODEL=bq|ipm python3 contour_rect.py STATE LAM HS TAG X_LO X_HI Y_LO Y_HI [K=16] [S_START=-30]"""
import numpy as np, sys, os, time
from bq_newton import full
if os.environ.get('MODEL', 'ipm') == 'ipm':
    from ipm_solver import IPM as M
    from ipm_stability import IPMStab as St, spectrum_c
else:
    from bq_solver import BQ as M
    from bq_stability import BQStab as St, spectrum_c
f, lam, hs, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
x_lo, x_hi, y_lo, y_hi = (float(v) for v in sys.argv[5:9])
k = int(sys.argv[9]) if len(sys.argv) > 9 else 16
ss = float(sys.argv[10]) if len(sys.argv) > 10 else -30.0
Y = np.load(f)
if ss != -20.0:
    Y = full(M(lam, hs=hs, Nb=32), Y)[M(lam, hs=hs, Nb=32, s_start=ss).i0:]
S = St(lam, Y, hs=hs, s_start=ss)
log = open(f"crect_{tag}.log", "w")


def out(s):
    print(s, flush=True); log.write(s + "\n"); log.flush()


out(f"# rectangle count {tag}: λ={lam} m={S.m:.10f} hs={hs} s_start={ss}; R=[{x_lo},{x_hi}]x[{y_lo},{y_hi}], k={k}")
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
        out(f"μ={mu.real:.5f}{mu.imag:+.5f}i ρ={np.abs(vals).max():.3f} min|1−ν|={np.abs(1 - vals).min():.3f} ({time.time() - t:.0f}s)")
    return cache[key]


unresolved = []


def seg(a, b, fa, fb, depth=0):
    d = np.angle(fb / fa)
    if abs(d) <= np.pi / 6 or depth >= 8:
        if abs(d) > np.pi / 6:
            unresolved.append((a, b, d)); out(f"  UNRESOLVED Δarg={d:.3f} between {a:.5f} and {b:.5f}")
        return d
    c = 0.5 * (a + b); fc = fval(c)[0]
    return seg(a, c, fa, fc, depth + 1) + seg(c, b, fc, fb, depth + 1)


corners = [complex(x_hi, y_lo), complex(x_hi, y_hi), complex(x_lo, y_hi), complex(x_lo, y_lo), complex(x_hi, y_lo)]
total = 0.0
for a, b in zip(corners[:-1], corners[1:]):
    n = max(2, int(np.ceil(abs(b - a) / 0.15)))
    pts = [a + (b - a) * t for t in np.linspace(0, 1, n + 1)]
    for p, q in zip(pts[:-1], pts[1:]):
        total += seg(p, q, fval(p)[0], fval(q)[0])
    out(f"edge {a} → {b}: accumulated Δarg/2π = {total / (2 * np.pi):.4f}")
rho = max(np.abs(v[1]).max() for v in cache.values())
out(f"RESULT {tag}: zeros in R = {total / (2 * np.pi):.4f} (nearest integer {int(round(total / (2 * np.pi)))}); "
    f"max ρ(T_μ) on ∂R = {rho:.3f}; unresolved jumps: {len(unresolved)}; {len(cache)} evaluations")
