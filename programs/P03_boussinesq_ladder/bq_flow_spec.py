"""2D Boussinesq: lower unstable spectrum at branch points between the smooth rungs (spectral flow along the continuous
branch). Same method as hl_deep_spec.py: parity flips of the count of real ν > 1 of T_μ on a grid in ν̂ = μ/(λ−1),
refined by bisection on the parity. The state is extended to the origin truncation S_START by its exact
constant-strain structure (as in stab_contour2.py), which moves the truncation artifact (≈ 0.6/|s_start|) down.
usage: python3 bq_flow_spec.py STATE LAM HS NUHAT_MAX STEP [S_START=-30] [NUHAT_MIN=0.2]"""
import numpy as np, sys, time
from bq_solver import BQ
from bq_newton import full
from bq_stability import BQStab, spectrum_c

f, lam, hs = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
nu_max, step = float(sys.argv[4]), float(sys.argv[5])
ss = float(sys.argv[6]) if len(sys.argv) > 6 else -30.0
nu_min = float(sys.argv[7]) if len(sys.argv) > 7 else 0.2
t0 = time.time()
Y = np.load(f)
if ss != -20.0:
    B20 = BQ(lam, hs=hs, Nb=32); Bn = BQ(lam, hs=hs, Nb=32, s_start=ss)
    Y = full(B20, Y)[Bn.i0:]
S = BQStab(lam, Y, hs=hs, s_start=ss)
d = lam - 1
last = {'v': None}


def count(mu):
    try:
        v, V = spectrum_c(S, complex(mu), k=16, v0=last['v'])
    except Exception:
        v, V = spectrum_c(S, complex(mu), k=16)
    last['v'] = V[:, 0]
    return int(np.sum(v[np.abs(v.imag) < 1e-8].real > 1)), int(np.sum(np.abs(v[np.abs(v.imag) >= 1e-8]) > 1)), float(np.abs(v).min())


grid = np.arange(nu_max, nu_min - 1e-9, -step)
rows = []
for g in grid:
    rows.append((g,) + count(g * d))
    print(f"  ν̂={g:.3f}: (N_real, N_cplx, |ν_16|) = {rows[-1][1:]}", flush=True)
roots = []
for i in range(len(rows) - 1):
    if (rows[i + 1][1] - rows[i][1]) % 2 == 1:
        lo, hi = rows[i + 1][0], rows[i][0]; plo = rows[i + 1][1] % 2
        for _ in range(9):
            mid = 0.5 * (lo + hi)
            if count(mid * d)[0] % 2 == plo:
                lo = mid
            else:
                hi = mid
        roots.append(0.5 * (lo + hi))
roots = sorted(roots)
print(f"2D λ={lam:.7f} z={1/d:.3f} m={S.m:.9f} hs={hs} s_start={ss:g} (artifact ≈ ν̂ {0.6/abs(ss)/d:.2f}): "
      f"ν̂ = μ/(λ−1) = {np.round(roots, 4).tolist()}; μ = {np.round(np.array(roots) * d, 6).tolist()}; "
      f"max |ν_16| = {max(r[3] for r in rows):.3f} ({time.time() - t0:.0f}s)", flush=True)
np.save(f"bqflow_lam{lam:.6f}_hs{hs}.npy", np.array(rows))
