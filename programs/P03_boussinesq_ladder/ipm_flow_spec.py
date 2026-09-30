"""2D IPM: lower unstable spectrum of a profile on the least-singular branch (a rung or a branch point between
rungs), by parity flips of the count of real ν > 1 of T_μ on a grid in ν̂ = μ/λ, refined by bisection on the
parity (as bq_flow_spec.py for Boussinesq; the transported-scalar exponent is λ here, so ν̂ = μ/λ).
The state is extended to the origin truncation S_START by its exact constant-strain structure.
usage: python3 ipm_flow_spec.py STATE LAM NUHAT_MAX STEP [S_START=-30] [NUHAT_MIN=0.1] [K=16] [HS=0.025]"""
import numpy as np, sys, time
from ipm_solver import IPM
from bq_newton import full
from ipm_stability import IPMStab, spectrum_c

f, lam = sys.argv[1], float(sys.argv[2])
nu_max, step = float(sys.argv[3]), float(sys.argv[4])
ss = float(sys.argv[5]) if len(sys.argv) > 5 else -30.0
nu_min = float(sys.argv[6]) if len(sys.argv) > 6 else 0.1
K = int(sys.argv[7]) if len(sys.argv) > 7 else 16
hs = float(sys.argv[8]) if len(sys.argv) > 8 else 0.025
t0 = time.time()
Y = np.load(f)
if ss != -20.0:
    B20 = IPM(lam, hs=hs, Nb=32); Bn = IPM(lam, hs=hs, Nb=32, s_start=ss)
    Y = full(B20, Y)[Bn.i0:]
S = IPMStab(lam, Y, hs=hs, s_start=ss)
d = lam
last = {'v': None}


def count(mu):
    try:
        v, V = spectrum_c(S, complex(mu), k=K, v0=last['v'])
    except Exception:
        v, V = spectrum_c(S, complex(mu), k=K)
    last['v'] = V[:, 0]
    par = 0 if np.real(np.prod(1 - v)) > 0 else 1           # parity of the real ν > 1 from the sign of Re Π(1 − ν)
    nre = int(np.sum(v[np.abs(v.imag) < 1e-8].real > 1))
    nre = 2 * (nre // 2) + par if (nre % 2) != par else nre
    return nre, int(np.sum(np.abs(v[np.abs(v.imag) >= 1e-8]) > 1)), float(np.abs(v).min())


grid = np.arange(nu_max, nu_min - 1e-9, -step)
rows = []
for g in grid:
    rows.append((g,) + count(g * d))
    print(f"  ν̂={g:.3f}: (N_real, N_cplx, |ν_K|) = {rows[-1][1:]}", flush=True)
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
print(f"IPM λ={lam:.7f} z={1/d:.4f} m={S.m:.9f} hs={hs} s_start={ss:g} (artifact ≈ ν̂ {0.6/abs(ss)/d:.2f}): "
      f"ν̂ = μ/λ = {np.round(roots, 4).tolist()}; μ = {np.round(np.array(roots) * d, 6).tolist()}; "
      f"max |ν_K| = {max(r[3] for r in rows):.3f} ({time.time() - t0:.0f}s)", flush=True)
np.save(f"ipmflow_lam{lam:.6f}_hs{hs}.npy", np.array(rows))
