"""Hou–Luo unstable spectrum along the continuous (least-singular) branch down to λ − 1 ≈ 0.025, beyond the smooth
rungs.  For z = 1/(λ−1) > 10 the branch profiles are smooth to |m − 2| < 1e-5, so their spectra differ from those of
nearby smooth profiles by a negligible amount, and the limit λ → 1 of μ_k/(λ − 1) can be followed directly.

The states hl_q_F_* (N = 65536 on η ∈ [−25, 75]) are extended to the left with the constant strain of the exact local
structure, at unchanged grid spacing, so that the perturbation march can start at η₀ = ETA0 (default −100). This
moves the origin-truncation artifact (≈ 0.6/|η₀|) below the smallest genuine eigenvalue c(λ − 1).

Real eigenvalues are located as parity flips of the count of real ν > 1 of T_μ on a grid in ν̂ = μ/(λ−1), then refined
by bisection on the parity.
usage: python3 hl_deep_spec.py ETA0 NUHAT_MAX STEP FILES...   (env: HLGRID="N0,L1,L2" grid of the states, default F; HLK eigenvalues per point, default 16; HLNUMIN lower end, default 0.2)"""
import numpy as np, sys, re, time, os
from hl_stability import HLStab

eta0, nu_max, step = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
KEIG = int(os.environ.get('HLK', '16')); NU_MIN = float(os.environ.get('HLNUMIN', '0.2'))
N0, L1_0, L2_0 = [t(x) for t, x in zip((int, float, float), os.environ.get('HLGRID', '65536,25,75').split(','))]
h = (L1_0 + L2_0) / N0
K = max(0, int(np.ceil((abs(eta0) + 2 - L1_0) / h / 4096.0)) * 4096)  # extra points on the left
N, L1 = N0 + K, L1_0 + K * h
for f in sys.argv[4:]:
    t0 = time.time()
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    q0 = np.load(f)
    q = np.concatenate([np.full(K, q0[0]), q0])
    St = HLStab(lam, q, N, L1=L1, L2=L2_0, eta_start=eta0)
    d = lam - 1

    def count(mu):
        v, _ = St.spectrum(mu, k=KEIG)
        # parity from the sign of det(I − T_μ) ≈ Re Π(1 − ν) over the computed set (complete while |ν_k| < 1): robust
        # where two real ν > 1 collide into a nearly real complex pair, which the |Im ν| test can split between the
        # "real" and "complex" counts
        par = 0 if np.real(np.prod(1 - v)) > 0 else 1
        nre = int(np.sum(v[np.abs(v.imag) < 1e-8].real > 1))
        return (2 * (nre // 2) + par if (nre % 2) != par else nre, int(np.sum(np.abs(v[np.abs(v.imag) >= 1e-8]) > 1)),
                float(np.abs(v).min()))       # the last entry must stay < 1 for the parity count to be complete

    grid = np.arange(nu_max, NU_MIN - 1e-9, -step)
    rows = [(g,) + count(g * d) for g in grid]
    roots = []
    for i in range(len(rows) - 1):
        if (rows[i + 1][1] - rows[i][1]) % 2 == 1:
            lo, hi = rows[i + 1][0], rows[i][0]; plo = rows[i + 1][1] % 2
            for _ in range(12):
                mid = 0.5 * (lo + hi)
                if count(mid * d)[0] % 2 == plo:
                    lo = mid
                else:
                    hi = mid
            roots.append(0.5 * (lo + hi))
    roots = sorted(roots)
    print(f"HL λ={lam:.6f} z={1/d:.3f} m={St.m:.9f} η0={eta0:g} (artifact ≈ ν̂ {0.6/abs(eta0)/d:.2f}): "
          f"ν̂ = μ/(λ−1) = {np.round(roots, 4).tolist()}; μ = {np.round(np.array(roots) * d, 6).tolist()}; "
          f"(N_real, N_cplx) at ν̂ = {rows[-1][0]:.2f}: {rows[-1][1:3]}; max |ν_{KEIG}| = {max(r[3] for r in rows):.3f} "
          f"({time.time() - t0:.0f}s)", flush=True)
    np.save(f"hldeep_lam{lam:.6f}_e{int(abs(eta0))}.npy", np.array(rows))
