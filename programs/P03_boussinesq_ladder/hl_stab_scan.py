"""Hou–Luo stability of the smooth profiles: fine real-μ scan (count of real ν > 1 of T_μ), Brent refinement of each
crossing, i.e. of each real unstable eigenvalue, and the ratios μ_k/(λ−1).
usage: python3 hl_stab_scan.py N ETA_START FILES..."""
import numpy as np, sys, re, time
from scipy.optimize import brentq
from hl_stability import HLStab
N, eta_start = int(sys.argv[1]), float(sys.argv[2])
for f in sys.argv[3:]:
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    St = HLStab(lam, np.load(f), N, eta_start=eta_start)
    def spec(mu):
        vals, _ = St.spectrum(mu, k=12)
        return vals
    def g_near(mu):
        v = spec(mu); re_ = v[np.abs(v.imag) < 1e-8].real
        return re_[np.argmin(np.abs(re_ - 1))] - 1
    mus = np.concatenate([np.arange(1.3, 0.2, -0.01), np.arange(0.2, 0.004, -0.002)])
    rows = []
    for mu in mus:
        v = spec(mu)
        rows.append((mu, int(np.sum(v[np.abs(v.imag) < 1e-8].real > 1)), int(np.sum(np.abs(v[np.abs(v.imag) >= 1e-8]) > 1))))
    rows = np.array(rows)
    def count(mu):
        v = spec(mu); return int(np.sum(v[np.abs(v.imag) < 1e-8].real > 1))
    # det(I − T_μ) is real on the real axis and its sign is (−1)^{#real ν > 1}: each eigenvalue μ (a real ν crossing 1)
    # flips the parity of the count; collisions of two real ν > 1 into a complex pair do not.
    roots = []
    for i in range(len(rows) - 1):
        if (rows[i + 1, 1] - rows[i, 1]) % 2 == 1:
            lo, hi = rows[i + 1, 0], rows[i, 0]; plo = rows[i + 1, 1] % 2
            for _ in range(14):
                mid = 0.5 * (lo + hi)
                if count(mid) % 2 == plo:
                    lo = mid
                else:
                    hi = mid
            roots.append(0.5 * (lo + hi))
    nontriv = sorted([r for r in roots if abs(r - 1) > 1e-4], reverse=True)
    low = rows[rows[:, 0] < 0.03]
    print(f"HL λ={lam:.8f} z={1/(lam-1):.4f} m={St.m:.10f} N={N} η0={eta_start}: {len(nontriv)} real crossings "
          f"μ={np.round(nontriv, 6).tolist()}; μ/(λ−1)={np.round(np.array(sorted(nontriv))/(lam-1), 4).tolist()}; "
          f"low-μ (N_real, N_cplx) at μ<0.03: {low[:, 1:].astype(int).tolist()[:4]}", flush=True)
    np.save(f"hlstab_scan_lam{lam:.8f}_N{N}_e{int(-eta_start)}.npy", rows)
