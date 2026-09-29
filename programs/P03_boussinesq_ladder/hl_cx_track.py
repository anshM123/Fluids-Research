"""Robustness of a complex Hou–Luo eigenvalue: secant on ν(μ) = 1 from a given start, for several origin truncations.
usage: python3 hl_cx_track.py STATE RE IM ETA0,ETA0,...   (μ/(λ−1) = RE + i IM; env HLGRID="N0,L1,L2", default F grid)"""
import numpy as np, sys, re, os
from hl_stability import HLStab
f = sys.argv[1]; st = complex(float(sys.argv[2]), float(sys.argv[3])); etas = [float(x) for x in sys.argv[4].split(',')]
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1)); d = lam - 1
N0, L1_0, L2_0 = [t(x) for t, x in zip((int, float, float), os.environ.get('HLGRID', '65536,25,75').split(','))]; h = (L1_0 + L2_0) / N0; q0 = np.load(f)
for eta0 in etas:
    K = max(0, int(np.ceil((abs(eta0) + 2 - L1_0) / h / 4096.0)) * 4096)
    St = HLStab(lam, np.concatenate([np.full(K, q0[0]), q0]), N0 + K, L1=L1_0 + K * h, L2=L2_0, eta_start=eta0)
    def g(mu):
        v, _ = St.spectrum(mu, k=24)
        return v[np.argmin(np.abs(v - 1))] - 1
    m0 = st * d; m1 = m0 * (1 + 1e-3); g0, g1 = g(m0), g(m1)
    for it in range(30):
        m2 = m1 - g1 * (m1 - m0) / (g1 - g0); m0, g0 = m1, g1; m1 = m2; g1 = g(m1)
        if abs(g1) < 1e-10:
            break
    print(f"z={1/d:.3f} η0={eta0:g}: μ/(λ−1) = {m1.real/d:.5f} {m1.imag/d:+.5f}i  (|ν−1| = {abs(g1):.1e}, {it+1} it)", flush=True)
