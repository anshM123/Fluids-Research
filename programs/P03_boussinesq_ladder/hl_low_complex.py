"""Complex eigenvalues of the Hou–Luo linearisation near the bottom of the unstable ladder (small μ/(λ−1)) at branch
points: secant iteration on ν(μ) = 1 from a grid of complex starts, as in hl_find_complex.py, on the extended grid
of hl_deep_spec.py.
usage: python3 hl_low_complex.py ETA0 STATE NUHAT_LO NUHAT_HI   (env HLGRID, HLK as in hl_deep_spec.py)"""
import numpy as np, sys, re, os
from hl_stability import HLStab
eta0, f = float(sys.argv[1]), sys.argv[2]; lo, hi = float(sys.argv[3]), float(sys.argv[4])
N0, L1_0, L2_0 = [t(x) for t, x in zip((int, float, float), os.environ.get('HLGRID', '65536,25,75').split(','))]
KEIG = int(os.environ.get('HLK', '20'))
h = (L1_0 + L2_0) / N0
K = max(0, int(np.ceil((abs(eta0) + 2 - L1_0) / h / 4096.0)) * 4096)
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1)); d = lam - 1
q0 = np.load(f); St = HLStab(lam, np.concatenate([np.full(K, q0[0]), q0]), N0 + K, L1=L1_0 + K * h, L2=L2_0, eta_start=eta0)
def g(mu):
    v, _ = St.spectrum(mu, k=KEIG)
    return v[np.argmin(np.abs(v - 1))] - 1
found = []
for re0 in np.arange(lo, hi + 1e-9, 0.25):
    for im0 in (0.1, 0.3, 0.6):
        m0 = complex(re0, im0) * d; m1 = m0 + 1e-3 * d
        g0, g1 = g(m0), g(m1)
        for it in range(40):
            if abs(g1 - g0) < 1e-300:
                break
            m2 = m1 - g1 * (m1 - m0) / (g1 - g0)
            m0, g0 = m1, g1; m1 = m2; g1 = g(m1)
            if abs(g1) < 1e-10:
                break
        if abs(g1) < 1e-8 and 0 < m1.real and all(abs(m1 - u) > 1e-5 * d for u in found):
            found.append(m1)
            print(f"  start {re0:.2f}+{im0:.1f}i → μ/(λ−1) = {m1.real/d:.4f} {m1.imag/d:+.4f}i  (|ν−1| = {abs(g1):.1e})", flush=True)
print(f"HL λ={lam:.6f} z={1/d:.3f} η0={eta0:g}: eigenvalues found in the window (units of λ−1): "
      f"{[f'{u.real/d:.4f}{u.imag/d:+.4f}i' for u in found]}", flush=True)
