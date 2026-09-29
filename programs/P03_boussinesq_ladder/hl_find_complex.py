"""Locate complex eigenvalues μ of the Hou–Luo linearisation: Newton/secant on ν(μ) = 1 in the complex plane, started
on a grid in the region where complex pairs of ν with |ν| > 1 occur on the real axis."""
import numpy as np, sys, re
from hl_stability import HLStab
f, N = sys.argv[1], int(sys.argv[2])
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
St = HLStab(lam, np.load(f), N)
def g(mu):
    v, _ = St.spectrum(mu, k=16)
    return v[np.argmin(np.abs(v - 1))] - 1
found = []
for re0 in np.arange(0.40, 0.72, 0.04):
    for im0 in (0.02, 0.06, 0.12):
        m0 = complex(re0, im0); m1 = m0 + 1e-3
        g0, g1 = g(m0), g(m1)
        for it in range(30):
            if abs(g1 - g0) < 1e-300:
                break
            m2 = m1 - g1 * (m1 - m0) / (g1 - g0)
            m0, g0 = m1, g1; m1 = m2; g1 = g(m1)
            if abs(g1) < 1e-10:
                break
        if abs(g1) < 1e-8 and m1.imag > 1e-6 and 0.03 < m1.real < 1.5 and all(abs(m1 - u) > 1e-5 for u in found):
            found.append(m1)
            print(f"complex eigenvalue μ = {m1.real:.6f} ± {m1.imag:.6f}i  (|ν−1| = {abs(g1):.1e})", flush=True)
print(f"HL λ={lam:.8f}: complex eigenvalues found: {[f'{u.real:.6f}±{u.imag:.6f}i' for u in found]}")
