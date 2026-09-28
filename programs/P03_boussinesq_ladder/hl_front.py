import numpy as np, sys, re
from hl_solver import HL
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    S = HL(lam); q = np.load(f); r = S.march(S.full(q))
    D, eps = r['D'], r['eps']
    eta = S.eta
    Dh = D / eps
    # front: where D̂ crosses 2 and 10 (first time after the dip)
    k2 = np.argmax((Dh > 2) & (eta > -8)); k10 = np.argmax((Dh > 10) & (eta > -8))
    kd = np.argmin(Dh[(eta > -8) & (eta < eta[k2])]) + np.argmax(eta > -8)
    grad = np.gradient(np.log(D), eta)
    print(f"λ={lam:.5f} z={1/(lam-1):6.2f} ε={eps:.4f}: dip D̂={Dh[kd]:.4f} at η={eta[kd]:+.3f}; D̂=2 at η={eta[k2]:+.4f}, D̂=10 at η={eta[k10]:+.4f} (width {eta[k10]-eta[k2]:.4f}); max d lnD/dη={grad[(eta>-8)&(eta<2)].max():.2f}")
