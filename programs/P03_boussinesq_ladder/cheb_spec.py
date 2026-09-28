import numpy as np, sys, re
from numpy.polynomial import chebyshev as C
from bq_solver import BQ
from bq_newton import full
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    mNb = re.search(r'Nb(\d+)_hs([0-9.]+?)_', f); Nb, hs = (int(mNb.group(1)), float(mNb.group(2))) if mNb else (32, 0.025)
    B = BQ(lam, Nb=Nb, hs=hs); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    x = np.cos(np.pi * np.arange(Nb + 1) / Nb)
    print(f"{f}  (Nb={Nb})")
    for sv in (-4, -1, 0, 1, 3, 6):
        i = int(np.argmin(np.abs(B.s - sv)))
        for name, F in (("Θ̂", r['Th'][i]), ("Ω̂", r['Om'][i]), ("X", X[i])):
            c = np.abs(C.chebfit(x, F, Nb)); c /= c.max()
            print(f"   s={sv:+d} {name}: |c_k|/max at k=8,16,24,{Nb}: {c[8]:.1e} {c[16]:.1e} {c[min(24,Nb)]:.1e} {c[Nb]:.1e}")
