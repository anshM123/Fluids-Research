"""Try the Hou–Luo WKB phase formula on 2D Boussinesq boundary data (sign of Ω flipped to the HL convention)."""
import numpy as np, sys, re
from bq_solver import BQ
from bq_newton import full
for f in sys.argv[1:]:
    lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
    B = BQ(lam); Y = np.load(f); X = full(B, Y)
    r = B.march(X / B.ea2[:, None], return_all=True)
    eps = 1 + lam - r['A']
    D = (1 + lam) + r['Ur'][:, 0]
    Om = -r['Omega'][:, 0]                         # 2D: Ω < 0 on y1 > 0 near the boundary; HL convention Ω > 0
    k = (-1 + np.sqrt(1 + 4j * Om)) / (2j * D)
    sel = slice(B.i0, None)
    Phi = np.trapezoid(k[sel], B.s[sel])
    print(f"λ={lam:.5f} z={1/(lam-1):.3f} ε={eps:.4f} m={r['m']:.6f}: εΦ = {eps*Phi.real:.5f}{eps*Phi.imag:+.5f}i  Δz_pred={np.pi/(r['m']*eps*Phi.real):.4f}")
