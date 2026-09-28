import numpy as np
from bq_local_eig import local_root, hl_root, shoot
for x in (0.05, 0.1, 0.2, 0.3, 0.5):
    Dh, ch, mu, G, Thy = 1.0, 2 * x, 4.0, -2 * x, 0.0
    g = hl_root(Dh, ch)
    k, ok = local_root(Dh, ch, mu, G, Thy, g)
    k2, ok2 = local_root(Dh, ch, mu, G, Thy, g * 0.8)
    print(f"x={x}: HL {g:.5f} -> 2D local root {k:.6f} ({ok}); from other guess {k2:.6f}; check |F|={abs(shoot(k, Dh, ch, mu, G, Thy)):.1e}")
