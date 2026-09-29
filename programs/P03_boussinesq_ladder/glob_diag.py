import numpy as np, sys
from bq_global import BQGlobal
f, lam0 = sys.argv[1], float(sys.argv[2])
G = BQGlobal(s_min=-8, s_max=30, hs=0.05, Nb=16)
U = G.initial_from_march(f, lam0)
R = G.residual(U); N = G.N
for name, Rp in (("Θ", R[:N]), ("Ω", R[N:2*N]), ("X", R[2*N:3*N])):
    A = np.abs(Rp).reshape(G.Ns, G.n)
    j, k = np.unravel_index(np.argmax(A), A.shape)
    print(f"{name}: max {A.max():.3e} at s={G.s[j]:.2f} β-index {k}; by s-range:", [f"{A[(G.s>=a)&(G.s<b)].max():.1e}" for a, b in ((-8,-7.5),(-7.5,-4),(-4,0),(0,4),(4,10),(10,20),(20,31))])
print("λ-eq:", R[-1])
