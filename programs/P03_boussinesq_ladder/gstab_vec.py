"""Eigenvector localisation of a global-discretisation eigenvalue near a shift σ (spurious boundary modes concentrate
at s_min or s_max; genuine modes live on the profile)."""
import numpy as np, scipy.sparse as sp, sys
from scipy.sparse.linalg import splu, LinearOperator, eigs
from bq_global import BQGlobal
fU = sys.argv[1]; smin, smax, hs, Nb = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
sig = complex(sys.argv[6].replace('i', 'j'))
G = BQGlobal(s_min=smin, s_max=smax, hs=hs, Nb=Nb); U = np.load(fU); N = G.N
J = G.jacobian(U)[:3 * N, :3 * N].tocsc()
mdiag = np.concatenate([np.where(G.inflow, 0.0, 1.0), np.where(G.inflow, 0.0, 1.0), np.zeros(N)])
M = sp.diags(mdiag).tocsc()
lu = splu((-J - sig * M).tocsc().astype(complex), permc_spec='COLAMD')
op = LinearOperator((3 * N, 3 * N), matvec=lambda x: lu.solve(M @ x), dtype=complex)
th, V = eigs(op, k=3, which='LM', tol=1e-10)
for i in range(3):
    mu = sig + 1 / th[i]; v = V[:, i]
    T = np.abs(v[:N]).reshape(G.Ns, G.n).max(1); O = np.abs(v[N:2 * N]).reshape(G.Ns, G.n).max(1)
    X = np.abs(v[2 * N:]).reshape(G.Ns, G.n).max(1)
    def where(A):
        j = np.argmax(A); return f"max at s={G.s[j]:.2f}; value at s_min {A[0]/A.max():.1e}, at s=0 {A[np.argmin(abs(G.s))]/A.max():.1e}, at s_max {A[-1]/A.max():.1e}"
    print(f"μ = {mu:.5f}\n  θ': {where(T)}\n  ω': {where(O)}\n  X': {where(X)}")
