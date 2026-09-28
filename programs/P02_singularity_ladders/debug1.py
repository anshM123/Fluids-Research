import numpy as np, time
from ccf_profile import CCFProfile
prof = CCFProfile(L1=30, L2=150, N=1024, c=0.75)
lam, p = 1.0, 2
Psi = prof.initial_guess(lam, p)
# finite-difference Jacobian check on a few columns
R0,G,W,D = prof.residual(Psi, lam, 0.0)
print("initial residual max", np.max(np.abs(R0)), " h1=", prof.h1(Psi), " G(left)=", G[3], "p_loc=", lam/(1+lam+prof.h1(Psi)))
# Jacobian (fixed lam) via same formulas
N=prof.N; c=prof.c
Dm = prof.circulant(prof.mD); Km = prof.circulant(prof.mK)
A = (-lam + (1+lam)*c)*np.eye(N) + (1+lam)*Dm + G[:,None]*(c*np.eye(N)+Dm) + (W*prof.E)[:,None]*Km
for j in [100, 500, 800]:
    e = np.zeros(N); e[j]=1e-6
    R1,_,_,_ = prof.residual(Psi+e, lam, 0.0)
    fd = (R1-R0)/1e-6
    print("col",j,"FD vs analytic max diff", np.max(np.abs(fd-A[:,j])), "scale", np.max(np.abs(A[:,j])))
# min(1+lam+G)
print("min(1+lam+G) =", np.min(1+lam+G))
# eigen/sing values of A
sv = np.linalg.svd(A, compute_uv=False)
print("smallest singular values of A:", sv[-5:], " largest", sv[0])
