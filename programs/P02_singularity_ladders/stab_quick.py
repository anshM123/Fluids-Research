import numpy as np, time, sys
from ccf_mapped import CCFMapped
from stability_mapped import MappedStability
from scipy.sparse.linalg import LinearOperator, eigs
d = np.load('mapped_l2_hs0.02_eps0.025.npy'); lam = 0.4713242277712
M = CCFMapped(30,120,hs=0.03,eta_d=-0.9977,eps=0.025,sigma=2.0,c=0.7)
ph, lam, ok = M.solve_p(M.from_other(d[0], d[1]), lam, 2.0)
st = MappedStability(M, ph, lam)
v = np.random.rand(M.N).astype(complex)
t=time.time(); 
for _ in range(10): st.T(v, 0.5)
print("T apply time", (time.time()-t)/10, flush=True)
cnt=[0]
def mv(x):
    cnt[0]+=1; return st.T(x, mu)
for mu in [0.3, 1.0]:
    cnt[0]=0; t=time.time()
    op = LinearOperator((M.N, M.N), matvec=mv, dtype=complex)
    vals = eigs(op, k=6, which='LM', tol=1e-8, ncv=40, maxiter=500, return_eigenvectors=False)
    print(f"mu={mu}: vals={np.round(vals,5)}  matvecs={cnt[0]} time={time.time()-t:.1f}s", flush=True)
