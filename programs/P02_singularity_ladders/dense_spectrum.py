"""Spectrum of the method-of-lines generator of the linearised self-similar CCF dynamics (same operator as the
linear evolution: 3rd-order upwind transport on the mapped grid, inflow δ=0, dense alternating-point Hilbert)."""
import numpy as np, sys, time
from block_evolution import setup
which = sys.argv[1]
M, phi, lam = setup(which)
Th = np.exp(phi + M.c*M.eta); d = 1 + lam + M.G(phi); Th_eta = lam*Th/d
pre, post = np.exp(-M.c*M.eta), np.exp((M.c-1)*M.eta)
v = d/(M.gp*M.hs); N = M.N
Ds = np.zeros((N, N))
for i in range(2, N-1):
    Ds[i, i+1] += 2/6; Ds[i, i] += 3/6; Ds[i, i-1] += -6/6; Ds[i, i-2] += 1/6
Ds[1,1], Ds[1,0] = 1, -1; Ds[N-1,N-1], Ds[N-1,N-2] = 1, -1
A = lam*np.eye(N) - v[:,None]*Ds - (post*Th_eta)[:,None]*(M.Hm*pre[None,:])
A[0,:] = 0.0; A[0,0] = -50.0          # inflow boundary: strongly damped (δ=0)
t0 = time.time()
ev = np.linalg.eigvals(A)
ev = ev[np.argsort(-ev.real)]
print(f"{which}: lam={lam:.10f} N={N}  top eigenvalues (Re desc): " + ", ".join(f"{z.real:.5f}{z.imag:+.5f}i" for z in ev[:10]) + f"  ({time.time()-t0:.0f}s)", flush=True)
