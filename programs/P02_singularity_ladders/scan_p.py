import numpy as np, sys
from ccf_picard import setup, G_of
from scipy.integrate import cumulative_trapezoid
def solve_from(lnP, lam, eta, h, mK, c, iters=3000, relax=0.3, tol=1e-12):
    i0 = np.argmin(np.abs(eta))
    for it in range(iters):
        G = G_of(np.exp(lnP), eta, mK, c)
        den = 1+lam+G
        if np.min(den) <= 0: return None, None
        phi = cumulative_trapezoid(lam/den - c, eta, initial=0.0)
        phi = phi - phi[i0] + lnP[i0]
        diff = np.max(np.abs(phi-lnP)[(eta>-20)&(eta<60)])
        lnP = (1-relax)*lnP + relax*phi
        if diff < tol: break
    h1 = -(2/np.pi)*h*np.sum(np.exp(lnP)*np.exp((c-1)*eta))
    return lnP, h1
L1, L2, N, c = 30, 150, 4096, 0.75
eta, h, mK = setup(L1, L2, N, c)
lam0 = 0.6
beta = lam0/(1+lam0); p0 = 2.0
Theta = np.exp(p0*eta - 0.5*(p0-beta)*np.logaddexp(0, 2*eta)); Psi = np.exp(-c*eta)*Theta
h1 = -(2/np.pi)*h*np.sum(Psi*np.exp((c-1)*eta)); Psi *= (lam0/p0-1-lam0)/h1
lnP = np.log(Psi)
res = []
for direction, lams in [(+1, np.arange(0.60, 0.99, 0.02)), (-1, np.arange(0.58, 0.19, -0.02))]:
    cur = lnP.copy() if direction > 0 else start_down.copy()
    for lam in lams:
        out, h1 = solve_from(cur, lam, eta, h, mK, c)
        if out is None:
            print(f"lam={lam:.3f}: sonic failure"); break
        cur = out
        p = lam/(1+lam+h1)
        res.append((lam, h1, p)); print(f"lam={lam:.3f}  h1={h1:.10f}  p={p:.8f}", flush=True)
        if direction > 0 and abs(lam-0.60) < 1e-9: start_down = out.copy()
np.save("scan_p_lambda.npy", np.array(sorted(res)))
