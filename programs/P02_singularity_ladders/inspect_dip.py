import numpy as np, time
from ccf_newton import CCFNewton
S = CCFNewton(30, 120, 4096, 0.68)
phi = S.guess(0.6, 2.0); phi, ok = S.solve_fixed(phi, 0.6)
for lam in [0.58, 0.56, 0.54, 0.53, 0.52]:
    phi, ok = S.solve_fixed(phi, lam)
    F, den = S.F(phi, lam, phi[S.i0])
    j = np.argmin(den)
    # curvature of den at min
    kap = (den[j+1]-2*den[j]+den[j-1])/S.h**2/2
    print(f"lam={lam} ok={ok} p={S.p_of(phi,lam):.8f} min den={den[j]:.5f} at eta={S.eta[j]:.3f}, kappa={kap:.3f}, width~{np.sqrt(den[j]/max(kap,1e-9)):.4f}", flush=True)
np.save("sol_lam052_N4096.npy", np.concatenate([[0.52], phi]))
xi = np.exp(S.eta); Th = np.exp(phi + S.c*S.eta)
m = (S.eta > -4) & (S.eta < 3)
for e in np.arange(-4, 3, 0.25):
    i = np.argmin(np.abs(S.eta-e)); print(f"  eta={S.eta[i]:+.2f} xi={xi[i]:.4f} Theta={Th[i]:.4e} den={den[i]:.4f} b={0.52/den[i]-S.c:+.3f}")
