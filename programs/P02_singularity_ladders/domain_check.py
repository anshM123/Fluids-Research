import numpy as np, time
from ccf_nk import CCFNK
from ccf_mapped import CCFMapped
S0 = CCFNK(30,120,16384,0.7)
d = np.load('ladder_F16k_up_cross2.npy'); lam0 = d[0]; phi0 = d[1:]
d2 = np.load('mapped_l2_hs0.02_eps0.025.npy'); eta2, phi2 = d2[0], d2[1]; lam2 = 0.4713242277712
for (L1, L2, c, hs) in [(30,120,0.7,0.03),(45,120,0.7,0.03),(30,170,0.7,0.03),(30,120,0.62,0.03),(30,120,0.78,0.03),(45,170,0.66,0.02)]:
    t0=time.time()
    M = CCFMapped(L1,L2,hs=hs,eta_d=-1.0,eps=1.0,sigma=2.0,c=c)
    ph = M.from_other(S0.eta, phi0 + (0.7 - c)*S0.eta)   # Psi weight change: phi = ln Theta - c eta
    ph, l0, ok0 = M.solve_p(ph, lam0, 2.0)
    M2 = CCFMapped(L1,L2,hs=0.02,eta_d=-0.9977,eps=0.025,sigma=2.0,c=c)
    ph2 = M2.from_other(eta2, phi2 + (0.7 - c)*eta2)
    ph2, l2, ok2 = M2.solve_p(ph2, lam2, 2.0)
    print(f"L1={L1} L2={L2} c={c}: lam0={l0:.13f} (beta0={l0/(1+l0):.10f}) ok={ok0} | lam2={l2:.12f} ok={ok2}  ({time.time()-t0:.0f}s)", flush=True)
