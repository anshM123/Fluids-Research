import numpy as np, time
from ccf_nk import CCFNK
from ccf_mapped import CCFMapped
S0 = CCFNK(30,120,16384,0.7)
d = np.load('ladder_F16k_up_cross1.npy'); lam1 = d[0]; phi1 = d[1:]
d = np.load('ladder_F16k_up_cross2.npy'); lam0 = d[0]; phi0_ = d[1:]
for hs in [0.06, 0.04, 0.03]:
    t0=time.time()
    M = CCFMapped(30,120,hs=hs,eta_d=-1.4,eps=1.0,sigma=3.0,c=0.7)
    ph, l, ok = M.solve_p(M.from_other(S0.eta, phi1), lam1, 2.0)
    ph0, l0, ok0 = M.solve_p(M.from_other(S0.eta, phi0_), lam0, 2.0)
    print(f"uniform-map hs={hs} N={M.N}: lam1={l:.13f} ok={ok}  lam0={l0:.13f} ok={ok0}  ({time.time()-t0:.1f}s)", flush=True)
print("FFT N=16k reference: lam1=%.13f lam0=%.13f" % (lam1, lam0))
