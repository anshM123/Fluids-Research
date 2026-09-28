import numpy as np, time, sys
from ccf_nk import CCFNK
from ccf_mapped import CCFMapped
S0 = CCFNK(30,120,65536,0.7)
d = np.load('ladder_F64k_dn_cross1.npy'); lam = d[0]; phi = d[1:]
F,den = S0.F(phi, lam, S0.normval(phi,lam)); etad = S0.eta[np.argmin(den)]
print("F64k lam2 =", lam, " dip at", etad, flush=True)
prev_eta, prev_phi = S0.eta, phi
for hs, eps, sig in [(0.04,0.2,2.0),(0.03,0.1,2.0),(0.03,0.05,2.0),(0.02,0.05,2.0),(0.02,0.025,2.0)]:
    t0=time.time()
    M = CCFMapped(30,120,hs=hs,eta_d=etad,eps=eps,sigma=sig,c=0.7)
    ph, l, ok = M.solve_p(M.from_other(prev_eta, prev_phi), lam, 2.0)
    Fm, denm = M.F(ph, l, M.normval(ph,l)); j=np.argmin(denm)
    print(f"hs={hs} eps={eps} N={M.N} h_min={hs*eps:.1e}: lam2={l:.13f} ok={ok} min_den={denm[j]:.6f} at {M.eta[j]:.4f} ({time.time()-t0:.0f}s)", flush=True)
    if ok:
        prev_eta, prev_phi, lam = M.eta, ph, l
        np.save(f"mapped_l2_hs{hs}_eps{eps}.npy", np.vstack([M.eta, ph]))
