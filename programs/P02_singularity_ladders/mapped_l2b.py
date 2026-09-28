import numpy as np, time
from ccf_nk import CCFNK
from ccf_mapped import CCFMapped
S0 = CCFNK(30,120,16384,0.7)
for f in ['ladder_N16k_cross1.npy', 'ladder_N32k_cross1.npy']:
    d = np.load(f); lam = d[0]; phi = d[1:]
    if len(phi) != S0.N:
        S1 = CCFNK(30,120,len(phi),0.7); phi = np.interp(S0.eta, S1.eta, phi)
    # this file used old normalisation; re-normalise implicitly by solving with new normval target
    F,den = S0.F(phi, lam, S0.normval(phi,lam)); etad = S0.eta[np.argmin(den)]
    print(f"{f}: start lam={lam:.10f} dip at {etad:.3f} min_den={den.min():.4f}", flush=True)
    prev_eta, prev_phi = S0.eta, phi
    for hs, eps in [(0.03,0.1),(0.03,0.05),(0.02,0.025)]:
        M = CCFMapped(30,120,hs=hs,eta_d=etad,eps=eps,sigma=2.0,c=0.7)
        ph, l, ok = M.solve_p(M.from_other(prev_eta, prev_phi), lam, 2.0)
        Fm, denm = M.F(ph, l, M.normval(ph,l))
        print(f"   hs={hs} eps={eps} N={M.N}: lam={l:.13f} ok={ok} min_den={denm.min():.6f}", flush=True)
        if ok: prev_eta, prev_phi, lam = M.eta, ph, l
