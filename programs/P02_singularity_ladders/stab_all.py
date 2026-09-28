import numpy as np, sys, time
from ccf_nk import CCFNK
from ccf_mapped import CCFMapped
from stability_mapped import MappedStability
S0 = CCFNK(30,120,16384,0.7)
cases = []
d = np.load('ladder_F16k_up_cross2.npy'); cases.append(('lambda0', d[0], S0.eta, d[1:], -1.0, 1.0))
d = np.load('ladder_F16k_up_cross1.npy'); cases.append(('lambda1', d[0], S0.eta, d[1:], -1.45, 0.3))
d = np.load('mapped_l2_hs0.02_eps0.025.npy'); cases.append(('lambda2', 0.4713242277712, d[0], d[1], -0.9977, 0.025))
which = sys.argv[1]
for name, lam, eta, phi, etad, eps in cases:
    if name != which: continue
    M = CCFMapped(30,120,hs=0.03,eta_d=etad,eps=eps,sigma=2.0,c=0.7)
    ph, lam, ok = M.solve_p(M.from_other(eta, phi), lam, 2.0)
    st = MappedStability(M, ph, lam)
    print(f"{name}: lam={lam:.12f} ok={ok} N={M.N}", flush=True)
    for mu, vr in st.crossings(np.round(np.arange(0.05, 1.61, 0.05), 3)):
        print(f"  mu={mu:.2f}: {np.round(vr[:6],4)}", flush=True)
