import numpy as np, sys, time
from ccf_nk import CCFNK
from ccf_mapped import CCFMapped
from stability_mapped import MappedStability
S0 = CCFNK(30,120,16384,0.7)
name = sys.argv[1]
if name == 'l0':
    d = np.load('ladder_F16k_up_cross2.npy'); lam, eta, phi, etad, eps = d[0], S0.eta, d[1:], -1.0, 1.0
elif name == 'l1':
    d = np.load('ladder_F16k_up_cross1.npy'); lam, eta, phi, etad, eps = d[0], S0.eta, d[1:], -1.45, 0.3
else:
    d = np.load('mapped_l2_hs0.02_eps0.025.npy'); lam, eta, phi, etad, eps = 0.4713242277712, d[0], d[1], -0.9977, 0.025
M = CCFMapped(30,120,hs=0.03,eta_d=etad,eps=eps,sigma=2.0,c=0.7)
ph, lam, ok = M.solve_p(M.from_other(eta, phi), lam, 2.0)
st = MappedStability(M, ph, lam)
print(f"{name}: lam={lam:.12f} ok={ok} N={M.N}", flush=True)
for mu in np.round(np.arange(0.05, 2.01, 0.05), 3):
    v = st.spectrum(mu, k=6)
    vr = np.sort(v[np.abs(v.imag) < 1e-6*np.maximum(1,abs(v))].real)[::-1]
    vc = v[np.abs(v.imag) >= 1e-6*np.maximum(1,abs(v))]
    print(f"  mu={mu:.2f}: real {np.round(vr,4)}  complex {np.round(vc[vc.imag>0],3)}", flush=True)
