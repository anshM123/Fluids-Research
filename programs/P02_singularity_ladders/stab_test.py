import numpy as np, time
from ccf_nk import CCFNK
from ccf_stability import CCFStability
S = CCFNK(30, 120, 16384, 0.7)
d = np.load('ladder_F16k_up_cross1.npy'); lam = d[0]; phi = d[1:]
st = CCFStability(S, phi, lam)
# check trivial time-translation mode: delta_T = lam*Theta*G/d  should satisfy T_1 delta = delta
G = st.d - 1 - lam
dT = lam * st.Theta * G / st.d
r = st.T(dT.astype(complex), 1.0)
m = np.abs(S.eta) < 20
print("lambda =", lam, " trivial-mode check |T1 dT - dT|/|dT| (|eta|<20):", np.max(np.abs(r - dT)[m]) / np.max(np.abs(dT)[m]))
for mu in [0.05, 0.2, 0.5, 1.0, 1.5]:
    t0 = time.time()
    vals, vecs = st.spectrum(mu, k=4)
    print(f"mu={mu}: top eigenvalues of T_mu: {np.round(vals, 5)}  ({time.time()-t0:.1f}s)", flush=True)
