import numpy as np, time, sys
from ccf_nk import CCFNK
from ccf_stability import CCFStability
N = int(sys.argv[1]); f = sys.argv[2]
S = CCFNK(30, 120, N, 0.7)
d = np.load(f); lam = d[0]; phi = d[1:]
if len(phi) != N:  # interpolate onto this grid
    S0 = CCFNK(30, 120, len(phi), 0.7); phi = np.interp(S.eta, S0.eta, phi)
    phi, lam, ok = S.solve_p(phi, lam, 2.0); print("re-solved on N", N, ok, lam)
st = CCFStability(S, phi, lam)
def nu_near1(mu, k=8):
    vals, _ = st.spectrum(mu, k=k)
    return vals
# scan mu on real axis and track eigenvalues of T_mu near 1
mus = np.arange(0.15, 2.01, 0.05)
prev = None
for mu in mus:
    vals = nu_near1(mu)
    vr = np.sort(vals[np.abs(vals.imag) < 1e-6].real)[::-1]
    print(f"mu={mu:.2f}: real eigenvalues of T_mu: {np.round(vr, 5)}", flush=True)
