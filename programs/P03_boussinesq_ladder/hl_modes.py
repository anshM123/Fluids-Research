"""Eigenfunctions of the lowest unstable Hou–Luo modes: θ'(η), ω'(η), q'(η) from the T_μ eigenvector with ν = 1,
and the local exponent p(η) = d ln|θ'|/dη inside the stalled layer, compared with the transport value
p_k = m(1 − μ_k/(λ−1)) = 2 − μ_k/ε (homogeneous solution of the temperature transport with speed εD̂ ≈ ε)."""
import numpy as np, sys, re
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from hl_stability import HLStab, _lin_gl2
f, N = sys.argv[1], int(sys.argv[2]); mus = [float(x) for x in sys.argv[3].split(',')]
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
St = HLStab(lam, np.load(f), N); S = St.S
eps = St.eps; eta = S.eta
D = (1 + lam) + St.q; Dh = D / eps
sel = (eta > -3) & (eta < 2); kd = np.argmax(sel) + np.argmin(Dh[sel]); kc = kd + np.argmax(Dh[kd:] > 3)
print(f"λ={lam:.6f} ε={eps:.5f}: dip at η={eta[kd]:.3f}, front (D/ε=3) at η={eta[kc]:.3f}")
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for mu in mus:
    vals, vecs = St.spectrum(mu, k=16)
    j = np.argmin(np.abs(vals - 1)); v = vecs[:, j]
    qp = np.empty(S.N, complex); qp[St.i0:] = v; qp[:St.i0] = v[0]
    q1 = S.st.offset(qp.real, S.c1) + 1j * S.st.offset(qp.imag, S.c1); q2 = S.st.offset(qp.real, S.c2) + 1j * S.st.offset(qp.imag, S.c2)
    th, om = _lin_gl2(S.h, St.i0, complex(mu), lam, eta, St.D1, St.D2, St.Tt1, St.Tt2, St.Ot1, St.Ot2, q1, q2)
    ph = th / (np.abs(th).max()); ph = ph * np.exp(-1j * np.angle(ph[np.argmax(np.abs(ph))]))
    r = ph.real
    p = np.gradient(np.log(np.abs(r) + 1e-300), eta)
    win = (eta > -4.5) & (eta < -2.0)
    zc = np.sum(np.diff(np.sign(r[(eta > -8) & (eta < eta[kc])])) != 0)
    print(f"  μ={mu:.6f} (ν={vals[j]:.6f}): μ/(λ−1)={mu/(lam-1):.4f}, transport exponent 2−μ/ε={2-mu/eps:.3f}; "
          f"median d ln|θ'|/dη on η∈[−4.5,−2] = {np.median(p[win]):.3f}; sign changes of θ' in the layer: {zc}")
    ax[0].plot(eta, r, label=f"μ={mu:.4f}"); ax[1].plot(eta, np.log10(np.abs(r) + 1e-30), label=f"μ={mu:.4f}")
for a in ax:
    a.axvline(eta[kd], color='k', ls=':'); a.axvline(eta[kc], color='r', ls=':'); a.set_xlim(-12, 3); a.set_xlabel('η = ln x')
ax[0].set_ylabel("θ' (normalised)"); ax[1].set_ylabel("log10|θ'|"); ax[0].legend(fontsize=7)
plt.tight_layout(); plt.savefig(f"hl_modes_lam{lam:.4f}.png", dpi=120)
