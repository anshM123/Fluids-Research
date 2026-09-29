"""Structure of the lowest Hou–Luo unstable modes along the deep branch, in the variables that make the transport
uniform: ζ = ln|Θ̄| (so that D ∂_η = (λ−1) ∂_ζ) and φ = θ'/Θ̄. In these variables the temperature equation is exactly
    φ_ζ + ν̂ φ = −q'/D,      ν̂ = μ/(λ−1),
so where φ is flat the perturbation strain must be q' = −ν̂ D φ. The script reports, for each mode:
  - φ in the stalled layer (flatness |φ_ζ|/|φ| away from the front), and the check q' ≈ −ν̂ D φ there;
  - the position and height of the peak of |φ| relative to the dip of D̂ = D/ε;
  - the number of sign changes of φ in the layer and across the front;
and saves φ(ζ), q'(ζ)/(ε ν̂), D̂(ζ) for collapse plots.
usage: python3 hl_deep_modes.py ETA0 STATE nuhat1,nuhat2,... [N0 L1 L2]   (defaults: the F states, 65536 25 75)"""
import numpy as np, sys, re
from hl_stability import HLStab, _lin_gl2

eta0 = float(sys.argv[1]); f = sys.argv[2]; nus = [float(x) for x in sys.argv[3].split(',')]
N0, L1_0, L2_0 = (int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])) if len(sys.argv) > 6 else (65536, 25.0, 75.0)
h = (L1_0 + L2_0) / N0
K = max(0, int(np.ceil((abs(eta0) + 2 - L1_0) / h / 4096.0)) * 4096)
lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
q0 = np.load(f); q = np.concatenate([np.full(K, q0[0]), q0])
St = HLStab(lam, q, N0 + K, L1=L1_0 + K * h, L2=L2_0, eta_start=eta0); S = St.S
d = lam - 1; eps = St.eps; eta = S.eta; i0 = St.i0
D = (1 + lam) + St.q; Dh = D / eps
Th = St.Th
zeta = np.log(np.abs(Th) + 1e-300)
sel = (eta > -3) & (eta < 2)
kd = np.argmax(sel) + np.argmin(np.where(sel, Dh, np.inf)[sel])       # dip of D̂
kc = kd + np.argmax(Dh[kd:] > 3)                                      # front: D̂ = 3 outward of the dip
print(f"HL λ={lam:.6f} z={1/d:.2f} ε={eps:.5f}: dip η={eta[kd]:.4f} (D̂_min={Dh[kd]:.4f}), front η={eta[kc]:.4f}, "
      f"ζ_dip={zeta[kd]:.3f}, ζ_front={zeta[kc]:.3f}")
out = {'eta': eta[i0:], 'zeta': zeta[i0:] - zeta[kd], 'Dh': Dh[i0:], 'lam': lam, 'eps': eps}
for nu in nus:
    mu = nu * d
    vals, vecs = St.spectrum(mu, k=16)
    j = np.argmin(np.abs(vals - 1)); v = vecs[:, j]
    qp = np.empty(S.N, complex); qp[i0:] = v; qp[:i0] = v[0]
    q1 = S.st.offset(qp.real, S.c1) + 1j * S.st.offset(qp.imag, S.c1)
    q2 = S.st.offset(qp.real, S.c2) + 1j * S.st.offset(qp.imag, S.c2)
    th, om = _lin_gl2(S.h, i0, complex(mu), lam, eta, St.D1, St.D2, St.Tt1, St.Tt2, St.Ot1, St.Ot2, q1, q2)
    ph = th / Th
    ph = ph * np.exp(-1j * np.angle(ph[kd]))                             # real, positive at the dip
    qq = qp * np.exp(-1j * np.angle(th[kd] / Th[kd]))
    r = ph.real; qr = qq.real
    lay = (eta > eta[kd] - 2.5) & (eta < eta[kd] - 1.0)                   # interior of the stalled layer
    dz = np.gradient(zeta, eta)
    flat = np.median(np.abs(np.gradient(r, eta)[lay] / dz[lay]) / np.abs(r[lay]))
    ratio = np.median(qr[lay] / (-nu * D[lay] * r[lay]))
    kp = np.argmax(np.abs(r[i0:kc + 200])) + i0
    zc_lay = np.sum(np.diff(np.sign(r[i0:kd])) != 0)
    zc_all = np.sum(np.diff(np.sign(r[i0:kc + int(1.0 / h)])) != 0)
    print(f"  ν̂={nu:.4f} (ν={vals[j].real:+.6f}{vals[j].imag:+.1e}i): layer φ={np.median(r[lay]):+.4f} "
          f"flatness |φ_ζ/φ|={flat:.2e}, q'/(−ν̂Dφ)={ratio:.4f}; peak |φ|={abs(r[kp]):.3f} at η−η_dip={eta[kp]-eta[kd]:+.4f}; "
          f"φ(front)={r[kc]:+.4f}, φ(η_dip+1)={r[kd + int(1.0 / h)]:+.4f}; sign changes: layer {zc_lay}, to front+1 {zc_all}",
          flush=True)
    out[f'phi_{nu:.4f}'] = r[i0:]; out[f'q_{nu:.4f}'] = qr[i0:] / (eps * nu)
np.savez(f"hldeepmodes_lam{lam:.6f}.npz", **out)
