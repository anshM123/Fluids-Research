"""Inner structure of the Hou–Luo dip/front in the transport-uniform coordinate ζ = ln|Θ̄| (D ∂_η = (λ−1) ∂_ζ).
For the deep branch states (hl_q_F_*, down to ε = 0.012) report, relative to the dip of D̂ = D/ε:
  - D̂_min, the dip-to-front distance in η and in ζ, the ζ-extent of D̂ < 2 D̂_min;
  - D̂/D̂_min sampled at fixed ζ − ζ_dip (collapse test: columns constant in ε if the dip has a limit shape in ζ);
  - Ω̄ ξ (local circulation density) at the dip."""
import numpy as np, glob, re
from hl_solver import HL
rows = []
for f in sorted(glob.glob("hl_q_F_lam1.0*.npy"), key=lambda s: -float(re.search(r'lam([0-9.]+?)\.npy', s).group(1))):
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    S = HL(lam, N=65536, L1=25, L2=75, eta_start=-20); q = S.full(np.load(f)); r = S.march(q)
    eps = r['eps']; eta = S.eta; D = (1 + lam) + q; Dh = D / eps; Th = r['Th']; zeta = np.log(np.abs(Th))
    sel = (eta > -3) & (eta < 2); kd = np.argmax(sel) + np.argmin(Dh[sel])
    kf = kd + np.argmax(Dh[kd:] > 1.0 / eps ** 0.5)        # front: D̂ = ε^{-1/2}, i.e. D = ε^{1/2}
    k2i = kd - np.argmax(Dh[kd::-1] > 2 * Dh[kd])           # inner side where D̂ = 2 D̂_min
    k2o = kd + np.argmax(Dh[kd:] > 2 * Dh[kd])
    zs = zeta - zeta[kd]
    lo_, hi_ = max(0, kd - 15000), kd + 15000
    samp = [np.interp(z0, zs[lo_:hi_], Dh[lo_:hi_]) / Dh[kd] for z0 in (-1.0, -0.5, -0.25, 0.25, 0.5)]
    rows.append((eps, Dh[kd], eta[kf] - eta[kd], zs[kf], zs[k2o] - zs[k2i], eta[k2o] - eta[k2i], *samp))
    print(f"ε={eps:.5f} D̂_min={Dh[kd]:.4f}  dip→front: Δη={eta[kf]-eta[kd]:.4f} Δζ={zs[kf]:.3f};  "
          f"width(D̂<2D̂min): Δη={eta[k2o]-eta[k2i]:.4f} Δζ={zs[k2o]-zs[k2i]:.3f};  D̂/D̂min at ζ−ζd=−1,−.5,−.25,+.25,+.5: "
          + " ".join(f"{v:.3f}" for v in samp), flush=True)
a = np.array(rows); le = np.log(a[:, 0])
for j, name in ((1, "D̂_min"), (2, "Δη dip→front"), (4, "Δζ width"), (5, "Δη width")):
    p = np.polyfit(le[-6:], np.log(np.abs(a[-6:, j])), 1)
    print(f"  fit over the 6 smallest ε: {name} ∝ ε^{p[0]:.3f}")
np.save("hl_inner_collapse.npy", a)
