"""Hou–Luo profiles approaching λ → 1: structure near the front, scaled by the front position x_c.
Reports x_c (maximum of Θ'/... slope; here the point where D̂ = D/ε first exceeds 1/ε^{1/2} beyond the dip), Θ(x_c),
the layer vorticity Ω·√(x_c − ξ) (constant if Ω ~ A/√(x_c − ξ)), the outer speed D/√(ξ − x_c) (constant if
D ~ k√(ξ − x_c)) and the outer vorticity near the front."""
import numpy as np, glob, re
from hl_solver import HL
files = sorted(glob.glob("hl_q_Fc_lam*.npy"), key=lambda s: -float(re.search(r'lam([0-9.]+?)\.npy', s).group(1)))
for f in files:
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    S = HL(lam, N=65536, L1=25, L2=75, eta_start=-20); q = S.full(np.load(f)); r = S.march(q)
    eps, eta = r['eps'], S.eta; xi = np.exp(eta); D = r['D']; Dh = D / eps; Th, Om = r['Th'], r['Om']
    sel = (eta > -3) & (eta < 2); kd = np.argmax(sel) + np.argmin(np.where(sel, Dh, np.inf)[sel])
    kf = kd + np.argmax(D[kd:] > 0.5 * np.sqrt(eps))       # front: D reaches √ε/2
    xc = xi[kf]
    out = [f"z={1/(lam-1):6.2f} ε={eps:.5f}: x_dip={xi[kd]:.5f} x_c={xc:.5f} Θ(x_c)={Th[kf]:.5f} Ω(x_c)={Om[kf]:.4f}"]
    for t in (0.5, 0.2, 0.1, 0.05, 0.02):
        k = np.argmin(np.abs(xi - xc * (1 - t)))
        out.append(f"  in: 1−ξ/x_c={t:5.2f}: Ω√(1−ξ/x_c)={Om[k]*np.sqrt(t):8.4f}  Θ={Th[k]:.4f}  D̂={Dh[k]:.4f}")
    for t in (0.02, 0.05, 0.1, 0.2, 0.5):
        k = np.argmin(np.abs(xi - xc * (1 + t)))
        out.append(f"  out: ξ/x_c−1={t:5.2f}: D/√t={D[k]/np.sqrt(t):7.4f}  Ω={Om[k]:8.4f}  Θ={Th[k]:.4f}")
    print("\n".join(out), flush=True)
