"""HL WKB phase (cut at the front, D/ε > Dcut) along the branch, interpolated to the crossings; checks ΔReΦ → π."""
import numpy as np, sys, re, glob
from hl_solver import HL
def phase_cut(f, N, L1, L2, eta_start, Dcut=3.0):
    lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
    q = np.load(f); S = HL(lam, N=N, L1=L1, L2=L2, eta_start=eta_start); r = S.march(S.full(q))
    D, Om, eps, eta = r['D'], r['Om'], r['eps'], S.eta
    Dh = D / eps
    sel = (eta > -3) & (eta < 2); kd = np.argmax(sel) + np.argmin(Dh[sel]); kc = kd + np.argmax(Dh[kd:] > Dcut)
    k = (-1 + np.sqrt(1 + 4j * Om)) / (2j * Dh)
    s = slice(S.i0, kc + 1)
    return lam, eps, r['m'], np.trapezoid(k[s], eta[s])
res = []
for f in sorted(glob.glob("hl_q_E_lam*.npy")):
    res.append(phase_cut(f, 32768, 40, 160, -30))
for f in sorted(glob.glob("hl_q_F_lam*.npy")):
    res.append(phase_cut(f, 65536, 25, 75, -20))
res.sort(key=lambda t: -t[0])
z = np.array([1 / (t[0] - 1) for t in res]); ePhi = np.array([t[3] for t in res])
for t, zz in zip(res, z):
    print(f"z={zz:7.3f} ε={t[1]:.5f} εΦ_cut={t[3].real:.5f}{t[3].imag:+.5f}i")
# crossings of HL (from scans E, F)
rows = np.vstack([np.load(f) for f in ("hl_scan_E.npy", "hl_scan_F.npy")])
lam = rows[:, 0]; zz = 1 / (lam - 1); fm = rows[:, 2] - 2; o = np.argsort(zz); zz, fm = zz[o], fm[o]
_, iu = np.unique(np.round(zz, 6), return_index=True); zz, fm = zz[iu], fm[iu]
cz = [zz[i] - fm[i] * (zz[i + 1] - zz[i]) / (fm[i + 1] - fm[i]) for i in range(len(zz) - 1) if fm[i] * fm[i + 1] < 0]
cz = np.array(cz)
ReP = np.interp(cz, z, ePhi.real) * 2 * cz
ImP = np.interp(cz, z, ePhi.imag) * 2 * cz
print("crossings z:", np.round(cz, 4))
print("Re Φ at crossings:", np.round(ReP, 4)); print("Δ Re Φ:", np.round(np.diff(ReP), 4))
print("Δ Im Φ:", np.round(np.diff(ImP), 4))
