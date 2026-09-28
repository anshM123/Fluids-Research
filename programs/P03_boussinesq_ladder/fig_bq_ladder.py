"""Main figure for P03: the 2D Boussinesq ladder and its Hou–Luo analogue."""
import numpy as np, glob, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator


def load_scans(pattern):
    rows = [np.load(f) for f in sorted(glob.glob(pattern))]
    D = np.vstack([r for r in rows if r.ndim == 2 and len(r)])
    D = D[np.argsort(D[:, 0])]
    _, iu = np.unique(np.round(D[:, 0], 9), return_index=True)
    return D[iu]


def crossings(D):
    lam, m = D[:, 0], D[:, 2]
    z = 1 / (lam - 1); f = m - 2; o = np.argsort(z); z, f = z[o], f[o]
    out = []
    for i in range(len(z) - 1):
        if f[i] * f[i + 1] < 0:
            out.append(z[i] - f[i] * (z[i + 1] - z[i]) / (f[i + 1] - f[i]))
    return np.array(out)


refined = {}
for f in glob.glob("cross_l*_nb32.log"):
    for line in open(f):
        mm = re.search(r"CROSSING Nb32_hs0.025_ss-20.0_sm100.0: λ = ([0-9.]+)", line)
        if mm:
            refined[int(re.search(r"cross_l(\d)", f).group(1))] = float(mm.group(1))
Bq = load_scans("scan2_*.npy")
H = load_scans("hl_scan_E.npy")
Hd = load_scans("hl_scan_D.npy")
fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))
a = ax[0, 0]
zB = 1 / (Bq[:, 0] - 1)
a.plot(zB, Bq[:, 2] - 2, 'k.-', ms=2.5, lw=0.8)
a.set_yscale('symlog', linthresh=1e-9); a.axhline(0, color='r', lw=0.6)
for n, l in sorted(refined.items()):
    a.plot(1 / (l - 1), 0, 'o', color='tab:blue', ms=6, zorder=5)
for n in range(9):
    a.axvline(1.0863 + 1.4187 * n, color='tab:gray', ls=':', lw=0.7)
a.set_xlabel('z = 1/(λ − 1)'); a.set_ylabel('m(λ) − 2   (symlog)')
a.set_title('(a) 2D Boussinesq: smoothness defect along the branch\n(blue: smooth profiles; dotted: two-point law of Wang et al.)', fontsize=9)
a.set_xlim(0.3, 13)
b = ax[0, 1]
Hall = np.vstack([Hd[1 / (Hd[:, 0] - 1) < 4], H])
zH = 1 / (Hall[:, 0] - 1); o = np.argsort(zH)
b.plot(zH[o], Hall[o, 2] - 2, 'k.-', ms=1.5, lw=0.6)
b.set_yscale('symlog', linthresh=1e-11); b.axhline(0, color='r', lw=0.6)
for zc in crossings(Hall):
    b.plot(zc, 0, 'o', color='tab:green', ms=4)
b.set_xlabel('z = 1/(λ − 1)'); b.set_ylabel('m(λ) − 2'); b.set_xlim(0.3, 16)
b.set_title('(b) Hou–Luo boundary model: same ladder', fontsize=9)
c = ax[1, 0]
from bq_solver import BQ
from bq_newton import full
for f, col in (("Y2_lam1.9200_Nb32_hs0.025.npy", 'tab:blue'), ("Y2_B_lam1.4000.npy", 'tab:orange'),
               ("Y2_B_lam1.2000.npy", 'tab:green'), ("Y2_B_lam1.1200.npy", 'tab:red'), ("Y2_B_lam1.0800.npy", 'tab:purple')):
    try:
        lam = float(re.search(r'lam([0-9.]+?)(?:_|\.npy)', f).group(1))
        B = BQ(lam); Y = np.load(f); X = full(B, Y)
        Ur = X @ B.Db.T
        A = -np.mean(Ur[B.i0:B.i0 + 40, 0]); eps = 1 + lam - A
        s = B.s; sel = (s > -4) & (s < 1.5)
        c.plot(np.exp(s[sel]), ((1 + lam) + Ur[sel, 0]) / eps, color=col, label=f"λ = {lam:.2f} (ε = {eps:.3f})")
    except Exception as e:
        print("skip", f, e)
c.set_ylim(0, 6); c.set_xscale('log'); c.set_xlabel('r (along the boundary)'); c.set_ylabel('(V_r/r)/ε  on the boundary')
c.legend(fontsize=7); c.set_title('(c) quasi-stagnant boundary region and front (x_c ≈ 0.5)', fontsize=9)
d = ax[1, 1]
zc_B = np.array([1 / (l - 1) for n, l in sorted(refined.items())])
nB = np.arange(len(zc_B))
d.plot(nB, zc_B, 'o-', color='tab:blue', label='2D Boussinesq (this work)')
zc_H = crossings(Hall)
d.plot(np.arange(len(zc_H)), zc_H, 's-', color='tab:green', ms=4, label='Hou–Luo (this work)')
nn = np.arange(0, 12)
d.plot(nn, 1.0863 + 1.4187 * nn, 'k:', lw=0.8, label='two-point law (Wang et al.)')
d.set_xlabel('n (instability index)'); d.set_ylabel('1/(λ_n − 1)'); d.legend(fontsize=8, loc='upper left')
d.set_title('(d) the ladder: 1/(λ_n − 1) vs n;  inset: WKB phase gained per profile / π', fontsize=9)
ins = d.inset_axes([0.58, 0.1, 0.38, 0.35])
dphi_B = np.array([2.743, 2.981, 3.030, 3.072, 3.110]) / np.pi
dphi_H = np.array([3.025, 3.052, 3.074, 3.069, 3.092, 3.094]) / np.pi
ins.plot(np.arange(1, len(dphi_B) + 1), dphi_B, 'o-', color='tab:blue', ms=4)
ins.plot(np.arange(3, 3 + len(dphi_H)), dphi_H, 's-', color='tab:green', ms=3)
ins.axhline(1, color='r', lw=0.7); ins.set_ylim(0.85, 1.03); ins.tick_params(labelsize=7)
ins.set_xlabel('interval n → n+1', fontsize=7); ins.set_ylabel('ΔReΦ/π', fontsize=7)
plt.tight_layout(); plt.savefig('fig_bq_ladder.png', dpi=150)
print("refined:", refined); print("HL crossings z:", np.round(zc_H, 4))
