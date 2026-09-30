"""Fig. 2: one branch, one phase.
(a) The smoothness defect |m − 2| along the continuous branch of least-singular profiles. Smooth profiles are its
    zeros: a resonance every half-turn of the phase (2D Boussinesq, Hou–Luo, IPM).
(b) The stalled layer that carries the phase: D/D(0) along the boundary of deep Hou–Luo profiles.
(c) The phase coefficient Re εΦ measured on the finite-ε profiles converges to the value a₀ = 1.2283 of the
    λ → 1 limit problem (ASYMPTOTICS §8).
(d) The rung spacing converges to the derived π/(2a₀) = 1.279."""
import numpy as np, glob, re, contextlib, io
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    from ladder_fit import bq_lam, hl_lam

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
COL = {'2D Boussinesq': '#2a78d6', 'Hou–Luo': '#eb6834', 'IPM': '#1f9e6e'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
A0 = 1.2283                                      # Re a₀ of the Ω_f → ∞ limit problem (hl_limit_table.txt)


def load(pattern, ipm=False):
    D = np.vstack([r for r in (np.load(f) for f in sorted(glob.glob(pattern))) if r.ndim == 2 and len(r)])
    z = 1 / D[:, 0] if ipm else 1 / (D[:, 0] - 1)
    o = np.argsort(z); z, f = z[o], (D[:, 2] - 2)[o]
    _, iu = np.unique(np.round(z, 7), return_index=True)
    return z[iu], f[iu]


hl_all = list(hl_lam) + [1.08113374, 1.07355524]
rungs = {"2D Boussinesq": [1 / (l - 1) for l in bq_lam], "Hou–Luo": [1 / (l - 1) for l in hl_all],
         "IPM": [1 / l for l in (1.0285722975, 0.4721297256, 0.3149621637, 0.2415660984, 0.1987299872, 0.1706180880, 0.15092)]}
fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.4))
a = ax[0, 0]
for name, pat, ipm, zmax, yr in (("2D Boussinesq", "scan2_[ABCF].npy", False, 9.0, 4e-10),
                                 ("Hou–Luo", "hl_scan_[DE].npy", False, 9.9, 2e-10),
                                 ("IPM", "ipm_branch_dn*.npy", True, 6.3, 1e-10)):
    z, f = load(pat, ipm)
    k = (z > 0.9) & (z < zmax) & (np.abs(f) > 0)
    a.semilogy(z[k], np.abs(f[k]), '-', color=COL[name], lw=1.3, label=name)
    a.plot(rungs[name], [yr * 1.6] * len(rungs[name]), '|', color=COL[name], ms=9, mew=1.6)
a.set_ylim(1e-10, 0.3); a.set_xlim(0.5, 14.2)
a.set_xlabel('1/ε along the branch  (ε = λ−1; IPM: λ)'); a.set_ylabel('|m − 2|  (smoothness defect)')
a.set_title('(a) smooth profiles are resonances of one branch', loc='left', fontsize=10)
a.legend(fontsize=8, frameon=False, loc='upper right')
a.text(0.6, 1.2e-9, 'ticks: smooth profiles', fontsize=7.5, color=MUTED)

b = ax[0, 1]
try:
    from hl_solver import HL
    for f, col in (("hl_q_Fc_lam1.09469700.npy", 0.35), ("hl_q_Fc_lam1.05537100.npy", 0.6), ("hl_q_Fc_lam1.03024804.npy", 0.9)):
        lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
        S = HL(lam, N=65536, L1=25, L2=75, eta_start=-20); r = S.march(S.full(np.load(f)))
        xi = np.exp(S.eta); sel = (xi > 1e-3) & (xi < 3)
        b.semilogy(xi[sel], (r['D'] / r['eps'])[sel], color=plt.cm.Oranges(col), lw=1.3,
                   label=f'Hou–Luo, 1/(λ−1) = {1 / (lam - 1):.1f}')
except Exception as e:
    b.text(0.5, 0.5, f'HL states unavailable ({e})', transform=b.transAxes, ha='center')
b.set_xscale('log'); b.set_xlabel('ξ along the boundary'); b.set_ylabel('D / D(0)   (radial self-similar speed)')
b.set_title('(b) the stalled layer: D = O(λ−1) up to the front', loc='left', fontsize=10)
b.legend(fontsize=8, frameon=False, loc='upper left')

c = ax[1, 0]
rows = [(float(m.group(1)), float(m.group(2)), float(m.group(3))) for line in open("hl_wkb2.out")
        for m in [re.search(r"z=\s*([\d.]+) ε=([\d.]+) εΦ_cut=([\d.]+)", line)] if m]
eps = np.array([1 / r[0] for r in rows]); ph = np.array([r[2] for r in rows])
c.plot(eps, ph, 's', color=COL['Hou–Luo'], ms=4, mec='white', mew=0.5, label='finite-ε profiles (measured)')
c.plot([0], [A0], '*', color=INK, ms=12, label=f'λ → 1 limit problem: a₀ = {A0}')
c.set_xlim(-0.006, 0.24); c.set_xlabel('λ − 1 = 1/z  along the branch'); c.set_ylabel('phase coefficient Re Φ/(2z)')
c.set_title('(c) the phase converges to the limit-problem value', loc='left', fontsize=10)
c.legend(fontsize=8, frameon=False, loc='lower left')

d = ax[1, 1]
z = np.array([1 / (l - 1) for l in hl_all]); dz = np.diff(z); zm = 0.5 * (z[1:] + z[:-1])
d.plot(1 / zm, dz, 's', color=COL['Hou–Luo'], ms=6, mec='white', mew=0.8, label='Hou–Luo rung spacings Δz_n')
pp = np.polyfit(1 / zm[2:], dz[2:], 1); xx = np.linspace(0, 1 / zm[2], 20)
d.plot(xx, np.polyval(pp, xx), '-', color=MUTED, lw=0.9, label=f'linear extrapolation → {pp[1]:.3f}')
d.plot([0], [np.pi / (2 * A0)], '*', color=INK, ms=12, label=f'derived π/(2a₀) = {np.pi / (2 * A0):.3f}')
d.set_xlim(-0.005, 0.37); d.set_ylim(1.245, 1.286)
d.set_xlabel('1/z at the interval midpoint  (z = 1/(λ−1))'); d.set_ylabel('spacing of consecutive smooth profiles')
d.set_title('(d) the spacing constant is derived, not fitted', loc='left', fontsize=10)
d.legend(fontsize=8, frameon=False, loc='lower left')
plt.tight_layout(); plt.savefig('fig_mechanism.png', dpi=160)
print('saved fig_mechanism.png')
