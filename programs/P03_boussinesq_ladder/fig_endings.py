"""Fig. 4: two ways a branch of self-similar profiles can end.
(a) Stalled-layer ending (ε → 0): the smoothness defect m − 2 oscillates about zero with an exponentially small
    amplitude, so the zeros (smooth profiles) never stop: 2D Boussinesq, Hou–Luo, IPM.
(b) Sonic-cusp ending (CCF, P02): the defect p − 2 oscillates log-periodically about p* − 2 = 0.0058 > 0 as the
    sonic depth δ → 0; after λ₂ it never returns to zero, and the ladder is finite.
(c) The stalled layer: D/ε along the boundary of deep Hou–Luo profiles (D/ε = O(1) up to the front at x_c).
(d) The sonic cusp: den = 1 + λ + HΘ/ξ of CCF profiles approaching the terminal cusp (den → 0 at an interior point)."""
import numpy as np, glob, re, sys, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
P02 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "P02_singularity_ladders")
sys.path.insert(0, P02)
from cusp_fit import parse

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
COL = {'2D Boussinesq': '#2a78d6', 'Hou–Luo': '#eb6834', 'IPM': '#1f9e6e', 'CCF': '#8a5cc2'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})


def load(pattern, ipm=False):
    rows = [np.load(f) for f in sorted(glob.glob(pattern))]
    D = np.vstack([r for r in rows if r.ndim == 2 and len(r)])
    lam, m = D[:, 0], D[:, 2]
    z = 1 / lam if ipm else 1 / (lam - 1)
    o = np.argsort(z); z, f = z[o], (m - 2)[o]
    _, iu = np.unique(np.round(z, 7), return_index=True)
    return z[iu], f[iu]


fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.4))
a = ax[0, 0]
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    from ladder_fit import bq_lam, hl_lam
rungs = {"2D Boussinesq": [1 / (l - 1) for l in bq_lam],
         "Hou–Luo": [1 / (l - 1) for l in list(hl_lam) + [1.08113374, 1.07355524]],
         "IPM": [1 / l for l in (1.0285722975, 0.4721297256, 0.3149621637, 0.2415660984, 0.1987299872, 0.1706180880, 0.15092)]}
# branch scans, restricted to where they are resolved (beyond, only the refined rungs are shown)
for name, pat, ipm, zmax, yr in (("2D Boussinesq", "scan2_[ABCF].npy", False, 9.0, 4e-10),
                                 ("Hou–Luo", "hl_scan_[DE].npy", False, 9.9, 2e-10),
                                 ("IPM", "ipm_branch_dn*.npy", True, 6.3, 1e-10)):
    z, f = load(pat, ipm)
    k = (z > 0.9) & (z < zmax) & (np.abs(f) > 0)
    a.semilogy(z[k], np.abs(f[k]), '-', color=COL[name], lw=1.3, label=name)
    a.plot(rungs[name], [yr * 1.6] * len(rungs[name]), '|', color=COL[name], ms=9, mew=1.6)
a.set_ylim(1e-10, 0.3); a.set_xlim(0.5, 14.2)
a.set_xlabel('1/ε along the branch'); a.set_ylabel('|m − 2|  (smoothness defect)')
a.set_title('(a) stalled-layer ending: zeros never stop', loc='left', fontsize=10)
a.legend(fontsize=8, frameon=False, loc='upper right')
a.text(0.6, 1.2e-9, 'ticks: smooth profiles (refined zeros of m − 2)', fontsize=7.5, color=MUTED)

b = ax[0, 1]
E = np.load(os.path.join(P02, "pdense_E.npy")); C = parse(os.path.join(P02, "pin_C.log"))
d = np.concatenate([E[:, 0], C[:, 0]]); p = np.concatenate([E[:, 2], C[:, 2]])
o = np.argsort(-d); d, p = d[o], p[o]
bm = np.load(os.path.join(P02, "branch_map.npy")); bm = bm[np.argsort(-bm[:, 0])]
b.semilogx(d, p - 2, 'o', color=COL['CCF'], ms=2.5, label='CCF branch near its end')
b.axhline(2.0057717 - 2, color=COL['CCF'], lw=1.0, ls='--', label='p* − 2 = 0.0058 (terminal cusp)')
b.axhline(0, color=MUTED, lw=0.9)
b.invert_xaxis(); b.set_xlabel('sonic depth δ  (→ 0 at the cusp)'); b.set_ylabel('p − 2  (smoothness defect)')
b.set_title('(b) sonic-cusp ending: oscillation about p* ≠ 2, zeros stop', loc='left', fontsize=10)
b.legend(fontsize=8, frameon=False, loc='lower left'); b.set_ylim(-0.001, 0.0075)

c = ax[1, 0]
try:
    from hl_solver import HL
    for f, col in (("hl_q_Fc_lam1.09469700.npy", 0.35), ("hl_q_Fc_lam1.05537100.npy", 0.6), ("hl_q_Fc_lam1.03024804.npy", 0.9)):
        lam = float(re.search(r'lam([0-9.]+?)\.npy', f).group(1))
        S = HL(lam, N=65536, L1=25, L2=75, eta_start=-20); q = S.full(np.load(f)); r = S.march(q)
        xi = np.exp(S.eta); Dh = r['D'] / r['eps']
        sel = (xi > 1e-3) & (xi < 3)
        c.semilogy(xi[sel], Dh[sel], color=plt.cm.Oranges(col), lw=1.3, label=f'Hou–Luo, 1/(λ−1) = {1 / (lam - 1):.1f}')
except Exception as e:
    c.text(0.5, 0.5, f'HL states unavailable ({e})', transform=c.transAxes, ha='center')
c.set_xscale('log'); c.set_xlabel('ξ (along the boundary)'); c.set_ylabel('D / D(0)   (D(0) = 1+λ−A = O(λ−1))')
c.set_title('(c) the stalled layer: D = O(ε) up to the front', loc='left', fontsize=10)
c.legend(fontsize=8, frameon=False, loc='upper left')

e = ax[1, 1]
cc = 0.7
for fn, lab, col in (("mapped_l2_hs0.02_eps0.025.npy", "λ₂ = 0.4713 (last smooth profile)", 0.45),
                     ("pin_C_final_state.npy", "near the cusp, λ = 0.45358", 0.9)):
    dd = np.load(os.path.join(P02, fn))
    eta, phi = dd[0], dd[1]
    lam = 0.471324227767 if 'l2' in fn else float(dd[2][0])
    slope = np.gradient(phi + cc * eta, eta)
    m = (eta > -5) & (eta < 5)
    e.semilogx(np.exp(eta[m]), lam / slope[m], color=plt.cm.Purples(col), lw=1.3, label=lab)
e.set_ylim(0, 1.6); e.set_xlabel('ξ'); e.set_ylabel('den = 1 + λ + HΘ/ξ  (characteristic speed)')
e.set_title('(d) the sonic cusp: den → 0 at an interior point', loc='left', fontsize=10)
e.legend(fontsize=8, frameon=False, loc='upper left')
plt.tight_layout(); plt.savefig('fig_endings.png', dpi=160)
print('saved fig_endings.png')
