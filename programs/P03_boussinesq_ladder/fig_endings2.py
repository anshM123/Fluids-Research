"""Fig. 5: how the stalled layer ends decides the hierarchy.
(a) Spacings of consecutive profiles in 1/ε: increasing towards π/C for a stalled layer of fixed extent (2D
    Boussinesq, Hou–Luo; π/C = 1.279 derived for Hou–Luo), contracting for the receding layer of IPM.
(b) IPM: anatomy of the phase coefficient I = D₀ Re Φ₀ = ∫Re κ ds with the repaired tracker
    (`ipm_phase_anatomy3.out`): the part before the dip core grows as the dip and front recede, the dip core
    shrinks and levels off, the front side vanishes.
(c) Front position x_f (D/D₀ = 2 beyond the dip) along the IPM branch (`ipm_travel_time.out`, to z = 9.5), against
    the converged fronts of Hou–Luo (x_c ≈ 0.61) and 2D Boussinesq (x_c ≈ 0.72).
(d) Sonic-cusp ending (CCF): the defect oscillates about p* − 2 = 0.0058 > 0, so the zeros stop."""
import numpy as np, glob, re, sys, os, contextlib, io
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
P02 = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "P02_singularity_ladders")
sys.path.insert(0, P02)
from cusp_fit import parse

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
COL = {'2D Boussinesq': '#2a78d6', 'Hou–Luo': '#eb6834', 'IPM': '#1f9e6e', 'CCF': '#8a5cc2'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.4))

# (a) spacings
a = ax[0, 0]
with contextlib.redirect_stdout(io.StringIO()):
    from ladder_fit import bq_lam, hl_lam
z2d = np.array([1 / (l - 1) for l in bq_lam])
zhl = np.array([1 / (l - 1) for l in list(hl_lam) + [1.08113374, 1.07355524]])
zip_ = 1 / np.array([1.0285722975, 0.4721297348, 0.3149618108, 0.2415663353, 0.1987224523, 0.1706180880, 0.15092])
for name, z in (("2D Boussinesq", z2d), ("Hou–Luo", zhl), ("IPM", zip_)):
    z = np.sort(z); n = np.arange(1, len(z))
    a.plot(n + 0.5, np.diff(z), 'o-', color=COL[name], ms=4, lw=1.3, label=name)
a.axhline(1.279, color=COL['Hou–Luo'], ls='--', lw=0.9); a.text(9.6, 1.29, 'π/C = 1.279 (derived)', fontsize=7.5, color=INK)
a.axhspan(1.476, 1.51, color=COL['2D Boussinesq'], alpha=0.12); a.text(0.7, 1.52, '2D asymptote 1.476–1.51', fontsize=7.5, color=INK)
a.set_xlabel('n (between profiles n−1 and n)'); a.set_ylabel('spacing of 1/ε'); a.set_ylim(0.6, 1.65)
a.legend(fontsize=8, frameon=False, loc='lower left')
a.set_title('a   Fixed layer: spacings grow; receding layer: they contract', loc='left', fontsize=10)

# (b) anatomy of the IPM phase (repaired tracker)
b = ax[0, 1]
rows = []
for line in open("ipm_phase_anatomy3.out"):
    p = line.split()
    if len(p) == 13 and re.match(r'^\d', p[1]):
        v = [float(x) for x in p[1:]]                       # z D0 I inner approach core front ...
        rows.append((v[0], v[2], v[3] + v[4], v[5], v[6]))
R = np.array(sorted(rows))
b.stackplot(R[:, 0], R[:, 2], R[:, 3], R[:, 4], colors=['#b9e2d1', COL['IPM'], '#0b4f37'], alpha=0.95,
            labels=['before the dip core (inner layer + approach)', 'dip core (below half depth)', 'front side (to D̂ = 2)'])
b.plot(R[:, 0], R[:, 1], 'k.', ms=4)
b.set_xlabel('z = 1/λ'); b.set_ylabel('I = D₀ Re Φ₀ = ∫ Re κ ds'); b.set_ylim(0, 2.2)
b.legend(fontsize=8, frameon=False, loc='upper left')
b.set_title('b   IPM: the phase grows as the layer recedes', loc='left', fontsize=10)

# (c) front position
c = ax[1, 0]
fr = []
for line in open("ipm_travel_time.out"):
    p = line.split()
    if len(p) == 13 and re.match(r'^\d', p[1]):
        fr.append((float(p[1]), float(p[12])))
F = np.array(sorted(fr))
c.plot(F[:, 0], F[:, 1], 'o-', color=COL['IPM'], ms=4, lw=1.3, label='IPM front (receding)')
c.axhline(0.61, color=COL['Hou–Luo'], ls='--', lw=1.0, label='Hou–Luo front, converged (x_c ≈ 0.61)')
c.axhline(0.72, color=COL['2D Boussinesq'], ls='--', lw=1.0, label='2D Boussinesq front, converged (x_c ≈ 0.72)')
c.set_xlabel('1/ε'); c.set_ylabel('front position x_f  (D/D₀ = 2)'); c.set_ylim(0.5, 1.2)
c.legend(fontsize=8, frameon=False, loc='upper left')
c.set_title('c   Where the stalled layer ends', loc='left', fontsize=10)

# (d) CCF sonic-cusp ending
d_ = ax[1, 1]
E = np.load(os.path.join(P02, "pdense_E.npy")); C = parse(os.path.join(P02, "pin_C.log"))
dd = np.concatenate([E[:, 0], C[:, 0]]); pp = np.concatenate([E[:, 2], C[:, 2]])
o = np.argsort(-dd); dd, pp = dd[o], pp[o]
d_.semilogx(dd, pp - 2, 'o', color=COL['CCF'], ms=2.5, label='CCF branch near its end')
d_.axhline(2.0057717 - 2, color=COL['CCF'], lw=1.0, ls='--', label='p* − 2 = 0.0058 (terminal cusp)')
d_.axhline(0, color=MUTED, lw=0.9)
d_.invert_xaxis(); d_.set_xlabel('sonic depth δ  (→ 0 at the cusp)'); d_.set_ylabel('p − 2  (smoothness defect)')
d_.legend(fontsize=8, frameon=False, loc='lower left'); d_.set_ylim(-0.001, 0.0075)
d_.set_title('d   Sonic cusp: oscillation about p* ≠ 2, zeros stop', loc='left', fontsize=10)
plt.tight_layout(); plt.savefig('fig_endings2.png', dpi=170); print('saved fig_endings2.png')
