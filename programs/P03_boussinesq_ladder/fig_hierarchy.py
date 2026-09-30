"""Fig. 1: the hierarchies. 1/ε_n against n, where ε is the exponent of the transported scalar (λ − 1 for 2D
Boussinesq and Hou–Luo, λ for IPM and CCF). Open symbols: profiles reported before this work; filled: new.
CCF ends at its sonic cusp (P02); the others continue (dashed: asymptotic spacing; Hou–Luo: derived π/C_∞)."""
import numpy as np, re, glob, contextlib, io
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    from ladder_fit import bq_lam, hl_lam

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
COL = {'2D Boussinesq': '#2a78d6', 'Hou–Luo': '#eb6834', 'IPM': '#1f9e6e', 'CCF': '#8a5cc2'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})


FINE = {5: 0.1706180880, 6: 0.15092}      # h_s = 0.0125 values supersede h_s = 0.025 (PREDICTIONS_IPM.md)


def ipm_rungs():
    out = {}
    for f in glob.glob("ipm_rung[0-9].out"):
        n = int(re.search(r"ipm_rung(\d)\.out", f).group(1))
        for line in open(f):
            m = re.search(r"CROSSING .*: λ = ([0-9.]+)", line)
            if m:
                out[n] = FINE.get(n, float(m.group(1)))
    return [out[n] for n in sorted(out)]


data = {
    '2D Boussinesq': ([1 / (l - 1) for l in bq_lam], 4),            # n < 4 reported (Wang et al. 2025)
    'Hou–Luo': ([1 / (l - 1) for l in list(hl_lam) + [1.08113374, 1.07355524]], 1),   # n = 0 (Chen–Hou)
    'IPM': ([1 / l for l in ipm_rungs()], 5),                       # n ≤ 3 (Wang et al. 2025), n = 4 (Wang–Léger–Lai–Buckmaster 2025)
    'CCF': ([1 / l for l in (1.180777662899, 0.6057337012, 0.471324227767)], 3),
}
fig, ax = plt.subplots(1, 1, figsize=(6.4, 4.8))
for name, (z, nknown) in data.items():
    n = np.arange(len(z)); c = COL[name]
    ax.plot(n[:nknown], z[:nknown], 'o', mfc='white', mec=c, mew=1.4, ms=7, zorder=3)
    if len(z) > nknown:
        ax.plot(n[nknown:], z[nknown:], 'o', color=c, mec='white', mew=0.8, ms=7, zorder=3)
    ax.plot(n, z, '-', color=c, lw=1.0, alpha=0.6, label=name)
    if name in ('2D Boussinesq', 'Hou–Luo') and len(z) > 3:
        slope = 1.279 if name == 'Hou–Luo' else np.diff(z)[-1]
        nn = np.array([n[-1], n[-1] + 2.5]); ax.plot(nn, z[-1] + slope * (nn - n[-1]), '--', color=c, lw=1.0)
zs = 1 / 0.4535843
ax.plot([2, 2.6], [data['CCF'][0][-1], zs], ':', color=COL['CCF'], lw=1.2)
ax.plot(2.6, zs, '*', color=COL['CCF'], ms=12, zorder=4)
ax.annotate('sonic cusp:\nCCF ladder ends', (2.6, zs), textcoords='offset points', xytext=(8, -26), fontsize=8,
            color=COL['CCF'])
ax.set_xlabel('instability index n  (the n-th profile has n unstable modes)')
ax.set_ylabel('1/ε_n   (ε = exponent of the transported scalar)')
ax.plot([], [], 'o', mfc='white', mec=MUTED, label='reported before'); ax.plot([], [], 'o', color=MUTED, label='this work')
ax.legend(fontsize=8, frameon=False, loc='upper left')
ax.set_xlim(-0.4, 12.9); ax.set_ylim(0, 16.5)
ins = ax.inset_axes([0.60, 0.08, 0.37, 0.34])
lam_ipm = np.array(ipm_rungs()); lc = 0.0370
zz = 1 / (lam_ipm - lc); nn = np.arange(len(zz))
ins.plot(nn, zz, 'o', color=COL['IPM'], ms=4.5, mec='white', mew=0.6)
pp = np.polyfit(nn[1:], zz[1:], 1); ins.plot(nn, np.polyval(pp, nn), '-', color=COL['IPM'], lw=0.9, alpha=0.7)
ins.set_title(f'IPM: 1/(λ − λ_c), λ_c = {lc}', fontsize=7.5, color=INK)
ins.tick_params(labelsize=7); ins.set_xlabel('n', fontsize=7.5)
ins.text(0.05, 0.8, f'spacing {pp[0]:.3f} (±0.3 %)', transform=ins.transAxes, fontsize=7, color=MUTED)
plt.tight_layout(); plt.savefig('fig_hierarchy.png', dpi=170)
for name, (z, k) in data.items():
    print(name, np.round(z, 4).tolist())
