"""The Hou–Luo λ → 1 limit problem and the derived spacing constant.
(a) limit solutions for several front vorticities Ω_f: layer vorticity Ω(ξ) (ξ < 1) and outer vorticity (ξ > 1);
(b) outer speed D = 2 + U/ξ, which vanishes on the stalled layer;
(c) phase constant as spacing π/C(Ω_f) against 1/Ω_f, with the Ω_f → ∞ extrapolation;
(d) measured rung spacings Δz_n = z_n − z_{n−1} against 1/z (z = 1/(λ − 1)) with the derived value at 1/z = 0."""
import numpy as np, glob, re, contextlib, io
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    from ladder_fit import hl_lam

INK, MUTED, GRID, ACC = '#0b0b0b', '#52514e', '#e4e3df', '#eb6834'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
J = 160; jj = np.arange(1, 2 * J, 2)


def layer_fields(b, th):
    x = np.cos(th)
    Om = (b[:, None] * np.cos(np.outer(jj, th))).sum(axis=0) / np.sin(th)
    Th = (b[:, None] * ((np.sin(jj * np.pi / 2) / jj)[:, None] - np.sin(np.outer(jj, th)) / jj[:, None])).sum(axis=0)
    return x, Om, Th


def phase(b, n=200000):
    u = (np.arange(n) + 0.5) / n * np.sqrt(np.pi / 2)
    th = u ** 2; dth = 2 * u * (u[1] - u[0])
    tot = 0
    for sl in np.array_split(np.arange(n), 20):
        x, Om, Th = layer_fields(b, th[sl])
        Dh = 2 * Th / (x * Om)
        tot += np.sum((-1 + np.sqrt(1 + 4j * Om)) / (2j * Dh) * np.tan(th[sl]) * dth[sl])
    return tot


rows = {}
for f in glob.glob("hl_limit_Omf*.npz"):
    rows[float(re.search(r"Omf([0-9.]+)\.npz", f).group(1))] = f
Omfs = np.array(sorted(rows))
C = np.array([2 * phase(np.load(rows[o])['b']).real for o in Omfs])
spacing = np.pi / C
# Ω_f → ∞: fit C = C_∞ + A Ω_f^{−q} on Ω_f ≥ 10
sel = Omfs >= 10
from scipy.optimize import curve_fit
p, _ = curve_fit(lambda x, Ci, A, q: Ci + A * x ** (-q), Omfs[sel], C[sel], p0=(2.456, 10, 3))
C_inf = p[0]; a_derived = np.pi / C_inf

fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
cols = plt.cm.viridis(np.linspace(0.05, 0.85, 4))
show = [o for o in (1.0, 3.0, 6.0, 18.0) if o in rows] or list(Omfs[::4])
for c, o in zip(cols, show):
    d = np.load(rows[o]); b = d['b']
    th = np.linspace(1e-4, np.pi / 2 - 1e-4, 800)
    x, Om, Th = layer_fields(b, th)
    yy, Oo, D = d['y'], d['Om_o'], d['D']
    k = yy < 3
    ax[0, 0].plot(x, Om, color=c, lw=1.4, label=f'Ω_f = {o:g}')
    ax[0, 0].plot(yy[k], Oo[k], color=c, lw=1.4, ls='--')
    ax[0, 1].plot(np.r_[x[::-1], yy[k]], np.r_[np.zeros_like(x), D[k]], color=c, lw=1.4, label=f'Ω_f = {o:g}')
xx = np.linspace(0, 0.999, 400)
ax[0, 0].plot(xx, 2 * xx / np.sqrt(1 - xx ** 2), color=MUTED, lw=1.0, ls=':', label='Ω_f = 0 (exact)')
ax[0, 0].set_xlim(0, 3); ax[0, 0].set_ylim(0, 8)
ax[0, 0].axvline(1, color=MUTED, lw=0.8)
ax[0, 0].text(0.5, 7.3, 'stalled layer\nΩ slaved to Θ′', ha='center', color=MUTED, fontsize=8)
ax[0, 0].text(2.0, 7.3, 'outer region\nΩ transported', ha='center', color=MUTED, fontsize=8)
ax[0, 0].set_xlabel('ξ / x_c'); ax[0, 0].set_ylabel('vorticity Ω')
ax[0, 0].set_title('(a) λ → 1 limit problem: vorticity', loc='left', fontsize=10)
ax[0, 0].legend(fontsize=7.5, frameon=False, loc='center right')
ax[0, 1].set_xlim(0, 3); ax[0, 1].axvline(1, color=MUTED, lw=0.8)
ax[0, 1].set_xlabel('ξ / x_c'); ax[0, 1].set_ylabel('self-similar speed D = 2 + U/ξ')
ax[0, 1].set_title('(b) D = 0 on the stalled layer, D ∝ √(ξ−1) past the front', loc='left', fontsize=10)
ax[0, 1].legend(fontsize=7.5, frameon=False)
# (c)
c = ax[1, 0]
c.plot(1 / Omfs, spacing, 's', color=ACC, ms=4.5, mec='white', mew=0.6, label='limit problem, π/C(Ω_f)')
xf = np.linspace(0, 1 / 10, 50); c.plot(xf, np.pi / (p[0] + p[1] * np.maximum(xf, 1e-9) ** p[2]), '-', color=ACC, lw=1.0)
c.plot([0], [a_derived], '*', color=INK, ms=11, label=f'Ω_f → ∞:  π/C_∞ = {a_derived:.3f}')
c.set_xlim(-0.01, 1.05); c.set_xlabel('1/Ω_f'); c.set_ylabel('predicted spacing π/C of z = 1/(λ−1)')
c.set_title('(c) the spacing constant from the limit problem', loc='left', fontsize=10)
c.legend(fontsize=8, frameon=False, loc='lower left')
# (d)
d_ = ax[1, 1]
lam = list(hl_lam) + [1.08113374, 1.07355524]
z = np.array([1 / (l - 1) for l in lam]); dz = np.diff(z); zm = 0.5 * (z[1:] + z[:-1])
d_.plot(1 / zm, dz, 's', color=ACC, ms=6, mec='white', mew=0.8, label='measured Δz_n (n = 1–10)')
for lo, ls in ((2, '-'), (4, '--')):
    pp = np.polyfit(1 / zm[lo:], dz[lo:], 1); xx = np.linspace(0, 1 / zm[lo], 20)
    d_.plot(xx, np.polyval(pp, xx), ls, color=MUTED, lw=0.9, label=f'linear fit n ≥ {lo + 1}: {pp[1]:.4f}')
d_.plot([0], [a_derived], '*', color=INK, ms=12, label=f'derived π/C_∞ = {a_derived:.3f}')
d_.set_xlim(-0.005, 0.37); d_.set_ylim(1.245, 1.286); d_.set_xlabel('1/z at the interval midpoint'); d_.set_ylabel('rung spacing Δz_n')
d_.set_title('(d) measured spacings extrapolate to the derived constant', loc='left', fontsize=10)
d_.legend(fontsize=8, frameon=False, loc='lower left')
plt.tight_layout(); plt.savefig('fig_limit.png', dpi=160)
print(f"C_∞ = {C_inf:.5f} (fit C = C_∞ + A Ω_f^−q, A = {p[1]:.3g}, q = {p[2]:.3f}); derived spacing π/C_∞ = {a_derived:.5f}")
for o, cc in zip(Omfs, C):
    print(f"  Ω_f = {o:7.3f}: C = {cc:.5f}, π/C = {np.pi / cc:.5f}")
