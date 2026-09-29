"""Fig. 2 of the manuscript: the instability ladder in 2D Boussinesq and Hou–Luo.
(a) unstable eigenvalues μ_k of the n-th profile; (b) right-half-plane counts; (c) μ_k/(λ_n−1) against λ_n − 1
with linear extrapolation (the arithmetic ladder); (d) Hou–Luo eigenfunctions localized at the front (n = 8)."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ladder_fit import bq_lam, bq, hl_lam, hl
C2D, CHL, GRID, INK, MUTED = '#2a78d6', '#eb6834', '#e4e3df', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
hl_cx = {8: [0.587182 + 0.022206j], 9: [0.606624 + 0.056719j], 10: [0.623297 + 0.083614j, 0.463219 + 0.020897j]}
hl_real_9_10 = {9: [0.800835, 0.469297, 0.404074, 0.313429, 0.228228, 0.142781, 0.056808],
                10: [0.812899, 0.358132, 0.283289, 0.205978, 0.128911, 0.050809]}
fig, ax = plt.subplots(2, 2, figsize=(11, 8.2))
# (a) spectra
a = ax[0, 0]
for n, v in bq.items():
    a.plot([n - 0.12] * len(v), v, 'o', color=C2D, ms=6, mec='white', mew=0.8, label='2D Boussinesq' if n == 1 else None)
hlall = dict(hl); hlall.update(hl_real_9_10)
for n, v in hlall.items():
    a.plot([n + 0.12] * len(v), v, 's', color=CHL, ms=5.5, mec='white', mew=0.8, label='Hou–Luo (real)' if n == 1 else None)
for n, cs in hl_cx.items():
    for c in cs:
        a.plot(n + 0.12, c.real, 's', mfc='white', mec=CHL, mew=1.4, ms=6, label='Hou–Luo complex pair (Re μ)' if (n == 8 and c == cs[0]) else None)
a.axhline(1.0, color=MUTED, lw=1.0, ls='--'); a.text(4.2, 1.015, 'trivial mode μ = 1 (time translation)', color=MUTED, fontsize=8)
a.set_xlabel('rung n'); a.set_ylabel('unstable eigenvalue μ  (perturbation ∝ e$^{μτ}$)')
a.set_title('(a) the n-th profile has n unstable eigenvalues', fontsize=10, loc='left'); a.set_ylim(0, 1.08)
a.set_xticks(range(0, 11)); a.legend(fontsize=8, loc='upper left', bbox_to_anchor=(0.0, 0.93), frameon=False)
# (b) counts
b = ax[0, 1]
c2d = {1: 1.997, 2: 3.000, 3: 4.000, 4: 5.000, 5: 5.952, 6: 6.942, 7: 7.921}
chl = {1: 1.991, 2: 3.000, 3: 4.000, 4: 5.000, 5: 6.000, 6: 7.000, 7: 8.000, 8: 9.000, 9: 9.959, 10: 11.000}
nn = np.arange(0, 11); b.plot(nn, nn, color=MUTED, lw=1.0, ls='--', label='count = n')
b.plot(np.array(list(c2d)) - 0.12, [v - 1 for v in c2d.values()], 'o', color=C2D, ms=7, mec='white', mew=0.8, label='2D Boussinesq')
b.plot(np.array(list(chl)) + 0.12, [v - 1 for v in chl.values()], 's', color=CHL, ms=6, mec='white', mew=0.8, label='Hou–Luo')
b.set_xlabel('rung n'); b.set_ylabel('zeros of det(I − T$_μ$) in the right half-plane box − 1')
b.set_title('(b) argument-principle counts (trivial zero removed)', fontsize=10, loc='left')
b.legend(fontsize=8, loc='upper left', frameon=False); b.set_xticks(range(0, 11))
# (c) ladder
c = ax[1, 0]
def series(lam, ev, k):
    ns = [n for n in sorted(ev) if len(ev[n]) > k and n >= 2]
    return np.array([lam[n] - 1 for n in ns]), np.array([sorted(ev[n])[k] / (lam[n] - 1) for n in ns])
for (lam, ev, col, mk, name) in ((bq_lam, bq, C2D, 'o', '2D'), (hl_lam, hl, CHL, 's', 'Hou–Luo')):
    for k in range(3):
        d, y = series(lam, ev, k)
        c.plot(d, y, mk, color=col, ms=5.5, mec='white', mew=0.8, label=name if k == 0 else None)
        p = np.polyfit(d[-3:], y[-3:], 1); dd = np.linspace(0, d[-3], 20)
        c.plot(dd, np.polyval(p, dd), '-', color=col, lw=1.2, alpha=0.8)
        c.plot(0, np.polyval(p, 0), mk, mfc='white', mec=col, mew=1.4, ms=6)
for k in range(3):
    c.text(0.235, k + 0.62, f'k = {k}', color=MUTED, fontsize=8)
c.set_xlim(-0.005, 0.28); c.set_xlabel('λ$_n$ − 1'); c.set_ylabel('μ$_{n,k}$ / (λ$_n$ − 1)')
c.set_title('(c) arithmetic ladder: spacing → λ$_n$ − 1 (open: extrapolated)', fontsize=10, loc='left')
c.legend(fontsize=8, loc='upper left', frameon=False)
# (d) localization of the Hou–Luo modes, n = 8
d_ = ax[1, 1]
try:
    from hl_stability import HLStab, _lin_gl2
    St = HLStab(1.09044036, np.load('hl_cross_E32768_lam1.09044036.npy'), 32768); S = St.S
    eta = S.eta; Dh = ((1 + St.lam) + St.q) / St.eps
    sel = (eta > -3) & (eta < 2); kd = np.argmax(sel) + np.argmin(Dh[sel])
    shades = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']
    for i, mu in enumerate((0.063722, 0.160042, 0.255826, 0.443727)):
        vals, vecs = St.spectrum(mu, k=16); j = np.argmin(np.abs(vals - 1)); v = vecs[:, j]
        qp = np.empty(S.N, complex); qp[St.i0:] = v; qp[:St.i0] = v[0]
        q1 = S.st.offset(qp.real, S.c1) + 1j * S.st.offset(qp.imag, S.c1); q2 = S.st.offset(qp.real, S.c2) + 1j * S.st.offset(qp.imag, S.c2)
        th, om = _lin_gl2(S.h, St.i0, complex(mu), St.lam, eta, St.D1, St.D2, St.Tt1, St.Tt2, St.Ot1, St.Ot2, q1, q2)
        r = th / th[np.argmax(np.abs(th))]
        d_.plot(np.exp(eta), np.abs(r), color=shades[i], lw=1.6, label=f'μ = {mu:.3f}')
    ax2 = d_.twinx(); ax2.plot(np.exp(eta), Dh, color=MUTED, lw=1.0, ls=':'); ax2.set_ylim(0, 4)
    ax2.set_ylabel('D/ε on the wall (dotted)', color=MUTED); ax2.grid(False)
    d_.set_xscale('log'); d_.set_xlim(0.05, 3); d_.set_xlabel('x (along the wall)'); d_.set_ylabel("|θ'| (normalized)")
    d_.axvline(np.exp(eta[kd]), color=MUTED, lw=0.8)
    d_.set_title('(d) Hou–Luo, n = 8: unstable modes sit at the dip/front', fontsize=10, loc='left')
    d_.legend(fontsize=8, loc='upper left', frameon=False)
except Exception as e:
    d_.text(0.1, 0.5, f'panel skipped: {e}')
plt.tight_layout(); plt.savefig('fig_instability_ladder.png', dpi=160)
print('saved')
