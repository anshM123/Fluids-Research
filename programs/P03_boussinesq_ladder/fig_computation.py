"""Fig. 4 (computational): the exponential wall and the phase observable that bypasses it (IPM).
(a) |m − 2| along the IPM branch on the coarse (h_s 0.025) and production (h_s 0.0125, s_start −20) solvers, with
    their noise floors; the envelope falls tenfold per half-period and meets the production floor two profiles below λ₆.
(b) The error audit at λ₆: effect on m of each numerical choice (production settings marked), against the defect
    amplitude at λ₆, λ₇, λ₈.
(c) The same choices seen by the phase: relative change of Re Φ (in units of π) versus relative change of the defect."""
import numpy as np, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
IPM = '#1f9e6e'; IPM_L = '#8fcfb6'; WARN = '#b4442a'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.3), gridspec_kw=dict(width_ratios=[1.35, 1.05, 1.0]))

# (a) defect along the branch
a = ax[0]
R = np.vstack([np.load(f) for f in sorted(glob.glob("ipm_branch_dn*.npy"))]); R = R[np.argsort(-R[:, 0])]
zc, dc = 1 / R[:, 0], R[:, 2] - 2
ok = zc < 6.0
a.semilogy(zc[ok], np.abs(dc[ok]), '-', color=IPM_L, lw=1.4, label='coarse solver (h_s = 0.025)')
a.semilogy(zc[~ok], np.abs(dc[~ok]), 'x', color=IPM_L, ms=4, label='coarse, below its floor')
fine = []
for tag in ("s20", "p7"):
    try:
        F = np.load(f"ipm_scan_{tag}.npy"); fine.append(F[:, [0, 3]])
    except FileNotFoundError:
        pass
if fine:
    F = np.vstack(fine); F = F[np.argsort(F[:, 0])]
    a.semilogy(F[:, 0], np.abs(F[:, 1]), 'o-', color=IPM, ms=3, lw=1.2, label='production solver (h_s = 0.0125)')
ext_z = np.array([1.32, 2.44, 3.48, 4.44, 5.30, 6.06, 6.87]); ext = np.array([8.05e-2, 7.70e-3, 7.71e-4, 7.55e-5, 7.04e-6, 5.45e-7, 5.07e-8])
a.semilogy(ext_z, ext, 'D--', color=MUTED, ms=3, lw=0.8, label='extrema: ÷10 per half-period (÷13 at λ₆)')
nxt = [(7.6, 5e-9), (8.25, 5e-10)]
a.semilogy([ext_z[-1]] + [p[0] for p in nxt], [ext[-1]] + [p[1] for p in nxt], ':', color=MUTED, lw=0.8)
a.annotate('λ₇, λ₈ extrapolated', (8.25, 5e-10), textcoords='offset points', xytext=(-70, 14), fontsize=7, color=MUTED)
a.axhspan(1e-7, 2e-6, color=WARN, alpha=0.08); a.text(1.1, 3e-7, 'coarse noise floor', color=WARN, fontsize=8)
a.axhspan(2e-10, 6e-10, color=WARN, alpha=0.16); a.text(1.1, 9e-11, 'production noise floor', color=WARN, fontsize=8)
for zr in (1/0.4721297348, 1/0.3149618108, 1/0.2415663353, 1/0.1987224523, 1/0.1706180880, 1/0.15092):
    a.axvline(zr, color=IPM, lw=0.5, alpha=0.5)
a.set_xlabel('z = 1/λ'); a.set_ylabel('|m − 2|  (smoothness defect)'); a.set_ylim(3e-11, 0.3); a.set_xlim(1, 8.3)
a.legend(loc='upper right', fontsize=7, frameon=False)
a.set_title('a   The defect falls tenfold per profile', loc='left', fontsize=10, color=INK)

# (b) error audit at λ6 (effect on m, relative to the production setting)
b = ax[1]
items = [("origin cut-off −14", 2.06e-7), ("origin cut-off −15", 2.1e-8), ("origin cut-off −16", 5.4e-9),
         ("origin cut-off −17", 1.4e-9), ("origin cut-off −19", 3e-10), ("Nb 48 at cut-off −16", 2.0e-8),
         ("Nb 48 at cut-off −20 (λ₄–λ₆)", 1.0e-9), ("Nb 64 vs 48", 9e-11), ("h_s 0.00625 vs 0.0125", 7e-10),
         ("s_max 130 vs 100", 2e-10), ("Newton floor (scatter)", 5e-10)]
y = np.arange(len(items))[::-1]
prod = {"origin cut-off −19", "Nb 48 at cut-off −20 (λ₄–λ₆)", "Nb 64 vs 48", "h_s 0.00625 vs 0.0125", "s_max 130 vs 100", "Newton floor (scatter)"}
b.barh(y, [v for _, v in items], color=[IPM if n in prod else IPM_L for n, _ in items], height=0.62)
b.set_xscale('log'); b.set_yticks(y); b.set_yticklabels([n for n, _ in items], fontsize=7.5)
for lab, v in (("defect at λ₆", 5.4e-7), ("at λ₇", 5e-8), ("at λ₈", 5e-9)):
    b.axvline(v, color=MUTED, lw=0.8, ls='--'); b.text(v * 1.1, len(items) - 0.4, lab, fontsize=7, color=MUTED, rotation=90, va='top')
b.set_xlim(3e-11, 1e-6); b.set_xlabel('|Δm| at λ₆'); b.grid(axis='y', visible=False)
b.set_title('b   Error audit (dark: production checks)', loc='left', fontsize=10, color=INK)

# (c) the same numerical changes, as the rung shift each method implies (fraction of a spacing):
#     from the defect, Δz/spacing = |Δm| / (π A_n) — grows tenfold per profile as A_n falls;
#     from the phase,  Δz/spacing = |Δ Re Φ| / π — independent of depth (Φ needs only the O(1) profile).
cax = ax[2]
var = [("Nb 48, cut-off −16", 20.8334, 1.50e-8), ("Nb 64, cut-off −16", 20.8333, 1.49e-8), ("cut-off −16 (Nb 32)", 20.8330, 5.4e-9)]
P0 = 20.8331
for n, A, mk, al in ((6, 5.4e-7, 'o', 1.0), (7, 5e-8, 's', 0.6), (8, 5e-9, '^', 0.35)):
    xs = [abs(dm) / (np.pi * A) for _, _, dm in var]; ys = [max(abs(P - P0) / np.pi, 2e-5) for _, P, _ in var]
    cax.scatter(xs, ys, s=26, marker=mk, color=IPM, alpha=al, zorder=3, label=f'at λ{chr(0x2080 + n)} (defect {A:.0e})')
cax.text(1.3e-3, 4e-4, 'variants: Nb 48, Nb 64 (cut-off −16),\ncut-off −16 (Nb 32); reference: production', fontsize=6.5, color=INK)
cax.plot([1e-5, 10], [1e-5, 10], ':', color=MUTED, lw=0.8); cax.text(3e-3, 6e-3, 'equal impact', rotation=33, fontsize=7, color=MUTED)
cax.axvspan(0.1, 10, color=WARN, alpha=0.07); cax.text(0.13, 3e-6, 'rung lost', fontsize=7, color=WARN)
cax.set_xscale('log'); cax.set_yscale('log'); cax.set_xlim(1e-3, 3); cax.set_ylim(1e-6, 3)
cax.set_xlabel('rung shift implied by Δm (fraction of a spacing)'); cax.set_ylabel('rung shift implied by Δ Re Φ')
cax.legend(loc='upper left', fontsize=7, frameon=False)
cax.set_title('c   Same errors: defect vs phase', loc='left', fontsize=10, color=INK)
fig.tight_layout(); fig.savefig("fig_computation.png", dpi=200); print("saved fig_computation.png")
