"""Figures for the Nature Computational Science manuscript (dossier/P03_MANUSCRIPT_NCS.md).
usage: python3 fig_ncs.py [1 2 3 4 5 6]  (default: all).  Every number is read from a logged output of this folder.
Palette (validated, light mode): defect #b4442a, geometry #2a78d6, phase #1f9e6e, original tracker #8a5cc2; grids
use one sequential blue ramp (light = coarse)."""
import numpy as np, re, glob, sys, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
DEF, GEO, PHA, OLD = '#b4442a', '#2a78d6', '#1f9e6e', '#8a5cc2'
HS = {0.025: '#9ecae1', 0.0125: '#4292c6', 0.00625: '#08519c'}
plt.rcParams.update({'font.size': 8.5, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6,
                     'axes.titlesize': 9.5, 'legend.frameon': False, 'legend.fontsize': 7.5})
RUNGS = [1.0285722975, 0.4721297348, 0.3149618108, 0.2415663353, 0.1987224523, 0.1706180880, 0.15092]
Z7, Z7_STAT, Z7_SYS = 7.3451, 0.0016, 0.08


def title(ax, s):
    ax.set_title(s, loc='left', color=INK)


def phase_table():
    """Repaired phase (wkb_phase3): z, Re Φ₀, Im Φ₀, I; rows with untracked front points excluded."""
    R = {}
    for f in ("ipm_phase3_rungs.out", "ipm_phase3_l7.out", "ipm_phase3_deep.out"):
        for line in open(f):
            m = re.search(r"z=\s*([0-9.]+) D0=([0-9.]+) D̂dip=([0-9.]+) s_cut=\s*(-?[0-9.]+) Φ=\s*([-0-9.]+)\s*([-+][0-9.]+)i I=([0-9.]+) front_fail=(\d+)", line)
            if m and int(m.group(8)) == 0:
                v = [float(x) for x in m.groups()]; R[round(v[0], 3)] = v
    A = np.array([R[k] for k in sorted(R)])
    return A[:, 0], A[:, 4], A[:, 5], A[:, 6]


def extrema():
    ext_z = np.array([1.32, 2.44, 3.48, 4.44, 5.30, 6.06, 6.87])
    ext = np.array([8.05e-2, 7.70e-3, 7.71e-4, 7.55e-5, 7.04e-6, 5.45e-7, 5.07e-8])
    return ext_z, ext


def crossing_slopes():
    """|dm/dz| at rungs 1–7: rungs 1–5 from the coarse branch (local cubic), λ₆ from the s20 scan (first two points),
    λ₇ from the shift-averaged h_s 0.00625 fit (ipm_lambda7_final.out)."""
    R = np.vstack([np.load(f) for f in sorted(glob.glob("ipm_branch_dn*.npy"))]); R = R[np.argsort(-R[:, 0])]
    z, d = 1 / R[:, 0], R[:, 2] - 2
    out = []
    for lr in RUNGS[1:6]:
        zr = 1 / lr; sel = np.abs(z - zr) < 0.15
        c = np.polyfit(z[sel] - zr, d[sel], min(3, sel.sum() - 1)); out.append((zr, abs(np.polyval(np.polyder(c), 0))))
    a = np.load("ipm_scan_s20.npy"); a = a[np.argsort(a[:, 0])]
    out.append((1 / RUNGS[6], abs((a[1, 3] - a[0, 3]) / (a[1, 0] - a[0, 0]))))
    sl7 = float(re.search(r"slope dm/dz = ([-0-9.e+]+)", open("ipm_lambda7_final.out").read()).group(1))
    out.append((Z7, abs(sl7)))
    return np.array(out)


def decay_pred(zs, Rs, Is):
    """−π dImΦ₀/dReΦ₀: differences between consecutive states where they are a rung apart (z < 7.3), local linear
    fits of Im against Re within ±0.45 in z on the dense deep table."""
    pz, pe = [], []
    for k in range(len(zs) - 1):
        if zs[k + 1] - zs[k] > 0.5:
            pz.append(0.5 * (zs[k] + zs[k + 1])); pe.append(-np.pi * (Is[k + 1] - Is[k]) / (Rs[k + 1] - Rs[k]))
    for zc in np.arange(7.75, 9.05, 0.15):
        w = np.abs(zs - zc) <= 0.45
        pz.append(zc); pe.append(-np.pi * np.polyfit(Rs[w], Is[w], 1)[0])
    return np.array(pz), np.array(pe)


# ---------------------------------------------------------------- Fig. 2: the exponential wall
def fig2():
    z, Re, Im, I = phase_table()
    ez, ea = extrema(); sl = crossing_slopes()
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.9))
    a = ax[0]
    R = np.vstack([np.load(f) for f in sorted(glob.glob("ipm_branch_dn*.npy"))]); R = R[np.argsort(-R[:, 0])]
    zc, dc = 1 / R[:, 0], np.abs(R[:, 2] - 2); ok = zc < 6.0
    a.semilogy(zc[ok], dc[ok], '-', color=MUTED, lw=1.0, label='|m − 2|, coarse branch (h_s 0.025)')
    for k, tg in enumerate(("s20", "h7")):
        F = np.load(f"ipm_scan_{tg}.npy"); F = F[np.argsort(F[:, 0])]
        a.semilogy(F[:, 0], np.abs(F[:, 3]), '-', color=INK, lw=1.2, label='|m − 2|, production scans (h_s 0.0125, 0.00625)' if k == 0 else None)
    a.semilogy(ez, ea, 'D', color=DEF, ms=5, label='extrema of the defect')
    for k, (x, y) in enumerate(zip(ez, ea)):
        if k in (0, 3, 6): a.annotate(f'{y:.0e}', (x, y), xytext=(5, 4), textcoords='offset points', fontsize=7, color=INK)
    for k, lr in enumerate(RUNGS[1:], 1):
        a.axvline(1 / lr, color=GRID, lw=1.0, zorder=0)
    a.axvline(Z7, color=GRID, lw=1.0, zorder=0)
    a.text(1 / RUNGS[1], 2e-9, 'λ₁', ha='center', fontsize=7.5, color=MUTED); a.text(Z7, 2e-9, 'λ₇', ha='center', fontsize=7.5, color=MUTED)
    a.axhspan(1e-9, 3.4e-9, color=DEF, alpha=0.10); a.text(1.1, 4.5e-9, 'production error floor at λ₇ (1–3 × 10⁻⁹)', fontsize=7, color=INK)
    a.set_xlim(1, 7.8); a.set_ylim(1e-9, 0.3); a.set_xlabel('z = 1/λ'); a.set_ylabel('|m − 2|')
    a.legend(loc='upper right'); title(a, 'a   The defect falls tenfold per half-period')

    b = ax[1]
    mz = 0.5 * (ez[1:] + ez[:-1]); meas = -np.log(ea[1:] / ea[:-1])
    o = np.argsort(z); zs, Rs, Is = z[o], Re[o], Im[o]
    pz, pe = decay_pred(zs, Rs, Is)
    b.plot(pz, pe, 'o-', color=PHA, ms=3.5, lw=1.4, label='predicted: −π dImΦ₀/dReΦ₀ (repaired phase)')
    b.plot(mz[1:], meas[1:], 'D', color=DEF, ms=5, label='measured: ln(A_k / A_{k+1}), consecutive extrema')
    for x, y in zip(mz[1:], meas[1:]):
        pr = np.interp(x, pz, pe)
        b.plot([x, x], [pr, y], '-', color=MUTED, lw=0.7)
    b.axhline(np.log(10), color=MUTED, lw=0.8, ls=':'); b.text(2.45, np.log(10) + 0.03, 'ln 10', fontsize=7, color=MUTED)
    b.text(2.45, 2.85, '(measured − predicted)/predicted = (0.7–1.4) λ\nfor z ≥ 3: the size of the O(λ) corrections\nleft out of the leading-order phase', fontsize=7, color=INK, va='top')
    b.set_xlim(2.3, 9.4); b.set_ylim(1.5, 3.0); b.set_xlabel('z = 1/λ'); b.set_ylabel('e-folds of the defect per half-period')
    b.legend(loc='lower right'); title(b, 'b   The imaginary part of the same phase sets the wall')

    c = ax[2]
    n = np.arange(1, 8); slope = sl[:, 1]
    sig_m = 1e-9
    c.semilogy(n, sig_m / slope, 'o-', color=DEF, ms=5, lw=1.3, label='defect: σ_m / |dm/dz|  (σ_m = 10⁻⁹)')
    c.fill_between(n, 0.5e-9 / slope, 3.4e-9 / slope, color=DEF, alpha=0.12, lw=0)
    # projected slopes for n = 8–10 from the predicted decay per half-period at the next rungs (Fig. b), × the λ₆→λ₇ excess
    o = np.argsort(z); zs, Rs = z[o], Re[o]
    pz, pe = decay_pred(zs, Rs, Im[o])
    zn = [8.009, 8.64, 9.22]; proj = [slope[-1]]
    for zz in zn:                                     # one half-period per rung; O(λ)-corrected rate (×(1 + λ))
        proj.append(proj[-1] * np.exp(-np.interp(zz - 0.3, pz, pe) * (1 + 1 / zz)))
    c.semilogy([7, 8, 9, 10], sig_m / np.array(proj), 'o--', color=DEF, mfc='white', ms=5, lw=1.0, label='defect, projected')
    dRe_n = np.interp(1 / np.array(RUNGS[1:] + [1 / Z7]), zs, np.gradient(Rs, zs))
    dlt = np.array([2.0331, 2.0227, 2.0320, 2.0360, 2.0301, 1.9742])
    sp = np.concatenate([np.hypot(0.005, 0.0105) / dRe_n[:5], [np.hypot(0.057, 0.0105) / dRe_n[5]], [0.007]])
    c.semilogy(n, sp, 's-', color=PHA, ms=5, lw=1.3, label='phase: offset scatter ⊕ grid, / dReΦ₀/dz')
    c.axhline(0.003, color=MUTED, lw=0.8, ls=':'); c.text(10.2, 0.0021, 'registered decisiveness\nthreshold (0.003)', fontsize=6.8, color=INK, ha='right', va='top')
    c.axhline(0.3, color=MUTED, lw=0.8, ls=':'); c.text(1.0, 0.42, 'half a spacing: rung unidentifiable', fontsize=6.8, color=INK)
    c.set_xlabel('profile n'); c.set_ylabel('uncertainty of z_n'); c.set_ylim(3e-9, 30); c.set_xlim(0.7, 10.3)
    c.set_xticks(range(1, 11)); c.legend(loc='lower right'); title(c, 'c   Direct location loses ×10–20 per profile')
    fig.tight_layout(); fig.savefig("fig_ncs2.png", dpi=200); plt.close(fig); print("saved fig_ncs2.png")


# ---------------------------------------------------------------- Fig. 5: the λ7 holdout
def fig5():
    fig, ax = plt.subplots(figsize=(8.6, 3.3))
    P = [("H1  shifted law, λ₂–λ₆", 7.3415, "c922751"), ("H1′ shifted law, λ₀–λ₆", 7.3404, "c922751"),
         ("H2  geometric contraction", 7.3343, "c922751"), ("H3  three-point law", 7.3342, "c922751"),
         ("H4  phase, coarse branch", 7.352, "273ae28")]
    for k, (lab, v, h) in enumerate(P):
        y = 4.6 - 0.42 * k
        ax.plot(v, y, 'v', color=INK, ms=6); ax.text(7.212, y, f"{lab}  ({h})", va='center', fontsize=7.5, color=INK)
        ax.plot([v, v], [y, 1.9], ':', color=GRID, lw=0.8)
    ax.text(7.212, 5.15, 'registered before computation', fontsize=8, color=MUTED, style='italic')
    y = 1.55
    ax.errorbar(Z7, y, xerr=Z7_SYS, fmt='none', ecolor=DEF, elinewidth=1.2, capsize=3, alpha=0.6)
    ax.errorbar(Z7, y, xerr=Z7_STAT, fmt='o', color=DEF, ecolor=DEF, elinewidth=3, ms=5)
    ax.text(7.212, y + 0.33, 'measured from the defect (h_s = 0.00625, four grid shifts)', fontsize=7.5, color=INK)
    ax.text(Z7 + 0.083, y, '± 0.08 systematic\n(angular resolution, origin\ntruncation, far field)', fontsize=6.8, color=INK, va='center')
    y = 0.6
    ax.plot([7.3476, 7.3602], [y, y], '-', color=PHA, lw=5, solid_capstyle='butt')
    ax.text(7.212, y + 0.33, 'phase reading, repaired tracker (not registered)', fontsize=7.5, color=INK)
    ax.text(7.366, y, 'z₇ = 7.354 ± 0.007', fontsize=7.2, color=INK, va='center')
    ax.axvspan(Z7 - 0.003, Z7 + 0.003, ymin=0.0, ymax=0.97, color=GRID, alpha=0.5, lw=0)
    ax.text(Z7, -0.25, '±0.003: registered\ndecisiveness threshold', ha='center', fontsize=6.5, color=MUTED, va='top')
    ax.set_xlim(7.21, 7.46); ax.set_ylim(-0.75, 5.5); ax.set_yticks([]); ax.grid(axis='y', visible=False)
    ax.set_xlabel('z₇ = 1/λ₇')
    title(ax, 'The seventh profile: a pre-registered holdout the defect cannot decide (rule f0b3f06)')
    fig.tight_layout(); fig.savefig("fig_ncs5.png", dpi=200); plt.close(fig); print("saved fig_ncs5.png")


# ---------------------------------------------------------------- Fig. 6: failure, repair, deep branch
def fig6():
    D = np.load("ipm_tracker_diag.npz")
    fig, ax = plt.subplots(2, 2, figsize=(10.6, 7.4))
    a = ax[0, 0]
    g1, K1, g3, K3 = D["z7.946_g1"], D["z7.946_K1"], D["z7.946_g3"], D["z7.946_K3"]
    s_dip, s_cut = -0.100, 0.024
    w1, w3 = (g1 > s_dip - 0.15) & (g1 <= s_cut + 1e-9), (g3 > s_dip - 0.15) & (g3 <= s_cut + 1e-9)
    a.plot(g1[w1], K1[w1].real, '-', color=OLD, lw=1.6, label='Re K, original tracker (holds κ)')
    a.plot(g3[w3], K3[w3].real, '-', color=PHA, lw=1.6, label='Re K, repaired (continues K = κD̂)')
    a.plot(g1[w1], -K1[w1].imag, '--', color=OLD, lw=1.0, label='−Im K, original')
    a.plot(g3[w3], -K3[w3].imag, '--', color=PHA, lw=1.0, label='−Im K, repaired')
    a.axvspan(s_dip, s_cut, color=GRID, alpha=0.6, lw=0); a.text(s_dip + 0.004, 4.9, 'front side:\nD̂ 0.31 → 2', fontsize=7, color=INK, va='top')
    a.axvline(s_dip, color=MUTED, lw=0.7); a.text(s_dip - 0.004, 5.3, 'dip', fontsize=7, color=MUTED, ha='right')
    a.set_xlabel('s = ln x'); a.set_ylabel('K = κ D̂'); a.set_ylim(0, 5.9)
    a.legend(loc='upper left', bbox_to_anchor=(0.0, 0.83)); title(a, 'a   z = 7.946: the original tracker inflates K')

    b = ax[0, 1]
    z, Re, Im, I = phase_table()
    old = []
    for f in ("ipm_wkb_e5_cut2.out", "ipm_wkb_fine_deep_cut2.out"):
        for line in open(f):
            m = re.search(r"λ=([0-9.]+).*I=([0-9.]+)[-+].*fail=(\d+)", line)
            if m: old.append((1 / float(m.group(1)), float(m.group(2)), int(m.group(3))))
    old = np.array(sorted(old))
    T = {}
    for line in open("ipm_travel_time.out"):
        p = line.split()
        if len(p) == 13 and re.match(r"^\d", p[1]): T[round(float(p[1]), 3)] = float(p[2])
    zt = np.array(sorted(T)); Tt = np.array([T[k] for k in zt])
    dd = z >= 6.5
    j = [np.argmin(abs(zt - x)) for x in z[dd]]
    cfit = np.polyfit(Tt[j], I[dd], 1)
    b.plot(zt[zt >= 6.5], np.polyval(cfit, Tt[zt >= 6.5]), '-', color=GEO, lw=3.0, alpha=0.35,
           label=f'geometry: {cfit[0]:.2f} T + const, T = ∫ds/D̂')
    b.plot(old[:, 0], old[:, 1], 'o-', color=OLD, ms=4, lw=1.0, label='I, original tracker')
    for x, y, nf in old:
        if nf: b.annotate(f'{int(nf)}', (x, y), xytext=(3, 3), textcoords='offset points', fontsize=6.5, color=OLD)
    b.plot(z[dd], I[dd], 's-', color=PHA, ms=3.5, lw=1.2, label='I, repaired tracker')
    b.text(9.35, 1.60, 'numbers: points the original tracker\ncounted as failed but carried on', fontsize=6.8, color=INK, ha='right')
    b.set_xlim(6.5, 9.4); b.set_xlabel('z = 1/λ'); b.set_ylabel('I = D₀ Re Φ₀')
    b.legend(loc='upper left'); title(b, 'b   Exposed by geometry, not by residuals')

    c = ax[1, 0]
    E = []
    for line in open("ipm_deep_e5.out"):
        m = re.search(r"z=([0-9.]+):.*D̂_dip=([0-9.]+) front s=([-0-9.]+) max∂ₓR=([0-9.]+) \|Ω_b\+∂ₓR\|/max=([0-9.e+-]+)", line)
        if m: E.append([float(x) for x in m.groups()])
    E = np.array(E)
    c.plot(E[:, 0], E[:, 3], 'o-', color=INK, ms=3.5, lw=1.2, label='max ∂ₓR on the front side (h_s 0.0125)')
    c.set_xlabel('z = 1/λ'); c.set_ylabel('max ∂ₓR'); c.set_yscale('log'); c.set_ylim(5, 30)
    c.set_yticks([5, 7, 10, 15, 20, 30]); c.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter()); c.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    c.legend(loc='upper left'); title(c, 'c   The deep front steepens')
    d_ = ax[1, 1]
    d_.plot(E[:, 0], 100 * E[:, 4], 'o-', color=HS[0.0125], ms=3.5, lw=1.2, label='h_s = 0.0125 (continuation)')
    fine = []
    for f in glob.glob("deepres_z*_hs0.00625.out"):
        m = re.search(r"z=([0-9.]+).*mis=([0-9.e+-]+)", open(f).read())
        if m: fine.append((float(m.group(1)), float(m.group(2))))
    if fine:
        fine = np.array(sorted(fine)); d_.plot(fine[:, 0], 100 * fine[:, 1], 's', color=HS[0.00625], ms=6, label='h_s = 0.00625 (spot checks)')
    d_.set_xlabel('z = 1/λ'); d_.set_ylabel('|Ω_b + ∂ₓR| / max ∂ₓR  (%)'); d_.set_ylim(0, 14)
    d_.legend(loc='upper left'); title(d_, 'd   Resolution of the exact identity Ω_b = −∂ₓR')
    fig.tight_layout(); fig.savefig("fig_ncs6.png", dpi=200); plt.close(fig); print("saved fig_ncs6.png")



# ---------------------------------------------------------------- Fig. 1: the object and the observable
def fig1():
    P = np.load("ipm_wall_profiles.npz", allow_pickle=True)
    fig, ax = plt.subplots(1, 3, figsize=(13.0, 3.9))
    ramp = ['#a8dcc5', '#6cc3a0', '#3aa57d', '#1f8a64', '#146c4d', '#0b4f37']
    a = ax[0]
    for k in range(6):
        a.plot(P[f"x{k}"], P[f"D{k}"], '-', color=ramp[k], lw=1.5, label=f"z = {float(P[f'z{k}']):.2f}")
    a.axhline(2.0, color=MUTED, lw=0.8, ls=':'); a.text(0.052, 1.86, 'front cut-off D̂ = 2', fontsize=7, color=INK)
    a.set_xscale('log'); a.set_xlim(0.05, 3); a.set_ylim(0, 3.2)
    a.set_xlabel('distance along the wall, x'); a.set_ylabel('wall speed D̂ = D / D₀')
    a.legend(loc='upper left', title='profile depth', title_fontsize=7.5, ncol=2, bbox_to_anchor=(0.0, 1.0))
    a.annotate('stalled layer', (0.15, 1.0), xytext=(0.07, 0.45), fontsize=7.5, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.7))
    a.annotate('dip', (0.995, 0.2), xytext=(0.45, 0.12), fontsize=7.5, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.7))
    a.annotate('front', (1.08, 1.6), xytext=(1.5, 0.9), fontsize=7.5, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.7))
    title(a, 'a   A stalled layer that deepens and recedes')

    b = ax[1]
    x, k = P["phase_x"], P["phase_k"].real
    b.fill_between(x, 0, k, color=PHA, alpha=0.25, lw=0); b.plot(x, k, '-', color=PHA, lw=1.4)
    parts = None
    for line in open("ipm_phase_anatomy3.out"):
        if line.startswith("ipm_h7x_z7.3460"):
            parts = [float(v) for v in line.split()[3:8]]
    b.set_xscale('log'); b.set_xlim(np.exp(-10), 1.3); b.set_ylim(0, max(k) * 1.12)
    b.axvline(0.5, color=MUTED, lw=0.7, ls=':')
    if parts:
        b.text(2e-4, max(k) * 0.40, f'inner layer, x < 0.5: {parts[1]:.3f}', fontsize=7.5, color=INK)
        b.text(2e-4, max(k) * 0.72, f'x ≥ 0.5:  approach {parts[2]:.3f},  dip {parts[3]:.3f},\n               front {parts[4]:.3f}', fontsize=7.5, color=INK, va='top')
    b.text(2e-4, max(k) * 0.92, f'I = ∫ Re κ ds = {parts[0]:.4f}  →  Re Φ₀ = I / D₀ = {parts[0] / float(P["phase_D0"]):.3f}' if parts else '', fontsize=7.5, color=INK)
    b.set_xlabel('x (log scale; area = contribution to I)'); b.set_ylabel('Re κ  (per unit s = ln x)')
    title(b, 'b   The phase is a quadrature over the layer (z = 7.346)')

    c = ax[2]
    z, Re, Im, I = phase_table()
    o = np.argsort(z)
    c.plot(z[o], Re[o], '-', color=PHA, lw=1.6, label='Re Φ₀(z), repaired phase')
    for n in range(1, 11):
        c.axhline(n * np.pi + 2.031, color=GRID, lw=0.9, zorder=0)
        c.text(9.45, n * np.pi + 2.031, f'{n}π + δ', fontsize=6.5, color=MUTED, va='center')
    zr = [1 / l for l in RUNGS[1:]]
    c.plot(zr, [np.interp(x, z[o], Re[o]) for x in zr], 'o', color=DEF, ms=5, label='profiles located by the defect (n = 1–6)')
    c.errorbar([Z7], [np.interp(Z7, z[o], Re[o])], xerr=[Z7_SYS], fmt='o', color=DEF, mfc='white', ms=5, elinewidth=1.0, capsize=2,
               label='n = 7 from the defect (± 0.08)')
    c.set_xlim(1.8, 10.2); c.set_ylim(3, 37); c.set_xlabel('z = 1/λ'); c.set_ylabel('Re Φ₀')
    c.legend(loc='upper left'); title(c, 'c   One profile per half-turn: Re Φ₀ = nπ + δ')
    fig.tight_layout(); fig.savefig("fig_ncs1.png", dpi=200); plt.close(fig); print("saved fig_ncs1.png")


# ---------------------------------------------------------------- Fig. 3: anatomy of the error floor
def rich_values():
    V = {}
    for f in glob.glob("rich7_*.out"):
        mm = re.match(r"rich7_([0-9.]+)_d([0-9.]+)_hs([0-9.]+)\.out", f)
        for line in open(f):
            m = re.search(r"m-2=([-+0-9.e]+)", line)
            if m and mm: V[(float(mm.group(1)), float(mm.group(2)), float(mm.group(3)))] = float(m.group(1))
    for tg, h in (("p7", 0.0125), ("h7", 0.00625)):           # unshifted values from the scans
        a = np.load(f"ipm_scan_{tg}.npy")
        for r in a:
            V.setdefault((round(r[0], 3), 0.0, h), r[3])
    return V


def fig3():
    fig, ax = plt.subplots(2, 2, figsize=(10.6, 7.4))
    a = ax[0, 0]
    V = rich_values(); spread = {}
    for (z, d, h), v in V.items():
        spread.setdefault((z, h), []).append(v)
    pts = {h: [] for h in HS}
    for (z, h), vals in spread.items():
        if len(vals) >= 3 and h in HS: pts[h].append(0.5 * (max(vals) - min(vals)))
    hh = sorted(pts)
    for h in hh:
        a.plot([h] * len(pts[h]), pts[h], 'o', color=HS[h], ms=6, mec=INK, mew=0.4)
    med = [np.exp(np.mean(np.log(pts[h]))) for h in hh]
    a.plot(hh, med, '-', color=MUTED, lw=1.0)
    p = np.polyfit(np.log(hh), np.log(med), 1)[0]
    a.text(0.0085, 2e-8, f'∝ h_s^{p:.1f}', fontsize=8, color=INK)
    a.set_xscale('log'); a.set_yscale('log'); a.set_xlim(0.005, 0.03); a.set_ylim(3e-11, 1e-5)
    a.set_xticks(hh); a.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%g')); a.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xlabel('radial grid spacing h_s'); a.set_ylabel('half-range of m over four grid shifts')
    a.axhspan(1e-9, 3.4e-9, color=DEF, alpha=0.10); a.text(0.0052, 4.2e-9, 'erratic floor at λ₇', fontsize=7, color=INK)
    title(a, 'a   Front–grid locking (z = 7.26–7.38)')

    b = ax[0, 1]
    rows = []
    for line in open("ipm_phase_variants.out"):
        m = re.match(r"(.+?)\s+([-+][0-9.]+e[-+]\d+)\s+([0-9.]+)\s+([0-9.]+)", line)
        if m: rows.append((m.group(1).strip(), float(m.group(2))))
    base = rows[0][1]
    lab_map = {"grid shift 1/4": None, "grid shift 1/2": None, "grid shift 3/4": None}
    items = [(l, abs(v - base)) for l, v in rows[1:] if l not in lab_map]
    sh = [v for l, v in rows if l.startswith("grid shift")] + [base]
    items.append(("grid shift (half-range of 4)", 0.5 * (max(sh) - min(sh))))
    items.append(("Newton tolerance 1e-11 → 1e-13 (λ₆)", 1.7e-11))
    y = np.arange(len(items))[::-1]
    b.barh(y, [v for _, v in items], color=DEF, alpha=0.75, height=0.6)
    b.set_yticks(y); b.set_yticklabels([l for l, _ in items], fontsize=7.2); b.set_xscale('log'); b.set_xlim(5e-12, 2e-8)
    sl = 2.54e-8
    for frac, lab in ((0.003, 'moves z₇ by 0.003'), (0.3, 'moves z₇ by 0.3')):
        b.axvline(sl * frac, color=MUTED, lw=0.8, ls='--'); b.text(sl * frac * 1.08, len(items) - 0.6, lab, fontsize=6.8, color=INK, rotation=90, va='top')
    b.set_xlabel('|Δ(m − 2)| at z = 7.346 relative to production'); b.grid(axis='y', visible=False)
    title(b, 'b   Audit at λ₇ (h_s = 0.00625)')

    c = ax[1, 0]
    S6 = {}
    for f in ("ipm_sstart_lam6.out", "ipm_noise2_ss16.out", "ipm_noise2_ss20.out"):
        for line in open(f):
            m = re.search(r"s(?:s|_start)=(-?[0-9.]+).*?λ=0\.150920000.*m-2=([-+0-9.e]+)|s_start=(-?[0-9.]+):.*m-2=([-+0-9.e]+)", line)
            if m:
                if m.group(1): S6[float(m.group(1))] = float(m.group(2))
                else: S6[float(m.group(3))] = float(m.group(4))
    ss = np.array(sorted(S6)); off = np.array([abs(S6[s] - S6[-20.0]) for s in ss])
    k = ss > -20
    c.semilogy(ss[k], off[k], 'o-', color=DEF, ms=5, lw=1.2, label='λ₆: |m(s_start) − m(−20)|')
    S5 = {}
    for line in open("ipm_sstart_lam5.out"):
        m = re.search(r"s_start=(-?[0-9.]+):.*m-2=([-+0-9.e]+)", line)
        if m: S5[float(m.group(1))] = float(m.group(2))
    c.semilogy([-16, -24], [abs(S5[-16.0] - S5[-20.0]), abs(S5[-24.0] - S5[-20.0])], 's', color=DEF, mfc='white', ms=6,
               label='λ₅: |m(s_start) − m(−20)|')
    sg = np.linspace(-24.5, -14.5, 10); c.semilogy(sg, off[0] * np.exp(sg - ss[0]), ':', color=MUTED, lw=0.9, label='∝ e^{s_start} (neglected local terms)')
    c.annotate('round-off offset,\nλ-independent', (-24, abs(S5[-24.0] - S5[-20.0])), xytext=(-23.6, 3e-8), fontsize=7, color=INK,
               arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.7))
    c.set_xlabel('origin truncation s_start'); c.set_ylabel('offset in m'); c.set_xlim(-25, -14); c.set_ylim(1e-11, 1e-6)
    c.legend(loc='upper left'); title(c, 'c   Origin truncation: deeper is not better')

    d = ax[1, 1]
    NF = []
    for f in ("ipm_sstart_lam6.out", "ipm_sstart_lam5.out", "ipm_noise2_ss16.out", "ipm_noise2_ss20.out"):
        for line in open(f):
            m = re.search(r"s(?:s|_start)=(-?[0-9.]+).*?\|R\|=([0-9.e+-]+)", line)
            if m: NF.append((float(m.group(1)), float(m.group(2))))
    NF.append((-30.0, 1e-11))
    NF = np.array(NF)
    d.semilogy(NF[:-1, 0], NF[:-1, 1], 'o', color=INK, ms=5, label='converged Newton residual')
    d.semilogy([-30], [1e-11], 'o', color=INK, mfc='white', ms=6, label='s_start = −30 (production log)')
    d.set_xlabel('origin truncation s_start'); d.set_ylabel('Newton residual floor |R|'); d.set_xlim(-31, -14); d.set_ylim(1e-14, 1e-10)
    d.legend(loc='upper right'); title(d, 'd   The Newton floor rises with depth of the start')
    fig.tight_layout(); fig.savefig("fig_ncs3.png", dpi=200); plt.close(fig); print("saved fig_ncs3.png")


# ---------------------------------------------------------------- Fig. 4: the phase is stable where the defect is not
def fig4():
    fig, ax = plt.subplots(2, 2, figsize=(10.6, 7.4))
    a = ax[0, 0]
    n = np.arange(1, 7); dl = np.array([2.0331, 2.0227, 2.0320, 2.0360, 2.0301, 1.9742])
    a.axhspan(2.031 - 0.005, 2.031 + 0.005, color=PHA, alpha=0.15, lw=0); a.axhline(2.031, color=PHA, lw=0.8)
    a.errorbar(n, dl, yerr=[0, 0, 0, 0, 0, 0.009], fmt='s', color=PHA, ms=6, capsize=3)
    a.text(1.0, 2.039, 'n = 1–5: 2.031 ± 0.005', fontsize=7.5, color=INK)
    a.annotate('n = 6: −0.057\n(λ₆ ± 0.00005 → ± 0.009)', (6, 1.9742), xytext=(4.3, 1.985), fontsize=7.2, color=INK,
               arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.7))
    a.set_xlim(0.5, 6.6); a.set_ylim(1.95, 2.06); a.set_xlabel('profile n'); a.set_ylabel('δ_n = Re Φ₀(λ_n) − nπ')
    title(a, 'a   The offset is constant to 0.2 % of a spacing (n ≤ 5)')

    b = ax[0, 1]
    rows = []
    for line in open("ipm_phase_variants.out"):
        m = re.match(r"(.+?)\s+([-+][0-9.]+e[-+]\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([-+][0-9.]+)\s+([-+][0-9.]+)", line)
        if m: rows.append((m.group(1).strip(), abs(float(m.group(5))), abs(float(m.group(6)))))
    OFF = {"Nb 48": (-14, 8), "Nb 64": (-6, 8), "Nb 48 + s_start −22 + s_max 130": (-60, -14), "s_start −22": (-30, 8),
           "s_start −24": (4, 6), "s_max 130": (2, -11), "s_max 160": (4, 6), "grid shift 1/4": (4, -9),
           "grid shift 1/2": (4, 4), "grid shift 3/4": (4, 4)}
    for lab, dd, dp in rows[1:]:
        dp = max(dp, 2e-5)
        b.plot(dd, dp, 'o', color=PHA if not lab.startswith('grid') else HS[0.00625], ms=6, mec=INK, mew=0.4)
        b.annotate(lab, (dd, dp), xytext=OFF.get(lab, (4, 3)), textcoords='offset points', fontsize=6.5, color=INK)
    b.plot([1e-4, 1], [1e-4, 1], ':', color=MUTED, lw=0.8); b.text(1.5e-3, 2.2e-3, 'equal shift', rotation=38, fontsize=7, color=MUTED)
    b.axvspan(0.003, 1, color=DEF, alpha=0.06, lw=0); b.axhline(2e-5, color=GRID, lw=0.8)
    b.text(1.4e-4, 2.4e-5, 'below 2 × 10⁻⁵ (plotted at the line)', fontsize=6.5, color=MUTED)
    b.set_xscale('log'); b.set_yscale('log'); b.set_xlim(1.2e-4, 0.4); b.set_ylim(1e-5, 0.4)
    b.set_xlabel('|Δz₇| implied by the defect, Δm / |dm/dz|'); b.set_ylabel('|Δz₇| implied by the phase, ΔRe Φ₀ / (dRe Φ₀/dz)')
    title(b, 'b   Same changes at λ₇: phase moves 10–1000× less')

    c = ax[1, 0]
    K = np.load("ipm_K_checks.npz")
    for tg, col, lab in (("l7", PHA, 'z = 7.346'), ("deep", '#0b4f37', 'z = 8.426')):
        G, KK = K[f"{tg}_G"], K[f"{tg}_K"]
        c.plot(np.abs(G), KK.real, '.', color=col, ms=3, label=f'Re K along the wall, {lab}')
    gg = np.linspace(0, 0.6, 10); c.plot(gg, gg, ':', color=MUTED, lw=0.9); c.text(0.25, 0.12, 'K = G', fontsize=7, color=MUTED)
    c.set_xscale('symlog', linthresh=0.1); c.set_xlim(0, 30); c.set_ylim(0, 1.3)
    c.set_xlabel('wall density gradient G = ∂ₓR'); c.set_ylabel('Re K = Re κ D̂')
    c.legend(loc='upper left'); title(c, 'c   K follows G, then saturates: Φ₀ is a travel time')

    d = ax[1, 1]
    inv = K["inv"]
    labs = {}
    for r in inv:
        labs.setdefault(int(r[0]), []).append(r)
    names = ["inner layer (x = 0.5)", "dip", "front side (D̂ ≈ 1)"]
    for (i, rr), nm, col in zip(sorted(labs.items()), names, ['#a8dcc5', PHA, '#0b4f37']):
        rr = np.array(rr); Kr = rr[:, 2] + 1j * rr[:, 3]; K0 = Kr[np.argmin(abs(rr[:, 1] - 1.0))]
        d.plot(rr[:, 1], np.maximum(np.abs(Kr - K0) / abs(K0), 1e-12), 'o-', color=col, ms=5, lw=1.0, label=nm)
    d.set_xscale('log'); d.set_yscale('log'); d.set_ylim(1e-12, 1e-2)
    d.set_xlabel('D̂ substituted in the local problem (ĉ, μ, G, R_y fixed)'); d.set_ylabel('|K(D̂) − K(1)| / |K(1)|')
    d.legend(loc='upper left'); title(d, 'd   Exact scaling κ = K/D̂ holds to solver precision')
    fig.tight_layout(); fig.savefig("fig_ncs4.png", dpi=200); plt.close(fig); print("saved fig_ncs4.png")


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6"]
    for w in which:
        globals()[f"fig{w}"]()
