"""Figures of the manuscript and the Supplementary Information, in Nature Portfolio format (183 mm or 89 mm wide,
sans-serif 6.5–7 pt, bold lowercase panel letters), saved as vector PDF and 600-dpi PNG.

usage:  python3 make_figures.py [names ...]     names: fig1 … fig6, sfig1 … sfig3 (default: all)
        env NCS_DATA (default ../../data/outputs) and NCS_FIGS (default ../../figures)

Every number is read from the logged outputs in data/outputs; nothing is typed in except the registered predictions
(Fig. 5, from PREDICTIONS_IPM.md) and the profile locations λ₀–λ₆ (data/outputs/ipm_rung*.out).
Colours: one categorical palette checked for colour-vision deficiency (defect #b4442a, geometry #2a78d6,
phase #1f9e6e, original tracker #8a5cc2) and one sequential ramp for grid spacing and profile depth."""
import os, re, sys, glob, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.environ.get("NCS_DATA", os.path.join(HERE, "..", "..", "data", "outputs")))
FIGS = os.path.abspath(os.environ.get("NCS_FIGS", os.path.join(HERE, "..", "..", "figures")))

MM = 1 / 25.4
W2, W1 = 183 * MM, 89 * MM
INK, MUTED, GRID = "#111111", "#55534f", "#e6e5e1"
DEF, GEO, PHA, OLD = "#b4442a", "#2a78d6", "#1f9e6e", "#8a5cc2"
HS = {0.025: "#9ecae1", 0.0125: "#4292c6", 0.00625: "#08519c"}
DEPTH = ["#a8dcc5", "#6cc3a0", "#3aa57d", "#1f8a64", "#146c4d", "#0b4f37"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "dejavusans", "font.size": 6.5, "axes.labelsize": 6.5, "xtick.labelsize": 6,
    "ytick.labelsize": 6, "legend.fontsize": 5.8, "legend.frameon": False, "legend.handlelength": 1.6,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5, "lines.linewidth": 0.9, "lines.markersize": 3, "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.4, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "axes.titlesize": 6.5})


def letter(ax, s, dx=-0.16, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"))
    plt.close(fig)
    print("saved", name)


def plain_log_ticks(axis, ticks):
    axis.set_ticks(ticks)
    axis.set_major_formatter(mt.FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_formatter(mt.NullFormatter())


# ------------------------------------------------------------------ data
def rungs():
    """λ₀–λ₆ from the rung logs (h_s 0.0125; λ₀ at 0.025; λ₆ from the fine scans, 0.15092 ± 0.00005)."""
    lam = {0: 1.0285722975, 6: 0.15092}
    for n in range(1, 6):
        txt = open(f"ipm_rung{n}_hs0125.out").read()
        lam[n] = float(re.findall(r"CROSSING .*?λ = ([0-9.]+)", txt)[-1])
    return [lam[n] for n in range(7)]


def phase_table():
    """Repaired phase (wkb_phase3): z, Re Φ₀, Im Φ₀, I; states with untracked front points excluded."""
    R = {}
    for f in ("ipm_phase3_rungs.out", "ipm_phase3_l7.out", "ipm_phase3_deep.out"):
        for line in open(f):
            m = re.search(r"z=\s*([0-9.]+) D0=([0-9.]+) D̂dip=([0-9.]+) s_cut=\s*(-?[0-9.]+) Φ=\s*([-0-9.]+)\s*"
                          r"([-+][0-9.]+)i I=([0-9.]+) front_fail=(\d+)", line)
            if m and int(m.group(8)) == 0:
                v = [float(x) for x in m.groups()]
                R[round(v[0], 3)] = v
    A = np.array([R[k] for k in sorted(R)])
    return A[:, 0], A[:, 4], A[:, 5], A[:, 6]


def z7_defect():
    txt = open("ipm_lambda7_final.out").read()
    z7 = float(re.search(r"z₇ = ([0-9.]+), slope", txt).group(1))
    sl = float(re.search(r"slope dm/dz = ([-0-9.e+]+)", txt).group(1))
    st = float(re.search(r"σ_stat = ([0-9.]+)", txt).group(1))
    return z7, sl, st


Z7_SYS = 0.08          # systematic range at the crossing (Nb, s_start, s_max variants; Table 1)


def coarse_branch():
    R = np.vstack([np.load(f) for f in sorted(glob.glob("ipm_branch_dn*.npy"))])
    R = R[np.argsort(-R[:, 0])]
    return 1 / R[:, 0], R[:, 2] - 2


def extrema(L):
    """Largest |m − 2| between consecutive profiles: coarse branch up to λ₆, production scan s20 beyond."""
    z, d = coarse_branch()
    edges = [1 / l for l in L]
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        s = (z > a) & (z < b)
        i = np.argmax(np.where(s, np.abs(d), -1))
        out.append((z[i], abs(d[i])))
    F = np.load("ipm_scan_s20.npy")
    i = np.argmax(np.abs(F[:, 3]))
    out.append((F[i, 0], abs(F[i, 3])))
    return np.array(out)


def crossing_slopes(L):
    """|dm/dz| at profiles 1–7: 1–5 from the coarse branch (local cubic), 6 from the s20 scan, 7 from the converged
    shift-averaged fit (ipm_lambda7_final.out)."""
    z, d = coarse_branch()
    out = []
    for lr in L[1:6]:
        zr = 1 / lr
        s = np.abs(z - zr) < 0.15
        c = np.polyfit(z[s] - zr, d[s], min(3, s.sum() - 1))
        out.append(abs(np.polyval(np.polyder(c), 0)))
    a = np.load("ipm_scan_s20.npy")
    a = a[np.argsort(a[:, 0])]
    out.append(abs((a[1, 3] - a[0, 3]) / (a[1, 0] - a[0, 0])))
    out.append(abs(z7_defect()[1]))
    return np.array(out)


def decay_pred(zs, Rs, Is):
    """−π dImΦ₀/dReΦ₀: differences between consecutive states a profile apart (z < 7.3), local linear fits of Im
    against Re within ±0.45 in z on the dense deep table."""
    pz, pe = [], []
    for k in range(len(zs) - 1):
        if zs[k + 1] - zs[k] > 0.5:
            pz.append(0.5 * (zs[k] + zs[k + 1]))
            pe.append(-np.pi * (Is[k + 1] - Is[k]) / (Rs[k + 1] - Rs[k]))
    for zc in np.arange(7.75, 9.05, 0.15):
        w = np.abs(zs - zc) <= 0.45
        pz.append(zc)
        pe.append(-np.pi * np.polyfit(Rs[w], Is[w], 1)[0])
    return np.array(pz), np.array(pe)


def two_solver():
    """|Δz_n| between the marching solver and the independent global solver, finest converged run per profile."""
    best = {1: ["l1"], 2: ["l2f", "l2"], 3: ["l3f", "l3"], 4: ["l4h", "l4f", "l4"], 5: ["l5h", "l5f"], 6: ["l6h", "l6f"]}
    out = []
    for n, tags in best.items():
        for tg in tags:
            if not os.path.exists(f"ipm_glob_{tg}.log"):
                continue
            txt = open(f"ipm_glob_{tg}.log").read()
            m = re.search(r"λ_∞ = ([0-9.]+)\s+\(march λ = ([0-9.]+)", txt)
            if m and "ok=False" not in txt:
                lg, lm = float(m.group(1)), float(m.group(2))
                out.append((n, abs(1 / lg - 1 / lm), abs(lg - lm)))
                break
    return np.array(out)


def deep_log():
    E = []
    for line in open("ipm_deep_e5.out"):
        m = re.search(r"z=([0-9.]+):.*D̂_dip=([0-9.]+) front s=([-0-9.]+) max∂ₓR=([0-9.]+) \|Ω_b\+∂ₓR\|/max=([0-9.e+-]+)",
                      line)
        if m:
            E.append([float(x) for x in m.groups()])
    return np.array(E)


def deepres():
    out = []
    for f in sorted(glob.glob("ipm_deepres_z*.out")):
        for line in open(f):
            m = re.search(r"z=([0-9.]+) hs=([0-9.]+): m-2=([-+0-9.e]+) D̂min=([0-9.]+).*max∂ₓR=([0-9.]+) mis=([0-9.e+-]+)"
                          r" T=([0-9.]+) ReΦ=([0-9.]+)", line)
            if m:
                out.append([float(x) for x in m.groups()])
    return np.array(out)


def travel_time():
    T = {}
    for line in open("ipm_travel_time.out"):
        p = line.split()
        if len(p) == 13 and re.match(r"^\d", p[1]):
            T[round(float(p[1]), 3)] = [float(x) for x in p[1:]]
    return np.array([T[k] for k in sorted(T)])


def rich_values():
    V = {}
    for f in glob.glob("rich7_*.out"):
        mm = re.match(r"rich7_([0-9.]+)_d([0-9.]+)_hs([0-9.]+)\.out", f)
        for line in open(f):
            m = re.search(r"m-2=([-+0-9.e]+)", line)
            if m and mm:
                V[(float(mm.group(1)), float(mm.group(2)), float(mm.group(3)))] = float(m.group(1))
    for tg, h in (("p7", 0.0125), ("h7", 0.00625)):
        for r in np.load(f"ipm_scan_{tg}.npy"):
            V.setdefault((round(r[0], 3), 0.0, h), r[3])
    return V


def variants():
    rows = []
    for line in open("ipm_phase_variants.out"):
        m = re.match(r"(.+?)\s+([-+][0-9.]+e[-+]\d+)\s+([0-9.]+)\s+([0-9.]+)\s+([-+][0-9.]+)\s+([-+][0-9.]+)", line)
        if m:
            rows.append((m.group(1).strip(), float(m.group(2)), float(m.group(3)), abs(float(m.group(5))),
                         abs(float(m.group(6)))))
    return rows


# ------------------------------------------------------------------ Fig. 1
def fig1():
    L = rungs()
    P = np.load("ipm_wall_profiles.npz", allow_pickle=True)
    fig, ax = plt.subplots(1, 3, figsize=(W2, 60 * MM), gridspec_kw=dict(width_ratios=[1, 1, 1]))
    a = ax[0]
    for k in range(6):
        a.plot(P[f"x{k}"], P[f"D{k}"], "-", color=DEPTH[k], lw=0.9, label=f"{float(P[f'z{k}']):.2f}")
    a.axhline(2.0, color=MUTED, lw=0.5, ls=":")
    a.text(0.055, 2.07, r"cut-off $\hat D=2$", fontsize=5.8, color=INK)
    a.set_xscale("log"); a.set_xlim(0.05, 3); a.set_ylim(0, 3.2)
    plain_log_ticks(a.xaxis, [0.1, 0.3, 1, 3])
    a.set_xlabel(r"distance along the wall, $x$"); a.set_ylabel(r"wall speed $\hat D = D/D_0$")
    a.legend(loc="upper left", title=r"$z=1/\lambda$", title_fontsize=5.8, ncol=2, columnspacing=0.8)
    a.annotate("stalled layer", (0.2, 1.0), xytext=(0.06, 0.45), fontsize=5.8, color=INK,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    a.annotate("dip", (0.995, 0.21), xytext=(0.45, 0.13), fontsize=5.8, color=INK,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    a.annotate("front", (1.08, 1.6), xytext=(1.45, 0.95), fontsize=5.8, color=INK,
               arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5))
    letter(a, "a")

    b = ax[1]
    x, k = P["phase_x"], P["phase_k"].real
    b.fill_between(x, 0, k, color=PHA, alpha=0.22, lw=0)
    b.plot(x, k, "-", color=PHA, lw=0.9)
    parts = None
    for line in open("ipm_phase_anatomy3.out"):
        if line.startswith("ipm_h7x_z7.3460"):
            parts = [float(v) for v in line.split()[3:8]]
    b.set_xscale("log"); b.set_xlim(np.exp(-10), 1.3); b.set_ylim(0, max(k) * 1.15)
    b.axvline(0.5, color=MUTED, lw=0.5, ls=":")
    ReP = parts[0] / float(P["phase_D0"])
    b.text(1.6e-4, max(k) * 1.0, rf"$I=\int \mathrm{{Re}}\,\kappa\,ds={parts[0]:.3f}$,  $\mathrm{{Re}}\,\Phi_0={ReP:.2f}$",
           fontsize=5.8, color=INK, va="top")
    b.text(1.6e-4, max(k) * 0.80, f"inner layer ($x<0.5$): {parts[1]:.3f}\napproach: {parts[2]:.3f}\n"
           f"dip core: {parts[3]:.3f}\nfront side: {parts[4]:.3f}", fontsize=5.8, color=INK, va="top")
    b.set_xlabel(r"$x$ (area under the curve = contribution to $I$)")
    b.set_ylabel(r"$\mathrm{Re}\,\kappa$ per unit $s=\ln x$")
    letter(b, "b")

    c = ax[2]
    z, Re, Im, I = phase_table()
    o = np.argsort(z)
    for n in range(1, 11):
        c.axhline(n * np.pi + 2.031, color=GRID, lw=0.6, zorder=0)
    c.plot(z[o], Re[o], "-", color=PHA, lw=1.0, label=r"$\mathrm{Re}\,\Phi_0(z)$")
    zr = [1 / l for l in L[1:]]
    c.plot(zr, [np.interp(v, z[o], Re[o]) for v in zr], "o", color=DEF, ms=3, label="profiles $n=1$–6 (defect zeros)")
    z7, _, _ = z7_defect()
    c.errorbar([z7], [np.interp(z7, z[o], Re[o])], xerr=[Z7_SYS], fmt="o", color=DEF, mfc="white", ms=3,
               elinewidth=0.6, capsize=1.5, label=r"$n=7$ from the defect ($\pm0.08$)")
    for n in (2, 4, 6, 8, 10):
        c.text(10.05, n * np.pi + 2.031, f"{n}π+δ", fontsize=5.2, color=MUTED, va="center")
    c.set_xlim(1.8, 10.0); c.set_ylim(3, 37)
    c.set_xlabel(r"$z=1/\lambda$"); c.set_ylabel(r"$\mathrm{Re}\,\Phi_0$")
    c.legend(loc="upper left")
    letter(c, "c")
    fig.tight_layout(w_pad=1.6)
    save(fig, "Fig1")


# ------------------------------------------------------------------ Fig. 2
def fig2():
    L = rungs()
    z, Re, Im, I = phase_table()
    o = np.argsort(z)
    zs, Rs, Is = z[o], Re[o], Im[o]
    ex = extrema(L)
    sl = crossing_slopes(L)
    fig, ax = plt.subplots(1, 3, figsize=(W2, 60 * MM))
    a = ax[0]
    zc, dc = coarse_branch()
    ok = zc < 6.0
    a.semilogy(zc[ok], np.abs(dc[ok]), "-", color=MUTED, lw=0.7, label=r"coarse branch ($h_s=0.025$)")
    for k, tg in enumerate(("s20", "h7")):
        F = np.load(f"ipm_scan_{tg}.npy"); F = F[np.argsort(F[:, 0])]
        a.semilogy(F[:, 0], np.abs(F[:, 3]), "-", color=INK, lw=0.9,
                   label=r"production ($h_s=0.0125$, $0.00625$)" if k == 0 else None)
    a.semilogy(ex[:, 0], ex[:, 1], "D", color=DEF, ms=3, label="extrema")
    for v in [1 / l for l in L[1:]] + [z7_defect()[0]]:
        a.axvline(v, color=GRID, lw=0.6, zorder=0)
    a.axhspan(1e-9, 3.4e-9, color=DEF, alpha=0.12, lw=0)
    a.text(1.08, 4.3e-9, r"error floor at $\lambda_7$", fontsize=5.6, color=INK)
    a.set_xlim(1, 7.8); a.set_ylim(1e-9, 0.5)
    a.set_xlabel(r"$z=1/\lambda$"); a.set_ylabel(r"$|m-2|$")
    a.legend(loc="lower left", bbox_to_anchor=(0.0, 0.1))
    letter(a, "a")

    b = ax[1]
    mz = 0.5 * (ex[1:, 0] + ex[:-1, 0]); meas = -np.log(ex[1:, 1] / ex[:-1, 1])
    pz, pe = decay_pred(zs, Rs, Is)
    b.plot(pz, pe, "o-", color=PHA, ms=2.2, lw=0.9, label=r"predicted, $-\pi\,d\,\mathrm{Im}\Phi_0/d\,\mathrm{Re}\Phi_0$")
    keep = mz > 2.5
    b.plot(mz[keep], meas[keep], "D", color=DEF, ms=3, label="measured between extrema")
    for xx, yy in zip(mz[keep], meas[keep]):
        b.plot([xx, xx], [np.interp(xx, pz, pe), yy], "-", color=MUTED, lw=0.5)
    b.axhline(np.log(10), color=MUTED, lw=0.5, ls=":")
    b.text(2.45, np.log(10) + 0.03, "ln 10", fontsize=5.6, color=MUTED)
    b.set_xlim(2.3, 9.3); b.set_ylim(1.5, 3.0)
    b.set_xlabel(r"$z=1/\lambda$"); b.set_ylabel("e-folds of the defect per half-period")
    b.legend(loc="upper left")
    letter(b, "b")

    c = ax[2]
    n = np.arange(1, 8)
    c.fill_between(n, 0.5e-9 / sl, 3.4e-9 / sl, color=DEF, alpha=0.14, lw=0)
    c.semilogy(n, 1e-9 / sl, "o-", color=DEF, ms=2.6, lw=0.9, label=r"defect, $\sigma_m/|dm/dz|$")
    zn = [7.999, 8.628, 9.201]
    proj = [sl[-1]]
    dprev = np.interp(7.3451, zs, np.gradient(Rs, zs))
    for zz in zn:
        dRe = np.interp(zz, zs, np.gradient(Rs, zs))
        proj.append(proj[-1] * np.exp(-np.interp(zz - 0.3, pz, pe) * (1 + 1 / zz)) * dRe / dprev)
        dprev = dRe
    c.semilogy([7, 8, 9, 10], 1e-9 / np.array(proj), "o--", color=DEF, mfc="white", ms=2.6, lw=0.7,
               label="defect, projected")
    dRe_n = np.interp(np.array([1 / l for l in L[1:]]), zs, np.gradient(Rs, zs))
    sp = np.concatenate([np.hypot(0.005, 0.0105) / dRe_n[:5], [np.hypot(0.057, 0.0105) / dRe_n[5]], [0.007]])
    c.semilogy(n, sp, "s-", color=PHA, ms=2.6, lw=0.9, label="phase")
    tw = two_solver()
    c.semilogy(tw[:, 0], tw[:, 1], "D", color=GEO, ms=3, label="two independent solvers")
    c.axhline(0.003, color=MUTED, lw=0.5, ls=":")
    c.text(10.25, 0.0042, "decisiveness\nthreshold", fontsize=5.4, color=INK, ha="right", va="bottom")
    c.axhline(0.3, color=MUTED, lw=0.5, ls=":")
    c.text(0.8, 0.42, "half a spacing", fontsize=5.4, color=INK)
    c.set_xlabel("profile $n$"); c.set_ylabel(r"uncertainty of $z_n$")
    c.set_ylim(3e-9, 30); c.set_xlim(0.6, 10.4); c.set_xticks(range(1, 11))
    c.legend(loc="lower right")
    letter(c, "c")
    fig.tight_layout(w_pad=1.6)
    save(fig, "Fig2")


# ------------------------------------------------------------------ Fig. 3
def fig3():
    fig, ax = plt.subplots(2, 2, figsize=(W2, 118 * MM))
    a = ax[0, 0]
    V = rich_values()
    spread = {}
    for (zz, d, h), v in V.items():
        spread.setdefault((zz, h), []).append(v)
    pts = {h: [] for h in HS}
    for (zz, h), vals in spread.items():
        if len(vals) >= 3 and h in HS:
            pts[h].append(0.5 * (max(vals) - min(vals)))
    hh = sorted(pts)
    for h in hh:
        a.plot([h] * len(pts[h]), pts[h], "o", color=HS[h], ms=3.5, mec=INK, mew=0.3)
    med = [np.exp(np.mean(np.log(pts[h]))) for h in hh]
    a.plot(hh, med, "-", color=MUTED, lw=0.6)
    p = np.polyfit(np.log(hh), np.log(med), 1)[0]
    a.text(0.0083, 1.5e-8, rf"$\propto h_s^{{{p:.1f}}}$", fontsize=6.5, color=INK)
    a.axhspan(1e-9, 3.4e-9, color=DEF, alpha=0.12, lw=0)
    a.text(0.0052, 4.2e-9, r"erratic floor at $\lambda_7$", fontsize=5.6, color=INK)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(0.005, 0.03); a.set_ylim(3e-11, 1e-5)
    plain_log_ticks(a.xaxis, hh)
    a.set_xlabel(r"radial grid spacing $h_s$"); a.set_ylabel("half-range of $m$ over four grid shifts")
    letter(a, "a", dx=-0.2)

    b = ax[0, 1]
    rows = variants()
    base = rows[0][1]
    items = [(l, abs(v - base)) for l, v, *_ in rows[1:] if not l.startswith("grid shift")]
    sh = [v for l, v, *_ in rows if l.startswith("grid shift")] + [base]
    items.append(("grid shift (half-range)", 0.5 * (max(sh) - min(sh))))
    items.append((r"Newton tol. $10^{-11}\to10^{-13}$", 1.7e-11))
    lab = {"Nb 48": r"$N_b$ 48", "Nb 64": r"$N_b$ 64", "s_start −22": r"$s_{\rm start}$ −22",
           "s_start −24": r"$s_{\rm start}$ −24", "s_max 130": r"$s_{\max}$ 130", "s_max 160": r"$s_{\max}$ 160",
           "Nb 48 + s_start −22 + s_max 130": "all three combined"}
    y = np.arange(len(items))[::-1]
    b.barh(y, [v for _, v in items], color=DEF, alpha=0.8, height=0.6)
    b.set_yticks(y); b.set_yticklabels([lab.get(l, l) for l, _ in items])
    b.set_xscale("log"); b.set_xlim(4e-12, 2e-8); b.set_ylim(-0.6, len(items) + 0.4)
    sl = abs(z7_defect()[1])
    for frac, t in ((0.003, r"$\Delta z_7=0.003$"), (0.3, r"$\Delta z_7=0.3$")):
        b.axvline(sl * frac, color=MUTED, lw=0.5, ls="--")
        b.text(sl * frac, len(items) - 0.35, t, fontsize=5.4, color=INK, ha="center", va="bottom")
    b.set_xlabel(r"$|\Delta(m-2)|$ at $z=7.346$"); b.grid(axis="y", visible=False)
    letter(b, "b", dx=-0.42)

    c = ax[1, 0]
    S6 = {}
    for f in ("ipm_sstart_lam6.out", "ipm_noise2_ss16.out", "ipm_noise2_ss20.out"):
        for line in open(f):
            m = re.search(r"s(?:s|_start)=(-?[0-9.]+).*?λ=0\.150920000.*m-2=([-+0-9.e]+)|s_start=(-?[0-9.]+):.*m-2=([-+0-9.e]+)",
                          line)
            if m:
                if m.group(1):
                    S6[float(m.group(1))] = float(m.group(2))
                else:
                    S6[float(m.group(3))] = float(m.group(4))
    ss = np.array(sorted(S6)); off = np.array([abs(S6[s] - S6[-20.0]) for s in ss]); k = ss > -20
    c.semilogy(ss[k], off[k], "o-", color=DEF, ms=3, lw=0.8, label=r"$\lambda_6$")
    S5 = {}
    for line in open("ipm_sstart_lam5.out"):
        m = re.search(r"s_start=(-?[0-9.]+):.*m-2=([-+0-9.e]+)", line)
        if m:
            S5[float(m.group(1))] = float(m.group(2))
    c.semilogy([-16, -24], [abs(S5[-16.0] - S5[-20.0]), abs(S5[-24.0] - S5[-20.0])], "s", color=DEF, mfc="white",
               ms=3.5, label=r"$\lambda_5$")
    deep6 = []
    for f in sorted(glob.glob("ipm_nfloor_ss*_lam6.out")):
        m = re.search(r"s_start=(-?[0-9.]+).*m-2=([-+0-9.e]+)", open(f).read())
        if m:
            deep6.append((float(m.group(1)), abs(float(m.group(2)) - S6[-20.0])))
    deep6 = np.array(sorted(deep6))
    c.semilogy(deep6[:, 0], deep6[:, 1], "o", color=DEF, mfc="white", ms=3.5, label=r"$\lambda_6$, deep start")
    sg = np.linspace(-21, -14.5, 10)
    c.semilogy(sg, off[k][-1] * np.exp(sg - ss[k][-1]), ":", color=MUTED, lw=0.6, label=r"$\propto e^{s_{\rm start}}$")
    c.text(-29.5, 4e-8, "round-off", fontsize=5.6, color=INK)
    c.set_xlabel(r"origin truncation $s_{\rm start}$"); c.set_ylabel(r"offset of $m$ from $s_{\rm start}=-20$")
    c.set_xlim(-31, -14); c.set_ylim(1e-11, 1e-6)
    c.legend(loc="upper right")
    letter(c, "c", dx=-0.2)

    d = ax[1, 1]
    NF = []
    for f in ("ipm_sstart_lam6.out", "ipm_sstart_lam5.out", "ipm_noise2_ss16.out", "ipm_noise2_ss20.out"):
        for line in open(f):
            m = re.search(r"s(?:s|_start)=(-?[0-9.]+).*?\|R\|=([0-9.e+-]+)", line)
            if m:
                NF.append((float(m.group(1)), float(m.group(2))))
    for f in sorted(glob.glob("ipm_nfloor_ss*_lam6.out")):
        m = re.search(r"s_start=(-?[0-9.]+).*\|R\|=([0-9.e+-]+)", open(f).read())
        if m:
            NF.append((float(m.group(1)), float(m.group(2))))
    NF = np.array(NF)
    d.semilogy(NF[:, 0], NF[:, 1], "o", color=INK, ms=3)
    d.set_xlabel(r"origin truncation $s_{\rm start}$"); d.set_ylabel(r"converged Newton residual $|R|$")
    d.set_xlim(-31, -14); d.set_ylim(1e-14, 1e-10)
    letter(d, "d", dx=-0.2)
    fig.tight_layout(h_pad=1.4, w_pad=1.6)
    save(fig, "Fig3")


# ------------------------------------------------------------------ Fig. 5 (phase stability)
def fig5():
    fig, ax = plt.subplots(2, 2, figsize=(W2, 118 * MM))
    a = ax[0, 0]
    L = rungs()
    z, Re, Im, I = phase_table()
    dl = np.array([np.interp(1 / L[n], z, Re) - n * np.pi for n in range(1, 7)])
    mean = dl[:5].mean()
    a.axhspan(mean - 0.005, mean + 0.005, color=PHA, alpha=0.15, lw=0)
    a.axhline(mean, color=PHA, lw=0.6)
    a.errorbar(np.arange(1, 7), dl, yerr=[0, 0, 0, 0, 0, 0.009], fmt="s", color=PHA, ms=3.5, capsize=1.5, elinewidth=0.6)
    a.text(1.0, mean + 0.008, rf"$n=1$–5: {mean:.3f} ± 0.005", fontsize=5.8, color=INK)
    a.set_xlim(0.5, 6.6); a.set_ylim(1.95, 2.06)
    a.set_xlabel("profile $n$"); a.set_ylabel(r"$\delta_n=\mathrm{Re}\,\Phi_0(\lambda_n)-n\pi$")
    letter(a, "a", dx=-0.2)

    b = ax[0, 1]
    rows = variants()
    lab = {"Nb 48": r"$N_b$ 48", "Nb 64": r"$N_b$ 64", "s_start −22": r"$s_{\rm start}$ −22",
           "s_start −24": r"$s_{\rm start}$ −24", "s_max 130": r"$s_{\max}$ 130", "s_max 160": r"$s_{\max}$ 160",
           "Nb 48 + s_start −22 + s_max 130": "combined", "grid shift 1/4": "shift ¼", "grid shift 1/2": "shift ½",
           "grid shift 3/4": "shift ¾"}
    off = {"Nb 48": (-10, 5), "Nb 64": (-4, 5), "Nb 48 + s_start −22 + s_max 130": (-14, -9), "s_start −22": (-22, 5),
           "s_start −24": (3, 4), "s_max 130": (2, -8), "s_max 160": (-6, 5), "grid shift 1/4": (3, -7),
           "grid shift 1/2": (3, 3), "grid shift 3/4": (3, 3)}
    for l, _, _, dd, dp in rows[1:]:
        dp = max(dp, 2e-5)
        b.plot(dd, dp, "o", color=HS[0.00625] if l.startswith("grid") else PHA, ms=3.5, mec=INK, mew=0.3)
        b.annotate(lab.get(l, l), (dd, dp), xytext=off.get(l, (3, 3)), textcoords="offset points", fontsize=5.4,
                   color=INK)
    b.plot([1e-4, 1], [1e-4, 1], ":", color=MUTED, lw=0.5)
    b.text(2.2e-3, 3.4e-3, "equal shift", rotation=33, fontsize=5.6, color=MUTED)
    b.axhline(2e-5, color=GRID, lw=0.6)
    b.set_xscale("log"); b.set_yscale("log"); b.set_xlim(1.2e-4, 0.4); b.set_ylim(1e-5, 0.4)
    b.set_xlabel(r"$|\Delta z_7|$ implied by the defect"); b.set_ylabel(r"$|\Delta z_7|$ implied by the phase")
    letter(b, "b", dx=-0.2)

    c = ax[1, 0]
    K = np.load("ipm_K_checks.npz")
    for tg, col, lb in (("l7", PHA, r"$z=7.346$"), ("deep", "#0b4f37", r"$z=8.426$")):
        c.plot(np.abs(K[f"{tg}_G"]), K[f"{tg}_K"].real, ".", color=col, ms=1.8, label=lb)
    gg = np.linspace(0, 0.6, 10)
    c.plot(gg, gg, ":", color=MUTED, lw=0.6)
    c.text(0.22, 0.1, r"$K=G$", fontsize=5.8, color=MUTED)
    c.set_xscale("symlog", linthresh=0.1); c.set_xlim(0, 30); c.set_ylim(0, 1.25)
    c.set_xlabel(r"wall density gradient $G=\partial_x R$"); c.set_ylabel(r"$\mathrm{Re}\,K=\mathrm{Re}\,\kappa\hat D$")
    c.legend(loc="upper left", markerscale=3)
    letter(c, "c", dx=-0.2)

    d = ax[1, 1]
    names = ["inner layer", "dip", "front side"]
    cols = ["#a8dcc5", PHA, "#0b4f37"]
    rowsK = [l for l in open("ipm_K_checks.out") if "D̂=0.012" in l]
    for line, nm, col in zip(rowsK, names, cols):
        pts = re.findall(r"D̂=([0-9.]+): ([0-9.e+-]+)", line)
        Dt = np.array([float(p[0]) for p in pts]); dv = np.array([float(p[1]) for p in pts])
        d.plot(Dt, dv, "o-", color=col, ms=3, lw=0.8, label=nm)
    d.set_xscale("log"); d.set_yscale("log"); d.set_ylim(1e-12, 1e-5)
    plain_log_ticks(d.xaxis, [0.01, 0.03, 0.1, 0.3, 1, 2])
    d.set_xlabel(r"$\hat D$ substituted in the local problem"); d.set_ylabel(r"relative change of $K=\kappa\hat D$")
    d.legend(loc="upper right")
    letter(d, "d", dx=-0.2)
    fig.tight_layout(h_pad=1.4, w_pad=1.6)
    save(fig, "Fig5")


# ------------------------------------------------------------------ Fig. 4 (holdout)
def fig4():
    z7, _, st = z7_defect()
    fig, ax = plt.subplots(figsize=(W2, 62 * MM))
    P = [("H1  shifted law, $\\lambda_2$–$\\lambda_6$", 7.3415, "c922751"),
         ("H1′ shifted law, $\\lambda_0$–$\\lambda_6$", 7.3404, "c922751"),
         ("H2  geometric contraction", 7.3343, "c922751"), ("H3  three-point law", 7.3342, "c922751"),
         ("H4  phase, coarse branch", 7.352, "273ae28")]
    ax.axvspan(z7 - 0.003, z7 + 0.003, color=GRID, alpha=0.8, lw=0)
    for k, (lb, v, h) in enumerate(P):
        yy = 4.55 - 0.45 * k
        ax.plot(v, yy, "v", color=INK, ms=3.5)
        ax.text(7.212, yy, f"{lb}   ({h})", va="center", fontsize=5.8, color=INK)
        ax.plot([v, v], [yy - 0.1, 1.85], ":", color=MUTED, lw=0.4)
    ax.text(7.212, 5.15, "registered before the computation", fontsize=6, color=MUTED, style="italic")
    yy = 1.45
    ax.errorbar(z7, yy, xerr=Z7_SYS, fmt="none", ecolor=DEF, elinewidth=0.8, capsize=2, alpha=0.7)
    ax.errorbar(z7, yy, xerr=st, fmt="o", color=DEF, ecolor=DEF, elinewidth=2.2, ms=3.5)
    ax.text(7.212, yy + 0.32, r"defect, $h_s=0.00625$, four grid shifts", fontsize=5.8, color=INK)
    ax.text(z7 + 0.084, yy, "±0.08 systematic", fontsize=5.6, color=INK, va="center")
    yy = 0.45
    ax.plot([7.3476, 7.3602], [yy, yy], "-", color=PHA, lw=3.5, solid_capstyle="butt")
    ax.text(7.212, yy + 0.32, "phase reading, repaired tracker (not registered)", fontsize=5.8, color=INK)
    ax.text(7.364, yy, r"$z_7=7.354\pm0.007$", fontsize=5.8, color=INK, va="center")
    ax.text(z7, -0.35, "±0.003\ndecisiveness\nthreshold", ha="center", fontsize=5.2, color=MUTED, va="top")
    ax.set_xlim(7.21, 7.46); ax.set_ylim(-1.15, 5.5); ax.set_yticks([]); ax.grid(axis="y", visible=False)
    ax.set_xlabel(r"$z_7=1/\lambda_7$")
    fig.tight_layout()
    save(fig, "Fig4")


# ------------------------------------------------------------------ Fig. 6
def fig6():
    D = np.load("ipm_tracker_diag.npz")
    fig, ax = plt.subplots(2, 2, figsize=(W2, 118 * MM))
    a = ax[0, 0]
    g1, K1, g3, K3 = D["z7.946_g1"], D["z7.946_K1"], D["z7.946_g3"], D["z7.946_K3"]
    s_dip, s_cut = -0.100, 0.024
    w1 = (g1 > s_dip - 0.15) & (g1 <= s_cut + 1e-9)
    w3 = (g3 > s_dip - 0.15) & (g3 <= s_cut + 1e-9)
    a.axvspan(s_dip, s_cut, color=GRID, alpha=0.7, lw=0)
    a.plot(g1[w1], K1[w1].real, "-", color=OLD, lw=1.0, label=r"original: holds $\kappa$")
    a.plot(g3[w3], K3[w3].real, "-", color=PHA, lw=1.0, label=r"repaired: continues $K=\kappa\hat D$")
    a.plot(g1[w1], -K1[w1].imag, "--", color=OLD, lw=0.7)
    a.plot(g3[w3], -K3[w3].imag, "--", color=PHA, lw=0.7)
    a.text(s_dip + 0.004, 5.0, "front side\n" + r"$\hat D$: 0.31 → 2", fontsize=5.6, color=INK, va="top")
    a.text(-0.235, 3.3, r"solid: $\mathrm{Re}\,K$" + "\n" + r"dashed: $-\mathrm{Im}\,K$", fontsize=5.6, color=MUTED)
    a.set_xlabel(r"$s=\ln x$"); a.set_ylabel(r"$K=\kappa\hat D$ at $z=7.946$"); a.set_ylim(0, 5.9)
    a.legend(loc="upper left")
    letter(a, "a", dx=-0.2)

    b = ax[0, 1]
    z, Re, Im, I = phase_table()
    old = []
    for f in ("ipm_wkb_e5_cut2.out", "ipm_wkb_fine_deep_cut2.out"):
        for line in open(f):
            m = re.search(r"λ=([0-9.]+).*I=([0-9.]+)[-+].*fail=(\d+)", line)
            if m:
                old.append((1 / float(m.group(1)), float(m.group(2)), int(m.group(3))))
    old = np.array(sorted(old))
    T = travel_time()
    zt, Tt = T[:, 0], T[:, 1]
    dd = z >= 6.5
    j = [np.argmin(abs(zt - v)) for v in z[dd]]
    cf = np.polyfit(Tt[j], I[dd], 1)
    b.plot(zt[zt >= 6.5], np.polyval(cf, Tt[zt >= 6.5]), "-", color=GEO, lw=2.4, alpha=0.3,
           label=rf"geometry: {cf[0]:.2f}$\,T$ + const")
    b.plot(old[:, 0], old[:, 1], "o-", color=OLD, ms=2.6, lw=0.7, label="original tracker")
    for xx, yy, nf in old:
        if nf:
            b.annotate(f"{int(nf)}", (xx, yy), xytext=(2, 2), textcoords="offset points", fontsize=5.2, color=OLD)
    b.plot(z[dd], I[dd], "s-", color=PHA, ms=2.2, lw=0.8, label="repaired tracker")
    b.set_xlim(6.5, 9.4)
    b.set_xlabel(r"$z=1/\lambda$"); b.set_ylabel(r"$I=D_0\,\mathrm{Re}\,\Phi_0$")
    b.legend(loc="upper left")
    letter(b, "b", dx=-0.2)

    E = deep_log()
    R2 = deepres()
    c = ax[1, 0]
    c.plot(E[:, 0], E[:, 3], "o-", color=HS[0.0125], ms=2.2, lw=0.8, label=r"$h_s=0.0125$")
    fine = R2[R2[:, 1] == 0.00625]
    c.plot(fine[:, 0], fine[:, 4], "s", color=HS[0.00625], ms=3.5, label=r"$h_s=0.00625$")
    c.set_yscale("log"); c.set_ylim(5, 30)
    plain_log_ticks(c.yaxis, [5, 7, 10, 15, 20, 30])
    c.set_xlabel(r"$z=1/\lambda$"); c.set_ylabel(r"steepest wall gradient, max $\partial_x R$")
    c.legend(loc="upper left")
    letter(c, "c", dx=-0.2)

    d = ax[1, 1]
    d.plot(E[:, 0], 100 * E[:, 4], "o-", color=HS[0.0125], ms=2.2, lw=0.8, label=r"$h_s=0.0125$")
    d.plot(fine[:, 0], 100 * fine[:, 5], "s", color=HS[0.00625], ms=3.5, label=r"$h_s=0.00625$")
    d.set_xlabel(r"$z=1/\lambda$"); d.set_ylabel(r"violation of $\Omega_b=-\partial_x R$ (%)"); d.set_ylim(0, 14)
    d.legend(loc="upper left")
    letter(d, "d", dx=-0.2)
    fig.tight_layout(h_pad=1.4, w_pad=1.6)
    save(fig, "Fig6")


# ------------------------------------------------------------------ Supplementary figures
def sfig2():
    """Geometry of the deep branch: phase anatomy, dip depth, dip and front positions."""
    fig, ax = plt.subplots(1, 3, figsize=(W2, 60 * MM))
    a = ax[0]
    rows = []
    for line in open("ipm_phase_anatomy3.out"):
        p = line.split()
        if len(p) == 13 and re.match(r"^\d", p[1]):
            v = [float(x) for x in p[1:]]
            rows.append((v[0], v[2], v[3] + v[4], v[5], v[6]))
    R = np.array(sorted(rows))
    a.stackplot(R[:, 0], R[:, 2], R[:, 3], R[:, 4], colors=["#b9e2d1", PHA, "#0b4f37"], alpha=0.95,
                labels=["before the dip core", "dip core", "front side"])
    a.plot(R[:, 0], R[:, 1], "k.", ms=2)
    a.set_xlabel(r"$z=1/\lambda$"); a.set_ylabel(r"$I=D_0\,\mathrm{Re}\,\Phi_0$"); a.set_ylim(0, 2.25)
    a.legend(loc="upper left")
    letter(a, "a")
    T = travel_time()
    E = deep_log()
    R2 = deepres()
    fine = R2[R2[:, 1] == 0.00625]
    b = ax[1]
    b.plot(T[:, 0], T[:, 6], "o-", color=INK, ms=2, lw=0.7, label="profiles along the branch")
    b.plot(E[E[:, 0] > T[:, 0].max(), 0], E[E[:, 0] > T[:, 0].max(), 1], "o-", color=MUTED, ms=2, lw=0.7,
           label="continuation beyond $z=9.5$")
    b.plot(fine[:, 0], fine[:, 3], "s", color=HS[0.00625], ms=3.5, label=r"$h_s=0.00625$")
    b.set_xlabel(r"$z=1/\lambda$"); b.set_ylabel(r"dip depth $\hat D_{\min}$"); b.set_ylim(0, 0.9)
    b.legend(loc="upper right")
    letter(b, "b")
    c = ax[2]
    c.plot(T[:, 0], T[:, 11], "o-", color=GEO, ms=2, lw=0.7, label=r"front, $\hat D=2$")
    c.plot(T[:, 0], T[:, 10], "s-", color=PHA, ms=2, lw=0.7, label=r"dip, $\hat D_{\min}$")
    c.set_xlabel(r"$z=1/\lambda$"); c.set_ylabel(r"position along the wall, $x$")
    c.legend(loc="lower right")
    letter(c, "c")
    fig.tight_layout(w_pad=1.6)
    save(fig, "SupplementaryFig2")


def sfig3():
    """Spacing of consecutive profiles: fixed stalled layer (2D Boussinesq, Hou–Luo) against receding (IPM)."""
    J = json.load(open("ladders.json"))
    fig, ax = plt.subplots(figsize=(W1, 62 * MM))
    for name, col in (("2D Boussinesq", GEO), ("Hou–Luo", "#c98a00"), ("IPM", PHA)):
        lam = np.array(J[name]["lambda"]); eps = lam - 1 if J[name]["eps"] == "lambda-1" else lam
        z = np.sort(1 / eps)
        ax.plot(np.arange(1, len(z)) + 0.5, np.diff(z), "o-", color=col, ms=2.6, lw=0.8, label=name)
    ax.axhline(1.279, color="#c98a00", lw=0.5, ls="--")
    ax.text(10.6, 1.30, r"$\pi/C=1.279$", fontsize=5.6, color=INK, ha="right")
    ax.set_xlabel(r"$n$ (between profiles $n-1$ and $n$)"); ax.set_ylabel(r"spacing of $1/\varepsilon$")
    ax.set_ylim(0.6, 1.65)
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "SupplementaryFig3")


def sfig1():
    """Raw grid-shift data: m − 2 on four shifted grids at each z, for three radial spacings."""
    V = rich_values()
    fig, ax = plt.subplots(1, 3, figsize=(W2, 58 * MM))
    for k, h in enumerate((0.025, 0.0125, 0.00625)):
        a = ax[k]
        for (zz, dl, hh), v in sorted(V.items()):
            if hh == h:
                a.plot(zz + 0.004 * (dl - 0.375), v, "o", color=HS[h], ms=3, mec=INK, mew=0.3)
        zz = sorted({key[0] for key in V if key[2] == h})
        means = [np.mean([V[key] for key in V if key[0] == q and key[2] == h]) for q in zz]
        a.plot(zz, means, "-", color=INK, lw=0.7)
        a.axhline(0, color=MUTED, lw=0.5)
        a.set_xlabel(r"$z=1/\lambda$"); a.set_ylabel(r"$m-2$")
        a.set_title(rf"$h_s={h}$", fontsize=6.5)
        a.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        letter(a, "abc"[k], dx=-0.22)
    fig.tight_layout(w_pad=1.6)
    save(fig, "SupplementaryFig1")


ALL = ["fig1", "fig2", "fig3", "fig4", "fig5", "fig6", "sfig1", "sfig2", "sfig3"]

if __name__ == "__main__":
    os.chdir(DATA)
    for name in (sys.argv[1:] or ALL):
        globals()[name]()
