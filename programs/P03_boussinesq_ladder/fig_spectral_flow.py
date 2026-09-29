"""Fig. 3: one phase for both ladders.
(a) Hou–Luo: lower unstable eigenvalues ν̂ = μ/(λ−1) along the continuous branch (filled: smooth rungs; open: branch
    points between them), with the lines predicted by the eigen-condition Re Φ(λ) − πμ/(λ−1) + θ₀ = (k + ½)π, i.e.
    ν̂ = s(z)·(c̄ + Ψ(z) − j), where Ψ(z) = n + f is the number of half-turns of the profile phase (linear between
    rungs), c̄ the rung offset and s(z) the local spacing.
(b) the same for 2D Boussinesq (when bq_flow.out exists).
(c) observed offset of the ladder against the prediction c_n + f (both models).
(d) local spacing s against 1/z: the prediction is s → 1."""
import numpy as np, re, glob, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    from ladder_fit import bq_lam, bq, hl_lam, hl

C2D, CHL, GRID, INK, MUTED = '#2a78d6', '#eb6834', '#e4e3df', '#0b0b0b', '#52514e'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})

# ---------------- rung data ----------------
hl_lam_all = list(hl_lam) + [1.08113374, 1.07355524]
hl_ev = dict(hl); hl_ev[9] = [0.056808, 0.142781, 0.228228, 0.313429, 0.404074, 0.469297]
hl_ev[10] = [0.050809, 0.128911, 0.205978, 0.283289, 0.358132]
zr_hl = {n: 1 / (l - 1) for n, l in enumerate(hl_lam_all)}
for n in range(11, 20):                                  # extrapolated rung positions (spacing growth ≈ 0.0015/rung)
    zr_hl[n] = zr_hl[n - 1] + (zr_hl[n - 1] - zr_hl[n - 2]) + 0.0015
zr_bq = {n: 1 / (l - 1) for n, l in enumerate(bq_lam)}


def rung_nuhat(lam, ev, n, kmax=4):
    v = sorted(x for x in ev[n] if x < 0.5)[:kmax]
    return [x / (lam[n] - 1) for x in v]


def ladder_params(nus):
    """spacing from the two lowest members and offset of the lowest member in spacing units"""
    s = nus[1] - nus[0]
    return s, nus[0] / s


def phase_count(z, zr):
    ns = sorted(zr)
    for n in ns[:-1]:
        if zr[n] <= z < zr[n + 1]:
            return n + (z - zr[n]) / (zr[n + 1] - zr[n]), n
    return np.nan, None


def parse(files, pat=r"λ=([0-9.]+) z=([0-9.]+).*?ν̂ = μ/\(λ−1\) = \[([^\]]*)\]"):
    out = []
    for fn in files:
        if not os.path.exists(fn):
            continue
        for line in open(fn):
            m = re.search(pat, line)
            if m and m.group(3).strip():
                out.append((float(m.group(1)), float(m.group(2)), [float(x) for x in m.group(3).split(',')]))
    return out


# states re-converged at their exact λ (reconverge.py); the earlier runs on file-name-rounded λ are superseded
hl_between = parse(["hl_saw_exact.out", "hl_deep_exact.out"])
bq_between = parse(["bq_flow.out"])
# de-duplicate (the k = 28 reruns supersede the k = 16 runs of the same state)
seen = {}
for lam, z, nus in hl_between:
    seen[round(lam, 6)] = (lam, z, nus)
hl_between = sorted(seen.values(), key=lambda t: t[1])


def classify(z, nus, zr, cbar):
    """split the roots into ladder members (near s(j + c̄ + f)) and off-ladder roots"""
    Psi, n = phase_count(z, zr)
    if n is None or len(nus) < 2:
        return None
    pred = cbar + Psi - n                      # offset predicted for the lowest member present at the rung n
    s = np.median(np.diff(sorted(nus)[-3:])) if len(nus) >= 3 else np.diff(sorted(nus))[-1]
    lad, off = [], []
    for v in sorted(nus):
        r = v / s - pred
        (lad if abs(r - round(r)) < 0.15 else off).append(v)
    return Psi, n, pred, s, lad, off


fig, ax = plt.subplots(2, 2, figsize=(11, 8.4))
rows_off = []                                      # (model, observed offset mod 1, predicted mod 1)
spacing_pts = []                                   # (model, 1/z, s)
for model, ax_, lam, ev, zr, between, col, mk, zlim in (
        ("Hou–Luo", ax[0, 0], hl_lam_all, hl_ev, zr_hl, hl_between, CHL, 's', (9.4, 16.4)),
        ("2D Boussinesq", ax[0, 1], bq_lam, bq, zr_bq, bq_between, C2D, 'o', (5.2, 11.6))):
    a = ax_
    ns = [n for n in sorted(ev) if n >= 3 and zlim[0] - 1.5 < zr[n] < zlim[1] + 1.5]
    cs = []
    for n in ns:
        nus = rung_nuhat(lam, ev, n)
        s, c = ladder_params(nus); cs.append(c)
        spacing_pts.append((model, 1 / zr[n], s))
        a.plot([zr[n]] * len(nus), nus, mk, color=col, ms=6, mec='white', mew=0.8,
               label='smooth profiles (rungs)' if n == ns[0] else None)
    cbar = float(np.mean(cs[-4:]))
    # predicted lines
    zz = np.linspace(zlim[0], zlim[1], 800)
    sfit = np.polyfit([p[1] for p in spacing_pts if p[0] == model], [p[2] for p in spacing_pts if p[0] == model], 1)
    for j in range(-2, 30):
        yy = []
        for zv in zz:
            Psi, n0 = phase_count(zv, zr)
            yy.append(np.polyval(sfit, 1 / zv) * (cbar + Psi - j) if n0 is not None else np.nan)
        yy = np.array(yy); yy[yy < 0] = np.nan
        a.plot(zz, yy, '-', color=col, lw=0.9, alpha=0.45, label='eigen-condition (one phase)' if j == 5 else None)
    first = True
    for lamb, z, nus in between:
        if not (zlim[0] < z < zlim[1]):
            continue
        res = classify(z, nus, zr, cbar)
        if res is None:
            continue
        Psi, n, pred, s, lad, off = res
        a.plot([z] * len(lad), lad, mk, mfc='white', mec=col, mew=1.3, ms=6,
               label='branch points between rungs' if first else None)
        if off:
            a.plot([z] * len(off), off, 'x', color=MUTED, ms=6, mew=1.2, label='off-ladder root' if first else None)
        first = False
        spacing_pts.append((model, 1 / z, s))
        o = lad[0] / s if lad else np.nan
        rows_off.append((model, o % 1, pred % 1, z))
    for n in sorted(zr):
        if zlim[0] < zr[n] < zlim[1]:
            a.axvline(zr[n], color=MUTED, lw=0.7, ls=':')
            a.text(zr[n], 4.05, f'n={n}', color=MUTED, fontsize=7.5, ha='center')
    a.set_xlim(*zlim); a.set_ylim(0, 4.2)
    a.set_xlabel('z = 1/(λ − 1)  along the continuous branch'); a.set_ylabel('μ / (λ − 1)')
    a.set_title(f"({'a' if model == 'Hou–Luo' else 'b'}) {model}: spectral flow along the branch", fontsize=10, loc='left')
    a.legend(fontsize=7.5, loc='lower right', frameon=True, framealpha=0.9)
    if model == "2D Boussinesq" and not between:
        a.text(0.5, 0.5, 'between-rung spectra running', transform=a.transAxes, ha='center', color=MUTED)
# (c) offsets
c = ax[1, 0]
c.plot([0, 1], [0, 1], color=MUTED, lw=1.0, ls='--', label='prediction')
for model, col, mk in (("Hou–Luo", CHL, 's'), ("2D Boussinesq", C2D, 'o')):
    pts = [(p, o) for (m, o, p, z) in rows_off if m == model and np.isfinite(o)]
    if pts:
        p, o = np.array(pts).T
        c.plot(p, o, mk, color=col, ms=6, mec='white', mew=0.8, label=model)
c.set_xlim(0, 1); c.set_ylim(0, 1)
c.set_xlabel('predicted offset  c$_n$ + f  (mod 1)'); c.set_ylabel('observed offset of the lowest member (mod 1)')
c.set_title('(c) the offset is locked to the profile phase', fontsize=10, loc='left'); c.legend(fontsize=8, frameon=False)
# (d) spacing
d = ax[1, 1]
for model, col, mk in (("Hou–Luo", CHL, 's'), ("2D Boussinesq", C2D, 'o')):
    pts = np.array([(x, s) for (m, x, s) in spacing_pts if m == model])
    if len(pts):
        d.plot(pts[:, 0], pts[:, 1], mk, color=col, ms=6, mec='white', mew=0.8, label=model)
        p = np.polyfit(pts[:, 0], pts[:, 1], 1); xx = np.linspace(0, pts[:, 0].max(), 20)
        d.plot(xx, np.polyval(p, xx), '-', color=col, lw=1.0, alpha=0.7)
        d.text(0.004, np.polyval(p, 0) + (0.012 if model == "Hou–Luo" else -0.02), f'{np.polyval(p, 0):.3f}', color=col, fontsize=8)
d.axhline(1.0, color=MUTED, lw=1.0, ls='--')
d.set_xlim(0, None); d.set_xlabel('1/z = λ − 1'); d.set_ylabel('local spacing of the ladder / (λ − 1)')
d.set_title('(d) spacing → λ − 1 (dashed: prediction)', fontsize=10, loc='left'); d.legend(fontsize=8, frameon=False)
plt.tight_layout(); plt.savefig('fig_spectral_flow.png', dpi=160)
for r in rows_off:
    print(f"{r[0]:14s} z={r[3]:7.3f} offset obs {r[1]:.3f} pred {r[2]:.3f}")
print('saved')
