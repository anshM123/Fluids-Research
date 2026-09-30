"""Fig. 3: the same phase quantizes the instabilities.
(a) The unstable spectrum of the n-th smooth profile, ν̂ = μ/ε (ε = λ−1 for 2D Boussinesq and Hou–Luo, λ for
    IPM), for all three models: n real members on a lattice of step ≈ 1, i.e. index n.
(b) Spectral flow along the continuous Hou–Luo branch: the lattice moves up by one step per rung interval (lines:
    eigen-condition Re Φ − πμ/ε + θ₀ = (k + ½)π with the phase linear between rungs), a new member entering at μ = 0.
(c) The IPM blind test: predicted (open) versus computed (filled) spectra at U₃ (stage 2a) and at three branch
    points between U₁ and U₂ (P5).
(d) The step of the lattice (gap between the two lowest members) tends to ε: extrapolations to ε → 0."""
import numpy as np, re, glob, os, contextlib, io
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    from ladder_fit import bq_lam, bq, hl_lam, hl

INK, MUTED, GRID = '#0b0b0b', '#52514e', '#e4e3df'
COL = {'2D Boussinesq': '#2a78d6', 'Hou–Luo': '#eb6834', 'IPM': '#1f9e6e'}
MK = {'2D Boussinesq': 'o', 'Hou–Luo': 's', 'IPM': 'D'}
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.6})


def ipm_result(fn):
    """(λ, non-trivial ν̂ list) from an ipm_flow_spec output (the trivial root ν̂ = 1/λ is dropped)"""
    if not os.path.exists(fn):
        return None
    for line in open(fn):
        m = re.search(r"IPM λ=([0-9.]+) .*ν̂ = μ/λ = \[([^\]]*)\]", line)
        if m:
            lam = float(m.group(1)); v = [float(x) for x in m.group(2).split(',') if x.strip()]
            return lam, [x for x in v if abs(x - 1 / lam) > 0.02]
    return None


# rung spectra (ν̂) per model
spec = {'2D Boussinesq': {n: sorted(np.array(bq[n]) / (bq_lam[n] - 1)) for n in bq},
        'Hou–Luo': {n: sorted(np.array(hl[n]) / (hl_lam[n] - 1)) for n in hl if n <= 7}}
ipm = {}
for n in range(0, 8):
    r = ipm_result(f"ipm_spec_U{n}.out")
    if r:
        ipm[n] = (r[0], sorted(r[1]))
spec['IPM'] = {n: v for n, (l, v) in ipm.items() if n >= 1}

fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.6))
a = ax[0, 0]
off = {'2D Boussinesq': -0.18, 'Hou–Luo': 0.0, 'IPM': 0.18}
for name, sp in spec.items():
    for n, v in sp.items():
        a.plot([n + off[name]] * len(v), v, MK[name], color=COL[name], ms=5.5, mec='white', mew=0.6,
               label=name if n == min(sp) else None)
for j in range(0, 9):
    a.axhline(0.75 + j, color=GRID, lw=0.8, zorder=0)
a.set_xlabel('n (smooth profile)'); a.set_ylabel('unstable growth rates  μ/ε')
a.set_title('(a) the n-th profile has n unstable modes, on a lattice of step ε', loc='left', fontsize=10)
a.legend(fontsize=8, frameon=False, loc='upper left'); a.set_xlim(0.4, 7.6); a.set_ylim(0, 9)

# (b) Hou–Luo spectral flow (rungs and branch points), as in fig_spectral_flow.py
b = ax[0, 1]


def parse(files, pat=r"λ=([0-9.]+) z=([0-9.]+).*?ν̂ = μ/\(λ−1\) = \[([^\]]*)\]"):
    out = []
    for fn in files:
        if os.path.exists(fn):
            for line in open(fn):
                m = re.search(pat, line)
                if m and m.group(3).strip():
                    out.append((float(m.group(1)), float(m.group(2)), [float(x) for x in m.group(3).split(',')]))
    return out


between = {round(t[0], 6): t for t in parse(["hl_saw_exact.out", "hl_deep_exact.out", "hl_G_check.out"]) if t[1] < 16.5}
hl_all = list(hl_lam) + [1.08113374, 1.07355524]
zr = [1 / (l - 1) for l in hl_all]
for n in range(3, 9):
    v = spec['Hou–Luo'].get(n) or sorted(np.array(hl[n]) / (hl_lam[n] - 1))
    b.plot([zr[n]] * len(v), v, 's', color=COL['Hou–Luo'], ms=6, mec='white', mew=0.8,
           label='smooth profiles' if n == 3 else None)
first = True
for lam_, z_, nus in between.values():
    b.plot([z_] * len(nus), nus, 's', mfc='white', mec=COL['Hou–Luo'], mew=1.2, ms=5.5,
           label='branch points between them' if first else None); first = False
for n in range(3, len(zr)):
    b.axvline(zr[n], color=MUTED, lw=0.6, ls=':')
b.set_xlim(zr[3] - 0.4, 16.4); b.set_ylim(0, 8.5)
b.set_xlabel('z = 1/(λ−1) along the Hou–Luo branch'); b.set_ylabel('μ/ε')
b.set_title('(b) the lattice moves with the phase: one mode enters per rung', loc='left', fontsize=10)
b.legend(fontsize=8, frameon=False, loc='upper left')

# (c) IPM blind test
c = ax[1, 0]
pred = {}
for fnp in ("ipm_stage2a_predictions.out", "ipm_stage2b_predictions.out"):
  if not os.path.exists(fnp):
    continue
  for line in open(fnp):
    m = re.search(r"PREDICTION ν̂ = μ/λ of U_(\d): \[([^\]]*)\]", line)
    if m:
        pred[('U', int(m.group(1)))] = [float(x) for x in m.group(2).split(',')]
    m = re.search(r"f = ([0-9.]+) \(z = ([0-9.]+), λ = ([0-9.]+)\): predicted ν̂ = \[([^\]]*)\]", line)
    if m and (fnp.endswith("2a_predictions.out") or float(m.group(1)) == 0.25):
        pred[('bp', fnp[-17:-16], float(m.group(1)))] = (float(m.group(2)), [float(x) for x in m.group(4).split(',') if x.strip()])
for n, (lam_, v) in ipm.items():
    if n >= 1:
        c.plot([1 / lam_] * len(v), v, 'D', color=COL['IPM'], ms=6, mec='white', mew=0.8,
               label='computed (rungs)' if n == 1 else None)
for key, val in pred.items():
    if key[0] == 'U':
        zz = [1 / l for n, (l, v) in ipm.items() if n == key[1]]
        zz = zz[0] if zz else None
        if zz:
            c.plot([zz] * len(val), val, 'D', mfc='none', mec=INK, mew=1.0, ms=10,
                   label='predicted before computing (rungs)' if key[1] == 3 else None)
    else:
        zf, v = val
        c.plot([zf] * len(v), v, 'o', mfc='none', mec=INK, mew=1.0, ms=9, label='predicted (branch points)' if key[2] == 0.5 else None)
for fr, tag in ((0.25, 'bp025'), (0.5, 'bp050'), (0.75, 'bp075'), (0.25, '23b025')):
    r = ipm_result(f"ipm_spec_{tag}.out")
    if r:
        c.plot([1 / r[0]] * len(r[1]), r[1], 'o', color=COL['IPM'], ms=5, mec='white', mew=0.6,
               label='computed (branch points)' if fr == 0.5 else None)
c.set_xlabel('1/λ along the IPM branch'); c.set_ylabel('μ/λ')
c.set_title('(c) IPM: predictions committed before the computation', loc='left', fontsize=10)
c.legend(fontsize=7.5, frameon=False, loc='upper left')

# (d) lattice step
d = ax[1, 1]
for name, sp in spec.items():
    pts = [(1 / (hl_lam[n] - 1) if name == 'Hou–Luo' else 1 / (bq_lam[n] - 1) if name == '2D Boussinesq' else 1 / ipm[n][0],
            v[1] - v[0]) for n, v in sp.items() if len(v) >= 3]
    if pts:
        x, y = np.array(pts).T
        d.plot(1 / x, y, MK[name], color=COL[name], ms=6, mec='white', mew=0.8, label=name)
        if len(x) >= 3:
            p = np.polyfit(1 / x, y, 1); xx = np.linspace(0, (1 / x).max(), 20)
            d.plot(xx, np.polyval(p, xx), '-', color=COL[name], lw=0.9, alpha=0.7)
d.axhline(1.0, color=MUTED, lw=1.0, ls='--')
d.set_xlim(-0.005, None); d.set_xlabel('ε (= λ − 1; IPM: λ)'); d.set_ylabel('gap between the two lowest members / ε')
d.set_title('(d) the step of the lattice → ε (dashed: prediction)', loc='left', fontsize=10)
d.legend(fontsize=8, frameon=False, loc='upper left')
plt.tight_layout(); plt.savefig('fig_spectrum_main.png', dpi=160)
print('saved fig_spectrum_main.png')
