"""Supplementary tables, generated from the logged outputs so that no number is transcribed by hand.
usage: python3 make_si_tables.py      (reads data/outputs, writes supplementary/tables/*.tex)"""
import os, re, sys, glob
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "figures"))
import make_figures as M                                   # parsers shared with the figures

DATA = M.DATA
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "supplementary", "tables"))
os.makedirs(OUT, exist_ok=True)
os.chdir(DATA)


def sci(v, d=1):
    if v == 0:
        return "0"
    e = int(np.floor(np.log10(abs(v))))
    m = v / 10 ** e
    return f"${m:.{d}f}\\times10^{{{e}}}$"


def write(name, header, rows, colspec, note=None):
    with open(os.path.join(OUT, name), "w") as f:
        f.write("\\begin{tabular}{" + colspec + "}\n\\toprule\n" + header + "\\\\\n\\midrule\n")
        for r in rows:
            f.write(" & ".join(r) + "\\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
        if note:
            f.write("\n\\smallskip\\par\\footnotesize " + note + "\n")
    print("wrote", name)


# ---- resonance phase along the branch
rows = []
for fn in ("ipm_phase3_rungs.out", "ipm_phase3_l7.out", "ipm_phase3_deep.out"):
    for line in open(fn):
        m = re.search(r"z=\s*([0-9.]+) D0=([0-9.]+) D̂dip=([0-9.]+) s_cut=\s*(-?[0-9.]+) Φ=\s*([-0-9.]+)\s*"
                      r"([-+][0-9.]+)i I=([0-9.]+) front_fail=(\d+)", line)
        if m:
            v = m.groups()
            rows.append((float(v[0]), [f"{float(v[0]):.4f}", v[1], v[2], f"{float(v[3]):+.3f}", v[4], v[5], v[6],
                                       v[7] + (" (excluded)" if int(v[7]) else "")]))
rows = [r for _, r in sorted(rows)]
write("phase_table.tex", "$z$ & $D_0$ & $\\hat D_{\\rm dip}$ & $s_{\\rm cut}$ & $\\mathrm{Re}\\,\\Phi_0$ & "
      "$\\mathrm{Im}\\,\\Phi_0$ & $I$ & untracked", rows, "rrrrrrrl",
      "Rungs $n=1$--5 and $\\lambda_6$ on $h_s=0.0125$; $z=7.316$--$7.376$ on $h_s=0.00625$; deep states on "
      "$h_s=0.0125$ (continuation). $I=D_0\\,\\mathrm{Re}\\,\\Phi_0$.")

# ---- offsets by cut-off
txt = open("ipm_cutoff_drift.out").read()
rows = []
for line in txt.splitlines():
    p = line.split()
    if len(p) == 9 and p[0].isdigit():
        rows.append([p[0], p[1], p[3], p[5], p[7]])
write("cutoff_table.tex", "$n$ & cut-off 1.5 & cut-off 2 & cut-off 3 & cut-off 5", rows, "rrrrr",
      "Offsets $\\delta_n=\\mathrm{Re}\\,\\Phi_0(\\lambda_n)-n\\pi$ with the repaired tracker for four front cut-offs "
      "$D/D_0$; no untracked points. Cut-off 2 is the registered convention.")

# ---- decay per half-period
L = M.rungs()
z, Re, Im, I = M.phase_table()
o = np.argsort(z)
pz, pe = M.decay_pred(z[o], Re[o], Im[o])
ex = M.extrema(L)
rows = []
for k in range(1, len(ex)):
    zm = 0.5 * (ex[k, 0] + ex[k - 1, 0])
    meas = -np.log(ex[k, 1] / ex[k - 1, 1])
    if zm < 2.5:
        continue
    pr = np.interp(zm, pz, pe)
    rows.append([f"{zm:.2f}", f"{1 / zm:.3f}", f"{pr:.2f}", f"{meas:.2f}", f"{(meas - pr) / pr * zm:.2f}"])
write("decay_table.tex", "$z$ (midpoint) & $\\lambda$ & predicted & measured & (meas.$-$pred.)/(pred.$\\,\\lambda$)", rows,
      "rrrrr", "E-folds of the defect per half-period: predicted $-\\pi\\,d\\,\\mathrm{Im}\\Phi_0/d\\,\\mathrm{Re}\\Phi_0$ "
      "from the repaired phase; measured from consecutive extrema of $|m-2|$.")

# ---- grid shifts
V = M.rich_values()
keys = sorted({(k[0], k[2]) for k in V})
rows = []
for zz, h in keys:
    vals = [V.get((zz, d, h)) for d in (0.0, 0.25, 0.5, 0.75)]
    if sum(v is not None for v in vals) < 3:
        continue
    have = [v for v in vals if v is not None]
    rows.append([f"{zz:.3f}", f"{h:g}"] + [sci(v, 2) if v is not None else "--" for v in vals] +
                [sci(0.5 * (max(have) - min(have)), 1)])
write("shift_table.tex", "$z$ & $h_s$ & $\\delta=0$ & $\\delta=\\tfrac14$ & $\\delta=\\tfrac12$ & $\\delta=\\tfrac34$ & "
      "half-range", rows, "rrrrrrr", "$m-2$ on radial grids shifted by $\\delta h_s$ ($s_{\\min}\\to-120+\\delta h_s$).")

# ---- two solvers
tw = M.two_solver()
sl = M.crossing_slopes(L)
rows = []
for n, dz, dl in tw:
    n = int(n)
    lam = L[n]
    dmdl = sl[n - 1] / lam ** 2
    rows.append([str(n), f"{lam:.10f}", sci(dl, 1), sci(dz, 1), sci(dmdl, 1), sci(dmdl * dl, 1)])
write("twosolver_table.tex", "$n$ & $\\lambda_n$ (marching) & $|\\Delta\\lambda|$ & $|\\Delta z|$ & $|dm/d\\lambda|$ & "
      "$|\\Delta m|$", rows, "rrrrrr",
      "Independent global solver ($h_s=0.0125$, $s_{\\min}=-16$, extrapolated in $s_{\\max}$; $n=1$ at $h_s=0.025$, "
      "$s_{\\min}=-12$; $n=4$ at $h_s=0.00625$) against the marching solver. $|\\Delta m|=|dm/d\\lambda|\\,|\\Delta\\lambda|$.")

# ---- origin truncation and Newton floor
rows = []
for fn, lab in (("ipm_sstart_lam6.out", "6"), ("ipm_sstart_lam5.out", "5")):
    for line in open(fn):
        m = re.search(r"s_start=(-?[0-9.]+):.*\|R\|=([0-9.e+-]+).*m-2=([-+0-9.e]+)", line)
        if m:
            rows.append((lab, float(m.group(1)), float(m.group(2)), float(m.group(3))))
for fn in ("ipm_noise2_ss16.out", "ipm_noise2_ss20.out"):
    for line in open(fn):
        m = re.search(r"ss=(-?[0-9.]+) λ=0\.150920000 \|R\|=([0-9.e+-]+) m-2=([-+0-9.e]+)", line)
        if m:
            rows.append(("6", float(m.group(1)), float(m.group(2)), float(m.group(3))))
for fn in sorted(glob.glob("ipm_nfloor_ss*_lam6.out")):
    m = re.search(r"s_start=(-?[0-9.]+).*\|R\|=([0-9.e+-]+) m-2=([-+0-9.e]+)", open(fn).read())
    if m:
        rows.append(("6", float(m.group(1)), float(m.group(2)), float(m.group(3))))
rows = sorted(set(rows), key=lambda r: (r[0], -r[1]))
write("truncation_table.tex", "profile & $s_{\\rm start}$ & Newton residual $|R|$ & $m-2$",
      [[f"$\\lambda_{r[0]}$", f"{r[1]:g}", sci(r[2], 1), sci(r[3], 2)] for r in rows], "rrrr",
      "$h_s=0.0125$, $N_b=32$; $\\lambda_6=0.15092$, $\\lambda_5=0.1706180880$.")

# ---- deep geometry
T = M.travel_time()
E = M.deep_log()
rows = []
for r in T:
    e = E[np.argmin(abs(E[:, 0] - r[0]))] if len(E) else None
    ext = [f"{e[3]:.2f}", f"{100 * e[4]:.1f}"] if e is not None and abs(e[0] - r[0]) < 2e-3 else ["--", "--"]
    rows.append([f"{r[0]:.3f}", f"{r[6]:.4f}", f"{r[10]:.4f}", f"{r[11]:.4f}", f"{r[1]:.4f}"] + ext)
for e in E[E[:, 0] > T[:, 0].max() + 1e-6]:
    rows.append([f"{e[0]:.3f}", f"{e[1]:.4f}", "--", "--", "--", f"{e[3]:.2f}", f"{100 * e[4]:.1f}"])
write("geometry_table.tex", "$z$ & $\\hat D_{\\min}$ & $x_{\\rm dip}$ & $x_{\\rm cut}$ & $T$ & max $\\partial_xR$ & "
      "indicator (\\%)", rows, "rrrrrrr",
      "$T=\\int_{-10}^{s_{\\rm cut}}ds/\\hat D$. Indicator: $\\max|\\Omega_b+\\partial_xR|/\\max\\partial_xR$ over the dip-to-front "
      "window (deep continuation, $h_s=0.0125$).")

R2 = M.deepres()
rows = [[f"{r[0]:.4f}", f"{r[1]:g}", sci(r[2], 2), f"{r[3]:.4f}", f"{r[4]:.2f}", f"{100 * r[5]:.1f}", f"{r[6]:.4f}",
         f"{r[7]:.4f}"] for r in R2]
write("deepres_table.tex", "$z$ & $h_s$ & $m-2$ & $\\hat D_{\\min}$ & max $\\partial_xR$ & indicator (\\%) & $T$ & "
      "$\\mathrm{Re}\\,\\Phi_0$", rows, "rrrrrrrr")

# ---- tracker
old = []
for fn in ("ipm_wkb_e5_cut2.out", "ipm_wkb_fine_deep_cut2.out"):
    for line in open(fn):
        m = re.search(r"λ=([0-9.]+).*I=([0-9.]+)[-+].*fail=(\d+)", line)
        if m:
            old.append((1 / float(m.group(1)), float(m.group(2)), int(m.group(3))))
rows = []
for zz, Io, nf in sorted(old):
    j = np.argmin(abs(z - zz))
    Ir = f"{I[j]:.5f}" if abs(z[j] - zz) < 2e-3 else "--"
    rows.append([f"{zz:.3f}", f"{Io:.5f}", str(nf), Ir])
write("tracker_table.tex", "$z$ & $I$, original tracker & failed points (counted) & $I$, repaired tracker", rows, "rrrr")

# ---- stage-4 predictions
txt = open("ipm_stage4_predictions.out").read()
rows = []
for line in txt.splitlines():
    m = re.match(r"\s*(\d+)\s+([0-9.]+) \[([0-9.]+), ([0-9.]+)\]\s+([0-9.]+) \[([0-9.]+), ([0-9.]+)\]\s+([0-9.]+)\s+"
                 r"([0-9.e+-]+)\s+([0-9.e+-]+)", line)
    if m:
        g = m.groups()
        slope = float(g[8])
        s003 = float(g[9])                                   # as registered (unrounded slope)
        rows.append([g[0], g[1], f"[{g[2]}, {g[3]}]", g[4], f"[{g[5]}, {g[6]}]", g[7], sci(slope, 1),
                     sci(s003 / 0.3, 1), sci(s003, 1)])
write("stage4_table.tex", "$n$ & $z_n$ & band & $\\lambda_n$ & $\\lambda$ band & e-folds/profile & "
      "$|dm/dz|$ & $\\sigma_m$ ($\\sigma_z=0.01$) & $\\sigma_m$ ($\\sigma_z=0.003$)", rows, "rrcrcrrrr")
