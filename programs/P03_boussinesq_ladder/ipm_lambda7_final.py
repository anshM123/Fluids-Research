"""λ₇ from the converged grid (h_s = 0.00625), averaged over grid shifts, with the systematic checks at the crossing.
Supplements ipm_lambda7_analysis.py (the pre-committed rule, which could not complete because the production scan
was paused before its crossing to test the front–grid locking; see ipm_lambda7_analysis_committed_rule.out).

Statistical part: at z = 7.256, 7.316, 7.346, 7.376 the value of m − 2 is the mean over four grid shifts
(δ = 0, ¼, ½, ¾); its error is the shift scatter / 2.  A cubic through these four means (weighted) gives z₇ and the
local slope; σ_stat = σ(mean at the crossing)/|slope|.
Systematic part: the change of m at z = 7.346 for Nb 48, s_start −22 and s_max 130 (sys7_*.out), each converted to
a shift of z₇ by dividing by the slope; σ_sys is their quadrature sum (the h_s error at 0.00625 is bounded by the
shift scatter, because the locking error falls 50–120× per halving of h_s)."""
import numpy as np, glob, re

vals = {}
for f in glob.glob("rich7_*hs0.00625.out"):
    for line in open(f):
        m = re.search(r"z=([0-9.]+) shift=([0-9.]+).*m-2=([-+0-9.e]+)", line)
        if m: vals.setdefault(float(m.group(1)), []).append(float(m.group(3)))
for line in open("ipm_scan_h7.out"):
    m = re.search(r"^z=([0-9.]+) .*m-2=([-+0-9.e]+)", line)
    if m:
        z = float(m.group(1))
        if z in vals or round(z, 3) in (7.256, 7.316, 7.376):
            vals.setdefault(round(z, 3), []).append(float(m.group(2)))
Z = np.array(sorted(k for k in vals if len(vals[k]) >= 4))
M = np.array([np.mean(vals[z]) for z in Z]); E = np.array([np.std(vals[z], ddof=1) / np.sqrt(len(vals[z])) for z in Z])
print("shift-averaged m − 2 (h_s = 0.00625):")
for z, m_, e in zip(Z, M, E):
    print(f"  z = {z:.3f}: {m_:+.4e} ± {e:.1e}  (n = {len(vals[z])}, values {np.array2string(np.array(vals[z]), precision=3)})")
c = np.polyfit(Z - 7.346, M, 3, w=1 / E)
r = np.roots(c).real + 7.346; r = r[(r > Z.min()) & (r < Z.max())]
z7 = r[np.argmin(abs(r - 7.346))]
slope = np.polyval(np.polyder(c), z7 - 7.346)
e_cross = np.interp(z7, Z, E)
s_stat = e_cross / abs(slope)
print(f"z₇ = {z7:.4f}, slope dm/dz = {slope:.3e}, σ_stat = {s_stat:.4f}")

base = np.mean([v for v in vals[7.346]]) if 7.346 in vals else None
base0 = None
for line in open(glob.glob("rich7_7.346_d0.0_hs0.00625.out")[0]):
    m = re.search(r"m-2=([-+0-9.e]+)", line)
    if m: base0 = float(m.group(1))
sys_shifts = []
for f, lab in (("sys7_nb48.out", "Nb 48 vs 32"), ("sys7_ss22.out", "s_start −22 vs −20"), ("sys7_sm130.out", "s_max 130 vs 100")):
    try:
        mm = [float(x) for x in re.findall(r"m-2=([-+0-9.e]+)", open(f).read())]
    except FileNotFoundError:
        mm = []
    if mm:
        d = mm[-1] - base0; sys_shifts.append(d / abs(slope))
        print(f"  {lab:20s}: Δm = {d:+.2e} → Δz₇ = {d / abs(slope):+.4f}")
    else:
        print(f"  {lab:20s}: not available")
s_sys = np.sqrt(np.sum(np.square(sys_shifts))) if sys_shifts else np.nan
s_tot = np.hypot(s_stat, s_sys) if np.isfinite(s_sys) else s_stat
print(f"σ_sys = {s_sys:.4f}; σ_total = {s_tot:.4f};  λ₇ = {1/z7:.6f}  (z₇ = {z7:.4f} ± {s_tot:.4f})")
PRED = {"H1 (shifted law, λ2–λ6)": 7.3415, "H1' (shifted law, λ0–λ6)": 7.3404, "H2 (geometric)": 7.3343,
        "H3 (3-point)": 7.3342, "H4 (phase, coarse branch)": 7.352}
print(f"comparison (registered rule: disfavoured if |Δ| > 2σ; decisive only if σ ≤ 0.003 → {'decisive' if s_tot <= 0.003 else 'NOT decisive'}):")
for k, v in PRED.items():
    print(f"  {k:28s} {v:.4f}: Δ = {v - z7:+.4f} = {(v - z7) / s_tot:+.1f}σ {'(disfavoured)' if abs(v - z7) > 2 * s_tot else ''}")
