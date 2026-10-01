"""Stage 4: λ₈, λ₉, λ₁₀ from the repaired phase (registered in PREDICTIONS_IPM.md before any test).
Rule: Re Φ₀(z_n) = nπ + δ_n with Re Φ₀(z) from the repaired deep table (wkb_phase3 on the h_s = 0.0125 continuation
e5; states with untracked front points excluded), smoothed by a least-squares quadratic in z over a window of ±0.5
around the target (the state-to-state scatter of Re Φ₀ is ≈ 0.04 from front–grid sampling of the dip).
Offset scenarios (the only freedom; fixed here, not refitted later):
  A  δ = δ₆ = 1.974 (the latest measured offset)                                  ← central value
  B  δ = 2.031 (mean of n = 1–5)
  C  δ drifting on with the n = 5 → 6 step: δ_n = 1.974 − 0.057 (n − 6)
Uncertainty band = [min, max] over A–C, widened by the grid correction of Re Φ₀ (h_s 0.0125 → 0.00625, measured at
z = 8.066 and 8.786 by ipm_deepres_check.py; its mean is applied as a shift and its spread added) and by the
smoothing-fit error.
For each rung the script also states what a decisive test needs: the expected defect amplitude near the rung (the
λ₇ crossing slope continued with the decay rate predicted by Im Φ₀, Fig. 2b, times 1 + λ) and the error in m that would
locate z_n to 0.003 (σ_m ≤ 0.003 |dm/dz|)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, re, glob

R = {}
for f in ("ipm_phase3_l7.out", "ipm_phase3_deep.out"):
    for line in open(f):
        m = re.search(r"z=\s*([0-9.]+) .*Φ=\s*([-0-9.]+)\s*([-+][0-9.]+)i I=([0-9.]+) front_fail=(\d+)", line)
        if m and int(m.group(5)) == 0:
            R[round(float(m.group(1)), 4)] = (float(m.group(2)), float(m.group(3)))
z = np.array(sorted(R)); Re = np.array([R[k][0] for k in z]); Im = np.array([R[k][1] for k in z])

# grid correction of Re Φ₀ from the fine-grid spot checks
corr = []
for f in sorted(glob.glob("ipm_deepres_z*.out")):
    v = {}
    for line in open(f):
        m = re.search(r"z=([0-9.]+) hs=([0-9.]+):.*ReΦ=([0-9.]+)", line)
        if m: v[float(m.group(2))] = (float(m.group(1)), float(m.group(3)))
    if 0.0125 in v and 0.00625 in v:
        corr.append((v[0.00625][0], v[0.00625][1] - v[0.0125][1]))
corr = np.array(corr) if corr else np.zeros((0, 2))
dgrid = float(np.mean(corr[:, 1])) if len(corr) else 0.0
sgrid = float(np.max(np.abs(corr[:, 1] - dgrid))) if len(corr) > 1 else abs(dgrid)
print("grid correction of Re Φ₀ (h_s 0.00625 − 0.0125):", ", ".join(f"z = {a:.3f}: {b:+.4f}" for a, b in corr),
      f"→ applied {dgrid:+.4f}, spread ± {sgrid:.4f}")


def phase_at(zz):
    w = np.abs(z - zz) <= 0.5
    c, cov = np.polyfit(z[w] - zz, Re[w], 2, cov=True)
    return c, float(np.sqrt(cov[2, 2]))


def solve(target, z0):
    zz = z0
    for _ in range(30):
        c, e = phase_at(zz)
        f0, d0 = np.polyval(c, 0.0) + dgrid - target, np.polyval(np.polyder(c), 0.0)
        step = -f0 / d0; zz += step
        if abs(step) < 1e-9: break
    c, e = phase_at(zz)
    return zz, np.polyval(np.polyder(c), 0.0), e


# defect scale: λ₇ crossing slope continued with the predicted decay per half-period (Fig. 2b) × (1 + λ)
sl7 = abs(float(re.search(r"slope dm/dz = ([-0-9.e+]+)", open("ipm_lambda7_final.out").read()).group(1)))
o = np.argsort(z)
print(f"\n{'n':>2} {'scenario':>9} {'δ':>6} {'z_n':>8} {'λ_n':>9}")
rows = []
slope = sl7; zprev = 7.3451; dprev = 4.54
for n in (8, 9, 10):
    zs = {}
    for lab, d in (("A", 1.974), ("B", 2.031), ("C", 1.974 - 0.057 * (n - 6))):
        zn, dRe, efit = solve(n * np.pi + d, 7.9 + 0.6 * (n - 8))
        zs[lab] = (zn, dRe, efit)
        print(f"{n:2d} {lab:>9} {d:6.3f} {zn:8.4f} {1 / zn:9.6f}")
    zA, dRe, efit = zs["A"]
    lo = min(v[0] for v in zs.values()) - (sgrid + efit) / dRe
    hi = max(v[0] for v in zs.values()) + (sgrid + efit) / dRe
    # decay between consecutive rungs: Re Φ₀ advances by π; Im Φ₀ changes by the local slope × π
    w = np.abs(z - 0.5 * (zprev + zA)) <= 0.45
    rate = -np.pi * np.polyfit(Re[w], Im[w], 1)[0] * (1 + 1 / zA)
    slope = slope * np.exp(-rate) * dRe / dprev; dprev = dRe
    rows.append((n, zA, lo, hi, 1 / zA, 1 / hi, 1 / lo, rate, slope, 0.003 * slope))
    zprev = zA
print("\nRegistered values (central = scenario A; band = scenarios A–C ⊕ grid spread ⊕ fit error):")
print(f"{'n':>2} {'z_n':>8} {'band':>19} {'λ_n':>9} {'λ band':>21} {'e-folds/rung':>12} {'|dm/dz| expected':>16} {'σ_m for ±0.003':>14}")
for n, zA, lo, hi, l, l1, l2, rate, sl, sm in rows:
    print(f"{n:2d} {zA:8.4f} [{lo:.4f}, {hi:.4f}] {l:9.6f} [{l1:.6f}, {l2:.6f}] {rate:12.2f} {sl:16.1e} {sm:14.1e}")
