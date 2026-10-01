"""IPM endpoint geometry and asymptotic class from the repaired phase (no use of m − 2).
Inputs: ipm_phase3_{rungs,l7,deep}.out (Re/Im Φ₀ and I = D₀ Re Φ₀ from wkb_phase3), ipm_travel_time.out (T = ∫ds/D̂
and its partition, dip depth/width/position, cut-off position), ipm_deep_e5.out (max ∂ₓR and the resolution
indicator |Ω_b + ∂ₓR|/max on the deep states), and the fine-grid spot checks deepres_*.out / ipm_phase3_deepres.out
when present.
Reports
 1. I(z) and T(z): local slopes and the effective wave number dI/dT (the phase growth is geometric if constant);
 2. the dip and front: D̂_min(z), x_dip(z), x_cut(z), max ∂ₓR(z) and the resolution indicator;
 3. three asymptotic classes for the phase, fitted to I(z) on z ≥ 4 with the shallow points weighted by their grid
    error and the deep points by their scatter:
      (A) polynomial growth, I = a + b z + c z²                 → λ_c = 0 (Re Φ₀ = 2zI grows like z² or z³)
      (B) pole,  I = a + b z + C/(z_c − z)                       → λ_c = 1/z_c > 0 (infinite accumulation)
      (C) log,   I = a + b z − C ln(1 − z/z_c)                    → λ_c = 1/z_c > 0 (infinite accumulation)
    For (B) and (C) the χ² profile in z_c gives the range of z_c the phase data allow (Δχ² ≤ 4); a profile that is
    flat out to z_c → ∞ means the data cannot detect a finite λ_c;
 4. finite termination (the branch ends with finite phase) would need D̂_min → 0 or max ∂ₓR → ∞ at finite z_e:
    linear extrapolations of D̂_min and 1/max ∂ₓR over the last 1.0 in z, against their local deceleration."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "lib"))  # solver library
import numpy as np, re, os

def phase_rows():
    R = {}
    for f in ("ipm_phase3_rungs.out", "ipm_phase3_l7.out", "ipm_phase3_deep.out"):
        for line in open(f):
            m = re.search(r"z=\s*([0-9.]+) D0=([0-9.]+) D̂dip=([0-9.]+) s_cut=\s*(-?[0-9.]+) Φ=\s*([-0-9.]+)\s*([-+][0-9.]+)i I=([0-9.]+) front_fail=(\d+)", line)
            if m:
                v = [float(x) for x in m.groups()]
                if v[7] > 0:
                    print(f"   (excluded: z = {v[0]:.3f}, {int(v[7])} front points not tracked)"); continue
                R[round(v[0], 3)] = v
    A = np.array([R[k] for k in sorted(R)])
    return A[:, 0], A[:, 4], A[:, 5], A[:, 6]

def tt_rows():
    R = {}
    for line in open("ipm_travel_time.out"):
        p = line.split()
        if len(p) == 13 and re.match(r"^\d", p[1]):
            v = [float(x) for x in p[1:]]; R[round(v[0], 3)] = v
    return np.array([R[k] for k in sorted(R)])

def e5_rows():
    out = []
    for line in open("ipm_deep_e5.out"):
        m = re.search(r"z=([0-9.]+):.*D̂_dip=([0-9.]+) front s=([-0-9.]+) max∂ₓR=([0-9.]+) \|Ω_b\+∂ₓR\|/max=([0-9.e+-]+)", line)
        if m: out.append([float(x) for x in m.groups()])
    return np.array(out)

z, Re, Im, I = phase_rows()
T = tt_rows()
E = e5_rows()
print("1. Phase coefficient I = D₀ Re Φ₀ and travel time T = ∫_{-10}^{cut} ds/D̂")
print("   z        I        T      dI/dz    dT/dz   dI/dT")
zt = T[:, 0]
for k in range(len(z)):
    j = np.argmin(abs(zt - z[k]))
    if abs(zt[j] - z[k]) > 2e-3: continue
    if k == 0:
        print(f"  {z[k]:6.3f}  {I[k]:.5f}  {T[j, 1]:.4f}"); continue
    jp = np.argmin(abs(zt - z[k - 1]))
    dz = z[k] - z[k - 1]
    if dz < 0.05: continue
    dI, dT = (I[k] - I[k - 1]) / dz, (T[j, 1] - T[jp, 1]) / dz
    print(f"  {z[k]:6.3f}  {I[k]:.5f}  {T[j, 1]:.4f}  {dI:7.4f}  {dT:7.4f}  {dI / dT:6.3f}")
deep = z >= 7.45
jj = [np.argmin(abs(zt - x)) for x in z[deep]]
cI = np.polyfit(z[deep], I[deep], 1); cT = np.polyfit(zt[jj], T[jj, 1], 1)
cIT = np.polyfit(T[jj, 1], I[deep], 1)
print(f"   deep (z ≥ 7.45): dI/dz = {cI[0]:.4f}, dT/dz = {cT[0]:.4f}, dI/dT = {cIT[0]:.3f} "
      f"(rms of I about the line in T: {np.std(I[deep] - np.polyval(cIT, T[jj, 1])):.4f}; about the line in z: {np.std(I[deep] - np.polyval(cI, z[deep])):.4f})")
sh = (z > 3) & (z < 7.4)
jj2 = [np.argmin(abs(zt - x)) for x in z[sh]]
print(f"   shallow (3 < z < 7.4): dI/dT = {np.polyfit(T[jj2, 1], I[sh], 1)[0]:.3f}")
print(f"   Im Φ₀ per half-period of Re Φ₀ (π dImΦ/dReΦ): " + ", ".join(
      f"{0.5*(z[k]+z[k+1]):.2f}: {np.pi*(Im[k+1]-Im[k])/(Re[k+1]-Re[k]):+.2f}" for k in range(len(z) - 1) if z[k+1] - z[k] > 0.5 or z[k] in (7.466, 8.426, 9.026)))

print("\n2. Geometry along the branch (travel-time table; deep states h_s = 0.0125 unless marked)")
print("   z       D̂min(grid) D̂min(parab) w_core  cells  x_dip   x_cut   T_approach T_core  T_front   max∂ₓR  |Ω_b+∂ₓR|/max")
for r in T:
    e = E[np.argmin(abs(E[:, 0] - r[0]))] if len(E) else None
    extra = f"  {e[3]:6.2f}   {e[4]:.1e}" if e is not None and abs(e[0] - r[0]) < 2e-3 else ""
    print(f"  {r[0]:6.3f}  {r[6]:.4f}     {r[7]:.4f}     {r[8]:.4f}  {int(r[9]):4d}   {r[10]:.4f}  {r[11]:.4f}  {r[3]:7.4f}   {r[4]:7.4f}  {r[5]:7.4f}{extra}")

print("\n3. Asymptotic class of the phase (fits to I(z), z ≥ 4)")
sel = z >= 4.0
zz, II = z[sel], I[sel]
sig = np.where(zz < 7.4, 6e-4, 2.2e-3)     # shallow: grid error of I (h_s 0.0125 vs 0.00625 at z = 7.316: 7e-4); deep: scatter
def chi2(model, p):
    return float(np.sum(((II - model(zz, *p)) / sig) ** 2))
X = np.vstack([np.ones_like(zz), zz, zz ** 2]).T / sig[:, None]
pA = np.linalg.lstsq(X, II / sig, rcond=None)[0]; chiA = float(np.sum((X @ pA - II / sig) ** 2))
print(f"   (A) a + bz + cz²: c = {pA[2]:+.5f}, χ² = {chiA:.1f} for {len(zz) - 3} dof")
X1 = np.vstack([np.ones_like(zz), zz]).T / sig[:, None]
p1 = np.linalg.lstsq(X1, II / sig, rcond=None)[0]; chi1 = float(np.sum((X1 @ p1 - II / sig) ** 2))
print(f"       linear only: b = {p1[1]:.4f}, χ² = {chi1:.1f} for {len(zz) - 2} dof")
for name, f in (("(B) pole  C/(z_c − z)", lambda x, zc: 1 / (zc - x)), ("(C) log  −ln(1 − z/z_c)", lambda x, zc: -np.log(1 - x / zc))):
    prof = []
    for zc in np.concatenate([np.linspace(9.8, 30, 203), [40, 60, 100, 300, 1000]]):
        Xc = np.vstack([np.ones_like(zz), zz, f(zz, zc)]).T / sig[:, None]
        p = np.linalg.lstsq(Xc, II / sig, rcond=None)[0]
        prof.append((zc, float(np.sum((Xc @ p - II / sig) ** 2)), p[2]))
    prof = np.array(prof); kmin = np.argmin(prof[:, 1])
    ok = prof[prof[:, 1] <= prof[kmin, 1] + 4]
    print(f"   {name}: best z_c = {prof[kmin, 0]:.1f} (χ² {prof[kmin, 1]:.1f}, C = {prof[kmin, 2]:+.4f}); "
          f"Δχ² ≤ 4 for z_c ∈ [{ok[:, 0].min():.1f}, {ok[:, 0].max():.0f}] → λ_c ∈ [{1/ok[:, 0].max():.4f}, {1/ok[:, 0].min():.4f}]; "
          f"χ² at z_c = 27 (λ_c = 0.037): {np.interp(27, prof[:, 0], prof[:, 1]):.1f}")

print("   deep range only (z ≥ 7.45; σ = rms scatter about a line):")
dd = z >= 7.45
zd, Id = z[dd], I[dd]
c1 = np.polyfit(zd, Id, 1); sd = np.std(Id - np.polyval(c1, zd), ddof=2)
c2, cov2 = np.polyfit(zd, Id, 2, cov=True)
print(f"     I: slope {c1[0]:.4f}, curvature c = {c2[0]:+.4f} ± {np.sqrt(cov2[0, 0]):.4f} (scatter {sd:.4f})")
jd = [np.argmin(abs(zt - x)) for x in zd]
Td = T[jd, 1]
t1 = np.polyfit(zd, Td, 1); t2, ct2 = np.polyfit(zd, Td, 2, cov=True)
print(f"     T: slope {t1[0]:.4f}, curvature c = {t2[0]:+.5f} ± {np.sqrt(ct2[0, 0]):.5f} (scatter {np.std(Td - np.polyval(t1, zd), ddof=2):.5f})")
for name, f in (("pole", lambda x, zc: 1 / (zc - x)), ("log", lambda x, zc: -np.log(1 - x / zc))):
    prof = []
    for zc in np.concatenate([np.linspace(9.5, 30, 206), [40, 60, 100, 300, 1000]]):
        Xc = np.vstack([np.ones_like(zd), zd, f(zd, zc)]).T
        p = np.linalg.lstsq(Xc, Id, rcond=None)[0]
        prof.append((zc, float(np.sum(((Xc @ p - Id) / sd) ** 2)), p[2]))
    prof = np.array(prof); kmin = np.argmin(prof[:, 1]); ok = prof[prof[:, 1] <= prof[kmin, 1] + 4]
    print(f"     {name}: best z_c = {prof[kmin, 0]:.1f} (C = {prof[kmin, 2]:+.4f}); Δχ² ≤ 4 for z_c ∈ [{ok[:, 0].min():.1f}, {ok[:, 0].max():.0f}]")

print("\n4. Finite termination indicators (deep states, h_s = 0.0125)")
if len(E):
    zE, dmin, Gm = E[:, 0], E[:, 1], E[:, 3]
    for lab, y in (("D̂_min (grid)", dmin), ("1/max ∂ₓR", 1 / Gm)):
        last = zE >= zE.max() - 1.0 - 1e-9; first = (zE >= 7.45) & (zE <= 8.5)
        s_last = np.polyfit(zE[last], y[last], 1); s_first = np.polyfit(zE[first], y[first], 1)
        print(f"   {lab:14s}: slope {s_first[0]:+.4f}/unit z on z ∈ [7.45, 8.5], {s_last[0]:+.4f} on the last unit; "
              f"linear zero of the last-unit fit at z_e = {-s_last[1] / s_last[0]:.1f}; "
              f"log-slope {np.polyfit(zE[first], np.log(y[first]), 1)[0]:+.3f} → {np.polyfit(zE[last], np.log(y[last]), 1)[0]:+.3f}")
for f in sorted(os.listdir(".")):
    if f.startswith("deepres_") and f.endswith(".out"):
        print("   fine-grid check:", open(f).read().strip())
