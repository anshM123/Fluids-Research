"""Stage-2 blind predictions for IPM (PREDICTIONS_IPM.md): from the computed rungs λ₀…λ_N and the unstable
spectra of U₁…U_N, predict λ_{N+1}, the index of U_{N+1} and its unstable eigenvalues, BEFORE U_{N+1} is computed.

Rung: z_n = 1/λ_n, Δ_n = z_n − z_{n−1}. P2 gives Δ_n = a − κ/z̄_n (z̄_n the interval midpoint). a and κ come
from the last two spacings, and Δ_{N+1} is solved self-consistently. Monotone spacing (P2) gives the hard bound
Δ_{N+1} > Δ_N, i.e. λ_{N+1} < 1/(z_N + Δ_N).

Spectrum: P3 gives index N+1. P4 gives ν̂ = μ/λ on a lattice with a displaced top member:
  - lowest member ν̂_0: linear in 1/z through the last two rungs;
  - lower gaps: s = 1 + β/z, with β from the lowest gap of U_N (P4: equal lower gaps);
  - top gap: linear in n through the last two rungs.
A parameter-free variant takes s = 1 and the offset of U_N (the pure asymptotic lattice).
usage: python3 ipm_stage2_predict.py   (reads ipm_stage2_data.py)"""
import numpy as np
from ipm_stage2_data import lam, nuhat        # lam: [λ0, …, λN]; nuhat: {n: sorted non-trivial ν̂ roots of U_n}

N = len(lam) - 1
z = 1 / np.array(lam)
dz = np.diff(z); zm = 0.5 * (z[1:] + z[:-1])
# rung
A = np.array([[1, -1 / zm[-2]], [1, -1 / zm[-1]]]); a, kap = np.linalg.solve(A, dz[-2:])
D = dz[-1]
for _ in range(50):
    D = a - kap / (z[-1] + D / 2)
z_next = z[-1] + D
lam_next = 1 / z_next
lam_bound = 1 / (z[-1] + dz[-1])
print(f"rungs used: λ = {np.round(lam, 10).tolist()}")
print(f"spacings Δ_n = {np.round(dz, 5).tolist()}; model Δ = a − κ/z̄ with a = {a:.5f}, κ = {kap:.5f}")
print(f"PREDICTION λ_{N+1} = {lam_next:.6f} (1/λ = {z_next:.5f}); hard bound from monotone spacing: λ_{N+1} < {lam_bound:.6f}")
print(f"   (for reference, the empirical law 1/λ = 1.1459 n + 0.9723 gives {1 / (1.1459 * (N + 1) + 0.9723):.6f})")
# spectrum
nus = {n: np.array(sorted(v)) for n, v in nuhat.items()}
lo = {n: v[0] for n, v in nus.items() if len(v)}
ns = sorted(lo)
p0 = np.polyfit([1 / z[n] for n in ns[-2:]], [lo[n] for n in ns[-2:]], 1)
nu0 = np.polyval(p0, 1 / z_next)
gl = nus[N][1] - nus[N][0] if len(nus[N]) >= 3 else None
# lower gap s = 1 + β/z: β from the lowest gap of U_N when U_N has one; otherwise the prior β = 0.83 ± 0.17
# (Hou–Luo 0.67, 2D Boussinesq 1.0: the only calibration taken from the other two models)
beta = (gl - 1) * z[N] if gl is not None else 0.83
s_next = 1 + beta / z_next
s_rng = (1 + (beta - 0.17) / z_next, 1 + (beta + 0.17) / z_next) if gl is None else (s_next, s_next)
tops = {n: v[-1] - v[-2] for n, v in nus.items() if len(v) >= 2}
nt = sorted(tops)
# top gap: linear in n through the last two rungs; with only one, grow it by 4 % ± 3 % per rung (HL 2 %, 2D 6 %)
gtop = np.polyval(np.polyfit(nt[-2:], [tops[n] for n in nt[-2:]], 1), N + 1) if len(nt) >= 2 else 1.04 * tops[nt[-1]]
pred = [nu0 + j * s_next for j in range(N)] + [nu0 + (N - 1) * s_next + gtop]
free = [nus[N][0] + j for j in range(N + 1)]
print(f"PREDICTION index(U_{N+1}) = {N + 1} (non-trivial unstable modes, all real)")
print(f"PREDICTION ν̂ = μ/λ of U_{N+1}: {np.round(pred, 3).tolist()}  (lowest {nu0:.3f}; lower gap s = {s_next:.3f}, "
      f"β = {beta:.3f}; top gap {gtop:.3f})")
print(f"           μ = {np.round(np.array(pred) * lam_next, 5).tolist()}")
print(f"           lower-gap range s ∈ [{s_rng[0]:.3f}, {s_rng[1]:.3f}]")
print(f"parameter-free lattice (s = 1, offset of U_{N}): ν̂ = {np.round(free, 3).tolist()}")
# P5: spectral flow inside the last computed interval [z_{N-1}, z_N], rigid lattice moving linearly in the fraction f
if N >= 2 and len(nus[N - 1]) >= 1 and len(nus[N]) >= 2:
    a_prev, b = nus[N - 1], nus[N]
    slope = b[-1] - a_prev[-1]                     # the member continuing from U_{N−1}'s top to U_N's top
    fstar = 1 - b[0] / slope
    print(f"P5 (interval z = {z[N - 1]:.4f} → {z[N]:.4f}): members move rigidly with slope {slope:.3f} per interval; "
          f"the new member enters μ = 0 at f* = {fstar:.3f} (z = {z[N - 1] + fstar * (z[N] - z[N - 1]):.4f})")
    for f in (0.25, 0.5, 0.75):
        zf = z[N - 1] + f * (z[N] - z[N - 1])
        pr = [v + slope * (f - 1) for v in b if v + slope * (f - 1) > 0]
        print(f"   f = {f:.2f} (z = {zf:.4f}, λ = {1 / zf:.7f}): predicted ν̂ = {np.round(pr, 3).tolist()}")
