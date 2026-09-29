"""Robust eigenvalues of the global-discretisation stability problem (glob_stab.py).
An eigenvalue counts if (i) it is found with residual < 1e-10 at ≥ 2 different shifts, and (ii) when the s_min = −20
run exists, it is also found there (to 2e-4). Truncation-induced modes (a line Re μ ≈ c/|s_min| with Im μ spaced by
2π/travel time, a non-normal pseudospectral cloud, and inflow-boundary modes localised at s_min) fail (i) or (ii)."""
import re, os, numpy as np
def parse(f, rmax=1e-10):
    shifts = []
    for line in open(f):
        if line.startswith('σ='):
            shifts.append(np.array([complex(m_ + 'j') for m_, r_ in
                                    re.findall(r"([-0-9.]+[+-][0-9.]+)i\(r=([0-9.e+-]+)\)", line) if float(r_) < rmax]))
    return shifts
def robust(n, tol=2e-4):
    f12, f20 = f"gstab_g{n}.log", f"gstab_g{n}s20.log"
    if not os.path.exists(f12):
        return None, False
    sh = parse(f12)
    allv = np.concatenate([s for s in sh if len(s)]) if sh else np.array([])
    cand = []
    for u in allv:
        if u.real <= 0 or any(abs(u - v) < tol for v in cand):
            continue
        if sum(np.min(np.abs(s - u)) < tol for s in sh if len(s)) >= 2:
            cand.append(u)
    both = os.path.exists(f20) and 'RESULT' in open(f20).read()
    if both:
        s20 = np.concatenate([s for s in parse(f20) if len(s)])
        cand = [u for u in cand if np.min(np.abs(s20 - u)) < tol]
    return sorted(cand, key=lambda c: -c.real), both
if __name__ == "__main__":
    for n in range(8):
        ev, both = robust(n)
        if ev is None:
            continue
        nontriv = [u for u in ev if abs(u - 1) > 1e-3]
        print(f"λ_{n} ({'s_min −12 & −20' if both else 's_min −12 only'}): " +
              ", ".join(f"{u.real:.5f}{u.imag:+.5f}i" for u in ev) + f"   → {len(nontriv)} non-trivial with Re μ > 0")
