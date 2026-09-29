"""Table 2: independent global-solver values of λ_n (s_max-extrapolated) at hs = 0.025 and 0.0125 vs the marching
solver (Table 1)."""
import re, os
march = [1.9205593, 1.3990961, 1.2523487, 1.1842533, 1.1449857, 1.1194738, 1.1015817, 1.0883384]
def lam_inf(f):
    if not os.path.exists(f):
        return None
    m = re.search(r"λ_∞ = ([0-9.]+)", open(f).read())
    return float(m.group(1)) if m else None
print("| n | march (hs 0.0125 from n=4) | global hs 0.025 | global hs 0.0125 | Δ(global − march) |")
print("|---|---|---|---|---|")
for n, lm in enumerate(march):
    a = lam_inf(f"glob_g{n}_h025_Nb24_s12.log"); b = lam_inf(f"glob_g{n}_h0125_Nb24_s12.log")
    best = b if b is not None else a
    fa = f"{a:.7f}" if a else "—"; fb = f"{b:.7f}" if b else "—"
    d = f"{best - lm:+.1e}" if best else "—"
    print(f"| {n} | {lm:.7f} | {fa} | {fb} | {d} |")
