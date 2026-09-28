"""Final analysis of the arc λ₂ → cusp: model-free extrema of p(δ) (local quadratic fits in ln δ), spacing ratios
vs the prediction e^{π/2τ} = 11.232, and asymptotic fits (fixed and free frequency) on the combined data."""
import numpy as np
from cusp_fit import parse, fit, OMEGA

G = np.load("pdense_G.npy"); E = np.load("pdense_E.npy")
C = parse("pin_C.log"); C = C[np.argsort(-C[:, 0])]
Cs = C[C[:, 0] < E[:, 0].min() * 0.99]
D = np.vstack([G[:, :3], E[:, :3], Cs[:, :3]])
D = D[np.argsort(D[:, 0])]
_, iu = np.unique(np.round(np.log(D[:, 0]), 9), return_index=True)
D = D[iu]
L, lam, p = np.log(D[:, 0]), D[:, 1], D[:, 2]
print(f"{len(D)} points, δ ∈ [{D[:, 0].min():.2e}, {D[:, 0].max():.2e}]")
ext = []
for i in range(1, len(D) - 1):
    if (p[i] - p[i - 1]) * (p[i + 1] - p[i]) < 0:
        j = slice(max(i - 2, 0), min(i + 3, len(D)))
        cf = np.polyfit(L[j] - L[i], p[j], 2)
        x0 = -cf[1] / (2 * cf[0])
        ext.append((np.exp(L[i] + x0), np.polyval(cf, x0), "max" if cf[0] < 0 else "min"))
print("extrema of p(δ):")
for d0, p0, kind in ext:
    print(f"   {kind}: δ = {d0:.4e}   p = {p0:.9f}")
print(f"predicted spacing factor e^(π/2τ) = {np.exp(np.pi / OMEGA):.4f}")
for k in range(len(ext) - 1):
    print(f"   ratio δ[{k+1}]/δ[{k}] = {ext[k+1][0]/ext[k][0]:.4f}")
for dmax in (3e-3, 1e-3):
    s = D[:, 0] <= dmax
    for name, y in (("p", p), ("λ", lam)):
        pf, rf = fit(D[s, 0], y[s], free_omega=False, use_log=False)
        pw, rw = fit(D[s, 0], y[s], free_omega=True, use_log=False)
        pl, rl = fit(D[s, 0], y[s], free_omega=True, use_log=True)
        print(f"δ ≤ {dmax:g} ({s.sum()} pts) {name}: fixed 2τ → {name}* = {pf[0]:.10f} (rms {rf:.1e}); "
              f"free ω = {pw[5]:.4f} → {name}* = {pw[0]:.10f}; free ω + δlnδ: ω = {pl[5]:.4f}, A = {pl[1]:+.4f}, {name}* = {pl[0]:.10f}")
