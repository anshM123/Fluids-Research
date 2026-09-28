"""Main figure: the complete CCF self-similar branch, its three smooth profiles and its terminal cusp.
Data: branch_map.npy (λ ≥ 0.4625, mapped grid), large_lambda_arc.log (λ0 … 8), pdense_G.npy (δ 0.0104→0.0215),
pdense_E.npy (δ 0.0104→8e−6, dense sampling, fully re-gridded every step), pin_C.log (δ → 2.6e−6)."""
import numpy as np, re, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cusp_fit import parse, fit, model, OMEGA

bm = np.load("branch_map.npy"); bm = bm[np.argsort(bm[:, 0])]
bm = bm[bm[:, 0] <= 1.19]
large = []
for line in open("large_lambda_arc.log"):
    m = re.search(r"lam=([\d.]+) p=([\d.]+) min_den=([\d.]+) eta_min=(-?[\d.]+)", line)
    if m and float(m.group(4)) < -29.9:
        large.append((float(m.group(1)), float(m.group(2))))
large = np.array(large)
G = np.load("pdense_G.npy")                      # (δ, λ, p, w, N), δ increasing
E = np.load("pdense_E.npy")                      # δ decreasing
C = parse("pin_C.log"); C = C[np.argsort(-C[:, 0])]
Cs = C[C[:, 0] < E[:, 0].min() * 0.99]
arc = np.vstack([G[::-1, :3], E[1:, :3], Cs[:, :3]])        # ordered along the branch (δ decreasing)
lam_s = [1.180777662899, 0.6057337012, 0.471324227767]
dall = np.concatenate([E[:, 0], Cs[:, 0]]); pall = np.concatenate([E[:, 2], Cs[:, 2]])
lall = np.concatenate([E[:, 1], Cs[:, 1]])
sel = dall <= 3e-3
parp, rmsp = fit(dall[sel], pall[sel], free_omega=False, use_log=False)
parl, rmsl = fit(dall[sel], lall[sel], free_omega=False, use_log=False)
parpf, _ = fit(dall[sel], pall[sel], free_omega=True, use_log=False)
parlf, _ = fit(dall[sel], lall[sel], free_omega=True, use_log=False)
pstar, lstar = parp[0], parl[0]

fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.8))
a = ax[0, 0]
for arr in (large, bm):
    a.plot(arr[:, 0], arr[:, 1], "k-", lw=1.4)
a.plot(arc[:, 1], arc[:, 2], "k-", lw=1.4)
a.axhline(2.0, color="tab:red", lw=0.8, ls="--")
labels = [("λ₀ = 1.18078\nstable", (8, -30)), ("λ₁ = 0.60573\n1 unstable mode", (-10, 14)),
          ("λ₂ = 0.47132\n2 unstable modes", (6, -34))]
for L, (lab, off) in zip(lam_s, labels):
    a.plot(L, 2.0, "o", color="tab:blue", ms=7, zorder=5)
    a.annotate(lab, (L, 2.0), textcoords="offset points", xytext=off, fontsize=8)
a.plot(lstar, pstar, "*", color="tab:orange", ms=13, zorder=6)
a.annotate("terminal cusp", (lstar, pstar), textcoords="offset points", xytext=(-8, 16), fontsize=8, color="tab:orange")
a.set_xscale("log"); a.set_xlabel("λ   (θ ~ (T−t)^λ)"); a.set_ylabel("local exponent p at the origin")
a.set_title("(a) the connected branch of CCF self-similar profiles", fontsize=10)
a.set_xlim(0.42, 8.5); a.set_ylim(1.05, 2.3)

b = ax[0, 1]
selb = bm[:, 0] < 0.5
b.plot(bm[selb, 0], bm[selb, 1], "k-", lw=1.2)
b.plot(arc[:, 1], arc[:, 2], "k-", lw=1.2)
b.axhline(2.0, color="tab:red", lw=0.8, ls="--")
b.plot(lam_s[2], 2.0, "o", color="tab:blue", ms=7)
b.annotate("λ₂", (lam_s[2], 2.0), textcoords="offset points", xytext=(6, 6), fontsize=9)
b.plot(lstar, pstar, "*", color="tab:orange", ms=13)
b.set_xlim(0.4530, 0.4745); b.set_ylim(1.9985, 2.0145)
b.set_xlabel("λ"); b.set_ylabel("p")
b.set_title("(b) beyond λ₂: no further p = 2 crossing", fontsize=10)
ins = b.inset_axes([0.08, 0.50, 0.40, 0.45])
ins.plot(arc[:, 1], arc[:, 2], "k.-", lw=0.8, ms=3)
ins.plot(lstar, pstar, "*", color="tab:orange", ms=10)
ins.set_xlim(lstar - 2e-5, lstar + 3.5e-4); ins.set_ylim(pstar - 9e-4, pstar + 2.5e-4)
ins.tick_params(labelsize=6); ins.set_title("log-periodic approach to the cusp (zoom)", fontsize=7)

c = ax[1, 0]
dd = np.logspace(np.log10(dall.min()), np.log10(3e-2), 500)
c.semilogx(np.concatenate([G[:, 0], dall]), np.concatenate([G[:, 2], pall]), "ko", ms=3.5,
           label="computed (layer-pinned continuation)")
c.semilogx(dd[dd <= 3e-3], model(parp, dd[dd <= 3e-3]), "tab:orange", lw=1.4,
           label=f"p* + δ[B + C cos(2τ ln δ) + D sin(2τ ln δ)],  2τ = {OMEGA:.4f} (no fit)")
c.axhline(2.0, color="tab:red", lw=0.8, ls="--")
c.axhline(pstar, color="tab:orange", lw=0.6, ls=":")
c.set_xlabel("sonic depth δ = min(1+λ+HΘ/ξ)"); c.set_ylabel("p")
c.set_ylim(1.999, 2.0145); c.legend(fontsize=7.5, loc="upper left")
c.set_title(f"(c) p* = {pstar:.6f} > 2;  free-ω fits: {parpf[5]:.3f} (p), {parlf[5]:.3f} (λ)", fontsize=10)

e = ax[1, 1]
prof = np.load("cusp_layer_profile.npy")
e.loglog(prof[:, 0], prof[:, 1] - 1, "b-", label="left flank")
e.loglog(prof[:, 0], prof[:, 2] - 1, "g-", label="right flank")
X = np.logspace(0.5, 3.8, 10)
e.loglog(X, 2.2 * X ** 0.5, "k--", lw=0.8, label="∝ X^{1/2}")
e.set_xlabel("X = |η − η_s| / w"); e.set_ylabel("den/δ − 1"); e.legend(fontsize=8)
e.set_title("(d) square-root cusp of the sonic factor; w ≈ 9.0 δ²", fontsize=10)
plt.tight_layout()
plt.savefig("fig_ladder.png", dpi=160)
print(f"p*={pstar:.10f} (rms {rmsp:.2e}), λ*={lstar:.10f} (rms {rmsl:.2e}); free ω: p {parpf[5]:.4f}, λ {parlf[5]:.4f}; "
      f"points: {len(dall)} (δ ≤ 3e−3: {sel.sum()})")
