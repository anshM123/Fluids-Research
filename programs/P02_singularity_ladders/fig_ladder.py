"""Main figure: the complete CCF self-similar branch, its three smooth profiles and its terminal cusp."""
import numpy as np, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cusp_fit import parse, fit, model, OMEGA

# ---------------- data ----------------
bm = np.load("branch_map.npy"); bm = bm[np.argsort(bm[:, 0])]
bm = bm[bm[:, 0] <= 1.19]                                  # mapped-grid natural continuation (reliable range)
large = []
for line in open("large_lambda_arc.log"):
    m = re.search(r"lam=([\d.]+) p=([\d.]+) min_den=([\d.]+) eta_min=(-?[\d.]+)", line)
    if m and float(m.group(4)) < -29.9:                   # sonic minimum at the origin: domain-converged part
        large.append((float(m.group(1)), float(m.group(2))))
large = np.array(large)
sac = np.load("sac_S2_branch.npy")                         # (λ, p, δ, layer, w, dλ/ds, N)
pin = parse("pin_C.log"); pin = pin[np.argsort(-pin[:, 0])]
lam_s = [1.180777662899, 0.6057337012, 0.471324227767]
d = pin[:, 0]
parp, rmsp = fit(d[d <= 1e-2], pin[d <= 1e-2, 2], free_omega=False, use_log=False)
parl, rmsl = fit(d[d <= 1e-2], pin[d <= 1e-2, 1], free_omega=False, use_log=False)
pstar, lstar = parp[0], parl[0]

fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
# (a) full branch
a = ax[0, 0]
a.plot(large[:, 0], large[:, 1], "k-", lw=1.4)
a.plot(bm[:, 0], bm[:, 1], "k-", lw=1.4)
a.plot(sac[:, 0], sac[:, 1], "k-", lw=1.4)
a.plot(pin[:, 1], pin[:, 2], "k-", lw=1.4)
a.axhline(2.0, color="tab:red", lw=0.8, ls="--")
for L, lab in zip(lam_s, ["λ₀ stable", "λ₁ (1 unstable)", "λ₂ (2 unstable)"]):
    a.plot(L, 2.0, "o", color="tab:blue", ms=7, zorder=5)
    a.annotate(lab, (L, 2.0), textcoords="offset points", xytext=(4, -16 if L > 1 else 8), fontsize=9)
a.plot(lstar, pstar, "*", color="tab:orange", ms=13, zorder=6)
a.annotate("terminal cusp", (lstar, pstar), textcoords="offset points", xytext=(6, 10), fontsize=9)
a.set_xscale("log"); a.set_xlabel("λ  (θ ~ (T−t)^λ)"); a.set_ylabel("local exponent p at the origin")
a.set_title("(a) the connected branch of CCF self-similar profiles")
a.set_xlim(0.44, 8.5); a.set_ylim(1.05, 2.25)
# (b) spiral end in the (λ, p) plane
b = ax[0, 1]
sel = bm[:, 0] < 0.5
b.plot(bm[sel, 0], bm[sel, 1], "k-", lw=1.2)
b.plot(sac[:, 0], sac[:, 1], "k-", lw=1.2)
b.plot(pin[:, 1], pin[:, 2], "k.-", lw=1.2, ms=5)
b.axhline(2.0, color="tab:red", lw=0.8, ls="--")
b.plot(lam_s[2], 2.0, "o", color="tab:blue", ms=7)
b.plot(lstar, pstar, "*", color="tab:orange", ms=13)
b.set_xlim(0.4530, 0.4745); b.set_ylim(1.9985, 2.0145)
b.set_xlabel("λ"); b.set_ylabel("p")
b.set_title("(b) no further p = 2 crossing: spiral into the cusp")
# (c) log-periodic approach
c = ax[1, 0]
dd = np.logspace(np.log10(d.min()), -2, 400)
c.semilogx(d, pin[:, 2], "ko", ms=5, label="computed (layer-pinned continuation)")
c.semilogx(dd, model(parp, dd), "tab:orange", lw=1.2,
           label=f"p* + δ[B + C cos(2τ ln δ) + D sin(2τ ln δ)],  2τ = {OMEGA:.4f}")
c.axhline(2.0, color="tab:red", lw=0.8, ls="--")
c.axhline(pstar, color="tab:orange", lw=0.6, ls=":")
c.set_xlabel("sonic depth δ = min(1+λ+HΘ/ξ)"); c.set_ylabel("p")
c.set_ylim(1.999, 2.013); c.legend(fontsize=8, loc="upper left")
c.set_title(f"(c) p* = {pstar:.6f} > 2,  λ* = {lstar:.6f}")
# (d) cusp structure of the sonic factor
e = ax[1, 1]
try:
    prof = np.load("cusp_layer_profile.npy")               # rows: X, den/δ (left), den/δ (right) at smallest δ
    e.loglog(prof[:, 0], prof[:, 1] - 1, "b-", label="left flank")
    e.loglog(prof[:, 0], prof[:, 2] - 1, "g-", label="right flank")
    X = np.logspace(0.5, 3.5, 10)
    e.loglog(X, 0.8 * X ** 0.5, "k--", lw=0.8, label="∝ X^{1/2}")
    e.set_xlabel("X = |η − η_s| / w"); e.set_ylabel("den/δ − 1"); e.legend(fontsize=8)
except FileNotFoundError:
    pass
e.set_title("(d) square-root cusp: den ≈ δ + k|η−η_s|^{1/2},  w ∝ δ²")
plt.tight_layout()
plt.savefig("fig_ladder.png", dpi=150)
print(f"p*={pstar:.10f} (rms {rmsp:.2e}), λ*={lstar:.10f} (rms {rmsl:.2e}); saved fig_ladder.png")
