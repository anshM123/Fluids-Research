"""Figure: linear stability of the eight smooth profiles.
(a) march-based method: number of real eigenvalues of T_μ above 1 (trivial time-translation mode removed) vs μ;
(b) unstable eigenvalues μ_k: brackets of (a) (grey) and the independent global-discretisation eigen-solver (red),
    counted only when found at ≥ 2 shifts with residual < 1e-10 (gstab_summary.robust) and inside a bracket of (a)."""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gstab_summary import robust
MU_ART = 0.045
lam = [1.9205593, 1.3990961, 1.2523487, 1.1842533, 1.1449857, 1.1194738, 1.1015817, 1.0883384]
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
a, b, c = ax
cols = plt.cm.viridis(np.linspace(0, 0.92, 8))
brackets = {}
for n in range(8):
    R = np.load(f"stab_l{n}.npy")
    keep = np.abs(R[:, 0] - 1.0) > 1e-6          # at μ = 1 the trivial eigenvalue is exactly 1 (±roundoff)
    mu, N = R[keep, 0], R[keep, 1] - (R[keep, 0] < 1.0)
    sel = mu >= MU_ART
    a.step(mu[sel], N[sel] + 0.04 * (n - 3.5), where='post', color=cols[n], label=f"λ$_{n}$ = {lam[n]:.4f}")
    brackets[n] = [(R[i + 1, 0], R[i, 0]) for i in range(len(R) - 1)
                   if R[i + 1, 1] > R[i, 1] and R[i + 1, 0] >= MU_ART and R[i, 0] < 0.99]
    for lo, hi in brackets[n]:
        b.plot([n, n], [lo, hi], color='0.75', lw=7, solid_capstyle='butt', zorder=1)
a.axvspan(0, MU_ART, color='0.92'); a.set_xlim(0, 1.05)
a.text(0.004, 6.6, 'origin-truncation\nartifacts', fontsize=7, color='0.4')
a.set_xlabel('μ  (perturbation ∝ e$^{μτ}$)'); a.set_ylabel('number of unstable eigenvalues > μ')
a.set_title('(a) march-based eigen-condition ν(μ) = 1', fontsize=9); a.legend(fontsize=6.5, ncol=2, loc='upper right')
for n in range(8):
    ev, both = robust(n)
    if not ev:
        continue
    ok = [e.real for e in ev if abs(e.imag) < 1e-6 and 0.05 < e.real < 0.99
          and any(lo - 0.01 <= e.real <= hi + 0.01 for lo, hi in brackets[n])]
    b.plot([n] * len(ok), ok, 'o', color='tab:red', ms=5, zorder=3, label='global eigen-solver' if n == 1 else None)
b.plot(range(8), [1.0] * 8, 'k_', ms=14, label='trivial mode μ = 1 (both methods)')
b.set_xlabel('profile n (continuation order)'); b.set_ylabel('unstable eigenvalue μ$_k$')
b.set_title('(b) unstable spectrum: brackets of (a) (grey), global eigen-solver (red)', fontsize=9)
b.set_ylim(0, 1.08); b.set_xticks(range(8)); b.legend(fontsize=7, loc='center left', bbox_to_anchor=(0.0, 0.66))
# (c) scaling of the lower eigenvalues with ε
for n in range(2, 8):
    ev, _ = robust(n)
    ok = sorted([e.real for e in (ev or []) if abs(e.imag) < 1e-6 and 0.05 < e.real < 0.99
                 and any(lo - 0.01 <= e.real <= hi + 0.01 for lo, hi in brackets[n])])
    if n == 7 and not ok:
        continue
    eps = (lam[n] - 1) / 2
    c.plot(np.arange(len(ok)), np.array(ok) / eps, 'o-', color=cols[n], ms=4, label=f"n = {n}")
c.set_xlabel('k (ordered from the smallest)'); c.set_ylabel('μ$_k$ / ε,   ε = (λ$_n$ − 1)/2')
c.set_title('(c) lower unstable eigenvalues scale with ε (observation)', fontsize=9); c.legend(fontsize=7)
plt.tight_layout(); plt.savefig('fig_bq_stability.png', dpi=150)
print("saved")
